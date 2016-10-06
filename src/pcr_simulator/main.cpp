#include <logger/logger.hpp>
#include <logger/log_writers.hpp>
#include <segfault_handler.hpp>
#include <perfcounter.hpp>
#include <boost/program_options.hpp>
#include <seqan/file.h>
#include <seqan/store.h>
#include <uchar.h>
#include <bitset>
#include "../ig_tools/utils/string_tools.hpp"

struct Options {
    struct PcrOptions {
        size_t cycles_count;
        double error_prob_first;
        double error_prob_last;
        double amplification_rate;
        double chimeras_rate;
        // 1 for barcode going with the left half, 2 for the right half, 3 for random choice
        size_t barcode_position;
    };

    std::string repertoire_file_path;
    std::string output_file_path;
    size_t barcode_length;
    PcrOptions pcr_options;
    size_t output_estimation_limit;
    std::string compress_repertoire;
};

Options parse_options(int argc, const char* const* argv) {
    Options options;
    Options::PcrOptions& pcrOptions = options.pcr_options;
    namespace po = boost::program_options;
    po::options_description cmdline_options("Command line options");
    cmdline_options.add_options()
            ("input-file,i", po::value<std::string>(&options.repertoire_file_path)->required(),
             "name of the input repertoire file (FASTA|FASTQ)")
            ("output-file,o", po::value<std::string>(&options.output_file_path)->required(),
             "output file for final repertoire")
            ("umi-length,l", po::value<size_t>(&options.barcode_length)->default_value(14),
             "length of generated barcodes (defaults to 14)")
            ("pcr-cycles,c", po::value<size_t>(&pcrOptions.cycles_count)->default_value(25),
             "number of PCR cycles to simulate (defaults to 25)")
            ("pcr-error1,e", po::value<double>(&pcrOptions.error_prob_first)->required(),
             "probability of PCR error on the first PCR cycle")
            ("pcr-error2,E", po::value<double>(&pcrOptions.error_prob_last)->required(),
             "probability of PCR error on the last PCR cycle (probability on the other cycles are interpolated)")
            ("pcr-rate,r", po::value<double>(&pcrOptions.amplification_rate)->default_value(0.1),
             "probability for each molecule to be amplified on each PCR cycle (defaults to 0.1)")
            ("chimeras-rate,k", po::value<double>(&pcrOptions.chimeras_rate)->default_value(0.001),
             "share of chimeric reads appering on each cycle (defaults to 0.001)")
            ("output-limit,m", po::value<size_t>(&options.output_estimation_limit)->default_value(static_cast<size_t>(1e8)),
             "the program will exit if expected number of reads exceeds this parameter (defaults to 100,000,000)")
            ("barcode-position,b", po::value<size_t>(&pcrOptions.barcode_position)->default_value(3),
             "indicator of barcode position in the read, used for chimeras simulation. "
             "1 for barcode going with the left half, 2 for the right half, 3 for random choice (defaults to 3)")
            ("compressed-path,w", po::value<std::string>(&options.compress_repertoire)->default_value(""),
             "for Ig Simulator output: report compressed repertoire either")
            ;
    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
    po::notify(vm);
    return options;
}

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

void read_repertoire(Options options, std::vector<seqan::CharString>& read_ids, std::vector<seqan::Dna5String>& reads, std::vector<size_t>& read_to_compressed) {
    seqan::SeqFileIn reads_file(options.repertoire_file_path.c_str());
    readRecords(read_ids, reads, reads_file);
    double exp_reads_count = static_cast<double>(reads.size()) *
                             pow(1.0 + options.pcr_options.amplification_rate + options.pcr_options.chimeras_rate, static_cast<double>(options.pcr_options.cycles_count));
    VERIFY(exp_reads_count <= options.output_estimation_limit);

    if (!options.compress_repertoire.empty()) {
        size_t current = 0;
        for (size_t i = 0; i < read_ids.size(); i ++) {
            auto id = seqan_string_to_string(read_ids[i]);
            if (boost::algorithm::ends_with(id, "_copy_1")) {
                current ++;
            }
            read_to_compressed.push_back(current - 1);
        }
    }
}

seqan::Dna5String generate_barcode(size_t length) {
    seqan::Dna5String barcode;
    for (size_t i = 0; i < length; i ++) {
        barcode += std::rand() % 4;
    }
    return barcode;
}

std::vector<seqan::Dna5String> generate_barcodes(size_t count, size_t barcode_length) {
    std::vector<seqan::Dna5String> barcodes(count);
    for (size_t i = 0; i < count; i ++) {
        barcodes[i] = generate_barcode(barcode_length);
    }
    return barcodes;
}

void report_average_error_rate(std::vector<size_t> error_count) {
    size_t total_errors = 0;
    size_t interesting_reads = 0;
    for (size_t errors : error_count) {
        if (errors < std::numeric_limits<size_t>::max()) {
            total_errors += errors;
            interesting_reads ++;
        }
    }
    INFO("Average amount of errors per read is " << ((double) total_errors / (double) interesting_reads) << " (total " << total_errors << " in " << interesting_reads << " of " << error_count.size() << " reads)");
}

void amplify(std::vector<seqan::Dna5String>& reads, std::vector<seqan::Dna5String>& barcodes, std::vector<seqan::CharString>& ids, std::vector<size_t>& error_count,
             double pcr_error_prob, Options::PcrOptions options, std::minstd_rand0& random_engine, std::vector<size_t>& read_to_compressed) {
    size_t size = reads.size();
    for (size_t read_idx = 0; read_idx < size; read_idx ++) {
        if (std::rand() <= static_cast<double>(RAND_MAX) * options.amplification_rate) {
            seqan::Dna5String read = reads[read_idx];
            seqan::Dna5String barcode = barcodes[read_idx];
            size_t errors = error_count[read_idx];
            for (size_t pos = 0; pos < length(barcode) + length(read); pos ++) {
                if (std::rand() <= static_cast<double>(RAND_MAX) * pcr_error_prob) {
                    auto current_value = pos < length(barcode) ? barcode[pos] : read[pos - length(barcode)];
                    int new_value_candidate = std::rand() % 3;
                    int new_value = new_value_candidate + (new_value_candidate >= current_value ? 1 : 0);
                    pos < length(barcode) ? barcode[pos] : read[pos - length(barcode)] = new_value;
                    if (pos >= length(barcode)) errors ++;
                }
            }
            barcodes.push_back(barcode);
            reads.push_back(read);
            ids.emplace_back(std::to_string(reads.size()) +
                                     (seqan_string_to_string(ids[read_idx]).find("chimera") == std::string::npos ? "" : "_chimera") +
                                     "_mutated_from_" +
                                     std::to_string(read_idx));
            error_count.push_back(errors);
            read_to_compressed.push_back(read_to_compressed[read_idx]);
        }
    }
    std::uniform_int_distribution<size_t> read_distribution(0, size - 1);
    std::uniform_int_distribution<size_t> barcode_position_distribution(1, std::bitset<2>(options.barcode_position).count());
    for (size_t chimera = 0; chimera < static_cast<double>(size) * options.chimeras_rate; chimera ++) {
        size_t left_idx = read_distribution(random_engine);
        size_t right_idx = read_distribution(random_engine);
        size_t barcode_idx = (options.barcode_position & 1) && barcode_position_distribution(random_engine) == 1 ? left_idx : right_idx;
        barcodes.push_back(barcodes[barcode_idx]);
        std::string left = seqan_string_to_string(reads[left_idx]);
        left = left.substr(0, left.length() / 2);
        std::string right = seqan_string_to_string(reads[right_idx]);
        right = right.substr(right.length() / 2);
        reads.push_back(seqan::Dna5String(left + right));
        error_count.push_back(std::numeric_limits<size_t>::max());
        read_to_compressed.push_back(std::numeric_limits<size_t>::max());
        ids.emplace_back(std::to_string(reads.size()) + "_chimera_from_" + std::to_string(left_idx) + "_" + std::to_string(right_idx));
    }
    VERIFY(reads.size() == barcodes.size() && reads.size() == ids.size() && reads.size() == read_to_compressed.size() && reads.size() == error_count.size());

//    report_average_error_rate(error_count);
}

void simulate_pcr(std::vector<seqan::Dna5String>& reads, std::vector<seqan::Dna5String>& barcodes, std::vector<seqan::CharString>& ids, Options::PcrOptions options,
                  std::minstd_rand0& random_engine, std::vector<size_t>& read_to_compressed) {
    std::vector<size_t> error_count(reads.size(), 0);
    for (size_t i = 0; i < options.cycles_count; i ++) {
        double pcr_error_prob = options.error_prob_first + (options.error_prob_last - options.error_prob_first) * static_cast<double>(i) / static_cast<double>(options.cycles_count - 1);
        amplify(reads, barcodes, ids, error_count, pcr_error_prob, options, random_engine, read_to_compressed);
    }

    report_average_error_rate(error_count);
}

void write_repertoire(const std::vector<seqan::Dna5String>& reads, const std::vector<seqan::Dna5String>& barcodes, const std::vector<seqan::CharString>& ids, std::vector<size_t> perm, std::string file_name) {
    seqan::SeqFileOut output_file(file_name.c_str());
    for (size_t i = 0; i < reads.size(); i ++) {
        std::stringstream sstr;
        sstr << ids[perm[i]] << "_UMI:" << barcodes[perm[i]];
        seqan::writeRecord(output_file, sstr.str(), reads[perm[i]]);
    }
}

void write_compressed(std::string file_name, std::vector<seqan::CharString>& read_ids, std::vector<seqan::Dna5String>& reads, std::vector<size_t>& read_to_compressed) {
    std::map<size_t, size_t> id_to_count;
    for (size_t id : read_to_compressed) {
        id_to_count[id] ++;
    }

    std::vector<seqan::CharString> compressed_ids;
    std::vector<seqan::Dna5String> compressed_reads;
    size_t current = 0;
    for (size_t i = 0; i < read_ids.size(); i ++) {
        read_to_compressed.push_back(current);
        auto id = seqan_string_to_string(read_ids[i]);
        if (boost::algorithm::ends_with(id, "_copy_1")) {
            current ++;
            compressed_ids.emplace_back((boost::format("cluster___%d___size___%d") % (current - 1) % id_to_count[current - 1]).str());
            compressed_reads.push_back(reads[i]);
        }
    }
    seqan::SeqFileOut compressed_file(file_name.c_str());
    seqan::writeRecords(compressed_file, compressed_ids, compressed_reads);
}

int main(int argc, const char* const* argv) {
    perf_counter pc;
    segfault_handler sh;
    create_console_logger();
    std::srand(495634);
    auto random_engine = std::default_random_engine(8356);

    const Options& options = parse_options(argc, argv);

    std::vector<seqan::CharString> original_ids;
    std::vector<seqan::Dna5String> original_reads;
    std::vector<size_t> read_to_compressed;
    read_repertoire(options, original_ids, original_reads, read_to_compressed);
    auto reads(original_reads);
    std::vector<seqan::CharString> ids(reads.size());
    for (size_t i = 0; i < reads.size(); i ++) {
        ids[i] = "original_" + std::to_string(i);
    }
    auto barcodes = generate_barcodes(reads.size(), options.barcode_length);

    simulate_pcr(reads, barcodes, ids, options.pcr_options, random_engine, read_to_compressed);

    std::vector<size_t> perm(reads.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), random_engine);

    write_repertoire(reads, barcodes, ids, perm, options.output_file_path);

    if (!options.compress_repertoire.empty()) {
        write_compressed(options.compress_repertoire, original_ids, original_reads, read_to_compressed);
    }
}