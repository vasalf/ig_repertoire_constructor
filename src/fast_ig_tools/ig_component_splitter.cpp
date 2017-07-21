#include <string>
#include <limits>
#include <cstddef>
#include <iostream>
#include <functional>
#include <memory>
#include <sstream>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <build_info.hpp>
#include "utils.hpp"
#include "fast_ig_tools.hpp"
#include "ig_final_alignment.hpp"

#include <seqan/seq_io.h>

struct IgComponentSplitterConfig {
    // Paths
    std::string reads_file;
    std::string rcm_file;
    std::string output_file;
    std::string output_rcm_file;
    std::string config_file;

    // Launch options
    int nthreads = 4;
    std::size_t max_votes = std::numeric_limits<std::size_t>::max() / 2;
    bool discard = false;
    bool recursive = false;
    bool flu = false;
    bool allow_unassigned = false;

    void Write() {
        INFO("reads_file=" << reads_file);
        INFO("rcm_file=" << rcm_file);
        INFO("output_file=" << output_file);
        INFO("output_rcm_file=" << output_rcm_file);
        INFO("config_file=" << config_file);
        INFO("nthreads=" << nthreads);
        INFO("max_votes=" << max_votes);
        INFO("discard=" << std::boolalpha << discard);
        INFO("recursive=" << std::boolalpha << recursive);
        INFO("flu=" << std::boolalpha << flu);
        INFO("allow_unassigned=" << std::boolalpha << allow_unassigned);
    }
    
    IgComponentSplitterConfig() {}
    
    void ReadOptions(int argc, char **argv) {
        namespace po = boost::program_options;

        // Declare a group of options that are
        // only allowed on a command line
        po::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version string")
            ("help,h", "produce help message")
            ("config,c", po::value<std::string>(&config_file)->default_value(config_file),
             "name of a file of a configuration")
            ("input-file,i", po::value<std::string>(&reads_file),
             "name of the input file (FASTA|FASTQ)")
            ("output-file,o", po::value<std::string>(&output_file)->default_value(output_file),
             "output file for final repertoire")
            ("rcm-file,R", po::value<std::string>(&rcm_file)->default_value(rcm_file),
             "input RCM-file")
            ("output-rcm-file,M", po::value<std::string>(&output_rcm_file)->default_value(output_rcm_file),
             "output RCM-file");

        // Declare a group of options that are
        // allowed both on command line and in
        // config file
        po::options_description config("Configuration");
        config.add_options()
            ("threads,t", po::value<int>(&nthreads)->default_value(nthreads),
             "the number of parallel threads")
            ("max-votes,V", po::value<std::size_t>(&max_votes)->default_value(max_votes),
             "max secondary votes threshold")
            ("discard,D", po::value<bool>(&discard)->default_value(discard),
             "whether to discard secondary votes")
            ("recursive,C", po::value<bool>(&recursive)->default_value(recursive),
             "whether to perform recursive splitting")
            ("flu,F", po::value<bool>(&flu)->default_value(flu),
             "Use FLU preset")
            ("allow-unassigned,A", po::value<bool>(&allow_unassigned)->default_value(allow_unassigned),
             "Allow unassigned reads in RCM");

        // Declare a group of options that are
        // allowed both on command line and in
        // config file but are not shown to user.
        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("help-hidden", "show all options, including developers' ones");
        
        po::options_description cmdline_options("All command line options");
        cmdline_options
            .add(generic)
            .add(config)
            .add(hidden);
        po::options_description config_file_options;
        config_file_options
            .add(config)
            .add(hidden);
        po::options_description visible_options("Allowed options");
        visible_options
            .add(generic)
            .add(config);

        po::positional_options_description positional_options;
        positional_options
            .add("input-file", 1);
        
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv)
                  .options(cmdline_options)
                  //.positional(positional_options)
                  .run(), vm);
        po::notify(vm);

        if (vm.count("help-hidden")) {
            std::cout << cmdline_options << std::endl;
            std::exit(0);
        }

        if (vm.count("help")) {
            std::cout << visible_options << std::endl;
            std::exit(0);
        }

        if (vm.count("version")) {
            std::cout << "IG Component Splitter, part of IgReC version "
                      << build_info::version
                      << "; git version: "
                      << build_info::git_hash7
                      << std::endl;
            std::exit(0);
        }

        if (vm.count("config-file")) {
            std::ifstream config(vm["config-file"].as<std::string>());
            if (!config) {
                std::cout << "could not open config file " << vm["config-file"].as<std::string>() << std::endl;
                std::exit(1);
            } else {
                po::store(po::parse_config_file(config, config_file_options), vm);
                po::notify(vm);
                // reparse cmd line again for update config defaults
                po::store(po::command_line_parser(argc, argv)
                          .options(cmdline_options)
                          .positional(positional_options)
                          .run(), vm);
                po::notify(vm);
            }
        }
    }
};

struct Read {
    seqan::Dna5String read;
    seqan::CharString id;
    Read(seqan::Dna5String& read_, seqan::CharString id_)
        : read(read_), id(id_) {}
};

template<class StringType>
class Component : public std::vector<std::shared_ptr<Read> > {
    StringType id_;
    seqan::Dna5String consensus_;
public:
    Component() {}
    
    Component(const StringType& id)
        : std::vector<std::shared_ptr<Read> >(),
        id_(id) {}

    virtual ~Component() = default;
    
    const StringType& GetId() {
        return id_;
    }

    bool operator<(const Component &other) const {
        return id_ < other.id_;
    }

    seqan::Dna5String & Consensus() {
        return consensus_;
    }
};

template<typename key, typename value>
class StdMapTag : public std::map<key, value> {};

template<typename key, typename value>
class StdUnorderedMapTag : public std::unordered_map<key, value> {};

template<template<typename,typename> class map_type>
class InputData {
    map_type<std::string, Component<std::string> > components_;
    map_type<std::string, std::string> rcm_;
    std::vector<std::shared_ptr<Read> > reads_;

    void read_rcm_file(std::string file_name) {
        std::ifstream rcm(file_name);
        if (!rcm) {
            // TODO: throw an exception
            ERROR("Could not open RCM file " << file_name);
            std::exit(1);
        }

        std::string id, target;
        std::string line;
        while (std::getline(rcm, line)) {
            std::vector<std::string> strs;
            boost::split(strs, line, boost::is_any_of("\t"));
            for (auto &s : strs) {
                boost::trim(s);
            }
            if (strs.size() == 2 && strs[0] != "" && strs[1] != "") {
                rcm_[strs[0]] = strs[1];
            }
        }
    }

    void add_read_record(std::string s_id, Read&& r) {
        reads_.push_back(std::make_shared<Read>(r));
        if (!components_.count(s_id))
            components_[s_id] = Component<std::string>(s_id);
        components_[s_id].push_back(reads_.back());
    }
    
public:
    InputData() {}

    void ReadData(const IgComponentSplitterConfig &cfg) {
        INFO("Input files: " << cfg.reads_file << ", " << cfg.rcm_file);

        std::vector<seqan::Dna5String> input_reads;
        std::vector<seqan::CharString> input_ids;
        // TODO: throw an exception on failure
        seqan::SeqFileIn seqinput(cfg.reads_file.c_str());

        INFO("Reading input reads starts");
        seqan::readRecords(input_ids, input_reads, seqinput);
        INFO(input_reads.size() << " reads were extracted from " << cfg.reads_file);

        INFO("Reading read-cluster map starts");
        read_rcm_file(cfg.rcm_file);

        std::size_t assigned_reads = 0;
        for (std::size_t i = 0; i != input_reads.size(); i++) {
            std::string id = seqan::toCString(input_ids[i]);
            if (rcm_.count(id))
                assigned_reads++;
            else
                TRACE("Read " << id <<  " not found in RCM " << cfg.rcm_file);
            add_read_record(rcm_[id], Read(input_reads[i], input_ids[i]));
        }

        INFO(assigned_reads << " over " << input_reads.size() << " reads assigned");
        if (assigned_reads < input_reads.size() && !cfg.allow_unassigned) {
            ERROR(input_reads.size() - assigned_reads << " unassigned reads in RCM " << cfg.rcm_file);
            ERROR("Unassigned reads are not allowed");
            ERROR("Pass option '--allow-unassigned=true' to allow");
            // TODO: exception
            std::exit(1);
        }

        INFO(components_.size() << " clusters were extracted from " << cfg.rcm_file);

        std::size_t max_component_size = 0;
        for (const auto& kv : components_) {
            max_component_size = std::max(max_component_size, kv.second.size());
        }
        INFO("Size of maximal cluster: " << max_component_size);
    }

    map_type<std::string, Component<std::string> > & Components() {
        return components_;
    }

    std::vector<Read> & Reads() {
        return reads_;
    }
};

template<template<typename,typename> class map_type>
class SplitterAlgorithm {
protected:
    IgComponentSplitterConfig &cfg_;
    InputData<map_type> &input_;
public:
    SplitterAlgorithm(IgComponentSplitterConfig &cfg,
                      InputData<map_type> &input)
        : cfg_(cfg), input_(input) {}

    virtual ~SplitterAlgorithm() = default;
    
    virtual std::vector<Component<std::string> > SplitComponent(Component<std::string> &component) = 0;
};

template<template<typename,typename> class map_type>
class DummySplitterAlgorithm : public SplitterAlgorithm<map_type> {
public:
    DummySplitterAlgorithm(IgComponentSplitterConfig &cfg,
                           InputData<map_type> &input)
        : SplitterAlgorithm<map_type>(cfg, input) {}

    virtual ~DummySplitterAlgorithm() {}

    virtual std::vector<Component<std::string> > SplitComponent(Component<std::string> &component) override {
        std::vector<seqan::Dna5String> reads;

        std::transform(component.begin(), component.end(),
                       std::back_inserter(reads),
                       [](std::shared_ptr<Read>& r) -> seqan::Dna5String {
                           return r->read;
                       });

        component.Consensus() = consensus(reads);
        return { component };
    }
};

// TODO: write a parallel processor here
template<template<typename,typename> class map_type>
class IgComponentSplitterProcessor {
    IgComponentSplitterConfig &cfg_;
    InputData<map_type> &input_;
    
public:
    IgComponentSplitterProcessor(IgComponentSplitterConfig &cfg,
                                 InputData<map_type> &input)
        : cfg_(cfg), input_(input) {}

    void Process() {
        std::vector<Component<std::string> > components;
        for (auto &p : input_.Components())
            components.push_back(p.second);
        std::sort(components.begin(), components.end());

        seqan::SeqFileOut seq_output(cfg_.output_file.c_str());
        std::ofstream out_rcm(cfg_.output_rcm_file.c_str());

        /*
        std::unique_ptr<SplitterAlgorithm<map_type> > splitter
            = std::make_unique<DummySplitterAlgorithm<map_type> >(cfg_, input_);
        */
        std::unique_ptr<SplitterAlgorithm<map_type> > splitter(new DummySplitterAlgorithm<map_type>(cfg_, input_));

        for (auto &component : components) {
            auto result = splitter->SplitComponent(component);
            for (std::size_t i = 0; i != result.size(); i++) {
                std::ostringstream cluster_id;
                cluster_id << component.GetId();
                if (result.size() > 1)
                    cluster_id << "X" << i;
                
                std::ostringstream id;
                id << "cluster___"
                   << cluster_id.str()
                   << "___size___"
                   << result[i].size();

                seqan::writeRecord(seq_output, id.str(), result[i].Consensus());
                for (auto& read : result[i]) {
                    std::string read_id = seqan::toCString(read->id);
                    out_rcm << read_id << "\t" << cluster_id.str() << "\n";
                }
            }
        }

        INFO("Final repertoire was written to " << cfg_.output_file);
        INFO("Final RCM was written to " << cfg_.output_rcm_file);
    }
};

int main(int argc, char *argv[]) {
    segfault_handler sh;
    perf_counter pc;
    create_console_logger("");

    INFO("Command line: " << join_cmd_line(argc, argv));

    IgComponentSplitterConfig cfg;

    try {
        cfg.ReadOptions(argc, argv);
    } catch(boost::program_options::error &e) {
        std::cout << "Parser error: " << e.what() << std::endl;
        return 1;
    } catch(std::exception &e) {
        std::cout << e.what() << std::endl;
        return 1;
    }
    cfg.Write();

    InputData<StdUnorderedMapTag> input;
    input.ReadData(cfg);

    IgComponentSplitterProcessor<StdUnorderedMapTag> processor(cfg, input);
    processor.Process();

    INFO("Running time: " << running_time_format(pc));
    
    return 0;
}
