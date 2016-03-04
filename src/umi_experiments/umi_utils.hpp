#pragma once

#include <verify.hpp>
#include <seqan/seq_io.h>
#include "../ig_tools/utils/string_tools.hpp"

class Umi {
public:
    explicit Umi(const seqan::Dna5String& umi) : umi_(umi) {}
    Umi() {}

    bool operator==(const Umi &other) const { return umi_ == other.umi_; }

    seqan::Dna5String GetString() const { return umi_; }

private:
    seqan::Dna5String umi_;
};

namespace std {
    template<>
    struct hash<seqan::Dna5String> {
        size_t operator()(const seqan::Dna5String& str) const {
            size_t h = 0;
            for (auto c : str) {
                h = h * 31 + seqan::ordValue(c);
            }
            return h;
        }
    };

    template<>
    struct hash<Umi> {
        size_t operator()(const Umi& umi) const {
            size_t h = hash<seqan::Dna5String>()(umi.GetString());
            return h;
        }
    };
}

void extract_barcodes_from_read_ids(const std::vector<seqan::CharString>& input_ids, std::vector<seqan::Dna5String>& umis,
                                    std::vector<seqan::DnaQString>& umi_quals);

void group_reads_by_umi(const std::vector<seqan::Dna5String>& umis, std::unordered_map<Umi, std::vector<size_t> >& umi_to_reads);
