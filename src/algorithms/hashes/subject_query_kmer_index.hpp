#pragma once

#include "polyhashes.hpp"
#include "kmer_index_primitives.hpp"

#include <unordered_map>
#include <memory>

namespace algorithms {
    // SubjectDatabase - type of sequence storage
    // StringType - type of sequences
    // QueryType - type of query object
    // default implementation works for any type of simple collections that support [] and .size(), e.g., std::vector
    // and StringType and QueryType are standard types of sequences, e.g., c-style, std::string and seqan sequences
    template<typename SubjectDatabase, typename StringType>
    class SubjectQueryKmerIndex {
    public:
        struct SubjectPosition {
            size_t subject_index;
            size_t position;
        };

    private:
        // input parameters
        const SubjectDatabase & db_;
        size_t k_;
        KmerIndexHelper<SubjectDatabase, StringType>& kmer_index_helper_;

        // inner structure
        std::map<size_t, std::unique_ptr<SubjectPosition[]>> kmer_query_pos_map_;
        std::map<size_t, int> kmer_query_pos_sizes_;
        std::map<size_t, int> kmer_query_pos_index_;

        void UpdateSize(const std::vector<size_t>& str_hashes, size_t str_length, size_t str_index) {
            for (size_t start = 0; start + k_ <= str_length; ++start) {
                kmer_query_pos_sizes_[str_hashes[start]]++;
            }
	    }

        void UpdateMap(const std::vector<size_t>& str_hashes, size_t str_length, size_t str_index) {
            for (size_t start = 0; start + k_ <= str_length; ++start) {
                int & i = kmer_query_pos_index_[str_hashes[start]];
                kmer_query_pos_map_[str_hashes[start]][i] = {str_index, start};
                i++;
            }
        }

        void Initialize() {
            for (size_t j = 0; j < kmer_index_helper_.GetDbSize(); ++j) {
                auto s = kmer_index_helper_.GetDbRecordByIndex(j);
                auto hashes = polyhashes(s, k_);
                UpdateSize(hashes, kmer_index_helper_.GetStringLength(s), j);
            }
            for (auto p : kmer_query_pos_sizes_) {
                size_t size = p.second;
                kmer_query_pos_map_[p.first] = std::unique_ptr<SubjectPosition[]>(new SubjectPosition[size]);
            }
            for (size_t j = 0; j < kmer_index_helper_.GetDbSize(); ++j) {
                auto s = kmer_index_helper_.GetDbRecordByIndex(j);
                auto hashes = polyhashes(s, k_);
                UpdateMap(hashes, kmer_index_helper_.GetStringLength(s), j);
            }
        }
    public:
        SubjectQueryKmerIndex(const SubjectDatabase &db,
                              size_t k,
                              KmerIndexHelper<SubjectDatabase, StringType>& kmer_index_helper) :
                db_(db),
                k_(k),
                kmer_index_helper_(kmer_index_helper) {
            Initialize();
        }

        const SubjectDatabase & Db() const { return db_; }

        bool SubjectsContainKmer(size_t kmer) const {
            return kmer_query_pos_map_.find(kmer) != kmer_query_pos_map_.end();
        }

        const std::unique_ptr<SubjectPosition[]>& GetSubjectPositions(size_t kmer) const {
            VERIFY_MSG(SubjectsContainKmer(kmer), "Subjects do not contain kmer " << kmer);
            return kmer_query_pos_map_.at(kmer);
        }

        size_t NumSubjects() const { return kmer_index_helper_.GetDbSize(); }

        size_t k() const { return k_; }

        SubjectKmerMatches GetSubjectKmerMatchesForQuery(const StringType &query_str) const {
            SubjectKmerMatches subj_kmer_matches(NumSubjects());
            auto query_hashes = polyhashes(query_str, k());
            size_t start = 0;
            size_t finish = kmer_index_helper_.GetStringLength(query_str);
            for(size_t j = start; j < finish - k() + 1; ++j) {
                auto kmer = query_hashes[j];
                if(!SubjectsContainKmer(kmer))
                    continue;
                const auto& subj_pos = GetSubjectPositions(kmer);
                //for (const auto &p : subj_it.second) {
                //for(auto it = subj_pos.begin(); it != subj_pos.end(); it++) {
                for (int i = 0; i < kmer_query_pos_sizes_.at(kmer); i++) {
                    auto it = &subj_pos[i];
                    size_t kmer_pos_in_query = j;
                    size_t kmer_pos_in_subject = it->position; // p.position;
                    subj_kmer_matches.Update(it->subject_index, {static_cast<int>(kmer_pos_in_subject),
                                                             static_cast<int>(kmer_pos_in_query)});
                }
            }
            return subj_kmer_matches;
        }
    };
}
