#include "match_threshold_alignment_quality_checker.hpp"

namespace vdj_labeler {

bool MatchThresholdAlignmentQualityChecker::AlignmentIsGood(
        const alignment_utils::ImmuneGeneReadAlignment &ig_gene_alignment) const
{
    int cnt_following_matches = 0;
    auto row1 = seqan::row(ig_gene_alignment.Alignment(), 0);
    auto row2 = seqan::row(ig_gene_alignment.Alignment(), 1);
    for (auto it_row1 = begin(row1), it_row2 = begin(row2);
         it_row1 != end(row1);
         ++it_row1, ++it_row2)
    {
        if (*it_row1 == *it_row2)
            cnt_following_matches++;
        else {
            cnt_following_matches = 0;
            continue;
        }
        if (cnt_following_matches == this->normalized_score_threshold_)
            return true;
    }
    return false;
}

} // End namespace vdj_labeler
