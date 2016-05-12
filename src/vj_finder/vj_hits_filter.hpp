#include "vj_finder_config.hpp"
#include "vj_alignment_structs.hpp"

namespace vj_finder {
    class VJHitsFilter {
    public:
        // returns true if VJ Hits are bad, otherwise returns false
        virtual bool Filter(const VJHits &vj_hits) const = 0;

        ~VJHitsFilter() { }
    };

    class EmptyHitsFilter : public VJHitsFilter {
    public:
        bool Filter(const VJHits &vj_hits) const {
            bool empty = vj_hits.NumJHits() == 0 or vj_hits.NumVHits() == 0;
            std::cout << "vj_hits.NumJHits() == 0 or vj_hits.NumVHits() == 0: " << empty << std::endl;
            return vj_hits.NumJHits() == 0 or vj_hits.NumVHits() == 0;
        }
    };

    class LeftCoverageFilter : public VJHitsFilter {
        size_t left_uncovered_limit_;

    public:
        LeftCoverageFilter(size_t left_uncovered_limit) :
                left_uncovered_limit_(left_uncovered_limit) { }

        bool Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.GetVHitByIndex(0).LeftUncovered() > left_uncovered_limit_;
            std::cout << "vj_hits.GetVHitByIndex(0).LeftUncovered() > left_uncovered_limit_: " << bad << std::endl;
            return vj_hits.GetVHitByIndex(0).LeftUncovered() > left_uncovered_limit_;
        }
    };

    class RightCoverageFilter : public VJHitsFilter {
        size_t right_uncovered_limit_;

    public:
        RightCoverageFilter(size_t right_uncovered_limit) :
                right_uncovered_limit_(right_uncovered_limit) { }

        bool Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.GetJHitByIndex(0).RightUncovered() > right_uncovered_limit_;
            std::cout << "vj_hits.GetJHitByIndex(0).RightUncovered() > right_uncovered_limit_: " << bad << std::endl;
            return vj_hits.GetJHitByIndex(0).RightUncovered() > right_uncovered_limit_;
        }
    };

    class VSegmentLengthFilter : public VJHitsFilter {
        size_t min_v_segment_length_;

    public:
        VSegmentLengthFilter(size_t min_v_segment_length) :
                min_v_segment_length_(min_v_segment_length) { }

        bool Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.GetVHitByIndex(0).SegmentLength() < min_v_segment_length_;
            std::cout << "vj_hits.GetVHitByIndex(0).SegmentLength() < min_v_segment_length_: " << bad << std::endl;
            return vj_hits.GetVHitByIndex(0).SegmentLength() < min_v_segment_length_;
        }
    };

    class JSegmentLengthFilter : public VJHitsFilter {
        size_t min_j_segment_length_;

    public:
        JSegmentLengthFilter(size_t min_j_segment_length) :
                min_j_segment_length_(min_j_segment_length) { }

        bool Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.GetJHitByIndex(0).SegmentLength() < min_j_segment_length_;
            std::cout << "vj_hits.GetJHitByIndex(0).SegmentLength() < min_j_segment_length_: " << bad << std::endl;
            return vj_hits.GetJHitByIndex(0).SegmentLength() < min_j_segment_length_;
        }
    };

    class AlignedSegmentLengthFilter : public VJHitsFilter {
        size_t min_aligned_segment_length_;

    public:
        AlignedSegmentLengthFilter(size_t min_aligned_segment_length) :
                min_aligned_segment_length_(min_aligned_segment_length) { }

        bool Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.AlignedSegmentLength() < min_aligned_segment_length_;
            std::cout << "vj_hits.AlignedSegmentLength (" << vj_hits.AlignedSegmentLength() << ") < min_aligned_segment_length_: " << bad << std::endl;
            return vj_hits.AlignedSegmentLength() < min_aligned_segment_length_;
        }
    };

    class CustomVjHitsFilter {
        std::vector<std::shared_ptr<VJHitsFilter> > vj_filters_;

    public:
        void Add(std::shared_ptr<VJHitsFilter> vj_filter);

        bool Filter(const VJHits &vj_hits) const;
    };

    class VersatileVjFilter {
        const vjf_config::AlgorithmParams::FilteringParams &filtering_params_;
        CustomVjHitsFilter custom_vj_filter_;

        void InitializeCustomVjFinder();

    public:
        VersatileVjFilter(const vjf_config::AlgorithmParams::FilteringParams &filtering_params) :
                filtering_params_(filtering_params) {
            InitializeCustomVjFinder();
        }

        bool Filter(const VJHits &vj_hits) const { return custom_vj_filter_.Filter(vj_hits); }
    };
}