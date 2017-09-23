#include "utils.hpp"
#include <fast_ig_tools.hpp>
#include <ig_final_alignment.hpp>
#include <ig_matcher.hpp>
#include "splitter_algorithm.hpp"

namespace ig_component_splitter {

template<class T>
struct SeqanLengthCompare {
    bool operator()(const T& first, const T& second) {
        return seqan::length(first) < seqan::length(second);
    }
};

struct PositionVote {
    size_t majory_votes;
    size_t majory_letter;
    size_t secondary_votes;
    size_t secondary_letter;
    size_t position;
    bool operator<(const PositionVote &other) const {
        return secondary_votes < other.secondary_votes;
    }
};


void LegacySplitterAlgorithm::split_component(Component &component,
                                              std::vector<Component> &res,
                                              size_t max_votes) {
    if (component.empty()) {
        return;
    } else if (component.size() == 1) {
        component.consensus = component[0]->read;
        res.push_back(component);
        return;
    }

    // What?
    if (max_votes == 0) {
        max_votes = std::numeric_limits<size_t>::max() / 2;
    }

    const std::vector<seq_string_t>& reads_in_component = component.GetReadSequences();
        
    profile_string_t profile;
    size_t len = seqan::length(*max_element(reads_in_component.begin(), reads_in_component.end(),
                                            SeqanLengthCompare<seqan::Dna5String>()));
    seqan::resize(profile, len);

    for (const seq_string_t& read : reads_in_component) {
        for (size_t j = 0; j != length(read); j++) {
            profile[j].count[seqan::ordValue(read[j])]++;
        }
    }

    size_t min_len = seqan::length(*min_element(reads_in_component.begin(), reads_in_component.end(),
                                                SeqanLengthCompare<seqan::Dna5String>()));
    std::vector<PositionVote> votes(min_len);
    for (size_t i = 0; i != min_len; i++) {
        std::array<std::pair<size_t, size_t>, seqan::ValueSize<char_t>::VALUE> letters;
        for (size_t k = 0; k < seqan::valueSize<char_t>(); k++) {
            letters[k] = std::make_pair(profile[i].count[k], k);
        }


        std::nth_element(letters.begin(), letters.begin() + 1, letters.end(),
                         std::greater<std::pair<size_t, size_t> >());

        PositionVote vote;
        std::tie(vote.majory_votes, vote.majory_letter) = letters[0];
        std::tie(vote.secondary_votes, vote.secondary_letter) = letters[1];
        vote.position = i;
        votes[i] = vote;
    }

    PositionVote maximal_mismatch = *std::max_element(votes.begin(), votes.end());
    VERIFY(maximal_mismatch.majory_votes >= maximal_mismatch.secondary_votes);

    TRACE("VOTES: " << maximal_mismatch.majory_votes
          << "/" << maximal_mismatch.secondary_votes
          << " POSITION: " << maximal_mismatch.position);

    bool do_split = false;
    size_t mmsv = maximal_mismatch.secondary_votes;

    // Uhm.
    if (this->cfg().flu) {
        const double component_size_coefficient = -0.0064174097073423269;
        const double secondary_votes_coefficient = 0.79633984048973583;
        const double weight_constant = 4.3364230321953841;

        double value = 0;
        value += component_size_coefficient * static_cast<double>(component.size());
        value += secondary_votes_coefficient * static_cast<double>(mmsv);
        do_split = value > weight_constant;
    } else {
        do_split = mmsv >= max_votes;
    }
    if (component.size() <= 5) {
        do_split = false;
    }
    // What??
    if (max_votes > component.size()) {
        do_split = false;
    }

    if (!do_split) {
        seq_string_t consensus;
        for (size_t i = 0; i != length(profile); i++) {
            size_t idx = seqan::getMaxIndex(profile[i]);
            if (idx < seqan::ValueSize<char_t>::VALUE) { // is not gap  TODO Check it!!
                seqan::appendValue(consensus, char_t(idx));
            }
        }
        component.consensus = consensus;
        res.push_back(component);
        return;
    }

    Component majory(component.GetId()), secondary(component.GetId());
    std::vector<ReadPtr> other;
    for (ReadPtr& read : component) {
        if (seqan::ordValue(read->read[maximal_mismatch.position]) == maximal_mismatch.majory_letter) {
            majory.push_back(read);
        } else if (seqan::ordValue(read->read[maximal_mismatch.position]) == maximal_mismatch.secondary_letter) {
            secondary.push_back(read);
        } else {
            other.push_back(read);
        }
    }

    VERIFY(majory.size() == maximal_mismatch.majory_votes);
    VERIFY(secondary.size() == maximal_mismatch.secondary_votes);

    auto majory_consensus = consensus_hamming(majory.GetReadSequences());
    auto secondary_consensus = consensus_hamming(secondary.GetReadSequences());

    for (ReadPtr& read : other) {
        auto dist_majory = hamming_rtrim(read->read, majory_consensus);
        auto dist_secondary = hamming_rtrim(read->read, secondary_consensus);

        if (dist_majory <= dist_secondary) {
            majory.push_back(read);
        } else {
            other.push_back(read);
        }
    }

    VERIFY(majory.size() + secondary.size() == component.size());

    INFO("Component splitted " << majory.size() << " + " << secondary.size());

    // What?
    if (!this->cfg().recursive) {
        max_votes = 0;
    }

    split_component(majory, res, max_votes);

    if (this->cfg().discard) {
        for (ReadPtr& read : secondary) {
            Component c(secondary.GetId());
            c.push_back(read);
            c.consensus = read->read;
            res.push_back(c);
        }
    } else {
        split_component(secondary, res, max_votes);
    }
}


std::vector<Component> LegacySplitterAlgorithm::SplitComponent(Component &component) {
    std::vector<Component> res;
    split_component(component, res, this->cfg().max_votes);
    return res;
}

}
