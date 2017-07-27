#pragma once

#include "config.hpp"
#include "component.hpp"
#include "input_data.hpp"

namespace ig_component_splitter {

// This is a base class for splitting algorithms
template<template<typename,typename> class map_type>
class SplitterAlgorithm {
    IgComponentSplitterConfig &cfg_;
    InputData<map_type> &input_;
protected:
    IgComponentSplitterConfig & cfg() {
        return cfg_;
    }
    InputData<map_type> & input() {
        return input_;
    }
public:
    SplitterAlgorithm(IgComponentSplitterConfig &cfg,
                      InputData<map_type> &input)
        : cfg_(cfg), input_(input) {}

    virtual ~SplitterAlgorithm() = default;
    
    virtual std::vector<Component> SplitComponent(Component &component) = 0;
};

template<template<typename,typename> class map_type>
class DummySplitterAlgorithm : public SplitterAlgorithm<map_type> {
public:
    DummySplitterAlgorithm(IgComponentSplitterConfig &cfg,
                           InputData<map_type> &input)
        : SplitterAlgorithm<map_type>(cfg, input) {}

    virtual ~DummySplitterAlgorithm() = default;

    virtual std::vector<Component> SplitComponent(Component &component) override {
        component.Consensus() = consensus(component.GetReadSequences());
        return { component };
    }
};

template<template<typename, typename> class map_type>
class LegacySplitterAlgorithm : public SplitterAlgorithm<map_type> {
    typedef seqan::String<seqan::ProfileChar<char_t> > profile_string_t;
    
    template<class T>
    struct SeqanLengthCompare {
        bool operator()(const T& first, const T& second) {
            return seqan::length(first) < seqan::length(second);
        }
    };

    struct PositionVote {
        std::size_t majory_votes;
        std::size_t majory_letter;
        std::size_t secondary_votes;
        std::size_t secondary_letter;
        std::size_t position;
        bool operator<(const PositionVote &other) const {
            return secondary_votes < other.secondary_votes;
        }
    };
    
    void split_component(Component component,
                         std::vector<Component> res,
                         std::size_t max_votes) {
        if (component.empty())
            return;
        else if (component.size() == 1) {
            component.Consensus() = component[0]->read;
            res.push_back(component);
            return;
        }

        // What?
        if (max_votes == 0)
            max_votes = std::numeric_limits<std::size_t>::max() / 2;

        std::vector<seq_string_t> reads_in_component = component.GetReadSequences();
        
        profile_string_t profile;
        std::size_t len = seqan::length(*max_element(reads_in_component.begin(), reads_in_component.end(),
                                                     SeqanLengthCompare<seqan::Dna5String>()));
        seqan::resize(profile, len);

        for (const seq_string_t& read : reads_in_component) {
            for (std::size_t j = 0; j != length(read); j++)
                profile[j].count[seqan::ordValue(read[j])]++;
        }

        std::size_t min_len = seqan::length(*min_element(reads_in_component.begin(), reads_in_component.end(),
                                                         SeqanLengthCompare<seqan::Dna5String>()));
        std::vector<PositionVote> votes(min_len);
        for (std::size_t i = 0; i != min_len; i++) {
            std::array<std::pair<std::size_t, std::size_t>, seqan::ValueSize<char_t>::VALUE> letters;
            for (std::size_t k = 0; k < seqan::valueSize<char_t>(); k++)
                letters[k] = std::make_pair(profile[i].count[k], k);


            std::nth_element(letters.begin(), letters.begin() + 1, letters.end(),
                             std::greater<std::pair<std::size_t, std::size_t> >());

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
        std::size_t mmsv = maximal_mismatch.secondary_votes;

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
        if (component.size() <= 5)
            do_split = false;
        // What??
        if (max_votes > component.size())
            do_split = false;

        if (!do_split) {
            component.Consensus() = consensus(reads_in_component);
            res.push_back(component);
            return;
        }

        Component majory(component.GetId()), secondary(component.GetId());
        std::vector<ReadPtr> other;
        for (ReadPtr& read : component) {
            if (read->read[maximal_mismatch.position] == maximal_mismatch.majory_letter)
                majory.push_back(read);
            else if (read->read[maximal_mismatch.position] == maximal_mismatch.secondary_letter)
                secondary.push_back(read);
            else
                other.push_back(read);
        }

        VERIFY(majory.size() == maximal_mismatch.majory_votes);
        VERIFY(secondary.size() == maximal_mismatch.secondary_votes);

        auto majory_consensus = consensus_hamming(majory.GetReadSequences());
        auto secondary_consensus = consensus_hamming(secondary.GetReadSequences());

        for (ReadPtr& read : other) {
            auto dist_majory = hamming_rtrim(read->read, majory_consensus);
            auto dist_secondary = hamming_rtrim(read->read, secondary_consensus);

            if (dist_majory <= dist_secondary)
                majory.push_back(read);
            else
                other.push_back(read);
        }

        VERIFY(majory.size() + secondary.size() == component.size());

        INFO("Component splitted " << majory.size() << " + " << secondary.size());

        // What?
        if (!this->cfg().recursive)
            max_votes = 0;

        split_component(majory, res, max_votes);

        if (this->cfg().discard) {
            for (ReadPtr& read : secondary) {
                Component c(secondary.GetId());
                c.push_back(read);
                c.Consensus() = read->read;
                res.push_back(c);
            }
        } else {
            split_component(secondary, res, max_votes);
        }
    }
    
public:
    LegacySplitterAlgorithm(IgComponentSplitterConfig &cfg,
                            InputData<map_type> &input)
        : SplitterAlgorithm<map_type>(cfg, input) {}

    virtual ~LegacySplitterAlgorithm() = default;

    virtual std::vector<Component> SplitComponent(Component &component) override {
        std::vector<Component> res;
        split_component(component, res, this->cfg().max_votes);
        return res;
    }
};

}
