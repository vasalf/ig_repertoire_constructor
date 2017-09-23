#pragma once

#include "config.hpp"
#include "component.hpp"
#include "input_data.hpp"

namespace ig_component_splitter {

// This is a base class for splitting algorithms
class SplitterAlgorithm {
    IgComponentSplitterConfig &cfg_;
    InputData &input_;
protected:
    IgComponentSplitterConfig & cfg() {
        return cfg_;
    }
    InputData & input() {
        return input_;
    }
public:
    SplitterAlgorithm(IgComponentSplitterConfig &cfg,
                      InputData &input)
        : cfg_(cfg), input_(input) {}

    virtual ~SplitterAlgorithm() = default;
    
    virtual std::vector<Component> SplitComponent(Component &component) = 0;
};

class DummySplitterAlgorithm : public SplitterAlgorithm {
public:
    DummySplitterAlgorithm(IgComponentSplitterConfig &cfg,
                           InputData &input)
        : SplitterAlgorithm(cfg, input) {}

    virtual ~DummySplitterAlgorithm() = default;

    virtual std::vector<Component> SplitComponent(Component &component) override {
        component.consensus = consensus(component.GetReadSequences());
        return { component };
    }
};

class LegacySplitterAlgorithm : public SplitterAlgorithm {
    typedef seqan::String<seqan::ProfileChar<char_t> > profile_string_t;
    
    void split_component(Component &component, std::vector<Component> &res, size_t max_votes);
    
public:
    LegacySplitterAlgorithm(IgComponentSplitterConfig &cfg,
                            InputData &input)
        : SplitterAlgorithm(cfg, input) {}

    virtual ~LegacySplitterAlgorithm() = default;

    virtual std::vector<Component> SplitComponent(Component &component) override;
};

}
