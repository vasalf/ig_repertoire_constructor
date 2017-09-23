#pragma once

#include <boost/make_unique.hpp>

#include <ig_final_alignment.hpp>
#include "config.hpp"
#include "component.hpp"
#include "input_data.hpp"
#include "splitter_algorithm.hpp"

#include <seqan/seq_io.h>

namespace ig_component_splitter {
    
class IgComponentSplitterProcessor {
    IgComponentSplitterConfig &cfg_;
    InputData &input_;
    
public:
    IgComponentSplitterProcessor(IgComponentSplitterConfig &cfg,
                                 InputData &input)
        : cfg_(cfg), input_(input) {}

    void Process();
};


class IgComponentSplitterParallelProcessor {
    IgComponentSplitterConfig &cfg_;
    InputData &input_;

    struct ResultRecord {
        Component comp;
        std::string cluster_id;
        std::string id;

        ResultRecord(Component &&comp_,
                     std::string cluster_id_,
                     std::string id_)
            : comp(comp_),
              cluster_id(cluster_id_),
              id(id_) {}
    };
    
public:
    IgComponentSplitterParallelProcessor(IgComponentSplitterConfig &cfg,
                                         InputData &input)
        : cfg_(cfg), input_(input) {}

    void Process();
};

}
