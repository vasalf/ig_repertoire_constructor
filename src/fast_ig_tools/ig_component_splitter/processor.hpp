#pragma once

#include "config.hpp"
#include "component.hpp"
#include "input_data.hpp"
#include "splitter_algorithm.hpp"

namespace ig_component_splitter {

    
template<template<typename,typename> class map_type>
class IgComponentSplitterProcessor {
    IgComponentSplitterConfig &cfg_;
    InputData<map_type> &input_;
    
public:
    IgComponentSplitterProcessor(IgComponentSplitterConfig &cfg,
                                 InputData<map_type> &input)
        : cfg_(cfg), input_(input) {}

    void Process() {
        std::vector<Component> components;
        for (auto &p : input_.Components())
            components.push_back(p.second);
        if (!map_type<std::string, Component>::is_sorted_container)
            std::sort(components.begin(), components.end());

        seqan::SeqFileOut seq_output(cfg_.output_file.c_str());
        std::ofstream out_rcm(cfg_.output_rcm_file.c_str());

        /*
        std::unique_ptr<SplitterAlgorithm<map_type> > splitter
            = std::make_unique<DummySplitterAlgorithm<map_type> >(cfg_, input_);
        */
        std::unique_ptr<SplitterAlgorithm<map_type> > splitter(new LegacySplitterAlgorithm<map_type>(cfg_, input_));

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


template<template<typename,typename> class map_type>
class IgComponentSplitterParallelProcessor {
    IgComponentSplitterConfig &cfg_;
    InputData<map_type> &input_;

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
                                         InputData<map_type> &input)
        : cfg_(cfg), input_(input) {}

    void Process() {
        std::vector<Component> components;
        for (auto &p: input_.Components())
            components.push_back(p.second);
        if (!map_type<std::string, Component>::is_sorted_container)
            std::sort(components.begin(), components.end());

        std::unique_ptr<SplitterAlgorithm<map_type> > splitter(new LegacySplitterAlgorithm<map_type>(cfg_, input_));

        std::vector<ResultRecord> output;
        
        omp_lock_t result_lock;
        omp_init_lock(&result_lock);
        
        omp_set_num_threads(cfg_.nthreads);
        INFO("ig_component_splitter starts splitting components in " << cfg_.nthreads << " threads");

        #pragma omp parallel for
        for (std::size_t comp_i = 0; comp_i < components.size(); comp_i++) {
            auto &component = components[comp_i];
            auto result = splitter->SplitComponent(component);
            
            omp_set_lock(&result_lock);

            for (std::size_t i = 0; i < result.size(); i++) {
                std::ostringstream cluster_id;
                cluster_id << component.GetId();
                if (result.size() > 1)
                    cluster_id << "X" << i;
                
                std::ostringstream id;
                id << "cluster___"
                   << cluster_id.str()
                   << "___size___"
                   << result[i].size();

                output.push_back(ResultRecord(std::move(result[i]), cluster_id.str(), id.str()));
            }
            
            omp_unset_lock(&result_lock);
        }

        omp_destroy_lock(&result_lock);

        seqan::SeqFileOut seq_output(cfg_.output_file.c_str());
        std::ofstream out_rcm(cfg_.output_rcm_file.c_str());
        
        for (auto &res : output) {
            seqan::writeRecord(seq_output, res.id, res.comp.Consensus());
            for (auto &read : res.comp) {
                std::string read_id = seqan::toCString(read->id);
                out_rcm << read_id << "\t" << res.cluster_id << "\n";
            }
        }

        INFO("Final repertoire was written to " << cfg_.output_file);
        INFO("Final RCM was written to " << cfg_.output_rcm_file);
    }
    
};

}
