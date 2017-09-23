#include "utils.hpp"
#include <fast_ig_tools.hpp>
#include "processor.hpp"

namespace ig_component_splitter {

void IgComponentSplitterProcessor::Process() {
    std::vector<Component> components;
    for (auto &p : input_.Components()) {
        components.push_back(p.second);
    }
    if (!map_type<std::string, Component>::is_sorted_container) {
        std::sort(components.begin(), components.end());
    }

    seqan::SeqFileOut seq_output(cfg_.output_file.c_str());
    std::ofstream out_rcm(cfg_.output_rcm_file.c_str());


    std::unique_ptr<SplitterAlgorithm > splitter
        = boost::make_unique<DummySplitterAlgorithm >(cfg_, input_);

    for (auto &component : components) {
        auto result = splitter->SplitComponent(component);
        for (size_t i = 0; i != result.size(); i++) {
            std::ostringstream cluster_id;
            cluster_id << component.GetId();
            if (result.size() > 1) {
                cluster_id << "X" << i;
            }
                
            std::ostringstream id;
            id << "cluster___"
               << cluster_id.str()
               << "___size___"
               << result[i].size();

            seqan::writeRecord(seq_output, id.str(), result[i].consensus);
            for (const auto& read : result[i]) {
                std::string read_id = seqan::toCString(read->id);
                out_rcm << read_id << "\t" << cluster_id.str() << "\n";
            }
        }
    }

    INFO("Final repertoire was written to " << cfg_.output_file);
    INFO("Final RCM was written to " << cfg_.output_rcm_file);
}

void IgComponentSplitterParallelProcessor::Process() {
    std::vector<Component> components; 
    for (const auto &p: input_.Components()) {
        components.push_back(p.second);
    }
    if (!map_type<std::string, Component>::is_sorted_container) {
        std::sort(components.begin(), components.end());
    }

    std::unique_ptr<SplitterAlgorithm > splitter(new LegacySplitterAlgorithm(cfg_, input_));

    std::vector<ResultRecord> output;
        
    omp_set_num_threads(cfg_.nthreads);
    INFO("ig_component_splitter starts splitting components in " << cfg_.nthreads << " threads");

#pragma omp parallel for
    for (size_t comp_i = 0; comp_i < components.size(); comp_i++) {
        auto &component = components[comp_i];
        auto result = splitter->SplitComponent(component);

#pragma omp critical
        {
            // Storing the result into the vector is a critical section
                
            for (size_t i = 0; i < result.size(); i++) {
                std::ostringstream cluster_id;
                cluster_id << component.GetId();
                if (result.size() > 1) {
                    cluster_id << "X" << i;
                }
                
                std::ostringstream id;
                id << "cluster___"
                   << cluster_id.str()
                   << "___size___"
                   << result[i].size();

                output.push_back(ResultRecord(std::move(result[i]), cluster_id.str(), id.str()));
            }
        }
            
    }

    seqan::SeqFileOut seq_output(cfg_.output_file.c_str());
    std::ofstream out_rcm(cfg_.output_rcm_file.c_str());
        
    for (auto &res : output) {
        seqan::writeRecord(seq_output, res.id, res.comp.consensus);
        for (auto &read : res.comp) {
            std::string read_id = seqan::toCString(read->id);
            out_rcm << read_id << "\t" << res.cluster_id << "\n";
        }
    }

    INFO("Final repertoire was written to " << cfg_.output_file);
    INFO("Final RCM was written to " << cfg_.output_rcm_file);

}

}
