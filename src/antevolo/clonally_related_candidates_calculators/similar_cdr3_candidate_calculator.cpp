#include <path_helper.hpp>
#include "similar_cdr3_candidate_calculator.hpp"
#include "../../graph_utils/graph_io.hpp"
#include "../../graph_utils/graph_splitter.hpp"
#include <convert.hpp>

namespace antevolo {
    void SimilarCDR3CandidateCalculator::Clear() {
        unique_cdr3s_.clear();
        unique_cdr3s_map_.clear();
    }

    void SimilarCDR3CandidateCalculator::CreateUniqueCDR3Map(
            core::DecompositionClass decomposition_class) {
        for(auto it = decomposition_class.begin(); it != decomposition_class.end(); it++) {
            if(clone_set_[*it].RegionIsEmpty(annotation_utils::StructuralRegion::CDR3))
                continue;
            auto cdr3 = core::dna5String_to_string(clone_set_[*it].CDR3());
            if(unique_cdr3s_map_.find(cdr3) == unique_cdr3s_map_.end())
                unique_cdr3s_map_[cdr3] = std::vector<size_t>();
            unique_cdr3s_map_[cdr3].push_back(*it);
        }
        for(auto it = unique_cdr3s_map_.begin(); it != unique_cdr3s_map_.end(); it++)
            unique_cdr3s_.push_back(it->first);
    }

    std::string SimilarCDR3CandidateCalculator::GetFastaFname(core::DecompositionClass decomposition_class) {
        std::stringstream ss;
        size_t key = *decomposition_class.begin();
        ss << "CDR3_sequences_key_" << key << ".fasta";
        return path::append_path(output_params_.cdr_graph_dir, ss.str());
    }

    std::string SimilarCDR3CandidateCalculator::GetGraphFname(core::DecompositionClass decomposition_class) {
        std::stringstream ss;
        size_t key = *decomposition_class.begin();
        ss << "CDR3_sequences_key_" << key << ".graph";
        return path::append_path(output_params_.cdr_graph_dir, ss.str());
    }

    std::string SimilarCDR3CandidateCalculator::WriteUniqueCDR3InFasta(core::DecompositionClass decomposition_class) {
        std::string output_fname = GetFastaFname(decomposition_class);
        std::ofstream out(output_fname);
        for(auto it = unique_cdr3s_.begin(); it != unique_cdr3s_.end(); it++)
            out << ">" << *it << std::endl << *it << std::endl;
        return output_fname;
    }

    // return connected components of Hamming graph on CDR3s
    std::vector<SparseGraphPtr> SimilarCDR3CandidateCalculator::ComputeCDR3HammingGraphs(std::string cdr_fasta,
                                                                                         std::string graph_fname) {
        std::string run_graph_constructor = "./build/release/bin/ig_swgraph_construct";
        std::stringstream ss;
        ss << run_graph_constructor << " -i " << cdr_fasta <<
                " -o " << graph_fname << " --tau " << num_mismatches_ << " -k 10 > " << output_params_.trash_output;
        int err_code = system(ss.str().c_str());
        VERIFY_MSG(err_code == 0, "Graph constructor finished abnormally, error code: " << err_code);
        auto sparse_cdr_graph_ = GraphReader(graph_fname).CreateGraph();
        TRACE("Hamming graph contains " << sparse_cdr_graph_->N() << " edges and " << sparse_cdr_graph_->NZ() << " edges");
        //if(sparse_cdr_graph_->NZ() != 0)
        //    std::cout << *sparse_cdr_graph_ << std::endl;
        auto connected_components = ConnectedComponentGraphSplitter(sparse_cdr_graph_).Split();
        graph_component_ = sparse_cdr_graph_->GetGraphComponentMap();
        return connected_components;
    }

    ClonallyRelatedCandidates SimilarCDR3CandidateCalculator::ComputeCandidatesForGraph(SparseGraphPtr hg_component,
                                                                                        size_t component_id) {
        ClonallyRelatedCandidates candidates(clone_set_);
        // adding edges between identical CDR3s
        for(size_t i = 0; i < hg_component->N(); i++) {
            size_t old_index = graph_component_.GetOldVertexByNewVertex(component_id, i);
            auto clones_sharing_cdr3 = unique_cdr3s_map_[unique_cdr3s_[old_index]];
            for(size_t it1 = 0; it1 < clones_sharing_cdr3.size(); it1++)
                for(size_t it2 = it1 + 1; it2 < clones_sharing_cdr3.size(); it2++)
                    candidates.AddCandidatePair(clones_sharing_cdr3[it1], clones_sharing_cdr3[it2]);
        }
        // adding edges between similar CDR3s
        for(size_t i = 0; i < hg_component->N(); i++)
            for(size_t j = hg_component->RowIndex()[i]; j < hg_component->RowIndex()[i + 1]; j++) {
                size_t old_index1 = graph_component_.GetOldVertexByNewVertex(component_id, i);
                size_t old_index2 = graph_component_.GetOldVertexByNewVertex(component_id, hg_component->Col()[j]);
                auto indices_1 = unique_cdr3s_map_[unique_cdr3s_[old_index1]];
                auto indices_2 = unique_cdr3s_map_[unique_cdr3s_[old_index2]];
                for(auto it1 = indices_1.begin(); it1!= indices_1.end(); it1++)
                    for(auto it2 = indices_2.begin(); it2!= indices_2.end(); it2++)
                        candidates.AddCandidatePair(*it1, *it2);
            }
        return candidates;
    }

    std::vector<ClonallyRelatedCandidates> SimilarCDR3CandidateCalculator::ComputeCandidates(
            core::DecompositionClass decomposition_class) {
        Clear();
        CreateUniqueCDR3Map(decomposition_class);
        std::string cdrs_fasta = WriteUniqueCDR3InFasta(decomposition_class);
        std::string graph_fname = GetGraphFname(decomposition_class);
        TRACE("--------------------------");
        TRACE("CDR3 fasta: "<< cdrs_fasta << ", CDR3 Hamming graph: " << graph_fname);
        auto connected_components = ComputeCDR3HammingGraphs(cdrs_fasta, graph_fname);
        TRACE("# connected components: " << connected_components.size());
        std::vector<ClonallyRelatedCandidates> vector_candidates;
        for(size_t i = 0; i < connected_components.size(); i++)
            vector_candidates.push_back(ComputeCandidatesForGraph(connected_components[i], i));
        return vector_candidates;
    }
}