#include "vdj_launch.hpp"

#include <read_archive.hpp>
#include <germline_db_generator.hpp>
#include <vj_parallel_processor.hpp>
#include <block_alignment/block_alignment_converter.hpp>
#include <vdj_alignments/vdj_hits.hpp>
#include <vdj_alignments/vdj_hits_storage.hpp>
#include <alignment_utils/alignment_positions.hpp>
#include <vdj_alignments/aligners/simple_d_aligner.hpp>
#include <model/recombination_model.hpp>
#include <germline_utils/chain_type.hpp>
#include <recombination_calculator/hc_model_based_recombination_calculator.hpp>
#include <recombination_utils/recombination_storage.hpp>

namespace vdj_labeler {

void VDJLabelerLaunch::TestRecombinationCalculator(const core::ReadArchive& reads_archive,
                                                   VDJHitsStoragePtr &hits_storage)
{
    size_t read_index = 3;
    core::ReadPtr read_3 = std::make_shared<core::Read>(reads_archive[read_index]);
    VDJHitsPtr hits_3 = (*hits_storage)[read_index];
    INFO("Read 3. #V: " << hits_3->VHitsNumber() <<
        ", #D: " << hits_3->DHitsNumber() <<
        ", #J: " << hits_3->JHitsNumber());

    auto v_alignment = hits_3->GetAlignmentByIndex(germline_utils::SegmentType::VariableSegment, 0);
    recombination_utils::CleavedIgGeneAlignment v_event_0(v_alignment, 0, 0, 0, 0);
    recombination_utils::CleavedIgGeneAlignment v_event_1(v_alignment, 0, -1, 0, 0);
    recombination_utils::CleavedIgGeneAlignment v_event_2(v_alignment, 0, -2, 0, 0);
    recombination_utils::CleavedIgGeneAlignment v_event_3(v_alignment, 0, -3, 0, 1);

    auto d_alignment = hits_3->GetAlignmentByIndex(germline_utils::SegmentType::DiversitySegment, 0);
    recombination_utils::CleavedIgGeneAlignment d_event_0(d_alignment, 1, 8, 0, 0);

    auto j_alignment = hits_3->GetAlignmentByIndex(germline_utils::SegmentType::JoinSegment, 0);
    recombination_utils::CleavedIgGeneAlignment j_event_0(j_alignment, 0, 0, 1, 0);
    recombination_utils::CleavedIgGeneAlignment j_event_1(j_alignment, 1, 0, 0, 0);

    recombination_utils::NongenomicInsertion vd_insertion_0(425, 441);
    recombination_utils::NongenomicInsertion vd_insertion_1(426, 441);
    recombination_utils::NongenomicInsertion vd_insertion_2(427, 441);
    recombination_utils::NongenomicInsertion vd_insertion_3(428, 441);

    recombination_utils::NongenomicInsertion dj_insertion_0(453, 452);
    recombination_utils::NongenomicInsertion dj_insertion_1(453, 453);

    recombination_utils::RecombinationStorage<recombination_utils::HCRecombination> recombination_storage(read_3);
    recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_0, d_event_0, j_event_0,
                                                           vd_insertion_0, dj_insertion_0));
    recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_1, d_event_0, j_event_0,
                                                           vd_insertion_1, dj_insertion_0));
    recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_2, d_event_0, j_event_0,
                                                           vd_insertion_2, dj_insertion_0));
    recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_3, d_event_0, j_event_0,
                                                           vd_insertion_3, dj_insertion_0));
    recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_0, d_event_0, j_event_1,
                                                           vd_insertion_0, dj_insertion_1));
    recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_1, d_event_0, j_event_1,
                                                           vd_insertion_1, dj_insertion_1));
    recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_2, d_event_0, j_event_1,
                                                           vd_insertion_2, dj_insertion_1));
    recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_3, d_event_0, j_event_1,
                                                           vd_insertion_3, dj_insertion_1));
    INFO(recombination_storage.size() << " recombinaions were generated");
}

void VDJLabelerLaunch::Launch() {
    INFO("VDJ labeler starts");
    std::string input_filename = config_.io_params.input_params.input_sequences;
    std::string v_germline_genes_fname = config_.io_params.input_params.germlines.igh_genes.variable_genes;
    std::string d_germline_genes_fname = config_.io_params.input_params.germlines.igh_genes.diversity_genes;
    std::string j_germline_genes_fname = config_.io_params.input_params.germlines.igh_genes.join_genes;

    core::ReadArchive read_archive(input_filename);
    read_archive.FixSpacesInHeaders();

    using namespace germline_utils;
    CustomGeneDatabase v_db(SegmentType::VariableSegment);
    CustomGeneDatabase d_db(SegmentType::DiversitySegment);
    CustomGeneDatabase j_db(SegmentType::JoinSegment);

    v_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::VariableSegment), v_germline_genes_fname);
    d_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::DiversitySegment), d_germline_genes_fname);
    j_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::JoinSegment), j_germline_genes_fname);

    germline_utils::ChainDatabase hc_db(germline_utils::ImmuneChainType::HeavyIgChain);
    hc_db.AddGenesFromFile(SegmentType::VariableSegment, v_germline_genes_fname);
    hc_db.AddGenesFromFile(SegmentType::DiversitySegment, d_germline_genes_fname);
    hc_db.AddGenesFromFile(SegmentType::JoinSegment, j_germline_genes_fname);

    INFO("Generation of DB for join segments...");
    INFO("Alignment against VJ germline segments");
    vj_finder::VJParallelProcessor processor(read_archive, config_.vj_finder_config.algorithm_params,
                                             v_db, j_db,
                                             config_.run_params.threads_count);
    vj_finder::VJAlignmentInfo alignment_info = processor.Process();
    INFO(alignment_info.NumVJHits() << " reads were aligned; " << alignment_info.NumFilteredReads() <<
        " reads were filtered out");

    // INFO(alignment_info.AlignmentRecords()[6].Read());
    // INFO(alignment_info.AlignmentRecords()[6].VHits()[0].ImmuneGene());
    // INFO(alignment_info.AlignmentRecords()[0].VHits()[0].Read());
    // INFO(alignment_info.AlignmentRecords()[0].JHits().size());
    // auto vdj_hits = VDJHits(alignment_info.AlignmentRecords()[6]);
    // INFO((*(vdj_hits.VHits().cbegin()))->Alignment());

    auto vdj_storage = VDJHitsStorage(alignment_info);
    for (auto it = vdj_storage.cbegin(); it != vdj_storage.cend(); ++it) {
        INFO("\nVgene");
        for (auto vhit_it = (*it)->VHits().cbegin(); vhit_it != (*it)->VHits().cend(); ++vhit_it) {
            INFO((*vhit_it)->EndQueryPosition());
        }
    }
    // INFO(*(vdj_storage[0]->Read()));

    alignment_utils::AlignmentPositions alignment_positions(std::make_pair<size_t, size_t>(0, read_archive[0].length() - 1),
                                                            std::make_pair<size_t, size_t>(0, 10));

    alignment_utils::ImmuneGeneAlignmentPositions immune_alignment_positions(alignment_positions,
                                                                             d_db[0],
                                                                             read_archive[0]);

    INFO(SimpleDAligner().ComputeAlignment(immune_alignment_positions)->Alignment());

    // Andy: Blank model "tested" here
    {
        std::ifstream in("src/vdj_labeler_0.3/test/blank_model.csv");
        // HCProbabilityRecombinationModel model(in, hc_db);
        // std::cout << model;
        IgGeneProbabilityModel model_V(in, hc_db.GetDb(germline_utils::SegmentType::VariableSegment));
        // std::cout << model_V;
        IgGeneProbabilityModel model_D(in, hc_db.GetDb(germline_utils::SegmentType::DiversitySegment));
        // std::cout << model_D;
        // IgGeneProbabilityModel model_J(in, hc_db.JoinGenes());
        // // cout << model_J;
        // NongenomicInsertionModel modelVD(in);
        // NongenomicInsertionModel modelDJ(in);
        // // cout << modelVD;
        // // cout << modelDJ;
        // PalindromeDeletionModel modelDelV(in, hc_db.VariableGenes());
        // // cout << modelDelV;
        // // cout << modelDelV.GetIgGeneDatabase() -> GetByIndex(0) -> name() << " ";
        // // cout << modelDelV.GetDeletionProbability(0, 0) << endl;
        // PalindromeDeletionModel modelDelJ(in, hc_db.JoinGenes());
        // // cout << modelDelJ.GetIgGeneDatabase() -> GetByIndex(1) -> name() << " ";
        // // cout << modelDelJ.GetDeletionProbability(1, -2) << endl;
        // PalindromeDeletionModel modelDelDL(in, hc_db.DiversityGenes());
        // PalindromeDeletionModel modelDelDR(in, hc_db.DiversityGenes());

        // HCModelBasedRecombinationCalculator recombination_calculator(model);
    }
}

} // End namespace vdj_labeler