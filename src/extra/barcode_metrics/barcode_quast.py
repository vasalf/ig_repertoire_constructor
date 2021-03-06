#!/usr/bin/env python

import os.path
import sys
import shutil
import logging
import argparse
import reading_utils
import metrics

def ParseCommandLine():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Bc", type=str, dest="barcode_clusters", help="CLUSTERS.FA or MiGEC file for assembled barcodes", required=True)
    parser.add_argument("--Br", type=str, dest="barcode_rcm", help="RCM file for assembled barcodes")
    parser.add_argument("-c", type=str, dest="data_clusters", help="CLUSTERS.FA for assembled barcodes", required=True)
    parser.add_argument("-r", type=str, dest="data_rcm", help="RCM for assembled barcodes")
    parser.add_argument("--tau", type=int, dest="tau", help="Tau parameter for algorithm", default=3)
    parser.add_argument("--rerun", action="store_true", dest='rerun', help="Run with match results from output directory")
    parser.add_argument("--threads-num", type=int, dest="threads_num", help="Number of threads to use", default=4)
    parser.add_argument("--out", type=str, dest="output_dir", help="Output directory", required=True)
    parser.add_argument("--rate", type=float, dest="rate", help="Rate for good and bad components", default=0.5)
    parser.add_argument("--deep-rcm-cmp", action="store_true", dest="deep_rcm_cmp", help="Compute alignment stats for good barcodes")
    return parser.parse_args()

def InitParams(params):
    params.ig_kplus_vj_finder = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        '../../build/release/bin/ig_kplus_vj_finder')
    params.ig_trie_compressor = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        '../../build/release/bin/ig_trie_compressor')
    params.ig_matcher = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        '../../build/release/bin/ig_matcher')
    params.germline_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        '../fast_ig_tools/germline')

def CreateLogger(params):
    log = logging.getLogger('barcode_quast')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    params.log_filename = os.path.join(params.output_dir, 'barcode_quast.log')
    if os.path.exists(params.log_filename):
        os.remove(params.log_filename)
    log_handler = logging.FileHandler(params.log_filename, mode='a')
    log.addHandler(log_handler)
    return log

def PrepareOutputDir(params):
    if not params.rerun and os.path.exists(params.output_dir):
        shutil.rmtree(params.output_dir)
    if not os.path.isdir(params.output_dir):
        os.makedirs(params.output_dir)

def ReadRepertoires(params, log):
    log.info(".. Reading repertoire from input files")
    params.barcode_repertoire, params.data_repertoire = \
        reading_utils.read_repertoires(params.barcode_clusters, params.barcode_rcm, 
                                       params.data_clusters, params.data_rcm)

def RunIgMatcherPreparations(params, log, rep, working_dir):
    vj_aln_file = os.path.join(working_dir, 'cleaned_reads.fa')
    command_line = params.ig_kplus_vj_finder + \
            ' -i ' + rep.clusters_filename + \
            ' -o ' + working_dir + \
            ' --db-directory ' + params.germline_dir + \
            ' | tee -a ' + params.log_filename
    if not params.rerun:
        exit_code = os.system(command_line)
        if exit_code != 0:
            log.info('ERROR: ig_kplus_vj_finder failed on ' + rep.clusters_filename + ' file')
            sys.exit(-1)
    if not os.path.exists(vj_aln_file):
        log.info('ERROR: cannot find VJ alignment file ' + vj_aln_file)
        sys.exit(-1)
    trie_compressed_file = os.path.join(working_dir, rep.name + '.fa')
    command_line = params.ig_trie_compressor + \
            ' -i ' + vj_aln_file + \
            ' -o ' + trie_compressed_file + \
            ' | tee -a ' + params.log_filename
    if not params.rerun:
        exit_code = os.system(command_line)
        if exit_code != 0:
            log.info('ERROR: ig_trie_compressor failed on ' + vj_aln_file + ' file')
            sys.exit(-1)
    if not os.path.exists(trie_compressed_file):
        log.info('ERROR: cannot find compressed trie file ' + trie_compressed_file)
        sys.exit(-1)
    cluster_num_to_ids = reading_utils.read_cluster_numbers_to_ids(trie_compressed_file)
    return trie_compressed_file, cluster_num_to_ids

def RunIgMatcher(params, log):
    log.info(".. Building neighbour clusters")

    barcode_aln_file, barcode_cluster_num_to_ids = RunIgMatcherPreparations(
        params, log, params.barcode_repertoire, 
        os.path.join(params.output_dir, params.barcode_repertoire.name))
    data_aln_file, data_cluster_num_to_ids = RunIgMatcherPreparations(
        params, log, params.data_repertoire, 
        os.path.join(params.output_dir, params.data_repertoire.name))

    params.compressed_barcode_repertoire = reading_utils.read_repertoire(barcode_aln_file, None, 'compressed_assembled_barcodes')
    params.compressed_data_repertoire = reading_utils.read_repertoire(data_aln_file, None, 'compressed_' + params.data_repertoire.name)

    barcode_matches_file = os.path.join(params.output_dir, params.barcode_repertoire.name + '.match')
    data_matches_file = os.path.join(params.output_dir, params.data_repertoire.name + '.match')
    command_line = params.ig_matcher + \
            ' -i ' + barcode_aln_file + \
            ' -I ' + data_aln_file + \
            ' --tau ' + str(params.tau) + ' --threads ' + str(params.threads_num) + \
            ' -o ' + barcode_matches_file + \
            ' -O ' + data_matches_file + ' | tee -a ' + params.log_filename

    if not params.rerun:
        exit_code = os.system(command_line)
        if exit_code != 0:
            log.info('ERROR: ig_matcher failed on ' + barcode_aln_file + ' and ' + data_aln_file + ' files')
            sys.exit(-1)
    if not os.path.exists(barcode_matches_file) or not os.path.exists(data_matches_file):
        log.info('ERROR: cannot find matches file') 
        sys.exit(-1)

    log.info(".. Matches between barcodes and data repertoire were computed")
    return (reading_utils.read_cluster_matches(barcode_matches_file, params.barcode_repertoire, 
                barcode_cluster_num_to_ids, data_cluster_num_to_ids, params.tau),
            reading_utils.read_cluster_matches(data_matches_file, params.data_repertoire,
                data_cluster_num_to_ids, barcode_cluster_num_to_ids, params.tau))

def Evaluate(params, log):
    ReadRepertoires(params, log)
    barcode_cluster_matches, data_cluster_matches = RunIgMatcher(params, log)
    barcode_metrics = metrics.BarcodeMetrics(
        params.barcode_repertoire, params.compressed_barcode_repertoire, barcode_cluster_matches, 
        params.data_repertoire, params.compressed_data_repertoire, data_cluster_matches, params.rate)
    barcode_metrics.evaluate(log, params.deep_rcm_cmp)
    metrics_file = os.path.join(params.output_dir, 'metrics.txt')
    barcode_metrics.write(metrics_file)
    # sizes_corr_filename = os.path.join(params.output_dir, 'sizes_corr.png')
    # barcode_metrics.draw_sizes_correlation_plot(sizes_corr_filename)
    if not os.path.exists(metrics_file):
        log.info("ERROR: barcode metrics were not found")
        sys.exit(-1)
    log.info("Barcode metrics were written to " + metrics_file)
    metrics_json_file = os.path.join(params.output_dir, 'metrics.json')
    barcode_metrics.write_json(metrics_json_file)
    if not os.path.exists(metrics_json_file):
        log.info("ERROR: barcode metrics in JSON were not found")
        sys.exit(-1)
    log.info("Barcode metrics in JSON were written to " + metrics_json_file)

    distr_filename, distr_nt_filename, isolated_distr_filename, isolated_distr_nt_filename = \
        barcode_metrics.draw_all_sizes_distributions(params.output_dir)
    if os.path.exists(distr_filename):
        log.info("Cluster sizes distribution were written to " + distr_filename)
    if os.path.exists(distr_nt_filename):
        log.info("Non-trivial cluster sizes distribution were written to " + distr_nt_filename)
    if os.path.exists(isolated_distr_filename):
        log.info("Isolated cluster sizes distribution were written to " + isolated_distr_filename)
    if os.path.exists(isolated_distr_nt_filename):
        log.info("Non-trivial isolated cluster sizes distribution were written to " + isolated_distr_nt_filename)
    length_distr_filename, isolated_length_distr_filename = \
        barcode_metrics.draw_all_lengths_distribution(params.output_dir)
    if os.path.exists(length_distr_filename):
        log.info("Cluster sequence's length distribution were written to " + length_distr_filename)
    if os.path.exists(isolated_length_distr_filename):
        log.info("Isolated cluster sequence's length distribution were written to " + isolated_length_distr_filename)

def main():
    params = ParseCommandLine()
    PrepareOutputDir(params)
    InitParams(params)
    log = CreateLogger(params)
    try:
        Evaluate(params, log)
    except KeyboardInterrupt:
        log.info("\nApplication was interrupted!")

if __name__ == "__main__":
    main()
