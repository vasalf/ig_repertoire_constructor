#include <iostream>

#include <boost/program_options.hpp>

#include "utils.hpp"
#include <logger/logger.hpp>
#include <build_info.hpp>

#include "config.hpp"

namespace ig_component_splitter {

void IgComponentSplitterConfig::Write() {
    INFO("reads_file=" << reads_file);
    INFO("rcm_file=" << rcm_file);
    INFO("output_file=" << output_file);
    INFO("output_rcm_file=" << output_rcm_file);
    INFO("config_file=" << config_file);
    INFO("nthreads=" << nthreads);
    INFO("max_votes=" << max_votes);
    INFO("discard=" << std::boolalpha << discard);
    INFO("recursive=" << std::boolalpha << recursive);
    INFO("flu=" << std::boolalpha << flu);
    INFO("allow_unassigned=" << std::boolalpha << allow_unassigned);
}

void IgComponentSplitterConfig::ReadOptions(int argc, char **argv) {
    namespace po = boost::program_options;

    // Declare a group of options that are
    // only allowed on a command line
    po::options_description generic("Generic options");
    generic.add_options()
        ("version,v", "print version string")
        ("help,h", "produce help message")
        ("config,c", po::value<std::string>(&config_file)->default_value(config_file),
         "name of a file of a configuration")
        ("input-file,i", po::value<std::string>(&reads_file),
         "name of the input file (FASTA|FASTQ)")
        ("output-file,o", po::value<std::string>(&output_file)->default_value(output_file),
         "output file for final repertoire")
        ("rcm-file,R", po::value<std::string>(&rcm_file)->default_value(rcm_file),
         "input RCM-file")
        ("output-rcm-file,M", po::value<std::string>(&output_rcm_file)->default_value(output_rcm_file),
         "output RCM-file");

    // Declare a group of options that are
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
        ("threads,t", po::value<int>(&nthreads)->default_value(nthreads),
         "the number of parallel threads")
        ("max-votes,V", po::value<std::size_t>(&max_votes)->default_value(max_votes),
         "max secondary votes threshold")
        ("discard,D", po::value<bool>(&discard)->default_value(discard),
         "whether to discard secondary votes")
        ("recursive,C", po::value<bool>(&recursive)->default_value(recursive),
         "whether to perform recursive splitting")
        ("flu,F", po::value<bool>(&flu)->default_value(flu),
         "Use FLU preset")
        ("allow-unassigned,A", po::value<bool>(&allow_unassigned)->default_value(allow_unassigned),
         "Allow unassigned reads in RCM");

    // Declare a group of options that are
    // allowed both on command line and in
    // config file but are not shown to user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("help-hidden", "show all options, including developers' ones");
        
    po::options_description cmdline_options("All command line options");
    cmdline_options
        .add(generic)
        .add(config)
        .add(hidden);
    po::options_description config_file_options;
    config_file_options
        .add(config)
        .add(hidden);
    po::options_description visible_options("Allowed options");
    visible_options
        .add(generic)
        .add(config);

    po::positional_options_description positional_options;
    positional_options
        .add("input-file", 1);
        
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
              .options(cmdline_options)
              //.positional(positional_options)
              .run(), vm);
    po::notify(vm);

    if (vm.count("help-hidden")) {
        std::cout << cmdline_options << std::endl;
        std::exit(0);
    }

    if (vm.count("help")) {
        std::cout << visible_options << std::endl;
        std::exit(0);
    }

    if (vm.count("version")) {
        std::cout << "IG Component Splitter, part of IgReC version "
                  << build_info::version
                  << "; git version: "
                  << build_info::git_hash7
                  << std::endl;
        std::exit(0);
    }

    if (vm.count("config-file")) {
        std::ifstream config(vm["config-file"].as<std::string>());
        if (!config) {
            std::cout << "could not open config file " << vm["config-file"].as<std::string>() << std::endl;
            std::exit(1);
        } else {
            po::store(po::parse_config_file(config, config_file_options), vm);
            po::notify(vm);
            // reparse cmd line again for update config defaults
            po::store(po::command_line_parser(argc, argv)
                      .options(cmdline_options)
                      .positional(positional_options)
                      .run(), vm);
            po::notify(vm);
        }
    }
}

}
