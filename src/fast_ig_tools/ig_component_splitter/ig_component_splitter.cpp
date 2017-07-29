#include <string>
#include <limits>
#include <cstddef>
#include <iostream>
#include <functional>
#include <memory>
#include <sstream>
#include <array>
#include <utility>

#include <chrono>
#include <atomic>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "utils.hpp"
#include <fast_ig_tools.hpp>
#include <build_info.hpp>
#include <ig_final_alignment.hpp>
#include <ig_matcher.hpp>

#include "config.hpp"
#include "component.hpp"
#include "input_data.hpp"
#include "splitter_algorithm.hpp"
#include "processor.hpp"

#include <seqan/seq_io.h>

using namespace ig_component_splitter;

// TODO: write a parallel processor here

int main(int argc, char *argv[]) {
    segfault_handler sh;
    perf_counter pc;
    create_console_logger("");

    INFO("Command line: " << join_cmd_line(argc, argv));

    IgComponentSplitterConfig cfg;

    try {
        cfg.ReadOptions(argc, argv);
    } catch(boost::program_options::error &e) {
        std::cout << "Parser error: " << e.what() << std::endl;
        return 1;
    } catch(std::exception &e) {
        std::cout << e.what() << std::endl;
        return 1;
    }
    cfg.Write();

    InputData<StdUnorderedMapTag> input;
    if (!input.ReadData(cfg))
        return 1;

    IgComponentSplitterParallelProcessor<StdUnorderedMapTag> processor(cfg, input);
    processor.Process();

    INFO("Running time: " << running_time_format(pc));
    
    return 0;
}
