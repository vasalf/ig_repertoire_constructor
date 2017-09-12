#pragma once

#include <string>
#include <cstddef>

namespace ig_component_splitter {

struct IgComponentSplitterConfig {
    // Paths
    std::string reads_file;
    std::string rcm_file;
    std::string output_file;
    std::string output_rcm_file;
    std::string config_file;

    // Launch options
    int nthreads = 4;
    size_t max_votes = std::numeric_limits<size_t>::max() / 2;
    bool discard = false;
    bool recursive = false;
    bool flu = false;
    bool allow_unassigned = false;

    void Write();
    
    IgComponentSplitterConfig() {}
    
    void ReadOptions(int argc, char **argv);
};

}
