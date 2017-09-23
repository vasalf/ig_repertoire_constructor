#pragma once

#include <map>
#include <unordered_map>
#include <fstream>
#include <string>

#include "component.hpp"

#include <boost/algorithm/string.hpp>

namespace ig_component_splitter {

template<typename key, typename value>
class StdMapTag : public std::map<key, value> {
public:
    static const bool is_sorted_container = false;
};

template<typename key, typename value>
class StdUnorderedMapTag : public std::unordered_map<key, value> {
public:
    static const bool is_sorted_container = true;
};

template<typename key, typename value>
using map_type = StdUnorderedMapTag<key, value>;
    
class InputData {
    map_type<std::string, Component> components_;
    map_type<std::string, std::string> rcm_;

    bool read_rcm_file(std::string filename);

    void add_read_record(std::string s_id, Read &&r);

public:
    InputData() {}

    bool ReadData(const IgComponentSplitterConfig &cfg);

    map_type<std::string, Component> &Components() {
        return components_;
    }
};
    
}
