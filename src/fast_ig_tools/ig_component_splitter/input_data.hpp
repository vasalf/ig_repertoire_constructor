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

template<template<typename, typename> class map_type>
class InputData {
    map_type<std::string, Component> components_;
    map_type<std::string, std::string> rcm_;

    bool read_rcm_file(std::string filename) {
        std::ifstream rcm(filename);
        if (!rcm) {
            INFO("Could not open RCM file " << filename);
            return false;
        }
        std::string id, target;
        std::string line;
        while (std::getline(rcm, line)) {
            std::vector<std::string> strs;
            boost::split(strs, line, boost::is_any_of("\t"));
            for (std::string &s : strs) {
                boost::trim(s);
            }
            if (strs.size() == 2 && strs[0] != "" && strs[1] != "") {
                rcm_[strs[0]] = strs[1];
            }
        }
        return true;
    }

    void add_read_record(std::string s_id, Read &&r) {
        auto it = components_.find(s_id);

        if (it == components_.end()) {
            components_[s_id] = Component(s_id);
            it = components_.find(s_id);
        }

        it->second.push_back(std::make_shared<Read>(r));
    }

public:
    InputData() {}

    bool ReadData(const IgComponentSplitterConfig &cfg) {
        INFO("Input files: " << cfg.reads_file << ", " << cfg.rcm_file);

        std::vector<seqan::Dna5String> input_reads;
        std::vector<seqan::CharString> input_ids;
        seqan::SeqFileIn seqinput(cfg.reads_file.c_str());

        INFO("Reading input reads starts");
        seqan::readRecords(input_ids, input_reads, seqinput);
        INFO(input_reads.size() << " reads were extracted from " << cfg.reads_file);

        INFO("Reading read-cluster map starts");
        if (!read_rcm_file(cfg.rcm_file)) {
            return false;
        }

        size_t assigned_reads = 0;
        for (size_t i = 0; i != input_reads.size(); i++) {
            std::string id = seqan::toCString(input_ids[i]);
            if (rcm_.count(id)) {
                assigned_reads++;
            } else {
                TRACE("Read " << id <<  " not found in RCM " << cfg.rcm_file);
            }
            add_read_record(rcm_[id], Read(input_reads[i], input_ids[i]));
        }

        INFO(assigned_reads << " over " << input_reads.size() << " reads assigned");
        if (assigned_reads < input_reads.size() && !cfg.allow_unassigned) {
            INFO(input_reads.size() - assigned_reads << " unassigned reads in RCM " << cfg.rcm_file);
            INFO("Unassigned reads are not allowed");
            INFO("Pass option '--allow-unassigned=true' to allow");
            return false;
        }

        INFO(components_.size() << " clusters were extracted from " << cfg.rcm_file);

        size_t max_component_size = 0;
        for (const auto& kv : components_) {
            max_component_size = std::max(max_component_size, kv.second.size());
        }
        INFO("Size of maximal cluster: " << max_component_size);
        return true;
    }

    map_type<std::string, Component> &Components() {
        return components_;
    }
};
    
}
