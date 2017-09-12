#include "component.hpp"

namespace ig_component_splitter {

std::vector<seq_string_t> Component::GetReadSequences() const {
    std::vector<seq_string_t> res;
    for (auto it = begin(); it != end(); it++) {
        res.push_back((*it)->read);
    }
    return res;
}
    
}
