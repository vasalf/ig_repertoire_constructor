#pragma once

#include <vector>
#include <memory>

#include <seqan/sequence.h>

namespace ig_component_splitter {
    typedef seqan::Dna5           char_t;
    typedef seqan::String<char_t> seq_string_t;
    typedef seqan::CharString     id_string_t;

    struct Read {
        seq_string_t read;
        id_string_t id;
        Read(const seq_string_t &read_, const id_string_t &id_)
            :read(read_), id(id_) {}
    };
    typedef std::shared_ptr<Read> ReadPtr;

    class Component : public std::vector<ReadPtr> {
        id_string_t id_;
    public:
        seq_string_t consensus;

        Component() {}

        Component(const id_string_t& id)
            : std::vector<ReadPtr>(),
            id_(id) {}

        virtual ~Component() = default;

        const id_string_t &GetId() {
            return id_;
        }

        bool operator<(const Component &other) const {
            return id_ < other.id_;
        }

        std::vector<seq_string_t> GetReadSequences() const;
    };
}
