#pragma once

#include <vector>
#include <limits>
#include <assert.h>
#include "verification.hpp"

template <typename T>
class StructPool {
public:
    typedef std::size_t Index;
    static constexpr Index INVALID_INDEX = std::numeric_limits<Index>::max();

protected:
    union ObjectOrIndex {
        Index index;
        T object;

        ~ObjectOrIndex() {}
    };

public:
    struct Chunk {
        VERIFICATION (
            bool active;
        )

        ObjectOrIndex data;

        Chunk() :
        VERIFICATION (
            active(false),
        )
            data()
        {}

        template <typename... T_Constructor_Args>
        void construct(T_Constructor_Args... args)
        {
            VERIFICATION (
                assert(!active);
            )
            VERIFICATION_EXIT (
                active = true;
            )

            // This constructs a new T object
            // inside the space already allocated
            // for this chunk.
            new(&data.object) T {args...};
        }

        void destruct()
        {
            VERIFICATION (
                assert(active);
            )
            VERIFICATION_EXIT (
                active = false;
            )

            data.object.~T();
        }

        inline T& object()
        {
            VERIFICATION (
                assert(active);
            )

            return data.object;
        }

        inline Index& next_available()
        {
            VERIFICATION (
                assert(!active);
            )

            return data.index;
        }
    };

protected:
    Index first_available_;
    std::vector<Chunk> chunks_;

public:

    StructPool()
        : first_available_(INVALID_INDEX), chunks_()
    { /* Done */ }

    StructPool(std::size_t reserve)
        : first_available_(INVALID_INDEX), chunks_(reserve)
    { /* Done */ }


    T& operator[](Index i)
    {
        VERIFICATION(
            assert(i < chunks_.size());
        )

        return chunks_[i].object();
    }

    template<typename... T_Constructor_Args>
    Index create(T_Constructor_Args... args)
    {
        Index i = first_available_;
        if (i != INVALID_INDEX) {
            assert(i < chunks_.size());

            first_available_ = chunks_[i].next_available();
        }
        else {
            i = chunks_.size();
            chunks_.push_back(Chunk());
        }

        chunks_[i].construct(args...);
        return i;
    }

    void destroy(Index i)
    {
        VERIFICATION (
            assert(i < chunks_.size());
        )

        Chunk& chunk = chunks_[i];

        // We assume that there really is an object here,
        // so destruct it.
        chunk.destruct();

        // Add this index onto the linkedlist of available
        // indices.
        chunk.next_available() = first_available_;
        first_available_ = i;
    }
};

