#pragma once

#include <vector>
#include <limits>
#include <bitset>
#include <iostream>
#include <typeinfo>
#include "verification.hpp"


struct StructPoolDefaultTag {};

template <
    typename T,
    typename SizeType = unsigned short,
    typename Tag = StructPoolDefaultTag
>
class StructPool {
public:
    class Index {
        friend class StructPool;

    private:
        SizeType value;

    public:
        constexpr Index() = default;
        explicit constexpr Index(SizeType value) : value(value) {}
        explicit constexpr operator SizeType() const {return value;}

        friend constexpr bool operator==(const Index lhs, const Index rhs)
        {
            return lhs.value == rhs.value;
        }

        friend constexpr bool operator!=(const Index lhs, const Index rhs)
        {
            return lhs.value != rhs.value;
        }

        friend std::ostream& operator<<(std::ostream& stream, const Index& index)
        {
            stream << index.value;
            return stream;
        }
    };

    static constexpr Index INVALID_INDEX
        {std::numeric_limits<SizeType>::max()};
    static constexpr SizeType MAX_SIZE
        {static_cast<SizeType>(INVALID_INDEX.value - 1)};

protected:
    union ObjectOrIndex {
        Index index;
        T object;

        ObjectOrIndex() {}
        ~ObjectOrIndex() {}
    };

public:
    // A Chunk is either active, and contains an object,
    // or it is inactive and contains an index pointing
    // to the next inactive chunk (or INVALID_INDEX if there are no more)
    // Thus, together the inactive chunks form a linked list,
    // while the active chunks store data.
    struct Chunk {
        ObjectOrIndex data;
        bool active;
        bool marked;

        Chunk() : data(), active(false), marked(false) {}

        ~Chunk()
        {
            VERIFY_EXIT(!active);
            if (active) {
                destruct();
            }
        }

        // In-place construction of a new object.
        template <typename... T_Constructor_Args>
        void construct(T_Constructor_Args... args)
        {
            VERIFY(!active);
            VERIFY_EXIT(active);

            // This constructs a new T object
            // inside the space already allocated
            // for this chunk.
            new(&data.object) T {args...};
            active = true;
        }

        // In-place destruction of current object.
        void destruct()
        {
            VERIFY(active);
            VERIFY_EXIT(!active);

            data.object.~T();
            active = false;
            marked = false;
        }

        // Get object held by active chunk
        inline T& object()
        {
            VERIFY(active);
            return data.object;
        }

        // Get index of next inactive chunk held by this inactive chunk
        inline Index& next_available()
        {
            VERIFY(!active);
            return data.index;
        }
    };

protected:
    Index first_available_;
    std::vector<Chunk> chunks_;

    // Helper functions used for sanity checks.
    VERIFICATION(
        // Ensures the StructPool is small enough for its SizeType
        // to reference all of its elements
        bool withinMaxSize() const
        {
            return chunks_.size() <= MAX_SIZE;
        }
    )

    SizeType toSizeType(std::size_t n) const
    {
        VERIFY(n < std::numeric_limits<SizeType>::max());
        return static_cast<SizeType>(n);
    }

    // This trick is used to create both const and nonconst versions
    // of the same function in one go.
    // The typename This will either be a pointer or a const pointer,
    // depending on whether the function calling this one is marked
    // const or not.
    template <typename This>
    static auto chunk_(This t, Index i) -> decltype(t->chunks_[i.value])
    {
    	VERIFY(i != INVALID_INDEX);
        VERIFY(i.value < t->chunks_.size());
        return t->chunks_[i.value];
    }
    Chunk& chunk(Index i) {return chunk_(this, i);}
    const Chunk& chunk(Index i) const {return chunk_(this, i);}

public:
    StructPool()
        : first_available_(INVALID_INDEX), chunks_()
    { /* Done */ }

    StructPool(std::size_t reserve)
        : first_available_(INVALID_INDEX), chunks_(reserve)
    { /* Done */ }

    SizeType size() const
    {
        return toSizeType(chunks_.size());
    }

    SizeType capacity() const
    {
        return toSizeType(chunks_.capacity());
    }

    SizeType computeOccupancy() const
    {
        SizeType occupancy = 0;
        for (auto it = begin(); it != end(); ++it) {
            occupancy += 1;
        }
        return occupancy;
    }

    bool active(Index i) const {return chunk(i).active;}
    bool inactive(Index i) const {return !chunk(i).active;}
    bool marked(Index i) const {return chunk(i).marked;}
    bool unmarked(Index i) const {return !chunk(i).marked;}

    void setActive(Index i, bool active) {chunk(i).active = active;}
    void setMarked(Index i, bool marked) {chunk(i).marked = marked;}

    inline T& operator[](Index i) {return chunk(i).object();}

    template<typename... T_Constructor_Args>
    Index create(T_Constructor_Args... args)
    {
        VERIFY_INVARIANT(withinMaxSize());

        Index i = first_available_;
        if (i != INVALID_INDEX) {
            first_available_ = chunk(i).next_available();
        }
        else {
            i = Index {static_cast<SizeType>(chunks_.size())};
            chunks_.push_back(Chunk());
        }

        chunk(i).construct(args...);
        return i;
    }

    void destroy(Index i)
    {
        VERIFY_INVARIANT(withinMaxSize());

        Chunk& chunk = this->chunk(i);

        // We assume that there really is an object here,
        // so destruct it.
        chunk.destruct();

        // Add this index onto the linkedlist of available
        // indices.
        chunk.next_available() = first_available_;
        first_available_ = i;
    }

    SizeType destroyUnmarkedAndResetFlags()
    {
        SizeType totalDestroyed = 0;
        for (auto it = begin(); it != end(); ++it) {
            if (!it.marked()) {
                destroy(it.getIndex());
                totalDestroyed += 1;
            } else {
                it.setMarked(false);
            }
        }
        return totalDestroyed;
    }


    // ==============
    // |  Iterator  |
    // ==============
    // Adapted from http://www.sj-vs.net/c-implementing-const_iterator-and-non-const-iterator-without-code-duplication/
    template<bool is_const_iterator>
    class maybe_const_iterator : public std::iterator<std::bidirectional_iterator_tag, T> {
        friend class StructPool<T>;
        friend class maybe_const_iterator<true>;

    public:
        using StructPoolPointer =   typename std::conditional<is_const_iterator,
                                        const StructPool<T>*,
                                        StructPool<T>*
                                    >::type;
        using ValueReferenceType =  typename std::conditional<is_const_iterator,
                                        const T&,
                                        T&
                                    >::type;
        using ValuePointerType =    typename std::conditional<is_const_iterator,
                                        const T*,
                                        T*
                                    >::type;

    protected:
        StructPoolPointer pool;
        Index index;

        maybe_const_iterator(StructPoolPointer pool)
            : pool(pool), index(0)
        {
            while (pool->inactive(index)) {
                ++index.value;
            }
        }

        maybe_const_iterator(StructPoolPointer pool, Index index)
            : pool(pool), index(index)
        { /* Done */ }


        bool active() const {return pool->active(index);}
        bool inactive() const {return pool->inactive(index);}
        void setActive(bool active) {pool->setActive(index, active);}

    public:
        Index getIndex() const {return index;}

        bool marked() const {return pool->marked(index);}
        bool unmarked() const {return pool->unmarked(index);}
        void setMarked(bool marked) {pool->setMarked(index, marked);}

        maybe_const_iterator(const maybe_const_iterator<false>& other)
            : pool(other.pool), index(other.index)
        { /* Done */ }

        bool operator== (const maybe_const_iterator& other) const
        {
            return index == other.index && pool == other.pool;
        }

        bool operator!= (const maybe_const_iterator& other) const
        {
            return !(*this == other);
        }

    public:

        // We can only produce a const version if it's a const iterator
        ValueReferenceType
        operator*() const {return pool->operator[](index);}

        ValueReferenceType
        operator*() {return pool->operator[](index);}



        // We can only produce const versions if it's a const iterator
        ValuePointerType
        operator->() const {return &pool->operator[](index);}

        ValuePointerType
        operator->() {return &pool->operator[](index);}



        // Prefix decrement
        maybe_const_iterator& operator--(){
            while (index.value-- > 0) {
                if (pool->active(index)) {
                    return *this;
                }
            }
            return *this;
        }

        // Postfix decrement
        maybe_const_iterator operator--(int){
            const maybe_const_iterator old(*this);
            --(*this);
            return old;
        }

        // Prefix increment
        maybe_const_iterator &operator++(){
            while (++index.value < pool->size()) {
                if (pool->active(index)) {
                    return *this;
                }
            }
            return *this;
        }

        // Postfix increment
        maybe_const_iterator operator++(int){
            const maybe_const_iterator old(*this);
            ++(*this);
            return old;
        }
    };

    using iterator = maybe_const_iterator<false>;
    using const_iterator = maybe_const_iterator<true>;

    iterator begin()
    {
        return iterator(this);
    }

    const_iterator begin() const
    {
        return const_iterator(this);
    }

    iterator end()
    {
        return iterator(this, Index {size()});
    }

    const_iterator end() const
    {
        return const_iterator(this, Index {size()});
    }

    // Auxiliary debug-type functions
    std::size_t get_memory_usage() {
        std::size_t memory = 0;
        memory += sizeof(Chunk) * chunks_.capacity();
        return memory;
    }
};

template <typename T, typename SizeType, typename Tag>
constexpr typename StructPool<T, SizeType, Tag>::Index
      StructPool<T, SizeType, Tag>::INVALID_INDEX;

template <typename T, typename SizeType, typename Tag>
constexpr SizeType StructPool<T, SizeType, Tag>::MAX_SIZE;
