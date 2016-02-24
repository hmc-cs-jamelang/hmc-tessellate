#pragma once

#include <vector>
#include <limits>
#include <bitset>
#include "verification.hpp"

template <typename T, typename SizeType = short>
class StructPool {
public:
    using Index = SizeType;
    static constexpr Index INVALID_INDEX = std::numeric_limits<Index>::max();

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

    // This trick is used to create both const and nonconst versions
    // of the same function in one go.
    // The typename This will either be a pointer or a const pointer,
    // depending on whether the function calling this one is marked
    // const or not.
    template <typename This>
    static auto _chunk(This t, Index i) -> decltype(t->chunks_[i])
    {
        VERIFY(i < t->chunks_.size());
        return t->chunks_[i];
    }
    Chunk& chunk(Index i) {return _chunk(this, i);}
    const Chunk& chunk(Index i) const {return _chunk(this, i);}

public:
    StructPool()
        : first_available_(INVALID_INDEX), chunks_()
    { /* Done */ }

    StructPool(std::size_t reserve)
        : first_available_(INVALID_INDEX), chunks_(reserve)
    { /* Done */ }

    Index size() const
    {
        return chunks_.size();
    }

    Index capacity() const
    {
        return chunks_.capacity();
    }

    Index computeOccupancy() const
    {
        Index occupancy = 0;
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
        Index i = first_available_;
        if (i != INVALID_INDEX) {
            first_available_ = chunk(i).next_available();
        }
        else {
            i = chunks_.size();
            chunks_.push_back(Chunk());
        }

        chunk(i).construct(args...);
        return i;
    }

    void destroy(Index i)
    {
        Chunk& chunk = this->chunk(i);

        // We assume that there really is an object here,
        // so destruct it.
        chunk.destruct();

        // Add this index onto the linkedlist of available
        // indices.
        chunk.next_available() = first_available_;
        first_available_ = i;
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
                ++index;
            }
        }

        maybe_const_iterator(StructPoolPointer pool, Index index)
            : pool(pool), index(index)
        { /* Done */ }

    public:
        Index getIndex() const {return index;}

        bool active() const {return pool->active(index);}
        bool inactive() const {return pool->inactive(index);}
        bool marked() const {return pool->marked(index);}
        bool unmarked() const {return pool->unmarked(index);}

        void setActive(bool active) {pool->setActive(index, active);}
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
        typename std::enable_if<is_const_iterator, ValueReferenceType>
        operator*() const {return pool[index];}

        typename std::enable_if<!is_const_iterator, ValuePointerType>
        operator*() {return pool[index];}



        // We can only produce const versions if it's a const iterator
        typename std::enable_if<is_const_iterator, ValuePointerType>
        operator->() const {return &pool[index];}

        typename std::enable_if<!is_const_iterator, ValuePointerType>
        operator->() {return &pool[index];}



        // Prefix decrement
        maybe_const_iterator& operator--(){
            while (index-- > 0) {
                if (pool.active(index)) {
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
            while (++index < pool->size()) {
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
        return iterator(this, size());
    }

    const_iterator end() const
    {
        return const_iterator(this, size());
    }
};

