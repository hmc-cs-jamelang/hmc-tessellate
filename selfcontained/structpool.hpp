#pragma once

#include <vector>
#include <limits>
#include <bitset>
#include <iostream>
#include <typeinfo>
#include <type_traits>

#include "verification.hpp"
#include "wrapper.hpp"


namespace hmc {
    namespace detail {
        struct StructPoolDefaultTag {};
    }

    template <
        typename T,
        typename SizeTypeParam = unsigned short,
        typename Tag = detail::StructPoolDefaultTag
    >
    class StructPool {
    public:
        using SizeType = SizeTypeParam;

        // The second argument is a tag, so that indices of different
        // kinds of StructPool are considered different types
        // so that they can't get mixed up
        using Index = Wrapper<SizeType, StructPool<T, SizeType, Tag>>;

        static constexpr Index INVALID_INDEX
            {std::numeric_limits<SizeType>::max()};
        static constexpr SizeType MAX_SIZE
            {static_cast<SizeType>(static_cast<SizeType>(INVALID_INDEX) - 1)};

    public:
        // A Chunk is either active, and contains an object,
        // or it is inactive and contains an index pointing
        // to the next inactive chunk (or INVALID_INDEX if there are no more)
        // Thus, together the inactive chunks form a linked list,
        // while the active chunks store data.
        struct Chunk {
            union ObjectOrIndex {
                Index index;
                T object;

                ObjectOrIndex() {}
                ~ObjectOrIndex() {}
            } data;
            bool active;

            Chunk() : data(), active(false) {}

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
        // ==================
        // |  Data Members  |
        // ==================
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

        void clear()
        {
            chunks_.clear();
            first_available_ = INVALID_INDEX;
        }

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

        void setActive(Index i, bool active) {chunk(i).active = active;}

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
                chunks_.emplace_back();
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
            StructPoolPointer pool_;
            Index index_;

            maybe_const_iterator(StructPoolPointer pool_)
                : pool_(pool_), index_(static_cast<SizeType>(0))
            {
                while (pool_->inactive(index_)) {
                    ++index_.value;
                }
            }

            maybe_const_iterator(StructPoolPointer pool_, Index index_)
                : pool_(pool_), index_(index_)
            { /* Done */ }


            bool active() const {return pool_->active(index_);}
            bool inactive() const {return pool_->inactive(index_);}
            void setActive(bool active) {pool_->setActive(index_, active);}

        public:
            Index index() const {return index_;}

            maybe_const_iterator(const maybe_const_iterator<false>& other)
                : pool_(other.pool_), index_(other.index_)
            { /* Done */ }

            bool operator== (const maybe_const_iterator& other) const
            {
                return index_ == other.index_ && pool_ == other.pool_;
            }

            bool operator!= (const maybe_const_iterator& other) const
            {
                return !(*this == other);
            }

        public:
            ValueReferenceType operator*() const {return (*pool_)[index_];}
            ValueReferenceType operator*() {return (*pool_)[index_];}

            ValuePointerType operator->() const {return &(*pool_)[index_];}
            ValuePointerType operator->() {return &(*pool_)[index_];}

            // Prefix decrement
            maybe_const_iterator& operator--(){
                while (index_.value-- > 0) {
                    if (pool_.active(index_)) {
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
                while (++index_.value < pool_->size()) {
                    if (pool_->active(index_)) {
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
}
