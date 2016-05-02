#pragma once

#include <vector>
#include <iterator>
#include <iostream>
#include <limits>

#include "vectormath.hpp"
#include "cellarray.hpp"

namespace hmc {
    template <typename iteratorOrInt>
    std::size_t iteratorOrIntDistance(iteratorOrInt a, iteratorOrInt b)
    {
        return std::distance(a, b);
    }

    template <>
    std::size_t iteratorOrIntDistance<std::size_t>(std::size_t a, std::size_t b)
    {
        return b - a;
    }

    template <typename PointHandle>
    class TrivialSDS {
    public:
        std::vector<PointHandle> points_;

    public:
        void clear()
        {
            points_.clear();
        }

        template <typename GetVectorFromHandle>
        void initialize(PointHandle begin, PointHandle end, GetVectorFromHandle)
        {
            points_.reserve(iteratorOrIntDistance(begin, end));
            for (; begin != end; ++begin) {
                points_.push_back(begin);
            }
        }

        const std::vector<PointHandle>& search(Vector3, double) const
        {
            return points_;
        }

        class ExpandingSearch {
            friend class TrivialSDS;

        protected:
            using Iterator = typename std::vector<PointHandle>::const_iterator;
            Iterator it_, end_;

            ExpandingSearch() = delete;
            ExpandingSearch(Iterator begin, Iterator end) : it_(begin), end_(end) {}

        public:
            Iterator begin() const { return it_; }
            Iterator end() const { return end_; }

            bool done() const { return it_ == end_; }
            void expandSearch(double) { it_ = end_; }
        };

        ExpandingSearch expandingSearch(Vector3) const
        {
            return ExpandingSearch {points_.begin(), points_.end()};
        }
    };

    template <typename PointHandle, typename SizeType_T = std::size_t>
    class CellArray {
    public:
        using SizeType = SizeType_T;
        spatial::Celery<SizeType> cellarray_;

    public:
        void clear()
        {
            // cellarray_.clear();
        }

        template <typename F>
        void initialize(PointHandle begin, PointHandle end, F getPointFromHandle)
        {
            cellarray_.initialize(begin, end, getPointFromHandle);
        }

        void search(Vector3 position, double searchRadius, std::vector<PointHandle>& result) const
        {
            cellarray_.findNeighborsInCellRadius(position.x, position.y, position.z, searchRadius, result);
        }

        class ExpandingSearch {
            friend class CellArray;

        protected:
            using Iterator = typename std::vector<PointHandle>::const_iterator;

            // const spatial::Celery<SizeType>& cellarray_;
            typename spatial::Celery<SizeType>::ExpandingSearch search_;
            // const Vector3 position_;
            std::vector<PointHandle> neighbors_;
            // unsigned shell_ = 0;
            bool done_ = false;

            // ExpandingSearch(const spatial::Celery<SizeType>& cellarray,
            //                 Vector3 position)
            //     //: search_(cellarray, position.x, position.y, position.z)
            // {
            //     startSearch(cellarray, position);
            //     //search_.expand(std::numeric_limits<double>::max(), neighbors_);
            // }

        public:
            ExpandingSearch() = default;

            Iterator begin() const { return neighbors_.begin(); }
            Iterator end() const { return neighbors_.end(); }

            bool done() const { return search_.done(); }

            void startSearch(const CellArray& cellarray, Vector3 position)
            {
                search_.initialize(cellarray.cellarray_, position.x, position.y, position.z);
                search_.expand(std::numeric_limits<double>::max(), neighbors_);
            }

            void expandSearch(double maxRadius)
            {
                neighbors_.clear();
                search_.expand(maxRadius, neighbors_);
            }
        };

        // ExpandingSearch expandingSearch(Vector3 position) const
        // {
        //     return ExpandingSearch {cellarray_, position};
        // }

        friend std::ostream& operator<<(std::ostream& out, const CellArray& c)
        {
            return out << c.cellarray_;
        }
    };

    template <typename PointHandle>
    class ShellArray {
    public:
        using SizeType = std::size_t;
        spatial::Celery<SizeType> cellarray_;

    public:
        void clear()
        {
            // cellarray_.clear();
        }

        template <typename F>
        void initialize(PointHandle begin, PointHandle end, F getPointFromHandle)
        {
            cellarray_.initialize(begin, end, getPointFromHandle);
        }

        const std::vector<PointHandle>& search(Vector3 position, double searchRadius) const
        {
            static std::vector<PointHandle> v;
            v.clear();
            cellarray_.findNeighborsInCellRadius(position.x, position.y, position.z, searchRadius, v);
            return v;
        }

        class ExpandingSearch {
            friend class ShellArray;

        protected:
            using Iterator = typename std::vector<PointHandle>::const_iterator;
            const spatial::Celery<SizeType>* cellarray_;
            Vector3 position_ = {0,0,0};
            static std::vector<PointHandle> neighbors_;
            unsigned shell_ = 0;
            bool done_ = false;

		public:
            ExpandingSearch() = default;

			void startSearch(const ShellArray<SizeType>& cellarray,
							 const Vector3 position)
            {
				cellarray_ = &cellarray.cellarray_;
				position_ = position;
				shell_ = 0;
				done_ = false;
                expandSearch(0);
            }

        public:
            Iterator begin() const { return neighbors_.begin(); }
            Iterator end() const { return neighbors_.end(); }

            bool done() const { return done_; }
            void expandSearch(double maxRadius)
            {
                neighbors_.clear();
                done_ = !cellarray_->findNeighborsInShell(
                            position_.x, position_.y, position_.z,
                            shell_, maxRadius, neighbors_
                        );
                ++shell_;
            }
        };

        ExpandingSearch expandingSearch(Vector3 position) const
        {
            return ExpandingSearch {cellarray_, position};
        }

        friend std::ostream& operator<<(std::ostream& out, const ShellArray& c)
        {
            return out << *c.cellarray_;
        }
    };

    template <typename PointHandle>
    std::vector<PointHandle> ShellArray<PointHandle>::ExpandingSearch::neighbors_ {};

}
