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

    template <typename PointHandle>
    class CellArray {
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
            // v.clear();
            // cellarray_.findNeighborsInCellRadius(position.x, position.y, position.z, searchRadius, v);
            return v;
        }

        class ExpandingSearch {
            friend class CellArray;

        protected:
            using Iterator = typename std::vector<PointHandle>::const_iterator;

            // const spatial::Celery<SizeType>& cellarray_;
            spatial::Celery<SizeType>::ExpandingSearch search_;
            // const Vector3 position_;
            std::vector<PointHandle> neighbors_;
            // unsigned shell_ = 0;
            bool done_ = false;

            ExpandingSearch() = delete;
            ExpandingSearch(const spatial::Celery<SizeType>& cellarray,
                            Vector3 position)
                : search_(cellarray, position.x, position.y, position.z)
            {
                search_.expand(std::numeric_limits<double>::max(), neighbors_);
            }

        public:
            Iterator begin() const { return neighbors_.begin(); }
            Iterator end() const { return neighbors_.end(); }

            bool done() const { return search_.done(); }
            void expandSearch(double maxRadius)
            {
                neighbors_.clear();
                search_.expand(maxRadius, neighbors_);
            }
        };

        ExpandingSearch expandingSearch(Vector3 position) const
        {
            return ExpandingSearch {cellarray_, position};
        }

        friend std::ostream& operator<<(std::ostream& out, const CellArray& c)
        {
            return out << c.cellarray_;
        }
    };

}
