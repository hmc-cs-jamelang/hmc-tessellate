#pragma once

#include <iostream>
#include <vector>
#include <limits>

#include "vectormath.hpp"
#include "polyhedron.hpp"
#include "verification.hpp"
#include "sds.hpp"

namespace hmc {
    using SizeType = std::size_t;
    using SDS = CellArray<SizeType>;

    class Diagram;
    class Cell;

    constexpr double NO_RADIUS = -1;

    struct Particle {
        int id;
        Vector3 position;

        Particle(int id, double x, double y, double z)
            : id(id), position(x, y, z)
        { /* Done */ }

        friend std::ostream& operator<<(std::ostream& out, const Particle& p)
        {
            return out << "{id = " << p.id << ", pos = " << p.position << "}";
        }
    };

    struct CellInfo {
        const Diagram* diagram;
        SizeType index;
        Vector3 position;
        double searchRadius;

        CellInfo(const Diagram& diagram, SizeType index, Vector3 position, double searchRadius = NO_RADIUS)
            : diagram(&diagram), index(index), position(position), searchRadius(searchRadius)
        { /* Done */ }
    };

    class Cell {
        friend class Diagram;

    public:
        const Diagram* diagram_;
        SizeType index_;
        Vector3 position_;
        Polyhedron poly_;
        bool polyComputed_;
        double searchRadius_;
        SDS::ExpandingSearch search_;

    public:
        Cell() : diagram_(nullptr), polyComputed_(false) {}

        Cell(const Diagram& diagram, SizeType index, Vector3 position, double searchRadius = NO_RADIUS)
            : diagram_(&diagram), index_(index),
              position_(position), polyComputed_(false), searchRadius_(searchRadius)
        { /* Done */ }

        Cell& operator=(const CellInfo& info)
        {
            clear();
            diagram_ = info.diagram;
            index_ = info.index;
            position_ = info.position;
            searchRadius_ = info.searchRadius;
            return *this;
        }

        void clear()
        {
            // diagram_ = nullptr;
            poly_.clear();
            polyComputed_ = false;
        }

        const Particle& getParticle() const;

        double computeVolume()
        {
            ensurePolyComputed();
            return poly_.computeVolume();
        }

        template <typename Collection>
        void computeNeighbors(Collection& result)
        {
            ensurePolyComputed();
            poly_.computeNeighbors(result);
        }

        void ensurePolyComputed()
        {
            VERIFY_EXIT(polyComputed_);

            if (!polyComputed_) {
                computeVoronoiCell();
                polyComputed_ = true;
            }
        }

    protected:
        void computeVoronoiCell();
    };

    class Diagram {
        friend class Cell;

        // using SDS = TrivialSDS<SizeType>;
        // using SDS = ShellArray<SizeType>;

    protected:
        Polyhedron containerShape_;
        std::vector<Particle> particles_;
        SDS spatialStructure_;

        Diagram() = default;

    public:
        Diagram(Polyhedron containerShape)
            : containerShape_(containerShape)
        { /* Done */ }

        static Diagram cube(double xmin, double xMAX,
                            double ymin, double yMAX,
                            double zmin, double zMAX)
        {
            Diagram diagram;
            diagram.containerShape_.buildCube(xmin, xMAX, ymin, yMAX, zmin, zMAX);
            return diagram;
        }

        SizeType size() const
        {
            return particles_.size();
        }

        void addParticle(int id, double x, double y, double z)
        {
            particles_.emplace_back(id, x, y, z);
        }

        void initialize()
        {
            spatialStructure_.initialize(0, particles_.size(),
                [this](SizeType index) -> Vector3 {
                    return particles_[index].position;
                }
            );
        }

        template <typename FillDiagramWithParticles>
        void initialize(FillDiagramWithParticles fill)
        {
            fill();
            initialize();
        }

        CellInfo getCell(SizeType index) const
        {
            VERIFY(index < particles_.size());
            return CellInfo(*this, index, particles_[index].position);
        }

        CellInfo getCell(SizeType index, double searchRadius) const
        {
            VERIFY(index < particles_.size());
            VERIFY(searchRadius > 0);
            return CellInfo(*this, index, particles_[index].position, searchRadius);
        }

        void computeVoronoiCell(SizeType particleIndex, Vector3 position, Polyhedron& poly,
                                SDS::ExpandingSearch& search, double searchRadius = NO_RADIUS) const
        {
            VERIFY(poly.isClear());
            poly = containerShape_;
            poly.translate(-position);

            // Use expanding search if no radius is given
            if (searchRadius == NO_RADIUS) {
    			searchRadius = std::numeric_limits<double>::max();
                for (search.startSearch(spatialStructure_, position);
                        !search.done();
    				    search.expandSearch(searchRadius = poly.maximumNeighborDistance()))
                {
                    for (SizeType index : search) {
                        Vector3 shiftedPosition = particles_[index].position - position;
                        if (index != particleIndex
                            && mag2(shiftedPosition) <= searchRadius)
                        {
                            poly.cutWithPlane(
                                index,
                                Plane::halfwayFromOriginTo(shiftedPosition)
                            );
                        }
                    }
                }
            }

            // Otherwise (given radius), do a basic search
            else {
                for (SizeType index : spatialStructure_.search(position, searchRadius)) {
                    if (index != particleIndex) {
                        poly.cutWithPlane(
                            index,
                            Plane::halfwayFromOriginTo(particles_[index].position - position)
                        );
                    }
                }
            }
        }
    };


    void Cell::computeVoronoiCell()
    {
        diagram_->computeVoronoiCell(index_, position_, poly_, search_, searchRadius_);
    }

    const Particle& Cell::getParticle() const
    {
        return diagram_->particles_[index_];
    }

}
