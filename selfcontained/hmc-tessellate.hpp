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

    inline void computeVoronoiCell(const Diagram&, SizeType, Vector3, Polyhedron&, SDS::ExpandingSearch&);
    inline void computeVoronoiCell(const Diagram&, SizeType, Vector3, Polyhedron&, SDS::ExpandingSearch&, double);

    inline const Particle& getParticle(const Diagram&, SizeType);

    struct CellInfo {
        const Diagram* diagram;
        SizeType index;
        Vector3 position;
        double searchRadius;

        CellInfo(const Diagram& diagram, SizeType index, Vector3 position, double searchRadius = 0)
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

        Cell(const Diagram& diagram, SizeType index, Vector3 position, double searchRadius = 0)
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
            diagram_ = nullptr;
            poly_.clear();
            polyComputed_ = false;
        }

        const Particle& getParticle() const
        {
            return ::hmc::getParticle(*diagram_, index_);
        }

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
                if (searchRadius_ == 0) {
                    computeVoronoiCell(*diagram_, index_, position_, poly_, search_);
                }
                else {
                    computeVoronoiCell(*diagram_, index_, position_, poly_, search_, searchRadius_);
                }
                polyComputed_ = true;
            }
        }
    };

    class Diagram {
        friend void computeVoronoiCell(const Diagram&, SizeType, Vector3, Polyhedron&, SDS::ExpandingSearch&);
        friend void computeVoronoiCell(const Diagram&, SizeType, Vector3, Polyhedron&, SDS::ExpandingSearch&, double);

        friend const Particle& getParticle(const Diagram&, SizeType);

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

            // std::cerr << spatialStructure_ << std::endl;
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

        void computeVoronoiCell(SizeType particleIndex, Vector3 position, Polyhedron& poly, SDS::ExpandingSearch& search) const
        {
            VERIFY(poly.isClear());
            poly = containerShape_;

            // auto furthestNeighborDistance = [&](){
            //     std::vector<SizeType> ns;
            //     poly.computeNeighbors(ns);
            //     double dist = 0;
            //     int is = -1;
            //     for (auto i : ns) if (i != -1) {
            //         double d = distance(particles_[i].position, position);
            //         // std::cerr << "  :: " << i << " -> " << d << std::endl;
            //         if (d > dist) {
            //             dist = d;
            //             is = i;
            //         }
            //     }
            //     // std::cerr << "[" << is << "] ";
            //     return dist;
            // };
			decltype(std::chrono::high_resolution_clock::now()) t1;

			double rad = std::numeric_limits<double>::max();
            for (search.startSearch(spatialStructure_, position);
                 !search.done();
				 t1 = std::chrono::high_resolution_clock::now(),
					 search.expandSearch(rad = poly.maximumNeighborDistance(position)),
					 expansionTime += std::chrono::high_resolution_clock::now() - t1)
            {
				++expansions;
                for (SizeType index : search) {
                    if (index != particleIndex && squaredDistance(position, particles_[index].position) <= rad) {
                        neededCuts += poly.cutWithPlane(
                            index,
                            Plane::between(position, particles_[index].position)
                        );
						++totalCuts;
                    }
                }
            }
            // std::cerr << particleIndex << " done at " << rad << std::endl;
            // std::cerr << "  " << particleIndex << " furthest neighbor "
            //     << furthestNeighborDistance() << " ("
            //     << 1.00001 * furthestNeighborDistance() << ")" << std::endl;
        }

        void computeVoronoiCell(SizeType particleIndex, Vector3 position,
                                Polyhedron& poly, SDS::ExpandingSearch& search,
                                double searchRadius) const
        {
            // std::cerr << particleIndex << " search radius is " << searchRadius << std::endl;
            constexpr bool tryExpanding = true;
            if (!tryExpanding) {
                VERIFY(poly.isClear());
                poly = containerShape_;

                // std::cerr << "  " << particleIndex
                //     // << " (" << spatialStructure_.cellarray_.getCells()[particleIndex]
                //     // << ") "
                //     << " ns:";
                for (SizeType index : spatialStructure_.search(position, searchRadius)) {
                    // std::cerr << " " << index;
                    if (index != particleIndex) {
                        poly.cutWithPlane(
                            index,
                            Plane::between(position, particles_[index].position)
                        );
                    }
                }
                // std::cerr << std::endl;
            }

            else {
                VERIFY(poly.isClear());
                poly = containerShape_;

                for (search.startSearch(spatialStructure_, position);
                     !search.done();
                     search.expandSearch(searchRadius))
                {
                    for (SizeType index : search) {
                        if (index != particleIndex) {
                            poly.cutWithPlane(
                                index,
                                Plane::between(position, particles_[index].position)
                            );
                        }
                    }
                }
            }

            // auto furthestNeighborDistance = [&](){
            //     std::vector<SizeType> ns;
            //     poly.computeNeighbors(ns);
            //     double dist = 0;
            //     int is = -1;
            //     for (auto i : ns) if (i != -1) {
            //         double d = distance(particles_[i].position, position);
            //         // std::cerr << "  :: " << i << " -> " << d << std::endl;
            //         if (d > dist) {
            //             dist = d;
            //             is = i;
            //         }
            //     }
            //     // std::cerr << "[" << is << "] ";
            //     return dist;
            // };

            // std::cerr << "  " << particleIndex << " furthest neighbor "
            //     << furthestNeighborDistance() << " ("
            //     << 1.00001 * furthestNeighborDistance() << ")" << std::endl;
        }
    };

    inline void computeVoronoiCell(const Diagram& diagram,
                                   SizeType index,
                                   Vector3 position,
                                   Polyhedron& poly,
                                   SDS::ExpandingSearch& search)
    {
        diagram.computeVoronoiCell(index, position, poly, search);
    }

    inline void computeVoronoiCell(const Diagram& diagram,
                                   SizeType index,
                                   Vector3 position,
                                   Polyhedron& poly,
                                   SDS::ExpandingSearch& search,
                                   double searchRadius)
    {
        diagram.computeVoronoiCell(index, position, poly, search, searchRadius);
    }

    inline const Particle& getParticle(const Diagram& diagram, SizeType index)
    {
        return diagram.particles_[index];
    }

}
