/**
 * \file hmc-tessellate.hpp
 *
 * \author 2015-2016 Sandia Clinic Team
 *
 * \brief
 *   Contains interface classes for the hmc-tessellate
 *   Voronoi tessellation package.
 */

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

	/// Value used to signify that no search radius is known for building a Cell.
    constexpr double NO_RADIUS = -1;

	/**
	 * \struct Particle
	 *
	 * \brief
	 *   A struct representing a single 3-D particle.
	 */
    struct Particle {
		/// An identifier for the particle.
        int id;

		/// A 3-D vector representing the position of the particle.
        Vector3 position;

		/// Constructor
        Particle(int id, double x, double y, double z)
            : id(id), position(x, y, z)
        { /* Done */ }

		/// Print information about the particle.
        friend std::ostream& operator<<(std::ostream& out, const Particle& p)
        {
            return out << "{id = " << p.id << ", pos = " << p.position << "}";
        }
    };

	/**
	 * \struct CellInfo
	 *
	 * \brief
	 *   Contains all the information for constructing a Cell in a Diagram.
	 */
    struct CellInfo {
		/// A pointer to the Diagram that contains the Cell.
        const Diagram* diagram;

		/// The index of the point in this Cell in the Diagram.
        SizeType index;

		/// The position of the point in this Cell.
        Vector3 position;

		/// The search radius for neighbors of the particle in the Cell. If no radius is known, set to NO_RADIUS.
        double searchRadius;

		/// Constructor
        CellInfo(const Diagram& diagram, SizeType index, Vector3 position, double searchRadius = NO_RADIUS)
            : diagram(&diagram), index(index), position(position), searchRadius(searchRadius)
        { /* Done */ }
    };

	/**
	 * \class Cell
	 *
	 * \brief
	 *   A Voronoi cell. Used to retrieve various parameters.
	 *
	 * /remarks
	 *   Each cell is only computed once, and the information is cached if it is needed later.
	 */
    class Cell {
        friend class Diagram;

    public:
		/// A pointer to the Diagram that contains this Cell.
        const Diagram* diagram_;

		/// The index of the point in this Cell in the Diagram.
        SizeType index_;

		/// The position of the point in this Cell.
        Vector3 position_;

		/// The Polyhedron representing the Voronoi cell.
        Polyhedron poly_;

		/// A flag that tells whether the Polyhedron has been computed.
        bool polyComputed_;

		/// The search radius for neighbors of the particle in the Cell. If no radius is known, set to NO_RADIUS.
        double searchRadius_;

		/// The ExpandingSearch object used to search for neighbors of the point in this Cell.
        SDS::ExpandingSearch search_;

    public:
		/// Default constructor
        Cell() : diagram_(nullptr), polyComputed_(false) {}

		/**
		 * \brief Constructor
		 *
		 * \param  diagram       A reference to the Diagram that contains this Cell.
		 * \param  index         The index of the point in this Cell in the Diagram.
		 * \param  position      The position of the point in this Cell.
		 * \param  searchRadius  The search radius for neighbors of the particle in the Cell. Optional.
		 *
		 * \remarks
		 *   If no search radius is given, neighbors will be found using an ExpaningSearch.
		 */
        Cell(const Diagram& diagram, SizeType index, Vector3 position, double searchRadius = NO_RADIUS)
            : diagram_(&diagram), index_(index),
              position_(position), polyComputed_(false), searchRadius_(searchRadius)
        { /* Done */ }

		/// Assignment operator
        Cell& operator=(const CellInfo& info)
        {
            clear();
            diagram_ = info.diagram;
            index_ = info.index;
            position_ = info.position;
            searchRadius_ = info.searchRadius;
            return *this;
        }

		/// Clear the Cell
        void clear()
        {
            // diagram_ = nullptr;
            poly_.clear();
            polyComputed_ = false;
        }

        const Particle& getParticle() const;

		/// Compute the volume of the Voronoi cell.
        double computeVolume()
        {
            ensurePolyComputed();
            return poly_.computeVolume();
        }

		/// Find the neighbors of the point in the Voronoi cell.
        template <typename Collection>
        void computeNeighbors(Collection& result)
        {
            ensurePolyComputed();
            poly_.computeNeighbors(result);
        }

		/**
		 * \brief
		 *   Checks whether the Polyhedron representing the Voronoi cell has been computed.
		 *   If it has not, performs the computation.
		 */
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

	/**
	 * \class Diagram
	 *
	 * \brief
	 *   A class that represents a Voronoi diagram.
	 *
	 * \details
	 *   Contains a vector of Particle objects, a spatial data structure, and a Polyhedron
	 *   representing the container.
	 */
    class Diagram {
        friend class Cell;

        // using SDS = TrivialSDS<SizeType>;
        // using SDS = ShellArray<SizeType>;

    protected:
		/// The Polyhedron that represents the container of the points.
        Polyhedron containerShape_;

		/// The vector of Particle objects in the Diagram.
        std::vector<Particle> particles_;

		/// The spatial data structure used to determine spatial relationships between points.
        SDS spatialStructure_;

		/// Default constructor
        Diagram() = default;

    public:
		/**
		 * \brief Constructor
		 *
		 * \param  containerShape  A Polyhedron represeting the shape of the container of points.
		 */
        Diagram(Polyhedron containerShape)
            : containerShape_(containerShape)
        { /* Done */ }

		/// Makes a Diagram in the shape of a cube.
        static Diagram cube(double xmin, double xMAX,
                            double ymin, double yMAX,
                            double zmin, double zMAX)
        {
            Diagram diagram;
            diagram.containerShape_.buildCube(xmin, xMAX, ymin, yMAX, zmin, zMAX);
            return diagram;
        }

		/// The number of Particle objects in the Diagram.
        SizeType size() const
        {
            return particles_.size();
        }

		/**
		 * \brief
		 *   Add a Particle to the Diagram.
		 *
		 * \param  id  An identifier for the Particle
		 * \param  x   The x-coordinate of the Particle
		 * \param  y   The y-coordinate of the Particle
		 * \param  z   The z-coordinate of the Particle
		 *
		 * \remarks
		 *   Should only be used before the initialize function is called.
		 */
        void addParticle(int id, double x, double y, double z)
        {
            particles_.emplace_back(id, x, y, z);
        }

		/**
		 * \brief
		 *   Initialize the Diagram.
		 *
		 * \remarks
		 *   Should only be called after all Particles have been added to the Diagram.
		 */
        void initialize()
        {
            spatialStructure_.initialize(0, particles_.size(),
                [this](SizeType index) -> Vector3 {
                    return particles_[index].position;
                }
            );
        }

		/**
		 * \brief
		 *   Initialize the Diagram.
		 *
		 * \param  fill  A lambda function that fills the Diagram with Particle objects
		 *
		 * \remarks
		 *   Should only be called after all Particles (other than those to be added by fill)
		 *   have been added to the Diagram.
		 */
        template <typename FillDiagramWithParticles>
        void initialize(FillDiagramWithParticles fill)
        {
            fill();
            initialize();
        }

		/**
		 * \brief
		 *   Get the CellInfo for a Particle.
		 *
		 * \param  index  The index of the Particle
		 */
        CellInfo getCell(SizeType index) const
        {
            VERIFY(index < particles_.size());
            return CellInfo(*this, index, particles_[index].position);
        }

		/**
		 * \brief
		 *   Get the CellInfo for a Particle using a specfic search radius.
		 *
		 * \param  index         The index of the Particle
		 * \param  searchRadius  The search radius to use for finding neighbors of the Particle
		 */
        CellInfo getCell(SizeType index, double searchRadius) const
        {
            VERIFY(index < particles_.size());
            VERIFY(searchRadius > 0);
            return CellInfo(*this, index, particles_[index].position, searchRadius);
        }

		/**
		 * \brief
		 *   Compute the Voronoi cell for a Particle.
		 *
		 * \param  particleIndex  The index of the Particle
		 * \param  position       The position of the Particle
		 * \param  poly           A reference to the Polyhedron to be cut to compute the cell
		 * \param  search         A reference to an ExpandingSearch for the Particle
		 * \param  searchRadius   The search radius to look for neighbors of the particle. Optional.
		 *
		 * \remarks
		 *   Uses the efficient expanding search strategy to find neighbors if searchRadius is not specified.
		 */
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


	/// A wrapper function for calling Diagram::ComputeVoronoiCell.
    void Cell::computeVoronoiCell()
    {
        diagram_->computeVoronoiCell(index_, position_, poly_, search_, searchRadius_);
    }

	/// A wrapper function for finding the Particle at an index of a Diagram.
    const Particle& Cell::getParticle() const
    {
        return diagram_->particles_[index_];
    }

}
