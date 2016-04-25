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

	using TargetGroup = std::unordered_set<int>;
	using SourceGroup = std::vector<int>;

    class Diagram;
    class Cell;

	/// Value used to signify that no group has been specified.
	constexpr int DEFAULT_GROUP = -1;

	/// Value used to signify that no search radius is known for building a Cell.
    constexpr double NO_RADIUS = -1;

	/// Value used to signify that a particle does not exist in a Diagram.
	constexpr size_t NO_INDEX = std::numeric_limits<size_t>::max();

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

		/// The target group of the Cell.
		TargetGroup targetGroup;

		/// Constructor
        CellInfo(const Diagram& diagram, SizeType index, Vector3 position, double searchRadius = NO_RADIUS)
            : diagram(&diagram), index(index), position(position), searchRadius(searchRadius)
        { /* Done */ }

		/// Constructor
        CellInfo(const Diagram& diagram, SizeType index, Vector3 position, TargetGroup targetGroup, double searchRadius = NO_RADIUS)
            : diagram(&diagram), index(index), position(position), targetGroup(targetGroup), searchRadius(searchRadius)
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

		// The target group of the Cell.
		TargetGroup target_group_;

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
		 * \param[in]  diagram       A reference to the Diagram that contains this Cell.
		 * \param[in]  index         The index of the point in this Cell in the Diagram.
		 * \param[in]  position      The position of the point in this Cell.
		 * \param[in]  searchRadius  The search radius for neighbors of the particle in the Cell. Optional.
		 *
		 * \remarks
		 *   If no search radius is given, neighbors will be found using an ExpaningSearch.
		 */
        Cell(const Diagram& diagram, SizeType index, Vector3 position, double searchRadius = NO_RADIUS)
            : diagram_(&diagram), index_(index),
              position_(position), polyComputed_(false), searchRadius_(searchRadius)
        { /* Done */ }

		/**
		 * \brief Constructor
		 *
		 * \param[in]  diagram       A reference to the Diagram that contains this Cell.
		 * \param[in]  index         The index of the point in this Cell in the Diagram.
		 * \param[in]  position      The position of the point in this Cell.
		 * \param[in]  searchRadius  The search radius for neighbors of the particle in the Cell. Optional.
		 * \param[in]  targetGroup   The target group of the Cell.
		 *
		 * \remarks
		 *   If no search radius is given, neighbors will be found using an ExpaningSearch.
		 */
        Cell(const Diagram& diagram, SizeType index, Vector3 position, TargetGroup targetGroup, double searchRadius = NO_RADIUS)
            : diagram_(&diagram), index_(index),
              position_(position), polyComputed_(false), searchRadius_(searchRadius),
			  target_group_(targetGroup)
        { /* Done */ }

		/**
		 * \brief Constructor
		 *
		 * \param[in]  info  A CellInfo object
		 *
		 * \remark
		 *   Creates the Cell specified by the CellInfo.
		 */
		Cell(const CellInfo info)
			: Cell()
		{
			*this = info;
		}
		
		/// Assignment operator
        Cell& operator=(const CellInfo info)
        {
            clear();
            diagram_ = info.diagram;
            index_ = info.index;
            position_ = info.position;
			target_group_ = info.targetGroup;
            searchRadius_ = info.searchRadius;
            return *this;
        }

		/// Clear the Cell
        void clear()
        {
            // diagram_ = nullptr;
            poly_.clear();
            polyComputed_ = false;
			target_group_.clear();
        }

        const Particle& getParticle() const;

		/**
		 * \brief
		 *   Compute the volume of the Voronoi cell.
		 *
		 * \return
		 *   The volume of the Voronoi cell.
		 */
        double computeVolume()
        {
            ensurePolyComputed();
            return poly_.computeVolume();
        }

		/**
		 * \brief
		 *   Find the neighbors of the point in the Voronoi cell.
		 *
		 * \param[out]  result  A container for storing the neighbors
		 */
        template <typename Collection>
        void computeNeighbors(Collection& result)
        {
            ensurePolyComputed();
            poly_.computeNeighbors(result);
        }

		/**
		 * \brief
		 *   Compute the vertices of the Voronoi cell.
		 *
		 * \param[out]  result  A container for storing the neighbors
		 */
		template<typename Collection>
		void computeVertices(Collection& result)
		{
			ensurePolyComputed();
			poly_.computeVertices(result);
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

		/// The vector of groups that contain each particle. The group at an index corresponds
		/// to the particle at that index.
		std::vector<int> groups_;

		/// The spatial data structure used to determine spatial relationships between points.
        SDS spatialStructure_;

		/// Default constructor
        Diagram() = default;

    public:
		/**
		 * \brief Constructor
		 *
		 * \param[in]  containerShape  A Polyhedron represeting the shape of the container of points.
		 */
        Diagram(Polyhedron containerShape)
            : containerShape_(containerShape)
        { /* Done */ }

		/**
		 * \brief
		 *   Makes a Diagram in the shape of a cube.
		 *
		 * \param[in]  xmin  The lower x boundary of the Diagram
		 * \param[in]  xMAX  The upper x boundary of the Diagram
		 * \param[in]  ymin  The lower y boundary of the Diagram
		 * \param[in]  yMAX  The upper y boundary of the Diagram
		 * \param[in]  zmin  The lower z boundary of the Diagram
		 * \param[in]  zMAX  The upper z boundary of the Diagram
		 *
		 * \return
		 *   A Diagram with the specified boundaries.
		 */
        static Diagram cube(double xmin, double xMAX,
                            double ymin, double yMAX,
                            double zmin, double zMAX)
        {
            Diagram diagram;
            diagram.containerShape_.buildCube(xmin, xMAX, ymin, yMAX, zmin, zMAX);
            return diagram;
        }

		/**
		 * \brief
		 *   Returns the number of Particle objects in the Diagram.
		 *
		 * \return
		 *   The number of Particle objects in the Diagram.
		 */
        SizeType size() const
        {
            return particles_.size();
        }

		/**
		 * \brief
		 *   Add a Particle to the Diagram.
		 *
		 * \param[in]  id     An identifier for the Particle
		 * \param[in]  x      The x-coordinate of the Particle
		 * \param[in]  y      The y-coordinate of the Particle
		 * \param[in]  z      The z-coordinate of the Particle
		 * \param[in]  group  The group containing the Particle. Optional.
		 *
		 * \remarks
		 *   Should only be used before the initialize function is called.
		 */
        void addParticle(int id, double x, double y, double z, int group = DEFAULT_GROUP)
        {
            particles_.emplace_back(id, x, y, z);
			groups_.push_back(group);
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
		 * \param[in]  fill  A lambda function that fills the Diagram with Particle objects
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
		 * \param[in]  index  The index of the Particle
		 *
		 * \return
		 *   The CellInfo for the particle.
		 */
        CellInfo getCell(SizeType index) const
        {
            VERIFY(index < particles_.size());
            return CellInfo(*this, index, particles_[index].position);
        }

		/**
		 * \brief
		 *   Get the CellInfo for a Particle.
		 *
		 * \param[in]  index        The index of the Particle
		 * \param[in]  targetGroup  The TargetGroup for cutting the Cell
		 *
		 * \return
		 *   The CellInfo for the particle.
		 */
        CellInfo getCell(SizeType index, TargetGroup targetGroup) const
        {
            VERIFY(index < particles_.size());
            return CellInfo(*this, index, particles_[index].position, targetGroup);
        }

		/**
		 * \brief
		 *   Get the CellInfo for a Particle using a specfic search radius.
		 *
		 * \param[in]  index         The index of the Particle
		 * \param[in]  searchRadius  The search radius to use for finding neighbors of the Particle
		 *
		 * \return
		 *   The CellInfo for the particle.
		 */
        CellInfo getCell(SizeType index, double searchRadius) const
        {
            VERIFY(index < particles_.size());
            VERIFY(searchRadius > 0);
            return CellInfo(*this, index, particles_[index].position, searchRadius);
        }

		/**
		 * \brief
		 *   Get the CellInfo for a Particle using a specfic search radius.
		 *
		 * \param[in]  index         The index of the Particle
		 * \param[in]  targetGroup   The TargetGroup for cutting the Cell
		 * \param[in]  searchRadius  The search radius to use for finding neighbors of the Particle
		 *
		 * \return
		 *   The CellInfo for the particle.
		 */
        CellInfo getCell(SizeType index, TargetGroup targetGroup, double searchRadius) const
        {
            VERIFY(index < particles_.size());
            VERIFY(searchRadius > 0);
            return CellInfo(*this, index, particles_[index].position, targetGroup, searchRadius);
        }

		/**
		 * \brief
		 *   Get the CellInfo for a hypothetical particle.
		 *
		 * \param[in]  x  The x-coordinate of the particle
		 * \param[in]  y  The y-coordinate of the particle
		 * \param[in]  z  The z-coordinate of the particle
		 *
		 * \return
		 *   The CellInfo for the particle.
		 *
		 * \remark
		 *   Used to find a cell for a particle that is not actually in the Diagram.
		 */
        CellInfo getCell(double x, double y, double z) const
        {
            VERIFY(index < particles_.size());
            return CellInfo(*this, NO_INDEX, {x, y, z});
        }

		/**
		 * \brief
		 *   Get the CellInfo for a Particle.
		 *
		 * \param[in]  x            The x-coordinate of the particle
		 * \param[in]  y            The y-coordinate of the particle
		 * \param[in]  z            The z-coordinate of the particle
		 * \param[in]  targetGroup  The TargetGroup for cutting the Cell
		 *
		 * \return
		 *   The CellInfo for the particle.
		 *
		 * \remark
		 *   Used to find a cell for a particle that is not actually in the Diagram.
		 */
        CellInfo getCell(double x, double y, double z, TargetGroup targetGroup) const
        {
            VERIFY(index < particles_.size());
            return CellInfo(*this, NO_INDEX, {x, y, z}, targetGroup);
        }

		/**
		 * \brief
		 *   Get the CellInfo for a Particle using a specfic search radius.
		 *
		 * \param[in]  x             The x-coordinate of the particle
		 * \param[in]  y             The y-coordinate of the particle
		 * \param[in]  z             The z-coordinate of the particle
		 * \param[in]  searchRadius  The search radius to use for finding neighbors of the Particle
		 *
		 * \return
		 *   The CellInfo for the particle.
		 *
		 * \remark
		 *   Used to find a cell for a particle that is not actually in the Diagram.
		 */
        CellInfo getCell(double x, double y, double z, double searchRadius) const
        {
            VERIFY(index < particles_.size());
            VERIFY(searchRadius > 0);
            return CellInfo(*this, NO_INDEX, {x, y, z}, searchRadius);
        }

		/**
		 * \brief
		 *   Get the CellInfo for a Particle using a specfic search radius.
		 *
		 * \param[in]  x             The x-coordinate of the particle
		 * \param[in]  y             The y-coordinate of the particle
		 * \param[in]  z             The z-coordinate of the particle
		 * \param[in]  targetGroup   The TargetGroup for cutting the Cell
		 * \param[in]  searchRadius  The search radius to use for finding neighbors of the Particle
		 *
		 * \return
		 *   The CellInfo for the particle.
		 *
		 * \remark
		 *   Used to find a cell for a particle that is not actually in the Diagram.
		 */
        CellInfo getCell(double x, double y, double z, TargetGroup targetGroup, double searchRadius) const
        {
            VERIFY(index < particles_.size());
            VERIFY(searchRadius > 0);
            return CellInfo(*this, NO_INDEX, {x, y, z}, targetGroup, searchRadius);
        }

		/**
		 * \brief
		 *   Create a TargetGroup containing the specified groups.
		 *
		 * \param[in]  groups  The groups to include in the TargetGroup.
		 *
		 * \return
		 *   The desired TargetGroup.
		 */
		template <typename... Groups>
		TargetGroup targetGroups(Groups... groups)
		{
			return TargetGroup {groups...};
		}

		// dani hacky constructor thing
		TargetGroup targetGroups(std::vector<int> groups)
		{
			TargetGroup targets = TargetGroup();
			for (int i = 0; i < groups.size(); ++i) {
				targets.insert(groups[i]);
			}
			return targets;
		}

		/**
		 * \brief
		 *   Create a SourceGroup containing particles from the specified groups. The SourceGroup
		 *   can be treated like a vector of indices of Particles in the Diagram.
		 *
		 * \param[in]  groups  The groups to be included.
		 *
		 * \return
		 *   The desired SourceGroup.
		 *
		 * \remark
		 *   The SourceGroup does not actually contain any group information once created, and
		 *   is simply included for convenience.
		 */
		template <typename... Groups>
		SourceGroup sourceGroups(Groups... groups)
		{
			std::unordered_set<int> sources = std::unordered_set<int> {groups...};
			SourceGroup src;

			int end = particles_.size();
			for (int i = 0; i < end; ++i) {
				if (sources.count(groups_[i])) {
					src.push_back(i);
				}
			}

			return src;
		}

		/**
		 * \brief
		 *   Compute the Voronoi cell for a Particle.
		 *
		 * \param[in]  particleIndex  The index of the Particle
		 * \param[in]  position       The position of the Particle
		 * \param[in]  poly           A reference to the Polyhedron to be cut to compute the cell
		 * \param[in]  search         A reference to an ExpandingSearch for the Particle
		 * \param[in]  targetGroup    If not empty, only cut the particle with the specified groups.
		 * \param[in]  searchRadius   The search radius to look for neighbors of the particle. Optional.
		 *
		 * \remarks
		 *   Uses the efficient expanding search strategy to find neighbors if searchRadius is not specified.
		 */
        void computeVoronoiCell(SizeType particleIndex, Vector3 position, Polyhedron& poly,
                                SDS::ExpandingSearch& search, TargetGroup& targetGroup,
								double searchRadius = NO_RADIUS) const
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
					if (targetGroup.empty()) {
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

					else {
						for (SizeType index : search) {
							Vector3 shiftedPosition = particles_[index].position - position;
							if (index != particleIndex
								&& mag2(shiftedPosition) <= searchRadius
								&& targetGroup.count(groups_[index]) != 0)
							{
								poly.cutWithPlane(
									index,
									Plane::halfwayFromOriginTo(shiftedPosition)
									);
							}
						}
					}
                }
            }

            // Otherwise (given radius), do a basic search
            else {
				if (targetGroup.empty()) {
					for (SizeType index : spatialStructure_.search(position, searchRadius)) {
						if (index != particleIndex) {
							poly.cutWithPlane(
								index,
								Plane::halfwayFromOriginTo(particles_[index].position - position)
								);
						}
					}
				}

				else {
					for (SizeType index : spatialStructure_.search(position, searchRadius)) {
						if (index != particleIndex && targetGroup.count(groups_[index]) != 0) {
							poly.cutWithPlane(
								index,
								Plane::halfwayFromOriginTo(particles_[index].position - position)
								);
						}
					}					
				}
            }
        }
    };


	/// A wrapper function for calling Diagram::ComputeVoronoiCell.
    void Cell::computeVoronoiCell()
    {
        diagram_->computeVoronoiCell(index_, position_, poly_, search_, target_group_, searchRadius_);
    }

	/**
	 * \brief
	 *   A wrapper function for finding the Particle at an index of a Diagram.
	 *
	 * \return
	 *   The desired Particle.
	 */
    const Particle& Cell::getParticle() const
    {
        return diagram_->particles_[index_];
    }

}
