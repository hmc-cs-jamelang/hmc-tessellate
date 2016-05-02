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
#include <assert.h>

#include "vectormath.hpp"
#include "polyhedron.hpp"
#include "verification.hpp"
#include "sds.hpp"

namespace hmc {
    using SizeType = unsigned int;
    using SDS = CellArray<SizeType, SizeType>;
    // using SDS = TrivialSDS<SizeType>;
    // using SDS = ShellArray<SizeType>;

	using TargetGroup = std::unordered_set<int>;
	using SourceGroup = std::vector<int>;

    class Diagram;
    class Cell;
	class VoronoiFace;

	/// Value used to signify that no group has been specified.
	constexpr int DEFAULT_GROUP = -1;

	/// Value used to signify that no search radius is known for building a Cell.
    constexpr double NO_RADIUS = -1;

	/// Value used to signify that a particle does not exist in a Diagram.
	constexpr SizeType NO_INDEX = std::numeric_limits<SizeType>::max();

	/**
	 * \struct CellConstructor
	 *
	 * \brief
	 *   Contains all the information for constructing a Cell in a Diagram.
	 */
    struct CellConstructor {
		/// A pointer to the Diagram that contains the Cell.
        const Diagram* diagram;

		/// The index of the point in this Cell in the Diagram.
        SizeType index;

		/// The position of the point in this Cell.
        Vector3 position;

		/// The search radius for neighbors of the particle in the Cell. If no radius is known, set to NO_RADIUS.
        double searchRadius;

		/// The target group of the Cell.
		const TargetGroup* targetGroup;

		/// Constructor
        CellConstructor(const Diagram& diagram, SizeType index, Vector3 position, double searchRadius = NO_RADIUS)
            : diagram(&diagram), index(index), position(position), searchRadius(searchRadius), targetGroup(nullptr)
        { /* Done */ }

		/// Constructor
        CellConstructor(const Diagram& diagram, SizeType index, Vector3 position, const TargetGroup& targetGroup, double searchRadius = NO_RADIUS)
            : diagram(&diagram), index(index), position(position), searchRadius(searchRadius), targetGroup(&targetGroup)
        { /* Done */ }

        CellConstructor& withSearchRadius(double searchRadius)
        {
            this->searchRadius = searchRadius;
            return *this;
        }

        CellConstructor& withTargetGroup(const TargetGroup& targetGroup)
        {
            this->targetGroup = &targetGroup;
            return *this;
        }
    };

	/**
	 * \class Cell
	 *
	 * \brief
	 *   A Voronoi cell. Used to retrieve various parameters.
	 *
	 * \remarks
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

		/// The target group of the Cell.
		const TargetGroup* targetGroup_;

		/// The ExpandingSearch object used to search for neighbors of the point in this Cell.
        SDS::ExpandingSearch search_;

    public:
		/// Default constructor
        Cell() : diagram_(nullptr), polyComputed_(false), targetGroup_(nullptr) {}

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
              position_(position), polyComputed_(false), searchRadius_(searchRadius), targetGroup_(nullptr)
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
        Cell(const Diagram& diagram, SizeType index, Vector3 position, const TargetGroup& targetGroup, double searchRadius = NO_RADIUS)
            : diagram_(&diagram), index_(index),
              position_(position), polyComputed_(false), searchRadius_(searchRadius),
			  targetGroup_(&targetGroup)
        { /* Done */ }

		/**
		 * \brief Constructor
		 *
		 * \param[in]  info  A CellConstructor object
		 *
		 * \remark
		 *   Creates the Cell specified by the CellConstructor.
		 */
		Cell(const CellConstructor info)
			: Cell()
		{
			*this = info;
		}

		/// Assignment operator
        Cell& operator=(const CellConstructor info)
        {
            clear();
            diagram_ = info.diagram;
            index_ = info.index;
            position_ = info.position;
			targetGroup_ = info.targetGroup;
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

        /**
         * \brief
         *   Get the position of the particle this
         *   Voronoi cell is constructed around.
         *
         * \return
         *   The position, as a Vector3.
         */
        inline const Vector3& getPosition() const;

        /**
         * \brief
         *   Gets the original index of the particle
         *   this cell is based on. The original index
         *   is the index position the particle was originally
         *   added to the diagram in.
         *   If the cell is not based on a particle in
         *   the Diagram, returns NO_INDEX instead
         *
         * \return
         *   The original index, or NO_INDEX.
         */
        inline SizeType getOriginalIndex() const;

        /**
         * \brief
         *   Gets the internal index of the particle
         *   this cell is based on. The internal index
         *   is actual index position the particle
         *   within the Diagram.
         *   If the cell is not based on a particle in
         *   the Diagram, returns NO_INDEX instead
         *
         * \return
         *   The internal index, or NO_INDEX.
         */
        inline SizeType getInternalIndex() const;

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
		 * \param[out]  result  A container of VoronoiFace objects for storing the neighbors
		 */
		template<typename Collection>
		void computeVertices(Collection& result)
		{
			ensurePolyComputed();
			poly_.computeVertices(result);
		}

		/**
		 * \brief
		 *   Compute the faces of the Voronoi cell.
		 *
		 * \param[out]  result  A container for storing the faces
		 */
		template<typename Collection>
		void computeFaces(Collection& result)
		{
			ensurePolyComputed();

			auto output = std::back_inserter(result);
			for (auto fit = poly_.faces_.begin(); fit != poly_.faces_.end(); ++fit) {
				*output++ = VoronoiFace(fit.index(), *this);
			}
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
        inline void computeVoronoiCell();
    };

	/**
	 * \class VoronoiFace
	 *
	 * \brief
	 *   Allows access to data for a single face of a Voronoi Cell.
	 */
	class VoronoiFace
	{
	private:

		// The FaceIndex corresponding to this face in the Polyhedron for the Voronoi Cell.
		FaceIndex fi_;

		// A pointer to the Voronoi Cell that contains this face.
		Cell* cell_;

	public:

		/**
		 * \brief
		 *   Constructor
		 *
		 * \param[in]  fi    The face index of this face in the cell Polyhedron
		 * \param[in]  cell  The cell that contains this face
		 */
		VoronoiFace(FaceIndex fi, Cell& cell)
			: fi_(fi), cell_(&cell)
		{}

		/**
		 * \brief
		 *   Get the vertices of the VoronoiFace.
		 *
		 * \param[out]  result  A container of doubles, stored sequentially as x, y, z coordinates of each vertex
		 */
		template<typename Collection>
		void computeVertices(Collection& result)
		{
			cell_->poly_.computeFaceVertices(fi_, result);
		}

		/**
		 * \brief
		 *   Get the area of the VoronoiFace.
		 *
		 * \return
		 *   The area of the VoronoiFace.
		 */
		double computeArea()
		{
			return 0.5 * mag(cell_->poly_.weightedNormal(fi_));
		}

		/**
		 * \brief
		 *   Get the neighbor of the Voronoi Cell which cut the polyhedron to make this face.
		 *
		 * \return
		 *  The index of the neighbor in the Diagram that contains the Voronoi Cell.
		 */
		SizeType getNeighbor()
		{
			return cell_->poly_.faces_[fi_].id;
		}
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
        /// Whether or not the diagram has been readied for Voronoi computations
        bool initialized_ = false;

		/// The Polyhedron that represents the container of the points.
        Polyhedron containerShape_;

        /// A bounding box enclosing all the points
        BoundingBox boundingBox_;

		/// The vector of Particle objects in the Diagram.
        std::vector<Vector3> particles_;

        /// A vector mapping internal indices to original indices
        std::vector<SizeType> originalIndices_;

        /// A vector mapping original indices to internal indices
        std::vector<SizeType> internalIndices_;

		/// The vector of groups that contain each particle. The group at an index corresponds
		/// to the particle at that index.
		std::vector<int> groups_;

		/// The spatial data structure used to determine spatial relationships between points.
        SDS spatialStructure_;

    public:
		/// Default constructor
        Diagram() = default;

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

        void clear()
        {
            initialized_ = false;
            boundingBox_.reset();

            particles_.clear();
            groups_.clear();
            originalIndices_.clear();
            internalIndices_.clear();
            spatialStructure_.clear();
        }

		/**
		 * \brief
		 *   Add a Particle to the Diagram.
		 *
		 * \param[in]  x      The x-coordinate of the Particle
		 * \param[in]  y      The y-coordinate of the Particle
		 * \param[in]  z      The z-coordinate of the Particle
		 * \param[in]  group  The group containing the Particle. Optional.
		 *
		 * \remarks
		 *   Should only be used before the initialize function is called.
		 */
        void addParticle(double x, double y, double z, int group = DEFAULT_GROUP)
        {
            addParticle(Vector3 {x, y, z}, group);
        }

		/**
		 * \brief
		 *   Add a Particle to the Diagram.
		 *
		 * \param[in]  position  A Vector3 containing the {x,y,z} coordinates of the particle.
		 * \param[in]  group     The group containing the Particle. Optional.
		 *
		 * \remarks
		 *   Should only be used before the initialize function is called.
		 */
        void addParticle(Vector3 position, int group = DEFAULT_GROUP)
        {
            particles_.push_back(position);
            boundingBox_.adjustToContain(position);
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
            // We are not certain whether or not padding
            // is a good or a bad idea.
            // It does not appear to be necessary, so,
            // we are saying that one should supply
            // their own bounding box if padding is desired.
            // However, here's the relevant call:
            //boundingBox_.pad(1e-6);

            assert(!initialized_ && "Diagram does not support multiple initializations");

            sortParticlesByMortonIndex();

            if (containerShape_.isClear()) {
                containerShape_.buildCube(boundingBox_);
            }

            spatialStructure_.initialize(0, particles_.size(),
                [this](SizeType index) -> Vector3 {
                    return particles_[index];
                }
            );

            initialized_ = true;
        }

		/**
		 * \brief
		 *   Sort particles in Morton order.
		 *
		 * \details
		 *   Morton order (or Z-order) is a method of sorting multidimensional data according
		 *   to a single dimension which preserves locality of data points.
		 */
        void sortParticlesByMortonIndex()
        {
            const BoundingBox& bb = boundingBox_;

            constexpr unsigned int numLevels = 10;
            constexpr double numCellsPerSide = 1 << numLevels;

            const Vector3 sideLengths = (bb.high - bb.low) / numCellsPerSide;

            using Index3 = Vector3_of<SizeType>;
            auto cellIndices = [&](SizeType i) -> Index3 {
                auto offset = particles_[i] - bb.low;
                return Index3(
                    offset.x / sideLengths.x,
                    offset.y / sideLengths.y,
                    offset.z / sideLengths.z
                );
            };

            auto mortonIndex = [&](SizeType i) -> SizeType {
                Index3 indices = cellIndices(i);

                SizeType morton = 0;
                unsigned int bit = 0;

                for (unsigned int level = 0; level < numLevels; ++level) {
                    morton |= ((indices.x >> level) & 1u) << bit++;
                    morton |= ((indices.y >> level) & 1u) << bit++;
                    morton |= ((indices.z >> level) & 1u) << bit++;
                }
                return morton;
            };


            // -----------------------------------------------------

            const SizeType n = particles_.size();

            std::vector<SizeType>& permutation = internalIndices_;
            std::vector<SizeType>& morton = originalIndices_;

            permutation.resize(n);
            morton.resize(n);

            for (SizeType i = 0; i < n; ++i) {
                permutation[i] = i;
                morton[i] = mortonIndex(i);
            }

            std::sort(permutation.begin(), permutation.end(),
                [&morton](SizeType i, SizeType j) -> bool {
                    return morton[i] < morton[j];
                });

            constexpr SizeType COMPLETED = std::numeric_limits<SizeType>::max();
            for (SizeType cycleStart = 0; cycleStart < n; ++cycleStart) {
                if (permutation[cycleStart] == COMPLETED) {
                    continue;
                }

                auto startParticle = particles_[cycleStart];
                auto startGroup = groups_[cycleStart];

                SizeType current = cycleStart;
                SizeType next = originalIndices_[current] = permutation[current];

                while (next != cycleStart) {
                    originalIndices_[current] = permutation[current];
                    permutation[current] = COMPLETED;

                    particles_[current] = particles_[next];
                    groups_[current] = groups_[next];

                    current = next;
                    next = permutation[next];
                }


                originalIndices_[current] = permutation[current];
                permutation[current] = COMPLETED;

                particles_[current] = startParticle;
                groups_[current] = startGroup;
            }

            for (SizeType i = 0; i < n; ++i) {
                internalIndices_[originalIndices_[i]] = i;
            }
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
		 *   Get the CellConstructor for a Particle.
		 *
		 * \param[in]  index  The index of the Particle
		 *
		 * \return
		 *   The CellConstructor for the particle.
		 */
        CellConstructor getCell(SizeType index) const
        {
            assert(initialized_ &&
                "Remember to call diagram.initialize() before doing Voronoi computations");
            VERIFY(index < particles_.size());
            auto x = CellConstructor(*this, index, particles_[index]);
            return x;
        }

		/**
		 * \brief
		 *   Get the CellConstructor for a Particle.
		 *
		 * \param[in]  index        The index of the Particle
		 * \param[in]  targetGroup  The TargetGroup for cutting the Cell
		 *
		 * \return
		 *   The CellConstructor for the particle.
		 */
        CellConstructor getCell(SizeType index, const TargetGroup& targetGroup) const
        {
            assert(initialized_ &&
                "Remember to call diagram.initialize() before doing Voronoi computations");
            VERIFY(index < particles_.size());
            return CellConstructor(*this, index, particles_[index], targetGroup);
        }

		/**
		 * \brief
		 *   Get the CellConstructor for a Particle using a specfic search radius.
		 *
		 * \param[in]  index         The index of the Particle
		 * \param[in]  searchRadius  The search radius to use for finding neighbors of the Particle
		 *
		 * \return
		 *   The CellConstructor for the particle.
		 */
        CellConstructor getCell(SizeType index, double searchRadius) const
        {
            assert(initialized_ &&
                "Remember to call diagram.initialize() before doing Voronoi computations");
            VERIFY(index < particles_.size());
            VERIFY(searchRadius > 0);
            return CellConstructor(*this, index, particles_[index], searchRadius);
        }

		/**
		 * \brief
		 *   Get the CellConstructor for a Particle using a specfic search radius.
		 *
		 * \param[in]  index         The index of the Particle
		 * \param[in]  targetGroup   The TargetGroup for cutting the Cell
		 * \param[in]  searchRadius  The search radius to use for finding neighbors of the Particle
		 *
		 * \return
		 *   The CellConstructor for the particle.
		 */
        CellConstructor getCell(SizeType index, const TargetGroup& targetGroup, double searchRadius) const
        {
            assert(initialized_ &&
                "Remember to call diagram.initialize() before doing Voronoi computations");
            VERIFY(index < particles_.size());
            VERIFY(searchRadius > 0);
            return CellConstructor(*this, index, particles_[index], targetGroup, searchRadius);
        }

		/**
		 * \brief
		 *   Get the CellConstructor for a hypothetical particle.
		 *
		 * \param[in]  x  The x-coordinate of the particle
		 * \param[in]  y  The y-coordinate of the particle
		 * \param[in]  z  The z-coordinate of the particle
		 *
		 * \return
		 *   The CellConstructor for the particle.
		 *
		 * \remark
		 *   Used to find a cell for a particle that is not actually in the Diagram.
		 */
        CellConstructor getCell(double x, double y, double z) const
        {
            assert(initialized_ &&
                "Remember to call diagram.initialize() before doing Voronoi computations");
            return CellConstructor(*this, NO_INDEX, {x, y, z});
        }

		/**
		 * \brief
		 *   Get the CellConstructor for a Particle.
		 *
		 * \param[in]  x            The x-coordinate of the particle
		 * \param[in]  y            The y-coordinate of the particle
		 * \param[in]  z            The z-coordinate of the particle
		 * \param[in]  targetGroup  The TargetGroup for cutting the Cell
		 *
		 * \return
		 *   The CellConstructor for the particle.
		 *
		 * \remark
		 *   Used to find a cell for a particle that is not actually in the Diagram.
		 */
        CellConstructor getCell(double x, double y, double z, const TargetGroup& targetGroup) const
        {
            assert(initialized_ &&
                "Remember to call diagram.initialize() before doing Voronoi computations");
            return CellConstructor(*this, NO_INDEX, {x, y, z}, targetGroup);
        }

		/**
		 * \brief
		 *   Get the CellConstructor for a Particle using a specfic search radius.
		 *
		 * \param[in]  x             The x-coordinate of the particle
		 * \param[in]  y             The y-coordinate of the particle
		 * \param[in]  z             The z-coordinate of the particle
		 * \param[in]  searchRadius  The search radius to use for finding neighbors of the Particle
		 *
		 * \return
		 *   The CellConstructor for the particle.
		 *
		 * \remark
		 *   Used to find a cell for a particle that is not actually in the Diagram.
		 */
        CellConstructor getCell(double x, double y, double z, double searchRadius) const
        {
            assert(initialized_ &&
                "Remember to call diagram.initialize() before doing Voronoi computations");
            VERIFY(searchRadius > 0);
            return CellConstructor(*this, NO_INDEX, {x, y, z}, searchRadius);
        }

		/**
		 * \brief
		 *   Get the CellConstructor for a Particle using a specfic search radius.
		 *
		 * \param[in]  x             The x-coordinate of the particle
		 * \param[in]  y             The y-coordinate of the particle
		 * \param[in]  z             The z-coordinate of the particle
		 * \param[in]  targetGroup   The TargetGroup for cutting the Cell
		 * \param[in]  searchRadius  The search radius to use for finding neighbors of the Particle
		 *
		 * \return
		 *   The CellConstructor for the particle.
		 *
		 * \remark
		 *   Used to find a cell for a particle that is not actually in the Diagram.
		 */
        CellConstructor getCell(double x, double y, double z, const TargetGroup& targetGroup, double searchRadius) const
        {
            assert(initialized_ &&
                "Remember to call diagram.initialize() before doing Voronoi computations");
            VERIFY(searchRadius > 0);
            return CellConstructor(*this, NO_INDEX, {x, y, z}, targetGroup, searchRadius);
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
			for (unsigned int i = 0; i < groups.size(); ++i) {
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
                                SDS::ExpandingSearch& search, const TargetGroup* targetGroup= nullptr,
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
					if (targetGroup == nullptr) {
						for (SizeType index : search) {
							Vector3 shiftedPosition = particles_[index] - position;

							if (index != particleIndex
								&& mag2(shiftedPosition) <= searchRadius)
							{
                    			poly.cutWithPlane(
									originalIndices_[index],
									Plane::halfwayFromOriginTo(shiftedPosition)
                                    );
							}
						}
					}

					else {
						for (SizeType index : search) {
							Vector3 shiftedPosition = particles_[index] - position;
							if (index != particleIndex
								&& mag2(shiftedPosition) <= searchRadius
								&& targetGroup->count(groups_[index]) != 0)
							{
								poly.cutWithPlane(
									originalIndices_[index],
									Plane::halfwayFromOriginTo(shiftedPosition)
									);
							}
						}
					}
                }
            }

            // Otherwise (given radius), do a basic search
            else {
				if (targetGroup == nullptr) {
					for (SizeType index : spatialStructure_.search(position, searchRadius)) {
						if (index != particleIndex) {
							poly.cutWithPlane(
								originalIndices_[index],
								Plane::halfwayFromOriginTo(particles_[index] - position)
								);
						}
					}
				}

				else {
					for (SizeType index : spatialStructure_.search(position, searchRadius)) {
						if (index != particleIndex && targetGroup->count(groups_[index]) != 0) {
							poly.cutWithPlane(
								originalIndices_[index],
								Plane::halfwayFromOriginTo(particles_[index] - position)
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
        diagram_->computeVoronoiCell(index_, position_, poly_, search_, targetGroup_, searchRadius_);
    }

	/**
	 * \brief
	 *   A wrapper function for finding the Particle at an index of a Diagram.
	 *
	 * \return
	 *   The desired Particle.
	 */
    const Vector3& Cell::getPosition() const
    {
        return diagram_->particles_[index_];
    }

    SizeType Cell::getOriginalIndex() const
    {
        return (index_ == NO_INDEX)
                ? NO_INDEX
                : diagram_->originalIndices_[index_];
    }

    SizeType Cell::getInternalIndex() const
    {
        // index_ might be NO_INDEX, but in that
        // event we do in fact want to return NO_INDEX
        return index_;
    }
}
