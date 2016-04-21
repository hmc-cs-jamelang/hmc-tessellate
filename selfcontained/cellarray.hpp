/**
 * \file cellarray.hpp
 *
 * \author 2015-2016 Sandia Clinic Team
 *
 * \brief Implementation of a basic 3-D cell array.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "vectormath.hpp"

namespace hmc {
namespace spatial {
	/**
	 * \class Celery
	 *
	 * \brief
	 *   A custom cell array for storing 3-D points.
	 *
	 * \details
	 *   This cell array is used for storing points in 3-D space to be Voronoi
	 *   tessellated.
	 *
	 * \remarks
	 *   The cell array is a very efficient spatial data structure for storing evenly
	 *   distributed points. However, it is far less effective when points are
	 *   concentrated in a specific area.
	 */
	template<typename PointType>
	class Celery
	{

	private:

		struct DistanceIndex;

		// The ideal number of particles in a cell. Used to determine the
		// size of a cell.
		static constexpr double cell_density_ = 5.0/4;

		// The domain of the cell array.
		double xmin_, xmax_, ymin_, ymax_, zmin_, zmax_;

		// The array of points.
		std::vector<PointType> points_;

		// The cell correspoding to the point with the same index.
		std::vector<std::size_t> cells_;

		// The delimiters for each cell.
		std::vector<std::size_t> delimiters_;

		// The size of a cell in each dimension.
		double cell_size_x_, cell_size_y_, cell_size_z_;

		// The inverse size of a cell in each dimension.
		double cell_size_inv_x_, cell_size_inv_y_, cell_size_inv_z_;

		// The number of cells in each dimension.
		std::size_t num_cells_dim_;

		// The order in which to search through cells
		std::vector<DistanceIndex> search_order_;

	private:

		// Compute the attributes of the cell array for a specific number of
		// cells. The cells are assumed to be evenly distributed.
		void computeCellData(std::size_t numPoints);

	public:

		// Allow default constructor for use with the initialize function.
		/// Default constructor
		Celery() = default;

		/**
		 * \brief
		 *   Creates a 3-D cell array from a set of points. Assumes the points
		 *   have x, y, and z coordinate data members.
		 *
		 * \param  xmin   The lower bounding plane of the region on the x-axis
		 * \param  xmax   The upper bounding plane of the region on the x-axis
		 * \param  ymin   The lower bounding plane of the region on the y-axis
		 * \param  ymax   The upper bounding plane of the region on the y-axis
		 * \param  zmin   The lower bounding plane of the region on the z-axis
		 * \param  zmax   The upper bounding plane of the region on the z-axis
		 * \param  begin  An iterator to the first point to add
		 * \param  end    An iterator to the end of the list of points
		 */
		template<typename PointIterator>
		Celery(double xmin, double xmax,
			   double ymin, double ymax,
			   double zmin, double zmax,
			   PointIterator begin, PointIterator end);

	public:

		/**
		 * \brief
		 *   Creates a 3-D cell array from a set of points, automatically fining
		 *   the boundaries of the region from the set of points. To be used
		 *   with the default constructor.
		 *
		 * \param  begin  An iterator to the first point to add
		 * \param  end    An iterator to the end of the list of points
		 *
		 * \remarks
		 *   Assumes that points have public x, y, and z fields.
		 */
		template<typename PointIterator>
		void initialize(PointIterator begin, PointIterator end);

		/**
		 * \brief
		 *   Creates a 3-D cell array from a set of points, automatically fining
		 *   the boundaries of the region from the set of points. To be used
		 *   with the default constructor.
		 *
		 * \param  begin  An iterator to the first point to add
		 * \param  end    An iterator to the end of the list of points
		 * \param  F      A lambda function for getting a point from the dereferenced iterator
		 *
		 * \remarks
		 *   Assumes that points have public x, y, and z fields.
		 */
		template<typename PointIterator, typename F>
		void initialize(PointIterator begin, PointIterator end, F getPoint);

		/**
		 * \brief
		 *   Creates a 3-D cell array from a set of points, automatically fining
		 *   the boundaries of the region from the set of points. To be used
		 *   with the default constructor.
		 *
		 * \param  begin  An index to the first point to add
		 * \param  end    An index to the end of the list of points
		 * \param  F      A lambda function for getting a point from the index
		 *
		 * \remarks
		 *   Assumes that points have public x, y, and z fields.
		 */
		template<typename F>
		void initialize(PointType begin, PointType end, F getPoint);

	private:

		/*
		 * computeBoundsFromPoints
		 *
		 * Automatically finds and sets the boundaries for a set of points.
		 */
		template<typename PointIterator, typename F>
		void computeBoundsFromPoints(PointIterator begin, PointIterator end, F getPoint);

		template<typename F>
		void computeBoundsFromPoints(PointType begin, PointType end, F getPoint);

		/*
		 * insert
		 *
		 * Inserts a single point into the cell array. Assumes the points have
		 * x, y, and z coordinate data members.
		 */
		template<typename F>
		void insert(PointType point, F getPoint);

		template<typename PointIterator, typename F>
		void insert(PointIterator begin, PointIterator end, F getPoint);

		template<typename F>
		void insert(PointType begin, PointType end, F getPoint);

		/*
		 * initializeCellArray
		 *
		 * Sort the array of points by cell, so that points are grouped by cell and
		 * arranged in order. Also populates the array of delimiters.
		 */
		void initializeCellArray();

	public:

		/// Get the number of cells in each dimension
		std::size_t getNumCellsDim() const
		{
			return num_cells_dim_;
		}

		/// Get the cell value that would be assigned to a point
		template <typename XYZPoint>
		std::size_t getCell(XYZPoint point) const
		{
			return getCell(point.x, point.y, point.z);
		}

		/// Get the cell value that would be assigned to a point with coordinates x, y, z.
		std::size_t getCell(double x, double y, double z) const
		{
			return getCellFromIndices(getXIndex(x), getYIndex(y), getZIndex(z));
		}

		/// Get the inverse size of a cell in the x dimension
		double getCellSizeInvX() const
		{
			return cell_size_inv_x_;
		}

		/// Get the inverse size of a cell in the y dimension
		double getCellSizeInvY() const
		{
			return cell_size_inv_y_;
		}

		/// Get the inverse size of a cell in the z dimension
		double getCellSizeInvZ() const
		{
			return cell_size_inv_z_;
		}

		/// Return the points_ vector
		const std::vector<PointType>& getPoints() const
		{
			return points_;
		}

		/// Return the cells_ vector
		const std::vector<std::size_t>& getCells() const
		{
			return cells_;
		}

		/// Return the delimiters_ vector
		const std::vector<std::size_t>& getDelimiters() const
		{
			return delimiters_;
		}

		/// Get lower x boundary
		double getXMin() const
		{
			return xmin_;
		}

		/// Get upper x boundary
		double getXMax() const
		{
			return xmax_;
		}

		/// Get lower y boundary
		double getYMin() const
		{
			return ymin_;
		}

		/// Get upper y boundary
		double getYMax() const
		{
			return ymax_;
		}

		/// Get lower z boundary
		double getZMin() const
		{
			return zmin_;
		}

		/// Get upper z boundary
		double getZMax() const
		{
			return zmax_;
		}

	private:

		// Helper for sorting points using cell information.
		struct MyComparator
		{
			const std::vector<std::size_t> & value_vector;

			MyComparator(const std::vector<std::size_t> & val_vect)
				: value_vector(val_vect)
			{};

			bool operator()(std::size_t i1, std::size_t i2)
			{
				return value_vector[i1] < value_vector[i2];
			}
		};

		// Helper function for computing the distance between two sets of 3-D coordinates
		double euclidDist(double x1, double y1, double z1, double x2, double y2, double z2) const
		{
			auto sq = [&](double n) -> double {
				return n * n;
			};

			return sqrt( sq(x1 - x2) + sq(y1 - y2) + sq(z1 - z2) );
		}

		// Get the x-index of a cell with given x-coordinate
		std::size_t getXIndex(double x) const
		{
			if (x >= xmax_) {return num_cells_dim_ - 1;}
			return std::min(std::size_t( (x - xmin_) * cell_size_inv_x_ ), num_cells_dim_ - 1);
		}

		// Get the y-index of a cell with given y-coordinate
		std::size_t getYIndex(double y) const
		{
			if (y >= ymax_) {return num_cells_dim_ - 1;}
			return std::min(std::size_t( (y - ymin_) * cell_size_inv_y_ ), num_cells_dim_ - 1);
		}

		// Get the z-index of a cell with given z-coordinate
		std::size_t getZIndex(double z) const
		{
			return std::min(std::size_t( (z - zmin_) * cell_size_inv_z_ ), num_cells_dim_ - 1);
		}

		// Get the cell from x-, y-, and z-indices in the cell array
		std::size_t getCellFromIndices(std::size_t xIndex, std::size_t yIndex, std::size_t zIndex) const
		{
			return xIndex * num_cells_dim_ * num_cells_dim_ + yIndex * num_cells_dim_ + zIndex;
		}

		// Check whether any part of a cell is within a certain distance of a point
		bool checkCellInRange(double x, double y, double z, double dist, std::size_t xIndex, std::size_t yIndex,
							  std::size_t zIndex) const;

	public:

		/**
		 * \brief
		 *   Given 3-D coordinates and a radius, finds all cells within that radius of the coordinates.
		 *   Stores pointers to the points in those cells in a vector.
		 *
		 * \param  x       The x-coordinate of the cell
		 * \param  y       The y-coordinate of the cell
		 * \param  z       The z-coordinate of the cell
		 * \param  radius  The maximum radius for searching
		 * \param  pts     A reference to the vector in which to store pointers to found points
		 */
		void findNeighborsInCellRadius(double x, double y, double z, double radius, std::vector<PointType*>& pts) const;

		/**
		 * \brief
		 *   Given 3-D coordinates and a radius, finds all cells within that radius of the coordinates.
		 *   Stores the points in those cells in a vector.
		 *
		 * \param  x       The x-coordinate of the cell
		 * \param  y       The y-coordinate of the cell
		 * \param  z       The z-coordinate of the cell
		 * \param  radius  The maximum radius for searching
		 * \param  pts     A reference to the vector in which to store found points
		 */
		void findNeighborsInCellRadius(double x, double y, double z, double radius, std::vector<PointType>& pts) const;

		/**
		 * \brief
		 *   Given 3-D coordinates and a radius, finds all points within that radius of the coordinates.
		 *   Requires validation of each point in addition to each cell, meaining that this function is
		 *   slower than findNeighborsInCellRadius. However, this function does not store any points
		 *   outside of the specified radius.
		 *
		 * \param  x       The x-coordinate of the cell
		 * \param  y       The y-coordinate of the cell
		 * \param  z       The z-coordinate of the cell
		 * \param  radius  The maximum radius for searching
		 * \param  pts     A reference to the vector in which to store pointers to found points
		 */
		void findNeighborsInRealRadius(double x, double y, double z, double radius, std::vector<PointType*>& pts) const;

		/**
		 * \brief
		 *   Find a shell of cells around a point, and store all of the points in those cells in a vector.
		 *
		 * \remarks
		 *   This function is meant for specific testing purposes and should probably not be used.
		 *
		 * \param  x          The x-coordinate of the cell
		 * \param  y          The y-coordinate of the cell
		 * \param  z          The z-coordinate of the cell
		 * \param  shell      The shell layer
		 * \param  maxRadius  The maximum radius for searching
		 * \param  pts         A reference to the vector in which to store found points
		 */
		bool findNeighborsInShell(double x, double y, double z, int shell, double maxRadius, std::vector<PointType>& pts) const;


		/**
		 * \brief
		 *   Necessary for compilation with some of our old files. Should not be used.
		 */
		void neighborQuery(const std::array<double, 3>& point, std::vector<PointType*>& NeighborParticles);

	private:

		// A struct that contains a combination of distance information and indices to a cell. Used for
		// sorting cells by distance when searching for neighbors.
		struct DistanceIndex {
			double dist;
			int i, j, k;

			DistanceIndex(double distance, int ind1, int ind2, int ind3)
				: dist(distance), i(ind1), j(ind2), k(ind3)
			{}

			friend bool operator<(const Celery<PointType>::DistanceIndex& di1,
								  const Celery<PointType>::DistanceIndex& di2)
			{
				return (di1.dist < di2.dist);
			}

			friend std::ostream& operator<<(std::ostream& out, const DistanceIndex& di)
			{
				return out << "{dist = " << di.dist << ", (" << di.i << ", " << di.j << ", " << di.k << ")}";
			}
		};

		/*
		 * createSearchArray
		 *
		 * Create a sorted array of cells to search in order of distance.
		 */
		void createSearchArray();

	public:

		/**
		 * class ExpandingSearch
		 *
		 * \brief
		 *   Efficiently searches outward from a point to find neighbors.
		 *
		 * \details
		 *   A class with the capability to automatically search outward from a point within a search
		 *   radius. Cells will never be searched more than once, since this class will only search
		 *   outward from the last searched cell.
		 */
		class ExpandingSearch
		{
		private:
			const Celery* stalk_;

			// Store the index of the last searched cell
			std::size_t last_search_index_;

			// Store the x, y, and z cell indices of the cell to search from
			int x_index_, y_index_, z_index_;

			bool done_;

		public:
			/// Default constructor
			ExpandingSearch() = default;

			/*
			 * Constructor
			 *
			 * Create an ExpandingSearch from a point at the given 3-D coordinates.
			 */
			// ExpandingSearch(const Celery& yourCrush, double x, double y, double z)
			// 	: stalk_(yourCrush)
			// {
			// 	x_index_ = stalk_->getXIndex(x);
			// 	y_index_ = stalk_->getYIndex(y);
			// 	z_index_ = stalk_->getZIndex(z);
			// }

			/**
			 * \brief
			 *   Creates an ExpandingSearch for a point in a Celery.
			 *
			 * \param  yourCrush  A reference to the Celery that contains the point
			 * \param  x          The x-coordinate of the point
			 * \param  y          The y-coordinate of the point
			 * \param  z          The z-coordinate of the point
			 *
			 * \remarks
			 *   The point does not need to exist in the Celery - hypothetical points
			 *   can be used.
			 */
			void initialize(const Celery& yourCrush, double x, double y, double z)
			{
				clear();
				stalk_ = &yourCrush;
				x_index_ = stalk_->getXIndex(x);
				y_index_ = stalk_->getYIndex(y);
				z_index_ = stalk_->getZIndex(z);
			}

			/// Reset the ExpandingSearch to start over.
			void clear()
			{
				last_search_index_ = 0;
				done_ = false;
			}

			/**
			 * \brief
			 *   Check whether the search for neighbors is complete.
			 *
			 * \returns
			 *   <CODE>true</CODE> if the search is complete, <CODE>false</CODE> otherwise.
			 */
			bool done() const
			{
				return done_;
			}

			/**
			 * \brief
			 *   Search outward, adding all points from previously unsearched cells within maxRadius of
			 *   the searching cell. The points are stored in searchPoints.
			 *
			 * \param  maxRadius     The maximum radius in which to search
			 * \param  searchPoints  A reference to the vector in which to store found points
			 *
			 * \remarks
			 *   If maxRadius is reached in an expansion step, the search ends and the expanding search
			 *   is considered complete.
			 */
			void expand(double maxRadius, std::vector<PointType>& searchPoints)
			{
				auto& searchOrder = stalk_->search_order_;
				auto searchOrderSize = searchOrder.size();

				auto invalid = [this](unsigned i) -> bool
				{
					return i >= stalk_->num_cells_dim_;
				};

				constexpr int searchIncrease = 1;
				if (searchIncrease == 1) {
					// std::cerr << "LSI: " << last_search_index_ << ", SOS: " << searchOrderSize << std::endl;
					int xToSearch;
					int yToSearch;
					int zToSearch;
					do {
						if (last_search_index_ >= searchOrderSize) {
							done_ = true;
							return;
						}
						if (searchOrder[last_search_index_].dist > maxRadius) {
							done_ = true;
							return;
						}

						xToSearch = x_index_ + stalk_->search_order_[last_search_index_].i;
						yToSearch = y_index_ + stalk_->search_order_[last_search_index_].j;
						zToSearch = z_index_ + stalk_->search_order_[last_search_index_].k;

						++last_search_index_;
					} while (invalid(xToSearch) || invalid(yToSearch) || invalid(zToSearch));

					std::size_t cellIndex = stalk_->getCellFromIndices(xToSearch, yToSearch, zToSearch);
					std::size_t cellBegin = stalk_->delimiters_[cellIndex];
					std::size_t cellEnd = stalk_->delimiters_[cellIndex + 1];

					for (std::size_t i = cellBegin; i < cellEnd; ++i) {
						searchPoints.push_back(stalk_->points_[i]);
					}
				}
				else {
					if (last_search_index_ >= searchOrderSize) {
						done_ = true;
						return;
					}

					if (searchOrder[last_search_index_].dist > maxRadius) {
						done_ = true;
						return;
					}

					const auto maxIndex = std::min(last_search_index_ + searchIncrease, stalk_->search_order_.size());

					auto index = last_search_index_;
					for (; index < maxIndex; ++index) {
						if (searchOrder[last_search_index_].dist > maxRadius) {
							return;
						}

						int xToSearch = x_index_ + stalk_->search_order_[index].i;
						int yToSearch = y_index_ + stalk_->search_order_[index].j;
						int zToSearch = z_index_ + stalk_->search_order_[index].k;

						if (invalid(xToSearch) || invalid(yToSearch) || invalid(zToSearch)) {
							continue;
						}

						std::size_t cellIndex = stalk_->getCellFromIndices(xToSearch, yToSearch, zToSearch);
						std::size_t cellBegin = stalk_->delimiters_[cellIndex];
						std::size_t cellEnd = stalk_->delimiters_[cellIndex + 1];

						for (std::size_t i = cellBegin; i < cellEnd; ++i) {
							searchPoints.push_back(stalk_->points_[i]);
						}
					}
					last_search_index_ = index;
				}
			}

			// void expand(double maxRadius, std::vector<PointType>& searchPoints)
			// {
			// 	std::size_t searchOrderSize = stalk_->search_order_.size();
			// 	std::size_t searchIndex = last_search_index_;
			// 	double finalDistance = stalk_->search_order_[searchIndex].dist;

			// 	if (searchIndex >= searchOrderSize || finalDistance > maxRadius) {
			// 		done_ = true;
			// 		return;
			// 	}

			// 	for (; searchIndex < searchOrderSize; ++searchIndex) {
			// 		last_search_index_ = searchIndex;

			// 		if (stalk_->search_order_[searchIndex].dist > finalDistance) {
			// 			return;
			// 		}

			// 		int xToSearch = x_index_ + stalk_->search_order_[searchIndex].i;
			// 		int yToSearch = y_index_ + stalk_->search_order_[searchIndex].j;
			// 		int zToSearch = z_index_ + stalk_->search_order_[searchIndex].k;

			// 		auto invalid = [this](int i) -> bool
			// 		{
			// 			return i < 0 || (unsigned) i >= stalk_->num_cells_dim_;
			// 		};

			// 		if (invalid(xToSearch) || invalid(yToSearch) || invalid(zToSearch)) {
			// 			continue;
			// 		}

			// 		std::size_t cellIndex = stalk_->getCellFromIndices(xToSearch, yToSearch, zToSearch);
			// 		std::size_t cellBegin = stalk_->delimiters_[cellIndex];
			// 		std::size_t cellEnd = stalk_->delimiters_[cellIndex + 1];

			// 		for (std::size_t i = cellBegin; i < cellEnd; ++i) {
			// 			searchPoints.push_back(stalk_->points_[i]);
			// 		}
			// 	}

			// 	// We broke out of the loop because we exceeded the cached searchOrder.
			// 	// We could be good little gnomes and do a more expensive search here,
			// 	// or we could, you know, not.
			// 	// static int times = 0;
			// 	// std::cerr << "Ohhhhh.... " << times++ << std::endl;
			// 	// std::cerr << searchIndex << " / " << searchOrderSize << "; " << finalDistance << " / " << stalk_->search_order_[searchIndex].dist << std::endl;
			// 	last_search_index_ = searchIndex;
			// }
		};
	};

	template<typename PointType>
	constexpr double Celery<PointType>::cell_density_;
}}

#include "cellarray-private.hpp"
