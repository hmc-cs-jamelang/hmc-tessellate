/*
 * cellarray.hpp
 *
 * Implementation of a basic 3-D cell array.
 *
 * The cell array is a very efficient spatial data structure for storing evenly
 * distributed points. However, it is far less effective when points are
 * concentrated in a specific area.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

namespace spatial
{

	template<typename PointType>
	class Celery
	{

	private:

		// The ideal number of particles in a cell. Used to determine the
		// size of a cell.
		static const std::size_t cell_density_ = 10;

		// The domain of the cell array.
		double xmin_, xmax_, ymin_, ymax_, zmin_, zmax_;

		// The array of points.
		std::vector<PointType> points_;

		// The cell correspoding to the point with the same index.
		std::vector<std::size_t> cells_;

		// The delimiters for each cell.
		std::vector<std::size_t> delimiters_;

		// The inverse size of a cell in each dimension.
		double cell_size_inv_x_, cell_size_inv_y_, cell_size_inv_z_;

		// The number of cells in each dimension.
		std::size_t num_cells_dim_;

	private:

		// Compute the attributes of the cell array for a specific number of
		// cells. The cells are assumed to be evenly distributed.
		void computeCellData(std::size_t numPoints);

	public:

		// Allow default constructor for use with the initialize function.
		Celery() = default;

		/*
		 * Constructor
		 *
		 * Creates a 3-D cell array from a set of points. Assumes the points
		 * have x, y, and z coordinate data members.
		 */
		template<typename PointIterator>
		Celery(double xmin, double xmax,
			   double ymin, double ymax,
			   double zmin, double zmax,
			   PointIterator begin, PointIterator end);

	public:

		/*
		 * initialize
		 *
		 * Creates a 3-D cell array from a set of points, automatically fining
		 * the boundaries of the region from the set of points. To be used
		 * with the default constructor.
		 */
		template<typename PointIterator>
		void initialize(PointIterator begin, PointIterator end);

		template<typename PointIterator, typename F>
		void initialize(PointIterator begin, PointIterator end, F getPoint);

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

		/*
		 * insert
		 *
		 * Inserts a set of points into the cell array. Assumes the points have
		 * x, y, and z coordinate data members.
		 */
		template<typename PointIterator, typename F>
		void insert(PointIterator begin, PointIterator end, F getPoint);

		/*
		 * insert
		 *
		 * Inserts a set of points into the cell array. Assumes the points have
		 * x, y, and z coordinate data members.
		 */
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
		void reserve(std::size_t size)
		{
			points_.reserve(size);
			cells_.reserve(size);
		}

		// Get the number of cells in each dimension
		std::size_t getNumCellsDim() const
		{
			return num_cells_dim_;
		}

		// Get the cell value that would be assigned to a point
		template <typename XYZPoint>
		std::size_t getCell(XYZPoint point) const
		{
			std::size_t xIndex = std::size_t( (point.x - xmin_) * cell_size_inv_x_ );
			std::size_t yIndex = std::size_t( (point.y - ymin_) * cell_size_inv_y_ );
			std::size_t zIndex = std::size_t( (point.z - zmin_) * cell_size_inv_z_ );

			return xIndex * pow(num_cells_dim_, 2) + yIndex * num_cells_dim_ + zIndex;
		}

		// Get the cell value that would be assigned to a point with coordinates
		// x, y, z.
		std::size_t getCell(double x, double y, double z) const
		{
			std::size_t xIndex = std::size_t( (x - xmin_) * cell_size_inv_x_ );
			std::size_t yIndex = std::size_t( (y - ymin_) * cell_size_inv_y_ );
			std::size_t zIndex = std::size_t( (z - zmin_) * cell_size_inv_z_ );

			return xIndex * pow(num_cells_dim_, 2) + yIndex * num_cells_dim_ + zIndex;
		}

		// Get the inverse size of a cell in the x dimension
		double getCellSizeInvX() const
		{
			return cell_size_inv_x_;
		}

		// Get the inverse size of a cell in the y dimension
		double getCellSizeInvY() const
		{
			return cell_size_inv_y_;
		}

		// Get the inverse size of a cell in the z dimension
		double getCellSizeInvZ() const
		{
			return cell_size_inv_z_;
		}

		// Return the points_ vector
		const std::vector<PointType>& getPoints() const
		{
			return points_;
		}

		// Return the cells_ vector
		const std::vector<std::size_t>& getCells() const
		{
			return cells_;
		}

		// Return the delimiters_ vector
		const std::vector<std::size_t>& getDelimiters() const
		{
			return delimiters_;
		}

		// Get lower x boundary
		double getXMin() const
		{
			return xmin_;
		}

		// Get upper x boundary
		double getXMax() const
		{
			return xmax_;
		}

		// Get lower y boundary
		double getYMin() const
		{
			return ymin_;
		}

		// Get upper y boundary
		double getYMax() const
		{
			return ymax_;
		}

		// Get lower z boundary
		double getZMin() const
		{
			return zmin_;
		}

		// Get upper z boundary
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

			bool operator()(std::size_t i1, std::size_t i2) const
			{
				return value_vector[i1] < value_vector[i2];
			}
		};

		// Helper function for computing the distance between two sets of 3-D coordinates
		double euclidDist(double x1, double y1, double z1, double x2, double y2, double z2) const
		{
			return sqrt( pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2) );
		}

		// Get the x-index of a cell with given x-coordinate
		std::size_t getXIndex(double x) const
		{
			return std::size_t( (x - xmin_) * cell_size_inv_x_ );
		}

		// Get the y-index of a cell with given y-coordinate
		std::size_t getYIndex(double y) const
		{
			return std::size_t( (y - ymin_) * cell_size_inv_y_ );
		}

		// Get the z-index of a cell with given z-coordinate
		std::size_t getZIndex(double z) const
		{
			return std::size_t( (z - zmin_) * cell_size_inv_z_ );
		}

		// Get the cell from x-, y-, and z-indices in the cell array
		std::size_t getCellFromIndices(std::size_t xIndex, std::size_t yIndex, std::size_t zIndex) const
		{
			return xIndex * pow(num_cells_dim_, 2) + yIndex * num_cells_dim_ + zIndex;
		}

		// Check whether any part of a cell is within a certain distance of a point
		bool checkCellInRange(double x, double y, double z, double dist, std::size_t xIndex, std::size_t yIndex,
							  std::size_t zIndex) const;

	public:

		/*
		 * findNeighborsInCellRadius
		 *
		 * Given 3-D coordinates and a radius, finds all cells within that radius of the coordinates.
		 * Stores pointers to the points in those cells in a vector.
		 */
		void findNeighborsInCellRadius(double x, double y, double z, double radius, std::vector<PointType*>& pts) const;

		void findNeighborsInCellRadius(double x, double y, double z, double radius, std::vector<PointType>& pts) const;


		/*
		 * findNeighborsInRealRadius
		 *
		 * Given 3-D coordinates and a radius, finds all points within that radius of the coordinates.
		 * Requires validation of each point in addition to each cell, meaining that this function is
		 * slower than findNeighborsInCellRadius. However, this function does not store any points
		 * outside of the specified radius.
		 */
		void findNeighborsInRealRadius(double x, double y, double z, double radius, std::vector<PointType*>& pts) const;

		/*
		 * neighborQuery
		 *
		 * Necessary for compilation. Does not do anything yet.
		 */
		void neighborQuery(const std::array<double, 3>& point, std::vector<PointType*>& NeighborParticles);
	};
}

#include "cellarray-private.hpp"
