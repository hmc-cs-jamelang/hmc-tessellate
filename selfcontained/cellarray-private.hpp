/**
 * \file cellarray-private.hpp
 *
 * \brief Implementation file for cellarray.hpp
 */

#include <limits>

namespace hmc { namespace spatial {

	template<typename PointType>
	template<typename PointIterator>
	Celery<PointType>::Celery(double xmin, double xmax,
							  double ymin, double ymax,
							  double zmin, double zmax,
							  PointIterator begin, PointIterator end) :
		xmin_(xmin), xmax_(xmax), ymin_(ymin), ymax_(ymax), zmin_(zmin), zmax_(zmax)
	{
		computeCellData(std::distance(begin, end));
		insert(begin, end, [](PointType& p){return p;});
		initializeCellArray();
		createSearchArray();
	}

	template<typename PointType>
	template<typename PointIterator>
	void Celery<PointType>::initialize(PointIterator begin, PointIterator end)
	{
		initialize(begin, end, [](PointType& p){return p;});
	}

	template<typename PointType>
	template<typename PointIterator, typename F>
	void Celery<PointType>::initialize(PointIterator begin, PointIterator end, F getPoint)
	{
		// Clear any existing point information
		points_.clear();
		cells_.clear();
		delimiters_.clear();
		search_order_.clear();

		auto size = std::distance(begin, end);
		points_.reserve(size);
		cells_.reserve(size);

		computeBoundsFromPoints(begin, end, getPoint);
		computeCellData(size);

		delimiters_.reserve(num_cells_dim_ * num_cells_dim_ * num_cells_dim_ + 1);

		insert(begin, end, getPoint);
		initializeCellArray();
		createSearchArray();
	}

	template<typename PointType>
	template<typename F>
	void Celery<PointType>::initialize(PointType begin, PointType end, F getPoint)
	{
		// Clear any existing point information
		points_.clear();
		cells_.clear();
		delimiters_.clear();
		search_order_.clear();

		auto size = end - begin;
		points_.reserve(size);
		cells_.reserve(size);

		computeBoundsFromPoints(begin, end, getPoint);
		computeCellData(size);

		delimiters_.reserve(num_cells_dim_ * num_cells_dim_ * num_cells_dim_ + 1);

		insert(begin, end, getPoint);
		initializeCellArray();
		createSearchArray();
	}

	template<typename PointType>
	template<typename PointIterator, typename F>
	void Celery<PointType>::computeBoundsFromPoints(PointIterator begin, PointIterator end, F getPoint)
	{
		xmin_ = ymin_ = zmin_ = std::numeric_limits<double>::max();
		xmax_ = ymax_ = zmax_ = std::numeric_limits<double>::min();

		auto adjust = [](double coord, double& min, double& max) {
			if (coord < min) { min = coord; }
			if (coord > max) { max = coord; }
		};

		for (auto i = begin; i < end; ++i) {
			auto&& p = getPoint(*i);
			adjust(p.x, xmin_, xmax_);
			adjust(p.y, ymin_, ymax_);
			adjust(p.z, zmin_, zmax_);
		}

		auto separate = [](double& /*min*/, double& /*max*/) {
			// if (min == max) {
			// 	max += std::numeric_limits<double>::epsilon();
			// }
		};

		separate(xmin_, xmax_);
		separate(ymin_, ymax_);
		separate(zmin_, zmax_);
	}

	template<typename PointType>
	template<typename F>
	void Celery<PointType>::computeBoundsFromPoints(PointType begin, PointType end, F getPoint)
	{
		xmin_ = ymin_ = zmin_ = std::numeric_limits<double>::max();
		xmax_ = ymax_ = zmax_ = std::numeric_limits<double>::min();

		auto adjust = [](double coord, double& min, double& max) {
			if (coord < min) { min = coord; }
			if (coord > max) { max = coord; }
		};

		for (auto i = begin; i < end; ++i) {
			auto&& p = getPoint(i);
			adjust(p.x, xmin_, xmax_);
			adjust(p.y, ymin_, ymax_);
			adjust(p.z, zmin_, zmax_);
		}

		auto separate = [](double& /*min*/, double& /*max*/) {
			// if (min == max) {
			// 	max += std::numeric_limits<double>::epsilon();
			// }
		};

		separate(xmin_, xmax_);
		separate(ymin_, ymax_);
		separate(zmin_, zmax_);
	}

	template<typename PointType>
	void Celery<PointType>::computeCellData(std::size_t numPoints)
	{
		// Round up the number of cells per dimension
		double cellsPerDim = (double) numPoints;
		double cellDensityPerDim = (double) cell_density_;

		// Compute the number of cells and inverse cell size in each dimension
		// Round up by adding 1 after truncation
		num_cells_dim_ = std::size_t( std::cbrt( cellsPerDim / cellDensityPerDim ) ) + 1;
		std::cerr << "num_cells_dim_: " << num_cells_dim_ << std::endl;

		cell_size_x_ = (xmax_ - xmin_) / num_cells_dim_;
		cell_size_y_ = (ymax_ - ymin_) / num_cells_dim_;
		cell_size_z_ = (zmax_ - zmin_) / num_cells_dim_;

		cell_size_inv_x_ = num_cells_dim_ / (xmax_ - xmin_);
		cell_size_inv_y_ = num_cells_dim_ / (ymax_ - ymin_);
		cell_size_inv_z_ = num_cells_dim_ / (zmax_ - zmin_);
	}

	template<typename PointType>
	template<typename F>
	void Celery<PointType>::insert(PointType point, F getPoint)
	{
		points_.push_back(point);
		cells_.push_back(getCell(getPoint(point)));
	}

	template<typename PointType>
	template<typename PointIterator, typename F>
	void Celery<PointType>::insert(PointIterator begin, PointIterator end, F getPoint)
	{
		for (auto i = begin; i != end; ++i) {
			insert(*i, getPoint);
		}
	}

	template<typename PointType>
	template<typename F>
	void Celery<PointType>::insert(PointType begin, PointType end, F getPoint)
	{
		for (auto i = begin; i != end; ++i) {
			insert(i, getPoint);
		}
	}

	template<typename PointType>
	void Celery<PointType>::initializeCellArray()
	{
		// Create a vector of indices for cells_ and points_
		std::vector<std::size_t> indices(cells_.size(), 0);
		for (std::size_t i = 0; i < indices.size(); ++i) {
			indices[i] = i;
		}

		// Sort points_ and cells_ together, both using cells_
		std::sort(indices.begin(), indices.end(), MyComparator(cells_));
		std::vector<PointType> pointsCopy = points_;
		for (std::size_t i = 0; i < points_.size(); ++i) {
			points_[i] = pointsCopy[indices[i]];
		}

		std::sort(cells_.begin(), cells_.end());

		// Build vector of delimiters
		delimiters_.push_back(0);
		std::size_t numPoints = cells_.size();

		std::size_t lastCell = num_cells_dim_ * num_cells_dim_ * num_cells_dim_;
		for (std::size_t i = 0; i < lastCell; ++i) {
			std::size_t offset = 0;
			std::size_t back = delimiters_.back();

			while (i == cells_[back + offset]) {
				++offset;

				// Avoid going out of bounds
				if (back + offset == numPoints) {
					for (; i < lastCell; ++i) {
						delimiters_.push_back(numPoints);
					}

					return;
				}
			}

			delimiters_.push_back(back + offset);
		}

		delimiters_.push_back(numPoints);
	}

	template<typename PointType>
	bool Celery<PointType>::checkCellInRange(double x, double y, double z, double dist,
											 std::size_t xIndex, std::size_t yIndex, std::size_t zIndex) const
	{
		auto sq = [&](double yomama) -> double {
			return yomama * yomama;
		};

		auto distance = [&](int i, int j, int k) -> double {
			return sq(i * cell_size_x_) + sq(j * cell_size_y_) + sq(k * cell_size_z_);
		};

		auto offset = [](int coord, int index) -> int {
			return std::max(0, std::abs(coord - index) - 1);
		};
		return distance(offset(getXIndex(x), xIndex), offset(getYIndex(y), yIndex), offset(getZIndex(z), zIndex)) <= dist*dist;
	}

	template<typename PointType>
	void Celery<PointType>::findNeighborsInCellRadius(double x, double y, double z, double radius,
													  std::vector<PointType*>& pts) const
	{
		double xlow = std::max(x - radius, xmin_);
		double ylow = std::max(y - radius, ymin_);
		double zlow = std::max(z - radius, zmin_);

		double xhigh = std::min(x + radius, xmax_);
		double yhigh = std::min(y + radius, ymax_);
		double zhigh = std::min(z + radius, zmax_);

		std::size_t xMaxIndex = getXIndex(xhigh);
		std::size_t yMaxIndex = getYIndex(yhigh);
		std::size_t zMaxIndex = getZIndex(zhigh);

		for (std::size_t i = getXIndex(xlow); i <= xMaxIndex; ++i) {
			for (std::size_t j = getYIndex(ylow); j <= yMaxIndex; ++j) {
				for (std::size_t k = getZIndex(zlow); k <= zMaxIndex; ++k) {

					// Since we search through cells in a grid, some cells (i.e. corners) might actually
					// be outside of the search radius. Therefore, we check for and skip these cells.
					if (checkCellInRange(x, y, z, radius, i, j, k)) {

						std::size_t cell = getCellFromIndices(i, j, k);
						for (std::size_t pt = delimiters_[cell]; pt < delimiters_[cell+1]; ++pt) {
							pts.push_back(const_cast<PointType*>(&points_[pt]));
						}
					}
				}
			}
		}
	}

	template<typename PointType>
	void Celery<PointType>::findNeighborsInCellRadius(double x, double y, double z, double radius,
													  std::vector<PointType>& pts) const
	{
		double xlow = std::max(x - radius, xmin_);
		double ylow = std::max(y - radius, ymin_);
		double zlow = std::max(z - radius, zmin_);

		double xhigh = std::min(x + radius, xmax_);
		double yhigh = std::min(y + radius, ymax_);
		double zhigh = std::min(z + radius, zmax_);

		std::size_t xMaxIndex = getXIndex(xhigh);
		std::size_t yMaxIndex = getYIndex(yhigh);
		std::size_t zMaxIndex = getZIndex(zhigh);

		for (std::size_t i = getXIndex(xlow); i <= xMaxIndex; ++i) {
			for (std::size_t j = getYIndex(ylow); j <= yMaxIndex; ++j) {
				for (std::size_t k = getZIndex(zlow); k <= zMaxIndex; ++k) {

					// Since we search through cells in a grid, some cells (i.e. corners) might actually
					// be outside of the search radius. Therefore, we check for and skip these cells.
					if (checkCellInRange(x, y, z, radius, i, j, k)) {

						std::size_t cell = getCellFromIndices(i, j, k);
						for (std::size_t pt = delimiters_[cell]; pt < delimiters_[cell+1]; ++pt) {
							pts.push_back(points_[pt]);
						}
					}
				}
			}
		}
	}


	template<typename PointType>
	void Celery<PointType>::findNeighborsInRealRadius(double x, double y, double z, double radius,
													  std::vector<PointType*>& pts) const
	{
		double xlow = std::max(x - radius, xmin_);
		double ylow = std::max(y - radius, ymin_);
		double zlow = std::max(z - radius, zmin_);

		double xhigh = std::min(x + radius, xmax_);
		double yhigh = std::min(y + radius, ymax_);
		double zhigh = std::min(z + radius, zmax_);

		std::size_t xMaxIndex = getXIndex(xhigh);
		std::size_t yMaxIndex = getYIndex(yhigh);
		std::size_t zMaxIndex = getZIndex(zhigh);

		for (std::size_t i = getXIndex(xlow); i <= xMaxIndex; ++i) {
			for (std::size_t j = getYIndex(ylow); j <= yMaxIndex; ++j) {
				for (std::size_t k = getZIndex(zlow); k <= zMaxIndex; ++k) {

					// Since we search through cells in a grid, some cells (i.e. corners) might actually
					// be outside of the search radius. Therefore, we check for and skip these cells.
					if (checkCellInRange(x, y, z, radius, i, j, k)) {

						std::size_t cell = getCellFromIndices(i, j, k);
						for (std::size_t pt = delimiters_[cell]; pt < delimiters_[cell+1]; ++pt) {

							// Some points in each cell may not be within the search radius, so we check for
							// and skip these points.
							if (euclidDist(x, y, z, points_[pt].x, points_[pt].y, points_[pt].z) <= radius) {
								pts.push_back(&points_[pt]);
							}
						}
					}
				}
			}
		}
	}

	template <typename PointType>
	void neighborQuery(std::array<double, 3> point, std::vector<PointType*> NeighborParticles) {

	}

	template<typename PointType>
	void Celery<PointType>::createSearchArray()
	{

		auto sq = [&](double yomama) -> double {
			return yomama * yomama;
		};

		auto distance = [&](int i, int j, int k) -> double {
			return sq(i * cell_size_x_) + sq(j * cell_size_y_) + sq(k * cell_size_z_);
		};

		int maxIndex = num_cells_dim_ - 1;

		auto cb = [](int n) { return n * n * n; };
		search_order_.reserve(cb(2 * maxIndex + 1));

		search_order_.emplace_back(-1, 0, 0, 0);

		for (int i = 0; i < maxIndex; ++i) {
			for (int j = 0; j < maxIndex; ++j) {
				for (int k = 0; k < maxIndex; ++k) {

					double dist = distance(i, j, k);
					search_order_.emplace_back(dist, i+1, j+1, k+1);
					search_order_.emplace_back(dist, i+1, j+1, -k-1);
					search_order_.emplace_back(dist, i+1, -j-1, k+1);
					search_order_.emplace_back(dist, -i-1, j+1, k+1);
					search_order_.emplace_back(dist, i+1, -j-1, -k-1);
					search_order_.emplace_back(dist, -i-1, j+1, -k-1);
					search_order_.emplace_back(dist, -i-1, -j-1, k+1);
					search_order_.emplace_back(dist, -i-1, -j-1, -k-1);
				}
			}
		}

		for (int i = 0; i < maxIndex; ++i) {
			for (int j = 0; j < maxIndex; ++j) {

				double dist = distance(i, j, 0);
				search_order_.emplace_back(dist, i+1, j+1, 0);
				search_order_.emplace_back(dist, i+1, -j-1, 0);
				search_order_.emplace_back(dist, -i-1, j+1, 0);
				search_order_.emplace_back(dist, -i-1, -j-1, 0);
			}
		}

		for (int i = 0; i < maxIndex; ++i) {
			for (int k = 0; k < maxIndex; ++k) {

				double dist = distance(i, 0, k);
				search_order_.emplace_back(dist, i+1, 0, k+1);
				search_order_.emplace_back(dist, i+1, 0, -k-1);
				search_order_.emplace_back(dist, -i-1, 0, k+1);
				search_order_.emplace_back(dist, -i-1, 0, -k-1);
			}
		}

		for (int j = 0; j < maxIndex; ++j) {
			for (int k = 0; k < maxIndex; ++k) {

				double dist = distance(0, j, k);
				search_order_.emplace_back(dist, 0, j+1, k+1);
				search_order_.emplace_back(dist, 0, j+1, -k-1);
				search_order_.emplace_back(dist, 0, -j-1, k+1);
				search_order_.emplace_back(dist, 0, -j-1, -k-1);
			}
		}

		for (int i = 0; i < maxIndex; ++i) {
			double dist = distance(i, 0, 0);
			search_order_.emplace_back(dist, i+1, 0, 0);
			search_order_.emplace_back(dist, -i-1, 0, 0);
		}

		for (int j = 0; j < maxIndex; ++j) {
			double dist = distance(0, j, 0);
			search_order_.emplace_back(dist, 0, j+1, 0);
			search_order_.emplace_back(dist, 0, -j-1, 0);
		}

		for (int k = 0; k < maxIndex; ++k) {
			double dist = distance(0, 0, k);
			search_order_.emplace_back(dist, 0, 0, k+1);
			search_order_.emplace_back(dist, 0, 0, -k-1);
		}

		std::sort(search_order_.begin(), search_order_.end());
	}

	template <typename PointType>
	bool Celery<PointType>::findNeighborsInShell(double x, double y, double z, int shell, double maxRadius, std::vector<PointType>& pts) const
	{
		auto addPoints = [&](unsigned c) -> void {
			++cellsSearched;
			for (unsigned pi = delimiters_[c]; pi < delimiters_[c+1]; ++pi) {
				pts.push_back(points_[pi]);
			}
		};

		auto valid = [&](unsigned index) -> bool {
			// By casting to unsigned, negative values
			// will wrap around to very large ones.
			// This way, we do not have to check
			// index >= 0 explicitly; that will be
			// caught by the index < num_cells_dim_
			// check.
			return index < num_cells_dim_;
		};

		if (shell == 0) {
			addPoints(getCell(x, y, z));
			return true;
		}

		int xi = getXIndex(x), yi = getYIndex(y), zi = getZIndex(z);
		auto inv = std::max(cell_size_inv_x_, std::max(cell_size_inv_y_, cell_size_inv_z_));
		decltype(shell) maxShell = maxRadius * inv + 1;

		if (shell > maxShell) { return false; }
		auto innerShell = shell - 1;


		for (int i : {xi - shell, xi + shell}) if (valid(i))
		for (int j = yi - shell; j <= yi + shell; ++j) if (valid(j))
		for (int k = zi - shell; k <= zi + shell; ++k) if (valid(k))
		{
			addPoints(getCellFromIndices(i, j, k));
		}

		for (int j : {yi - shell, yi + shell}) if (valid(j))
		for (int i = xi - innerShell; i <= xi + innerShell; ++i) if (valid(i))
		for (int k = zi - shell; k <= zi + shell; ++k) if (valid(k))
		{
			addPoints(getCellFromIndices(i, j, k));
		}

		for (int k : {zi - shell, zi + shell}) if (valid(k))
		for (int i = xi - innerShell; i <= xi + innerShell; ++i) if (valid(i))
		for (int j = yi - innerShell; j <= yi + innerShell; ++j) if (valid(j))
		{
			addPoints(getCellFromIndices(i, j, k));
		}

		return true;
	}

}}
