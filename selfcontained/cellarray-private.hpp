/*
 * cellarray-private.hpp
 *
 * Implementation file for cellarray.hpp
 */

namespace spatial
{

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

		computeBoundsFromPoints(begin, end, getPoint);
		computeCellData(std::distance(begin, end));
		insert(begin, end, getPoint);
		initializeCellArray();
	}

	template<typename PointType>
	template<typename F>
	void Celery<PointType>::initialize(PointType begin, PointType end, F getPoint)
	{
		// Clear any existing point information
		points_.clear();
		cells_.clear();
		delimiters_.clear();

		computeBoundsFromPoints(begin, end, getPoint);
		computeCellData(end - begin);
		insert(begin, end, getPoint);
		initializeCellArray();
	}

	template<typename PointType>
	template<typename PointIterator, typename F>
	void Celery<PointType>::computeBoundsFromPoints(PointIterator begin, PointIterator end, F getPoint)
	{
		std::vector<double> xvals, yvals, zvals;
		for (auto i = begin; i < end; ++i) {
			auto&& p = getPoint(*i);
			xvals.push_back(p.x);
			yvals.push_back(p.y);
			zvals.push_back(p.z);
		}

		std::sort(xvals.begin(), xvals.end());
		std::sort(yvals.begin(), yvals.end());
		std::sort(zvals.begin(), zvals.end());

		xmin_ = xvals.front();
		xmax_ = xvals.back();
		ymin_ = yvals.front();
		ymax_ = yvals.back();
		zmin_ = zvals.front();
		zmax_ = zvals.back();

		// Make the bounding box slightly larger than the extreme values of the points
		double xscl = 0.01 * (xmax_ - xmin_);
		double yscl = 0.01 * (ymax_ - ymin_);
		double zscl = 0.01 * (zmax_ - zmin_);

		xmin_ -= xscl;
		xmax_ += xscl;
		ymin_ -= zscl;
		ymax_ += zscl;
		zmin_ -= zscl;
		zmax_ += zscl;
	}

	template<typename PointType>
	template<typename F>
	void Celery<PointType>::computeBoundsFromPoints(PointType begin, PointType end, F getPoint)
	{
		std::vector<double> xvals, yvals, zvals;
		for (auto i = begin; i < end; ++i) {
			auto&& p = getPoint(i);
			xvals.push_back(p.x);
			yvals.push_back(p.y);
			zvals.push_back(p.z);
		}

		std::sort(xvals.begin(), xvals.end());
		std::sort(yvals.begin(), yvals.end());
		std::sort(zvals.begin(), zvals.end());

		xmin_ = xvals.front();
		xmax_ = xvals.back();
		ymin_ = yvals.front();
		ymax_ = yvals.back();
		zmin_ = zvals.front();
		zmax_ = zvals.back();

		// Make the bounding box slightly larger than the extreme values of the points
		double xscl = 0.01 * (xmax_ - xmin_);
		double yscl = 0.01 * (ymax_ - ymin_);
		double zscl = 0.01 * (zmax_ - zmin_);

		xmin_ -= xscl;
		xmax_ += xscl;
		ymin_ -= yscl;
		ymax_ += yscl;
		zmin_ -= zscl;
		zmax_ += zscl;
	}

	template<typename PointType>
	void Celery<PointType>::computeCellData(std::size_t numPoints)
	{
		// Round up the number of cells per dimension
		std::size_t cellsPerDim = std::size_t( cbrt((double) numPoints) ) + 1;
		std::size_t cellDensityPerDim = std::size_t( cbrt((double) cell_density_) ) + 1;

		// Compute the number of cells and inverse cell size in each dimension
		// Round up by adding 1 after truncation
		num_cells_dim_ = std::size_t( (double) cellsPerDim / cellDensityPerDim ) + 1;

		cell_size_inv_x_ = (double) num_cells_dim_ / (xmax_ - xmin_);
		cell_size_inv_y_ = (double) num_cells_dim_ / (ymax_ - zmin_);
		cell_size_inv_z_ = (double) num_cells_dim_ / (ymax_ - zmin_);
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
		reserve(std::distance(begin, end));
		for (auto i = begin; i != end; ++i) {
			insert(*i, getPoint);
		}
	}

	template<typename PointType>
	template<typename F>
	void Celery<PointType>::insert(PointType begin, PointType end, F getPoint)
	{
		//reserve(std::distance(begin, end));
		reserve(end - begin);
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
		std::sort(indices.begin(), indices.end(),
			[this](std::size_t i, std::size_t j) {
				return cells_[i] < cells_[j];
			});
			// MyComparator(cells_));
		std::vector<PointType> pointsCopy = points_;
		for (std::size_t i = 0; i < points_.size(); ++i) {
			points_[i] = pointsCopy[indices[i]];
		}

		std::sort(cells_.begin(), cells_.end());

		// Build vector of delimiters
		delimiters_.push_back(0);
		std::size_t numPoints = cells_.size();

		std::size_t lastCell = pow(num_cells_dim_, 3);
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
											 std::size_t xIndex, std::size_t yIndex, std::size_t zIndex)
	const {
		double cellx, celly, cellz;

		if (getXIndex(x) > xIndex) {
			cellx = (xIndex + 1) / cell_size_inv_x_;
		}
		else {
			cellx = xIndex / cell_size_inv_x_;
		}

		if (getYIndex(y) > yIndex) {
			celly = (yIndex + 1) / cell_size_inv_y_;
		}
		else {
			celly = yIndex / cell_size_inv_y_;
		}

		if (getZIndex(z) > zIndex) {
			cellz = (zIndex + 1) / cell_size_inv_z_;
		}
		else {
			cellz = zIndex / cell_size_inv_z_;
		}

		return (euclidDist(x, y, z, cellx, celly, cellz) <= dist);
	}

	template<typename PointType>
	void Celery<PointType>::findNeighborsInCellRadius(double x, double y, double z, double radius,
													  std::vector<PointType*>& pts)
	const {
		double xlow = std::max(x - radius, xmin_);
		double ylow = std::max(y - radius, ymin_);
		double zlow = std::max(z - radius, zmin_);

		double xhigh = std::min(x + radius, xmax_);
		double yhigh = std::min(y + radius, ymax_);
		double zhigh = std::min(z + radius, zmax_);

		std::size_t xMaxIndex = getXIndex(xhigh);
		std::size_t yMaxIndex = getYIndex(yhigh);
		std::size_t zMaxIndex = getZIndex(zhigh);

		// std::cerr << "Loop, x: 0 -> " << xMaxIndex << ", "
		//           << "y: 0 -> " << yMaxIndex << ", "
		//           << "z: 0 -> " << zMaxIndex << std::endl;

		for (std::size_t i = getXIndex(xlow); i <= xMaxIndex; ++i) {
			for (std::size_t j = getYIndex(ylow); j <= yMaxIndex; ++j) {
				for (std::size_t k = getZIndex(zlow); k <= zMaxIndex; ++k) {

					// Since we search through cells in a grid, some cells (i.e. corners) might actually
					// be outside of the search radius. Therefore, we check for and skip these cells.
					if (checkCellInRange(x, y, z, radius, i, j, k)) {

						std::size_t cell = getCellFromIndices(i, j, k);
						for (std::size_t pt = delimiters_[cell]; pt < delimiters_[cell+1]; ++pt) {
							pts.push_back(&points_[pt]);
						}
					}
				}
			}
		}
	}

	template<typename PointType>
	void Celery<PointType>::findNeighborsInCellRadius(double x, double y, double z, double radius,
													  std::vector<PointType>& pts)
	const {
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
					//if (checkCellInRange(x, y, z, radius, i, j, k)) {

						std::size_t cell = getCellFromIndices(i, j, k);
						// if (cell + 1 >= delimiters_.size()) {
						// 	std::cerr << "Error. Indices " << i << ", " << j << ", " << k
						// 	          << ", calculated cell " << cell << " / " << delimiters_.size() -1 << std::endl;
						// 	continue;
						// }
						for (std::size_t pt = delimiters_[cell]; pt < delimiters_[cell+1]; ++pt) {
							pts.push_back(points_[pt]);
						}
					//}
				}
			}
		}
	}

	template<typename PointType>
	void Celery<PointType>::findNeighborsInRealRadius(double x, double y, double z, double radius,
													  std::vector<PointType*>& pts)
	const {
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
	bool Celery<PointType>::findNeighborsInShell(double x, double y, double z, unsigned shell, double maxRadius, std::vector<PointType>& pts) const
	{
		auto addPoints = [&](unsigned c) {
			for (unsigned pi = delimiters_[c]; pi < delimiters_[c+1]; ++pi) {
				pts.push_back(points_[pi]);
			}
		};

		if (shell == 0) {
			addPoints(getCell(x, y, z));
			return true;
		}

		auto xi = getXIndex(x), yi = getYIndex(y), zi = getZIndex(z);
		auto inv = std::min(cell_size_inv_x_, std::min(cell_size_inv_y_, cell_size_inv_z_));
		auto dist = (shell - 1) / inv;
		if (dist > maxRadius) { return false; }

		for (auto xs : {xi - shell, xi, xi + shell}) if (xs >= 0 && xs < num_cells_dim_) {
			for (auto ys : {yi - shell, yi, yi + shell}) if (ys >= 0 && ys < num_cells_dim_) {
				for (auto zs : {zi - shell, zi, zi + shell}) if (zs >= 0 && zs < num_cells_dim_) {
					if (xs == xi && ys == yi && zs == zi) {continue;}
					addPoints(getCellFromIndices(xs, ys, zs));
				}
			}
		}

		return true;
	}

	template <typename PointType>
	void neighborQuery(std::array<double, 3> point, std::vector<PointType*> NeighborParticles) {

	}

}
