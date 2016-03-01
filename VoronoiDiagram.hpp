#ifndef VORONOI_DIAGRAM_HPP
#define VORONOI_DIAGRAM_HPP

#pragma once
#include <iterator>
#include <vector>
#include <unordered_set>
#include "vorocode.hpp"

namespace voro {

	using Polyhedron = voronoiCell;
	using our_size_t = int;

	using TargetGroup = std::unordered_set<int>;
	using SourceGroup = std::vector<int>;

	const int DEFAULT_GROUP = -1;
	const our_size_t DEFAULT_INDEX = -1;

	class Diagram;

	inline Polyhedron* computePolyhedron(Diagram*, our_size_t);
	inline Polyhedron* computePolyhedron(Diagram*, double, double, double);
	inline Polyhedron* computePolyhedron(Diagram*, our_size_t, TargetGroup);
	inline Polyhedron* computePolyhedron(Diagram*, double, double, double, TargetGroup);
	inline int getParticleIdFromIndex(Diagram*, our_size_t);

	class Cell {
	public:
		Cell() = default;

		Cell(Diagram* diagram, our_size_t point)
			// Needed to add no-argument constructor for voronoiCell to get this to work
			: diagram_(diagram), point_(point), poly_computed_(false)
		{
			// Nothing going on here
		}

		Cell(Diagram* diagram, our_size_t point, TargetGroup targetGroup)
			: diagram_(diagram), point_(point), poly_computed_(false), target_group_(targetGroup)
		{

		}

		Cell(Diagram* diagram, double x, double y, double z)
			// Needed to add no-argument constructor for voronoiCell to get this to work
			: diagram_(diagram), point_(DEFAULT_INDEX), x_(x), y_(y), z_(z), poly_computed_(false)
		{
			// Nothing going on here
		}

		Cell(Diagram* diagram, double x, double y, double z, TargetGroup targetGroup)
			: diagram_(diagram), point_(DEFAULT_INDEX), x_(x), y_(y), z_(z), poly_computed_(false), target_group_(targetGroup)
		{

		}

		int id() {
			return getParticleIdFromIndex(diagram_, point_);
		}

		double computeVolume() {
			ensurePolyComputed();
			return poly_->volume();
		}

		void computeVertices(std::vector<double> &v) {
			ensurePolyComputed();
			poly_->vertices(v);
		}
		void computeVertices(double x, double y, double z, std::vector<double> &v) {
			return computeVertices(v);
		}


		void computeFaceVertices(std::vector<int> &v) {
			//ensurePolyComputed();
			//poly_.face_vertices(v);
		}
		virtual void computeNeighbors(std::vector<int> &v) {
			ensurePolyComputed();
			poly_->neighbors(v);
		}
		void computeFaceAreas(std::vector<double> &v) {
			ensurePolyComputed();
			poly_->face_areas(v);
		}

		Cell& operator=(const Cell& cell) = default;

	private:
		void ensurePolyComputed() {
			if (!poly_computed_) {
				if (!target_group_.empty()) {
					delete poly_;
					if (point_ != DEFAULT_INDEX) {
						poly_ = computePolyhedron(diagram_, point_, target_group_);
					}
					else {
						poly_ = computePolyhedron(diagram_, x_, y_, z_, target_group_);
					}
				}
				else {
					delete poly_;
					if (point_ != DEFAULT_INDEX) {
						poly_ = computePolyhedron(diagram_, point_);
					}
					else {
						poly_ = computePolyhedron(diagram_, x_, y_, z_);
					}
				}
			}
		}

		Polyhedron* poly_ = nullptr;
		Diagram* diagram_;
		our_size_t point_;
		bool poly_computed_;
		TargetGroup target_group_;
		double x_, y_, z_;
	};

	class Diagram {
	public:
		Diagram() = default;
		
		Diagram(double xmin, double xMAX, double ymin, double yMAX, double zmin, double zMAX)
			: xmin_(xmin), xMAX_(xMAX), ymin_(ymin), yMAX_(yMAX), zmin_(zmin), zMAX_(zMAX)
		{
			// Nothing happening
		}

		template<class c_class>
		Diagram(c_class &con)
		{
			// Nothing to do here
		}

		our_size_t size()
		{
			return particles_.size();
		}

		our_size_t addParticle(double x, double y, double z, int id)
		{
			our_size_t index = size();
			particles_.push_back(Particle(id, x, y, z));
			groups_.push_back(DEFAULT_GROUP);
			return index;
		}

		our_size_t addParticle(double x, double y, double z, int id, int group)
		{
			our_size_t index = size();
			particles_.push_back(Particle(id, x, y, z));
			groups_.push_back(group);
			return index;
		}

		template <typename... Groups>
		TargetGroup targetGroups(Groups... groups)
		{
			return TargetGroup {groups...};
		}

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

		Cell operator[] (our_size_t i)
		{
			return Cell(this, i);
		}

		Cell getCell(our_size_t i)
		{
			return Cell(this, i);
		}

		Cell getCell(our_size_t i, TargetGroup targetGroup)
		{
			return Cell(this, i, targetGroup);
		}

		Cell getCell(double x, double y, double z)
		{
			return Cell(this, x, y, z);
		}

		Cell getCell(double x, double y, double z, TargetGroup targetGroup)
		{
			return Cell(this, x, y, z, targetGroup);
		}

		our_size_t numCells()
		{
			return 0;
		}

		Polyhedron* computePolyhedron(our_size_t i)
		{
			Polyhedron* poly = new Polyhedron("cube", 100000000, particles_[i], 100000000, xmin_, xMAX_, ymin_, yMAX_, zmin_, zMAX_);

			our_size_t end = particles_.size();
			for (our_size_t index = 0; index < end; ++index) {
				if (index != i) {
					poly->cutCell(particles_[index]);
				}
			}

			return poly;
		}

		Polyhedron* computePolyhedron(double x, double y, double z)
		{
			Particle p = Particle(DEFAULT_INDEX, x, y, z);
			Polyhedron* poly = new Polyhedron("cube", 10000000, p, 10000000, xmin_, xMAX_, ymin_, yMAX_, zmin_, zMAX_);

			for (auto particle : particles_) {
				poly->cutCell(particle);
			}

			return poly;
		}

		Polyhedron* computePolyhedron(our_size_t i, TargetGroup targetGroup)
		{
			Polyhedron* poly = new Polyhedron("cube", 10000000, particles_[i], 10000000, xmin_, xMAX_, ymin_, yMAX_, zmin_, zMAX_);

			our_size_t end = particles_.size();
			for (our_size_t index = 0; index < end; ++index) {
				if (targetGroup.count(groups_[index]) && index != i) {
					poly->cutCell(particles_[index]);
				}
			}

			return poly;
		}

		Polyhedron* computePolyhedron(double x, double y, double z, TargetGroup targetGroup)
		{
			Particle p = Particle(DEFAULT_INDEX, x, y, z);
			Polyhedron* poly = new Polyhedron("cube", 10000000, p, 10000000, xmin_, xMAX_, ymin_, yMAX_, zmin_, zMAX_);

			our_size_t end = particles_.size();
			for (our_size_t index = 0; index < end; ++index) {
				if (targetGroup.count(groups_[index])) {
					poly->cutCell(particles_[index]);
				}
			}

			return poly;
		}

		void draw_particles(const char *filename) {}
		void draw_cells_gnuplot(const char *filename) {}
		
//private:
		std::vector<Particle> particles_;
		std::vector<int> groups_;
		double xmin_, xMAX_, ymin_, yMAX_, zmin_, zMAX_;
	};

	inline Polyhedron* computePolyhedron(Diagram* diagram, our_size_t i) {
		return diagram->computePolyhedron(i);
	}

	inline Polyhedron* computePolyhedron(Diagram* diagram, double x, double y, double z) {
		return diagram->computePolyhedron(x, y, z);
	}

	inline Polyhedron* computePolyhedron(Diagram* diagram, our_size_t i, TargetGroup targetGroup) {
		return diagram->computePolyhedron(i, targetGroup);
	}

	inline Polyhedron* computePolyhedron(Diagram* diagram, double x, double y, double z, TargetGroup targetGroup) {
		return diagram->computePolyhedron(x, y, z, targetGroup);
	}

	inline int getParticleIdFromIndex(Diagram* diagram, our_size_t i) {
		return diagram->particles_[i].id;
	}

}

		

#endif
