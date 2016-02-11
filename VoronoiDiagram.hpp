#ifndef VORONOI_DIAGRAM_HPP
#define VORONOI_DIAGRAM_HPP

#pragma once
#include <iterator>
#include <vector>
#include "vorocode.hpp"

namespace voro {

	typedef voronoiCell Polyhedron;

	typedef int our_size_t;

	// Already in vorocode

	// typedef struct Particle {
	// 	CGALInterface::Point point;
	// 	our_size_t index;
	// 	int id;
	// 	int group;

	// 	Particle(double x, double y, double z, our_size_t index, int id, int group = 0)
	// 		: point(x,y,z), index(index), id(id), group(group)
	// 	{
	// 		// Nothing left to do
	// 	}
	// } Particle;

	class Blob;

	inline void computePolyhedron(Blob*, our_size_t, Polyhedron&);
	inline void computePolyhedron(Blob*, our_size_t, Polyhedron&, int);
	inline int getParticleIdFromIndex(Blob*, our_size_t);

	class FatCell {
	public:
		FatCell() = default;

		FatCell(Blob* diagram, our_size_t point)
			// Needed to add no-argument constructor for voronoiCell to get this to work
			: poly_(), diagram_(diagram), point_(point), poly_computed_(false)
		{

		}

		FatCell(Blob* diagram, our_size_t point, int targetGroup)
			: poly_(), diagram_(diagram), point_(point), poly_computed_(false), target_group_(targetGroup)
		{

		}

		int id() {
			return getParticleIdFromIndex(diagram_, point_);
		}

		double volume() {
			ensurePolyComputed();
			return poly_.volume();
		}

		void vertices(std::vector<double> &v) {
			ensurePolyComputed();
			poly_.vertices(v);
		}
		void vertices(double x, double y, double z, std::vector<double> &v) {
			return vertices(v);
		}


		void face_vertices(std::vector<int> &v) {
			//ensurePolyComputed();
			//poly_.face_vertices(v);
		}
		virtual void neighbors(std::vector<int> &v) {
			ensurePolyComputed();
			poly_.neighbors(v);
		}
		void face_areas(std::vector<double> &v) {
			ensurePolyComputed();
			poly_.face_areas(v);
		}

		// So apparently the compiler won't make a default assignment operator
		// WHYYYYYYYYYYYY
		FatCell& operator=(const FatCell& cell)
		{
			return *this;
		}

	private:
		void ensurePolyComputed() {
			if (!poly_computed_) {
				if (target_group_) {
					computePolyhedron(diagram_, point_, poly_, target_group_);
				}
				else {
					computePolyhedron(diagram_, point_, poly_);
				}
			}
		}

		Polyhedron poly_;
		Blob* diagram_;
		our_size_t point_;
		bool poly_computed_;
		int target_group_ = 0;
	};

	class Blob {
	public:
		Blob() = delete;
		
		Blob(double xmin, double xMAX, double ymin, double yMAX, double zmin, double zMAX)
			: xmin_(xmin), xMAX_(xMAX), ymin_(ymin), yMAX_(yMAX), zmin_(zmin), zMAX_(zMAX)
		{
			// Nothing left to do
		}

		template<class c_class>
		Blob(c_class &con)
		{
			// Do nothing
		}

		our_size_t size()
		{
			return particles_.size();
		}

		our_size_t put(int id, double x, double y, double z)
		{
			our_size_t index = size();
			// Needed to add new Particle constructor to get this to work
			particles_.push_back(Particle(x, y, z, index, id));
			return index;
		}

		our_size_t put(std::vector<int> group, int id, double x, double y, double z)
		{
			our_size_t index = put(id, x, y, z);
			group.push_back(index);
			return index;
		}

		our_size_t put(int group, int id, double x, double y, double z)
		{
			our_size_t index = size();
			particles_.push_back(Particle(x, y, z, index, id, group));
			return index;
		}

		FatCell operator[] (our_size_t i)
		{
			return FatCell(this, i);
		}

		FatCell getCell(our_size_t i, int group)
		{
			return FatCell(this, i, group);
		}

		our_size_t numCells()
		{
			return 0;
		}

		void computePolyhedron(our_size_t i, Polyhedron& poly)
		{
			// Needed to add clear() to voronoiCell
			poly.clear();
			// Needed to add initialize to voronoiCell
			poly.initialize(particles_[i].position, xmin_, xMAX_, ymin_, yMAX_, zmin_, zMAX_);
			for (auto particle : particles_) {
				poly.cutCell(particle);
			}
		}

		void computePolyhedron(our_size_t i, Polyhedron& poly, int targetGroup)
		{
			poly.clear();
			poly.initialize(particles_[i].position, xmin_, xMAX_, ymin_, yMAX_, zmin_, zMAX_);
			for (auto particle : particles_) {
				// Needed to add group member to Particle
				if (particle.group == targetGroup) {
					poly.cutCell(particle);
				}
			}
		}

		void draw_particles(const char *filename) {}
		void draw_cells_gnuplot(const char *filename) {}
		
//private:
		std::vector<Particle> particles_;
		double xmin_, xMAX_, ymin_, yMAX_, zmin_, zMAX_;
	};

	inline void computePolyhedron(Blob* diagram, our_size_t i, Polyhedron& poly) {
		diagram->computePolyhedron(i, poly);
	}

	inline void computePolyhedron(Blob* diagram, our_size_t i, Polyhedron& poly, int targetGroup) {
		diagram->computePolyhedron(i, poly, targetGroup);
	}

	inline int getParticleIdFromIndex(Blob* diagram, our_size_t i) {
		return diagram->particles_[i].id;
	}

}


#endif
