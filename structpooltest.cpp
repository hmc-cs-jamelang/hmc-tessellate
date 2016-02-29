#include <iostream>
#include <assert.h>
#include <bitset>
#include "structpool.hpp"


struct Vertex {
	double x, y, z;

	Vertex(double x, double y, double z) : x(x), y(y), z(z) {}
};

std::ostream& operator<<(std::ostream& os, const Vertex& v)
{
	return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

struct HalfEdge;

typedef StructPool<HalfEdge>::Index EdgeIndex;
typedef StructPool<Vertex>::Index VertexIndex;
constexpr EdgeIndex INVALID_EDGE = StructPool<HalfEdge>::INVALID_INDEX;
constexpr VertexIndex INVALID_VERTEX = StructPool<Vertex>::INVALID_INDEX;

struct HalfEdge {
	EdgeIndex flip;
	EdgeIndex next;
	VertexIndex target;

	HalfEdge(EdgeIndex flip, EdgeIndex next, VertexIndex target)
		: flip(flip), next(next), target(target)
	{}

	// HalfEdge() : flip(INVALID_EDGE), next(INVALID_EDGE), target(INVALID_VERTEX) {}
};

union union_thingy {
  	EdgeIndex index;
  	HalfEdge obj;
};

int main()
{
	std::cout << "HalfEdge size: " << sizeof(HalfEdge) << std::endl;
	std::cout << "StructPool<HalfEdge> union size: "
			  << sizeof(union_thingy)
			  << std::endl;
	std::cout << "Bitset size: "
			  << sizeof(std::bitset<2>)
			  << std::endl;
	std::cout << "StructPool<HalfEdge> chunk size: "
			  << sizeof(StructPool<HalfEdge>::Chunk)
			  << std::endl;

	StructPool<HalfEdge> edges {13};
	StructPool<Vertex> vertices;

	auto v1 = vertices.create(11.0,12.0,13.0);
	auto v2 = vertices.create(21.0,22.0,23.0);
	auto v3 = vertices.create(31.0,32.0,33.0);
	auto v4 = vertices.create(41.0,42.0,43.0);
	auto v5 = vertices.create(51.0,52.0,53.0);

	auto e1 = edges.create(INVALID_EDGE, INVALID_EDGE, v1);
	auto e2 = edges.create(e1, INVALID_EDGE, v2);
	auto e3 = edges.create(e2, e1, v3);
	auto e4 = edges.create(e3, e3, v4);
	auto e5 = edges.create(e1, e2, v5);

	std::cout << "Invalid: " << INVALID_EDGE << std::endl;
	std::cout << "Invalid as double: " << *((double*) &INVALID_EDGE) << std::endl;

	EdgeIndex uninit;
	std::cout << "Default constructed: " << uninit << std::endl;

	std::cout << v2 << ": " << vertices[v2] << std::endl;
	std::cout << vertices[edges[e1].target] << std::endl;

	vertices.destroy(v1);
	vertices.destroy(v3);
	// This should error out in debug mode
	// (i.e., without the flag -DNVERIFY)
	// std::cout << vertices[edges[e1].target] << std::endl;
	assert(vertices.inactive(v1));

	vertices.create(0.0,0.0,0.0);

	std::cout << "Vertices: Occupancy: " << vertices.computeOccupancy()
			  << " Size: " << vertices.size()
			  << " Capacity: " << vertices.capacity()
			  << std::endl;

	return 0;
}
