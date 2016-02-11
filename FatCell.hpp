// This is a dummy interface for our FatCell class.
// Currently, it exists only to allow us to compile MOAB using our
// interface design. However, if the real interface follows this design,
// we should not need to make any significant changes in MOAB.

// #pragma once
// #include <vector>

// namespace voro {

// 	class FatCell {
// 	public:
// 		// For now, use the default constructor
// 		double volume() {
// 			return 0.0;
// 		}

// 		void vertices(std::vector<double> &v) {}
// 		void vertices(double x, double y, double z, std::vector<double> &v) {}
// 		void face_vertices(std::vector<int> &v) {}
// 		virtual void neighbors(std::vector<int> &v) {}
// 		void face_areas(std::vector<double> &v) {}
// 	};

// }
#pragma once
#include "VoronoiDiagram.hpp"
