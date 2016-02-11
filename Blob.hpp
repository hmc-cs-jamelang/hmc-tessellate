// This is a dummy interface for the blob class we have been discussing.
// Currently, it exists only to allow us to compile MOAB using our
// interface design. However, if the real interface follows this design,
// we should not need to make any significant changes in MOAB.

// #pragma once
// #include "FatCell.hpp"
// #include <iterator>

// namespace voro {

// 	class Blob {
// 	public:
// 		template<class c_class>
// 		Blob(c_class &con)
// 		{
// 			// Do nothing
// 		}

// 		void put(int n, double x, double y, double z) {}

// 		FatCell operator[] (int i)
// 		{
// 			return FatCell();
// 		}

// 		unsigned numCells()
// 		{
// 			return 0;
// 		}
// 	};

// }
#pragma once
#include "VoronoiDiagram.hpp"
