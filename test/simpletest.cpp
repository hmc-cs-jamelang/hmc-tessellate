#include <iostream>

#include "../src/utilities.hpp"
#include "../src/hmc-tessellate.hpp"


using namespace hmc;
using namespace hmc::utilities;

int main()
{
    Polyhedron p;
    std::cerr << p << std::endl;
    Diagram d = Diagram::cube(0,1,0,1,0,1);
    d.initialize([&d](){
        for (unsigned i = 0; i < 25; ++i) {
            d.addParticle(uniform(0,1), uniform(0,1), uniform(0,1));
        }
    });

    double volume = 0;
    Cell c;
    for (unsigned i = 0; i < d.size(); ++i) {
        c = d.getCell(i);
        volume += c.computeVolume();
    }

    std::cout << "Volume: " << volume << std::endl;
    return 0;
}
