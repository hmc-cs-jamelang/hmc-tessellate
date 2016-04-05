#include <iostream>

#include "utilities.hpp"
#include "hmc-tessellate.hpp"


using namespace hmc;

int main()
{
    Polyhedron p; p.fullOutput();
    Diagram d = Diagram::cube(0,1,0,1,0,1);
    d.initialize([&d](){
        for (unsigned i = 0; i < 25; ++i) {
            d.addParticle(i, uniform(0,1), uniform(0,1), uniform(0,1));
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
