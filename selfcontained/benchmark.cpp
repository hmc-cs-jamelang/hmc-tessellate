#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <assert.h>
#include <chrono>
#include <limits>
#include <set>
#include <cstdlib>

std::size_t totalCuts = 0,
	neededCuts = 0,
	expansions = 0,
	cellsSearched = 0;

auto dont_use_this_name = std::chrono::high_resolution_clock::now();
auto expansionTime = dont_use_this_name - dont_use_this_name;

#include "utilities.hpp"
#include "vectormath.hpp"
#include "hmc-tessellate.hpp"

#include <voro++.hh>

#if defined(COUNT_MALLOCS)
    #include "malloc_count-0.7/malloc_count.h"

    void outputMallocs(char const* preamble = "")
    {
        std::cerr << preamble
                  << "Malloc called: "
                  << getNumberOfTimesMallocHasBeenCalled()
                  << " times."
                  << std::endl;
    }
#else
    void resetNumberOfTimesMallocHasBeenCalled() {}
    void outputMallocs(char const* = "") {}
#endif

// struct Particle {
//     int id;
//     hmc::Vector3 position;

//     Particle(int id, double x, double y, double z)
//         : id(id), position(x, y, z)
//     {}
// };

using Particle = hmc::Particle;

using ParticleList = std::vector<Particle>;

const std::vector<double> defaultDists {};

enum Tessellator {voro_pp, hmc_tessellate};

template <Tessellator package, bool useDistances = false, typename CheckType>
void runTrial(const double boxLength,
              const ParticleList& particles,
              CheckType& record,
              const std::vector<double>& distances = defaultDists)
{
    const std::size_t numPoints = particles.size();
    const double bl2 = boxLength/2;

    if (package == voro_pp) {
        std::cerr << "Voro++, ";
    }
    else /* hmc_tessellate */ {
        std::cerr << "HMC";
        if (useDistances) {
            std::cerr << " with shrinkwrap, ";
        }
        else {
            std::cerr << " auto-finding, ";
        }
    }
    std::cerr << numPoints << " particles, "
              << "checking " << CheckType::checkedData() << std::endl;

    resetNumberOfTimesMallocHasBeenCalled();

    double vvol = 0;
    const double expectedVol = boxLength*boxLength*boxLength;

    auto startTime = std::chrono::high_resolution_clock::now();

    if (package == voro_pp) {
        const int n_x=6, n_y=6, n_z=6;
        voro::container con {-bl2, bl2, -bl2, bl2, -bl2, bl2,
                            n_x, n_y, n_z, false, false, false, 8};

        for (auto& p : particles) {
            con.put(p.id, p.position.x, p.position.y, p.position.z);
        }

        voro::voronoicell_neighbor c;
        voro::c_loop_all vl(con);
        if (vl.start()) do {
            if (con.compute_cell(c, vl)) {
                double vol = c.volume();
                vvol += vol;
                record(vl.pid(), c);
            }
        } while (vl.inc());
    }

    else /* hmc tessellate */ {
        hmc::Diagram diagram = hmc::Diagram::cube(-bl2, bl2, -bl2, bl2, -bl2, bl2);
        diagram.initialize([&](){
            for (auto& p : particles) {
                diagram.addParticle(p.id, p.position.x, p.position.y, p.position.z);
            }
        });

        hmc::Cell c;
        for (std::size_t i = 0; i < diagram.size(); ++i) {
            if (useDistances) {
                c = diagram.getCell(i, distances[i]);
            }
            else {
                c = diagram.getCell(i);
            }
            double vol = c.computeVolume();

            vvol += vol;
            record(c);
        }
    }

    outputMallocs();

    auto totalTime = std::chrono::high_resolution_clock::now() - startTime;
    std::cerr << "Completed in "
              << std::chrono::duration<double, std::milli> {totalTime}.count()
              << " ms" << std::endl;


    using dbl = std::numeric_limits<double>;
    if (vvol == expectedVol) {
        std::cerr << "Volume correct (" << vvol << ")" << std::endl;
    }
    else if (std::abs(vvol - expectedVol) < 1e-8) {
        auto oldPrecision = std::cerr.precision();
        std::cerr << std::setprecision(dbl::max_digits10)
                  << "Volume nearly correct:"
                  << std::endl
                  << "  Got " << vvol << ", expected " << expectedVol
                  << std::endl
                  << "  Off by " << (vvol - expectedVol)
                  << std::endl
                  << "    = " << (vvol - expectedVol) / dbl::epsilon()
                  << " * machine epsilon"
                  << std::endl
                  << "  (Machine epsilon: " << dbl::epsilon() << ")"
                  << std::endl
                  << std::setprecision(oldPrecision);
    }
    else {
        auto oldPrecision = std::cerr.precision();
        std::cerr << "Volume incorrect ("
                  << std::setprecision(dbl::max_digits10)
                  << vvol << " vs " << expectedVol
                  << ")" << std:: endl
                  << std::setprecision(oldPrecision);
    }

    std::cerr << std::endl;
}


template <typename Data, int verbose = 1>
struct Check {
    static constexpr char const* checkedData() {return Data::checkedData();}
    std::vector<Data> voropp;
    std::vector<Data> hmc;

    Check(std::size_t numPoints) {
        voropp.resize(numPoints);
        hmc.resize(numPoints);
    }

    void operator() (int id, voro::voronoicell_neighbor& c) {
        voropp[id] = Data {c};
    }

    void operator() (hmc::Cell& c) {
        hmc[c.getParticle().id] = Data {c};
    }

    void check(const ParticleList& particles) {
        assert(hmc.size() == voropp.size());
        bool allMatch = true;
        for (std::size_t i = 0; i < hmc.size(); ++i) {
            if (verbose == 0) {
                if (!(hmc[i] == voropp[i])) {
                    allMatch = false;
                }
            }
            else {
                if (!(hmc[i] == voropp[i])) {
                    allMatch = false;
                    std::cerr << "!! Mismatched data for particle "
                              << i << ": " << particles[i] << std::endl;

                    std::cerr << "  Voro++:" << std::endl;
                    std::cerr << "    " << voropp[i] << std::endl;

                    std::cerr << "  HMC:" << std::endl;
                    std::cerr << "    " << hmc[i] << std::endl;
                }
                else if (verbose > 1) {
                    std::cerr << "** Matching data for particle "
                              << i << ": "<< particles[i] << std::endl;

                    std::cerr << "  Voro++:" << std::endl;
                    std::cerr << "    " << voropp[i] << std::endl;

                    std::cerr << "  HMC:" << std::endl;
                    std::cerr << "    " << hmc[i] << std::endl;
                }
            }
        }

        if (allMatch) {
            std::cerr << std::endl
                      << "Success!"
                      << std::endl
                      << "Voro++ and HMC agree on " << checkedData()
                      << std::endl;
        }
        else {
            std::cerr << std::endl
                      << "!! FAILURE !!"
                      << std::endl
                      << "Voro++ and HMC do NOT agree on " << checkedData()
                      << std::endl;
        }
    }
};

template <int verbose>
struct Check<void, verbose> {
    static constexpr char const* checkedData() {return "nothing";}
    Check(std::size_t) {}
    void operator() (int, voro::voronoicell_neighbor&) {}
    void operator() (hmc::Cell&) {}
    void check(const ParticleList&) {
        if (verbose) {
            std::cerr << "No verbose output: not checking any data." << std::endl;
        }
        std::cerr << "No data checked."
                  << " There is no guarantee that the results are correct."
                  << std::endl;
    }
};

struct Volume {
    static constexpr char const* checkedData() {return "volumes";}
    static constexpr double TOLERANCE = 1e-11;
    double volume;

    Volume() = default;
    Volume(voro::voronoicell_neighbor& c) : volume(c.volume()) {}
    Volume(hmc::Cell& c) : volume(c.computeVolume()) {}

    friend bool operator==(const Volume& a, const Volume& b) {
        return std::abs(a.volume - b.volume) < Volume::TOLERANCE;
    }
    friend std::ostream& operator<<(std::ostream& out, const Volume& v) {
        return out << v.volume;
    }
};

struct Neighbors {
    static constexpr char const* checkedData() {return "neighbors";}
    std::set<int> neighbors;

    Neighbors() = default;
    Neighbors(voro::voronoicell_neighbor& c) {
        std::vector<int> n;
        c.neighbors(n);
        for (auto i : n) {
            if (i >= 0) {neighbors.insert(i);}
        }
    }
    Neighbors(hmc::Cell& c) {
        std::vector<int> n;
        c.computeNeighbors(n);
        for (auto i : n) {
            if (i >= 0) {neighbors.insert(i);}
        }
    }

    friend bool operator==(const Neighbors& a, const Neighbors& b) {
        return a.neighbors == b.neighbors;
    }
    friend std::ostream& operator<<(std::ostream& out, const Neighbors& n) {
        out << "[";
        char* sep = (char*) "";
        for (auto a : n.neighbors) {
            out << sep << a;
            sep = (char*) ", ";
        }
        out << "]";
        return out;
    }
};

struct AllData {
    static constexpr char const* checkedData() {return "volumes and neighbors";}
    Volume v;
    Neighbors n;

    AllData() = default;
    AllData(voro::voronoicell_neighbor& c) : v(c), n(c) {}
    AllData(hmc::Cell& c) : v(c), n(c) {}

    friend bool operator==(const AllData& a, const AllData& b) {
        return a.v == b.v && a.n == b.n;
    }
    friend std::ostream& operator<<(std::ostream& out, const AllData& a) {
        return out << "Volume: " << a.v << std::endl
                   << "Neighbors: " << a.n;
    }
};

























































using CheckType = Check<void>;
constexpr std::size_t DEFAULT_NUM_POINTS = 100;
constexpr bool shrinkwrap = false;
constexpr double shrinkwrapPadding = 1.00001;


























































int main(int argc, char* argv[]) {
    std::size_t numPoints = DEFAULT_NUM_POINTS;

    if (argc > 1) {
        std::size_t n = atoi(argv[1]);
        if (n != 0) { numPoints = n; }
    }

    if (argc > 2) {
        hmc::utilities::set_random_state(argv[2]);
    }
    else {
        hmc::utilities::set_random_state_true_random();
    }



    ParticleList particles;
    double boxLength = 1;

    auto runDoubleTrial = [&](std::size_t numPoints) -> void {
        std::cerr << "Generating " << numPoints << " particles." << std::endl
                  << "Random state: \"" << hmc::utilities::get_random_state()
                  << "\"" << std::endl;
        particles.clear();
        const double bl2 = boxLength/2;
        for (std::size_t i = 0; i < numPoints; ++i) {
            particles.emplace_back(i,
                hmc::utilities::uniform(-bl2, bl2),
                hmc::utilities::uniform(-bl2, bl2),
                hmc::utilities::uniform(-bl2, bl2)
            );
        }

        std::vector<double> ds;
        Check<Neighbors> ns {numPoints};

        if (shrinkwrap) {
            runTrial<voro_pp>(boxLength, particles, ns);

            for (unsigned n = 0; n < ns.voropp.size(); ++n) {
                double d = 0;
                for (auto& i : ns.voropp[n].neighbors) {
                    double d2 = distance(particles[i].position, particles[n].position);
                    if (d2 > d) { d = d2; }
                }
                ds.push_back(shrinkwrapPadding * d);
            }
        }
        // else {
        //     Check<void> q {numPoints};
        //     runTrial<voro_pp>(boxLength, particles, q);
        // }

        CheckType check {numPoints};


        runTrial<voro_pp>(boxLength, particles, check);
        runTrial<hmc_tessellate, shrinkwrap>(boxLength, particles, check, ds);

        check.check(particles);

		std::cout << "Total number of cuts: " << totalCuts << std::endl;
		std::cout << "Needed cuts: " << neededCuts << std::endl;
		std::cout << "Number of expansions: " << expansions << std::endl;
		std::cout << "Number of cells searched: " << cellsSearched << std::endl;
		std::cout << "Time spent expanding: " << std::chrono::duration<double, std::milli> {expansionTime}.count() << " ms" << std::endl;
    };

    runDoubleTrial(numPoints);
}
