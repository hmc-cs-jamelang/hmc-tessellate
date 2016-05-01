#include <iostream>
#include <sstream>
#include <fstream>
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
	cellsSearched = 0,
    attemptedDestroyedVertices = 0,
    destroyedVertices = 0,
    attemptedDestroyedEdges = 0,
    destroyedEdges = 0,
    destroyedMSVertices = 0,
    destroyedMSEdges = 0;

auto dont_use_this_name = std::chrono::high_resolution_clock::now();
auto expansionTime = dont_use_this_name - dont_use_this_name;

#if defined(COUNT_MALLOCS)
    #include <malloc_count.h>

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

#include "../src/utilities.hpp"
#include "../src/vectormath.hpp"
#include "../src/hmc-tessellate.hpp"

#include <voro++.hh>

struct Particle {
    int id;
    hmc::Vector3 position;

    Particle(int id, double x, double y, double z)
        : id(id), position(x, y, z)
    { /* Done */ }

    friend std::ostream& operator<<(std::ostream& out, const Particle& p)
    {
        return out << "Particle {id = " << p.id
                    << ", position = " << p.position
                    << "}";
    }
};

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
#define USE_PRE_CONTAINER
#ifdef USE_PRE_CONTAINER
        voro::pre_container preContainer(-bl2, bl2, -bl2, bl2, -bl2, bl2,
                                         false, false, false);

        for (auto& p : particles) {
            preContainer.put(p.id, p.position.x, p.position.y, p.position.z);
        }
        // grr, signed...
        std::array<int, 3> numberOfBoxes;
        // ugly
        preContainer.guess_optimal(numberOfBoxes[0], numberOfBoxes[1],
                                   numberOfBoxes[2]);

        voro::container con(-bl2, bl2, -bl2, bl2, -bl2, bl2,
                            numberOfBoxes[0], numberOfBoxes[1], numberOfBoxes[2],
                            false, false, false, 8);
        preContainer.setup(con);
#else
        const int n_x=6, n_y=6, n_z=6;
        voro::container con {-bl2, bl2, -bl2, bl2, -bl2, bl2,
                            n_x, n_y, n_z, false, false, false, 8};

        for (auto& p : particles) {
            con.put(p.id, p.position.x, p.position.y, p.position.z);
        }
#endif

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
                diagram.addParticle(p.position.x, p.position.y, p.position.z);
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
    std::cerr << "##############################" << std::endl
		      << "Completed in " << std::endl
              << "        " << std::chrono::duration<double, std::milli> {totalTime}.count()
              << " ms" << std::endl
			  << "##############################" << std::endl;


    using dbl = std::numeric_limits<double>;
    if (vvol == expectedVol) {
        std::cerr << "Volume correct (" << vvol << ")" << std::endl;
    }
    else if (std::abs((vvol - expectedVol)/expectedVol) < 1e-8) {
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
        hmc[c.getOriginalIndex()] = Data {c};
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
        auto oldPrecision = out.precision();
        return out << std::setprecision(std::numeric_limits<double>::max_digits10)
                   << v.volume
                   << std::setprecision(oldPrecision);
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
constexpr int DEFAULT_NUM_POINTS = 100;
constexpr bool shrinkwrap = false;
constexpr double shrinkwrapPadding = 1 + 1e-5;















































int main(int argc, char* argv[]) {
    int np = DEFAULT_NUM_POINTS;

    if (argc > 1) {
        std::size_t n = atoi(argv[1]);
        if (n != 0) { np = n; }
    }

    if (argc > 2) {
        hmc::utilities::set_random_state(argv[2]);
    }
    else {
        hmc::utilities::set_random_state_true_random();
    }



    ParticleList particles;
    double boxLength = 1;

    if (np == -1) {
        boxLength = 2*57.8;
    }

    auto runDoubleTrial = [&](int np) -> void {
        if (np == -1) {
            int id;
            double x, y, z;

            std::cerr << "Reading points from 'particles'" << std::endl;
            std::ifstream inputFile;
            inputFile.open("particles");

            while (inputFile >> id >> x >> y >> z) {
                particles.emplace_back(id, x, y, z);
            }

            inputFile.close();
        }

        else {
            std::cerr << "Generating " << np << " particles." << std::endl
                      << "Random state: \"" << hmc::utilities::get_random_state()
                      << "\"" << std::endl;
            particles.clear();
            const double bl2 = boxLength/2;
            for (int i = 0; i < np; ++i) {
                particles.emplace_back(i,
                    hmc::utilities::uniform(-bl2, bl2),
                    hmc::utilities::uniform(-bl2, bl2),
                    hmc::utilities::uniform(-bl2, bl2)
                );
            }
        }

        std::size_t numPoints = particles.size();
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
    };

    runDoubleTrial(np);
}
