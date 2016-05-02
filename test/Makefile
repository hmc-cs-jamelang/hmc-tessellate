CC_INCLUDE =
LD_FLAGS =
FLAGS =
CC_FLAGS += -O3 -Wall -Wextra -std=c++11 $(FLAGS)
VORO_INCLUDE = -I../../../source/voro++-0.4.6/include/voro++ -L../../../source/voro++-0.4.6/lib -lvoro++ -Wno-unused-parameter
#  VORO_INCLUDE = -I../../voro++-0.4.6/build/include/voro++ -L../../voro++-0.4.6/build/lib -lvoro++ -Wno-unused-parameter
MALLOC_COUNT_DIRECTORY = ../../malloc_count-0.7
COUNT_MALLOCS = -I$(MALLOC_COUNT_DIRECTORY) -DCOUNT_MALLOCS malloc_count.o -ldl
#CC = clang++-3.5
CC = g++
HMC_DEP = ../src/sds.hpp ../src/hmc-tessellate.hpp ../src/polyhedron.hpp ../src/structpool.hpp ../src/verification.hpp ../src/wrapper.hpp ../src/vectormath.hpp ../src/cellarray.hpp ../src/cellarray-private.hpp

TARGETS = benchmark mallocs
SUBTARGETS = celerytest

# Number of points used in "make test"
# If left alone, defaults to the default in
# benchmark.cpp.
# Change this by running
#     make test n=...
n = 0

# Similarly, the random seed used in "make test"
s =

.PHONY: all clean again test

all: $(TARGETS)

malloc_count.o:
	$(CC) $(MALLOC_COUNT_DIRECTORY)/malloc_count.c -o malloc_count.o -c

benchmark: benchmark.cpp $(HMC_DEP)
	$(CC) $< -o $@ $(CC_FLAGS) $(CC_INCLUDE) $(LD_FLAGS) $(VORO_INCLUDE)

mallocs: benchmark.cpp malloc_count.o $(HMC_DEP)
	$(CC) $< -o $@ $(CC_FLAGS) $(CC_INCLUDE) $(LD_FLAGS) $(VORO_INCLUDE) $(COUNT_MALLOCS)

celerytest: celerytest.cpp cellarray.hpp cellarray-private.hpp
	$(CC) $< -o $@ $(CC_FLAGS) $(CC_INCLUDE) $(LD_FLAGS)

diagramtest: diagram-test.cpp $(HMC_DEP)
	$(CC) $< -o $@ -O3 -std=c++11 $(CC_INCLUDE) $(LD_FLAGS)
	./diagramtest

test: benchmark
	./benchmark $(n) $(s)

clean:
	rm -f $(TARGETS) $(SUBTARGETS)

again: clean $(TARGETS)
