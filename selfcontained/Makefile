CC_INCLUDE =
LD_FLAGS =
FLAGS =
CC_FLAGS += -O3 -Wall -Wextra -std=c++11 $(FLAGS)
VORO_INCLUDE = -Wno-unused-parameter -I../../../source/voro++-0.4.6/include/voro++ -L../../../source/voro++-0.4.6/lib -lvoro++
COUNT_MALLOCS = -DCOUNT_MALLOCS malloc_count.o -ldl
#CC = clang++-3.5
CC = g++
HMC_DEP = sds.hpp hmc-tessellate.hpp polyhedron.hpp structpool.hpp verification.hpp wrapper.hpp vectormath.hpp cellarray.hpp cellarray-private.hpp

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
	$(CC) malloc_count-0.7/malloc_count.c -o malloc_count.o -c

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
