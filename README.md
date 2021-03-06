# hmc-tessellate

Given a set of particles, a Voronoi tessellation divides a space into Voronoi cells. Each Voronoi
cell contains exactly one particle and all of the space closer to that particle than any other.
`hmc-tessellate` aims to compute 3D Voronoi tessellations.

![alt text](https://github.com/hmc-cs-jamelang/hmc-tessellate/blob/master/images/voronoithing.png)

## Initializing a Diagram

To create a Voronoi diagram, first you fill it with particles
by using the `put` function. You only _need_ to enter particles
which will be used as targets at some point.

After you are done entering particles, you must call
`initializeDiagram()`, which prepares the diagram for use.


###   Particles

Each particle has an x, y, z position, as well as an id.
Optionally, you may provide a group for the

The id is simply an arbitrary integer tag for your own bookkeeping.
There are no requirements on the id, uniqueness or otherwise,
and you are free to just make all the ids 0 or -1 if you'd like.

Groups are also arbitrary integers, but they must be distinct
from each other, and a group cannot be -1.
-1 is reserved to mean "no group."


###   Examples

####       Without Groups
```cpp
         Diagram diagram;
         for (std::size_t i = 0; i < numPoints; ++i) {
             diagram.put(xs[i], ys[i], zs[i], ids[i]);
         }
         diagram.initializeDiagram();
```

####       With Groups
```cpp
         int Red = 2;
         int Blue = 235;
         int Green = 17;

         Diagram diagram;
         for (std::size_t i = 0; i < numPoints; ++i) {
             diagram.put(xs[i], ys[i], zs[i], ids[i], colors[i]);
         }
         diagram.intializeDiagram();
```

## Constructing Voronoi Cells

Once the diagram has been initialized, you can construct Voronoi
cells and then query their properties. By doing each cell on its
own, you can reuse the same Cell object and avoid multiple expensive
memory allocations.


###   Groups

As far as source and target groups: Conceptually, you iterate over
the cells in the source group, and construct them using the particles
in the target group.

Thus `diagram.sourceGroups(...)` creates an object
that you can iterate over, while `diagram.targetGroups(...)` just
creates a simple flag specifying the groups.

Ideally, you should store the source and target groups in a variable
rather than reconstructing them each time through the loop.
This is especially important for source groups, which can actually
be expensive to create.


###   Virtual Source Groups

Since only the target group and the position of the
current cell is used for the actual construction, the particles
in the source group don't actually have to be inside the diagram.

In order to iterate over a source group that isn't actually contained
in the diagram, you simply need to provide an x, y, z position to
the `getCell` function rather than an index. You can find
an example of virtual source groups in `Examples` below.


###   Examples

####       Without Groups

We construct a Cell object outside the for loop.
Then we iterate over all the particles in the diagram
and create a cell for each one in turn.
```cpp
             Cell c;
             for (std::size_t i = 0; i < diagram.size(); ++i) {
                 c = diagram.getCell(i);

                 double volume = c.computeVolume();

                 std::vector<std::size_t> neighbors = c.computeNeighbors();
                 // Or if you already have a vector you want to be filled:
                 //     c.computeNeighbors(std::back_inserter(neighbors));

                 // ... do some things with the data ...
             }
```

####       With Target Groups

We create a target group object outside of the loop.
When we create each cell using `getCell`, we also pass
the target group object in.
```cpp
             auto targetGroups = diagram.targetGroups(Red, Blue);

             Cell c;
             for (std::size_t i = 0; i < diagram.size(); ++i) {
                 c = diagram.getCell(i, targetGroups);

                 // Using the cell is exactly the same otherwise
                 double volume = c.computeVolume();

                 // etc.
             }
```

####       With Source Groups

We create a source group object, and iterate over it instead
of over the whole diagram.
```cpp
             auto sourceGroups = diagram.sourceGroups(Green, Blue);

             Cell c;
             for (std::size_t i = 0; i < sourceGroups.size(); ++i) {
                 c = sourceGroups.getCell(i);

                 // Using the cell is exactly the same otherwise
             }
```
If you are looking at the ith particle in the source group,
and you want to know what index that would be inside the
overall diagram, you can do that with
```cpp
             sourceGroups.getTrueIndex(i);
```


####       With Source Groups and Target Groups

Both techniques can of course be used together.
```cpp
             auto sourceGroups = diagram.sourceGroups(Green, Blue);
             auto targetGroups = diagram.targetGroups(Red, Blue);

             Cell c;
             for (std::size_t i = 0; i < sourceGroups.size(); ++i) {
                 c = sourceGroups.getCell(i, targetGroups);

                 // Using the cell is exactly the same otherwise
             }
```

####       With "Virtual" Source Groups

Instead of iterating over indices in a source group,
you can also directly supply a position to `getCell`.
This allows you to construct cells for particles that
aren't actually contained in the diagram.

Supposing that you have some positions stored in the
arrays xs, ys and zs:
```cpp
             Cell c;
             for (std::size_t i = 0; i < numPoints; ++i) {
                 c = diagram.getCell(xs[i], ys[i], zs[i]);

                 // Using the cell is exactly the same otherwise
             }
```
Of course, this can be combined with target groups
by using
```cpp
             c = diagram.getCell(x, y, z, targetGroups);
```


####       Parallelization with OMP

All of the for loops in the previous examples,
including those iterating over source groups,
are designed so that they can be easily parallelized by OMP
by adding a single line:
```cpp
             Cell c;
             #pragma omp parallel for private(c) // <-- New
             for (std::size_t i = 0; i < diagram.size(); ++i) {
                 // etc.
             }
```
