struct FaceZipper {
    Polyhedron& poly;
    VertexIndex backwardVertex;
    EdgeIndex forwardEdge;
    EdgeIndex backwardEdge;

    FaceZipper(Polyhedron& poly, VertexIndex startVertex)
        : poly(poly), backwardVertex(vi), forwardEdge(), backwardEdge() {}

    inline void zipForward(VertexIndex newVertex)
    {
        EdgePointer newForwardEdge = poly.createEdge(newVertex);
        poly.edgePointer(forwardEdge)->next = newForwardEdge;
        forwardEdge = newForwardEdge;
    }

    inline void zipBackward(VertexIndex newVertex)
    {
        backwardEdge = poly.createEdge(backwardVertex, backwardEdge);
        backwardVertex = newVertex;
    }

    inline void zip(FaceZipper& backward, VertexIndex newVertex)
    {
        zipForward(newVertex);
        backward.zipBackward(newVertex);
        poly.edgePointer(edge)->setFlips(flip);
    }

    inline void done()
    {
        assert(poly.edgePointer(forwardEdge)->target == backwardVertex);
        poly.edgePointer(forwardEdge)->next = backwardEdge;
    }
};

// inline void zipFaces(FaceZipper& forward, FaceZipper& backward, const VertexIndex vertex)
// {
//     setFlips(
//         zipNewHalfEdge_Forward(forward, vertex),
//         zipNewHalfEdge_Backward(backward, vertex)
//     );
// }

// void createCube()
// {
//     FaceZipper down {LDF};
//     FaceZipper front = zipNewFace(down, LDF, RDF);
//     FaceZipper right = zipNewFace(down, RDF, RDB);
//     FaceZipper back = zipNewFace(down, RDB, LDB);
//     down.done();

//     zipFaces(right, front, RDF, RTF);
//     zipFaces(back, right, RDB, RTB);
//     zipFaces(left, back, LDB, LTB);
//     zipFaces(front, left, LDF, LTF);

//     FaceZipper top = zipNewFace(front, RTF);
//     front.done();

//     zipFaces(top, left, LTB);
//     zipFaces()
// }

// contract(class) : Class Invariants
// contract(public_method) : Checks class invariants
// contract(private_method) : Does not check class invariants
// contract(constructor) : Checks class invariants at the beginning


// If I'm going to use functions, then, might as well make them
// part of the polyhedron ...
//
// Sounds good. That's the way to do it. I should have paid attention
// to younger me sooner.


