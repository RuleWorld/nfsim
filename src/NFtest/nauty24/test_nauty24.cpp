#include "test_nauty24.hh"
#include "../../nauty24/nausparse.h"
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;

void testNautyBasicSparse() {
    cout << "  Testing nauty basic sparse graph..." << endl;

    // Define a simple square graph
    int nv = 4;
    int nde = 8; // directed edges: (0,1), (1,0), (1,2), (2,1), (2,3), (3,2), (3,0), (0,3)
    int m = (nv + WORDSIZE - 1) / WORDSIZE;

    nauty_check(WORDSIZE, m, nv, NAUTYVERSIONID);

    SG_DECL(sg);
    SG_DECL(cg);
    SG_ALLOC(sg, nv, nde, "malloc");
    sg.nv = nv;
    sg.nde = nde;

    // Define edges
    sg.v[0] = 0; sg.d[0] = 2; sg.e[0] = 1; sg.e[1] = 3;
    sg.v[1] = 2; sg.d[1] = 2; sg.e[2] = 0; sg.e[3] = 2;
    sg.v[2] = 4; sg.d[2] = 2; sg.e[4] = 1; sg.e[5] = 3;
    sg.v[3] = 6; sg.d[3] = 2; sg.e[6] = 2; sg.e[7] = 0;

    int* lab = new int[nv];
    int* ptn = new int[nv];
    int* orbits = new int[nv];
    setword* workspace = new setword[10 * m];

    // Initial partitions
    for (int i = 0; i < nv; i++) {
        lab[i] = i;
        ptn[i] = 1;
    }
    ptn[nv - 1] = 0; // last value must be 0

    static DEFAULTOPTIONS_SPARSEGRAPH(options);
    statsblk stats;

    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    options.digraph = FALSE;

    nauty((graph*)&sg, lab, ptn, NULL, orbits, &options, &stats,
          workspace, 10 * m, m, nv, (graph*)&cg);

    // Verify some expected canonical property, for instance all nodes are in same orbit
    if (orbits[0] != orbits[1] || orbits[1] != orbits[2] || orbits[2] != orbits[3]) {
        throw std::runtime_error("testNautyBasicSparse failed: Expected all nodes to be in the same orbit.");
    }

    SG_FREE(sg);
    SG_FREE(cg);
    delete[] lab;
    delete[] ptn;
    delete[] orbits;
    delete[] workspace;
}

void NFtest_nauty24::run() {
    cout << "Running nauty24 tests..." << endl;
    testNautyBasicSparse();
    cout << "  Passed." << endl;
}
