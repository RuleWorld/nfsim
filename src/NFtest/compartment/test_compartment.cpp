#include "test_compartment.hh"
#include <iostream>
#include <stdexcept>
#include "../../NFcore/compartment.hh"

using namespace std;
using namespace NFcore;

void NFtest_compartment::run()
{
    cout << "Running Compartment tests..." << endl;

    Compartment compA("A", 3, 10.0);
    Compartment compB("B", 3, 5.0, &compA);
    Compartment compC("C", 3, 2.0, &compB);
    Compartment compD("D", 3, 5.0, &compA); // Sibling of B

    cout << "  Testing isInside direct nesting..." << endl;
    if (!compB.isInside(&compA)) {
        throw std::runtime_error("isInside failed: B should be inside A");
    }
    if (!compC.isInside(&compB)) {
        throw std::runtime_error("isInside failed: C should be inside B");
    }

    cout << "  Testing isInside transitive nesting..." << endl;
    if (!compC.isInside(&compA)) {
        throw std::runtime_error("isInside failed: C should be inside A (transitively)");
    }

    cout << "  Testing isInside non-nesting / reverse nesting..." << endl;
    if (compA.isInside(&compB)) {
        throw std::runtime_error("isInside failed: A should NOT be inside B");
    }
    if (compB.isInside(&compC)) {
        throw std::runtime_error("isInside failed: B should NOT be inside C");
    }

    cout << "  Testing isInside sibling..." << endl;
    if (compB.isInside(&compD)) {
        throw std::runtime_error("isInside failed: B should NOT be inside D (sibling)");
    }
    if (compD.isInside(&compB)) {
        throw std::runtime_error("isInside failed: D should NOT be inside B (sibling)");
    }

    cout << "  Testing isInside NULL..." << endl;
    if (compA.isInside(NULL)) {
        throw std::runtime_error("isInside failed: A should NOT be inside NULL");
    }

    cout << "  Compartment tests passed!" << endl;
}
