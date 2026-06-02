#include "test_compartment.hh"
#include "../../NFcore/compartment.hh"
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace NFcore;

void NFtest_compartment::run()
{
	cout << "Running Compartment tests..." << endl;

	cout << "  Testing Compartment Constructors..." << endl;

    Compartment* root = new Compartment("root", 3, 100);
    Compartment* child1 = new Compartment("child1", 3, 50, root);
    Compartment* child2 = new Compartment("child2", 3, 50, root);
    Compartment* grandchild1 = new Compartment("grandchild1", 3, 20, child1);

	cout << "  Testing Compartment::isInside..." << endl;

    if (root->isInside(nullptr) != false) {
        throw std::runtime_error("isInside(nullptr) did not return false");
    }

    if (!root->isInside(root)) {
        throw std::runtime_error("isInside(this) did not return true");
    }

    if (!grandchild1->isInside(child1)) {
        throw std::runtime_error("isInside(parent) did not return true");
    }

    if (!grandchild1->isInside(root)) {
        throw std::runtime_error("isInside(grandparent) did not return true");
    }

    if (child1->isInside(grandchild1)) {
        throw std::runtime_error("parent isInside(child) returned true, expected false");
    }

    if (root->isInside(grandchild1)) {
        throw std::runtime_error("grandparent isInside(grandchild) returned true, expected false");
    }

    if (child1->isInside(child2)) {
        throw std::runtime_error("sibling isInside(sibling) returned true, expected false");
    }

    delete grandchild1;
    delete child2;
    delete child1;
    delete root;

	cout << "Compartment tests completed successfully." << endl;
}
