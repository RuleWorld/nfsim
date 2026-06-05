#include "test_compartment.hh"
#include "../../NFcore/compartment.hh"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

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

    // Test false return paths
    if (root->isInside(nullptr) != false) {
        throw std::runtime_error("isInside(nullptr) did not return false");
    }

    // Test early return for identity
    if (!root->isInside(root)) {
        throw std::runtime_error("isInside(this) did not return true");
    }

    // Test early return for identity on a non-root compartment
    if (!grandchild1->isInside(grandchild1)) {
        throw std::runtime_error("isInside(this) did not return true");
    }

    // Test pointer traversal logic (is inside parent)
    if (!grandchild1->isInside(child1)) {
        throw std::runtime_error("isInside(parent) did not return true");
    }

    // Test pointer traversal logic (is inside grandparent)
    if (!grandchild1->isInside(root)) {
        throw std::runtime_error("isInside(grandparent) did not return true");
    }

    // Test false return paths: passing a child to check if parent is inside
    if (child1->isInside(grandchild1)) {
        throw std::runtime_error("parent isInside(child) returned true, expected false");
    }

    if (root->isInside(grandchild1)) {
        throw std::runtime_error("grandparent isInside(grandchild) returned true, expected false");
    }

    // Test false return paths: checking siblings
    if (child1->isInside(child2)) {
        throw std::runtime_error("sibling isInside(sibling) returned true, expected false");
    }

    delete grandchild1;
    delete child2;
    delete child1;
    delete root;

	cout << "  Compartment::isInside tests passed!" << endl;

	// Test printDetails with parent
	{
		Compartment parent("cytoplasm", 3, 100.0);
		Compartment child("nucleus", 3, 20.0, &parent);

		// Redirect cout to a stringstream
		stringstream buffer;
		streambuf* old_cout = cout.rdbuf(buffer.rdbuf());

		child.printDetails();

		// Restore cout
		cout.rdbuf(old_cout);

		string expected = "Compartment 'nucleus': 3D, size=20, parent=cytoplasm\n";
		if (buffer.str() != expected) {
			throw runtime_error("Compartment::printDetails() did not match expected output with parent.");
		}
	}

	// Test printDetails without parent
	{
		Compartment c("membrane", 2, 50.0);

		// Redirect cout to a stringstream
		stringstream buffer;
		streambuf* old_cout = cout.rdbuf(buffer.rdbuf());

		c.printDetails();

		// Restore cout
		cout.rdbuf(old_cout);

		string expected = "Compartment 'membrane': 2D, size=50\n";
		if (buffer.str() != expected) {
			throw runtime_error("Compartment::printDetails() did not match expected output without parent.");
		}
	}

	cout << "  Compartment::printDetails tests passed!" << endl;
	cout << "Compartment tests completed successfully." << endl;
}
