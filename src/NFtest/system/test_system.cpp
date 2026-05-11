#include "test_system.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

using namespace std;
using namespace NFcore;

void NFtest_system::run()
{
	cout << "Running System tests..." << endl;

	cout << "  Testing System::addMoleculeType..." << endl;

	// Instantiate a System
	System* sys = new System("TestSystem");

	// Initial checks
	if (sys->getNumOfMoleculeTypes() != 0) {
		throw std::runtime_error("System did not initialize with 0 MoleculeTypes.");
	}

	// Create first MoleculeType
	vector<string> compNameA;
	compNameA.push_back("x");
	MoleculeType* mtA = new MoleculeType("MolA", compNameA, sys);

	// Verify first MoleculeType addition
	if (mtA->getTypeID() != 0) {
		throw std::runtime_error("First MoleculeType type ID is not 0.");
	}
	if (sys->getNumOfMoleculeTypes() != 1) {
		throw std::runtime_error("System number of MoleculeTypes is not 1 after adding first MoleculeType.");
	}
	if (sys->getMoleculeType(0) != mtA) {
		throw std::runtime_error("System::getMoleculeType(0) did not return the correct MoleculeType pointer.");
	}

	// Create second MoleculeType
	vector<string> compNameB;
	compNameB.push_back("y");
	compNameB.push_back("z");
	MoleculeType* mtB = new MoleculeType("MolB", compNameB, sys);

	// Verify second MoleculeType addition
	if (mtB->getTypeID() != 1) {
		throw std::runtime_error("Second MoleculeType type ID is not 1.");
	}
	if (sys->getNumOfMoleculeTypes() != 2) {
		throw std::runtime_error("System number of MoleculeTypes is not 2 after adding second MoleculeType.");
	}
	if (sys->getMoleculeType(1) != mtB) {
		throw std::runtime_error("System::getMoleculeType(1) did not return the correct MoleculeType pointer.");
	}

	// Verify that the first MoleculeType is still accessible correctly
	if (sys->getMoleculeType(0) != mtA) {
		throw std::runtime_error("System::getMoleculeType(0) did not return the correct MoleculeType pointer after adding second MoleculeType.");
	}

	cout << "  System::addMoleculeType tests passed!" << endl;

	cout << "  Testing System::getMoleculeTypeByName..." << endl;

	// Test happy path
	if (sys->getMoleculeTypeByName("MolA") != mtA) {
		throw std::runtime_error("System::getMoleculeTypeByName did not return correct MoleculeType for 'MolA'.");
	}
	if (sys->getMoleculeTypeByName("MolB") != mtB) {
		throw std::runtime_error("System::getMoleculeTypeByName did not return correct MoleculeType for 'MolB'.");
	}

	// Test error path
	bool caughtError = false;
	try {
		sys->getMoleculeTypeByName("NonExistentMol");
	} catch (const std::runtime_error& e) {
		caughtError = true;
	}

	if (!caughtError) {
		throw std::runtime_error("System::getMoleculeTypeByName did not throw std::runtime_error for non-existent molecule type.");
	}

	cout << "  System::getMoleculeTypeByName tests passed!" << endl;

	// Cleanup
	// Note: System destructor deletes MoleculeTypes added to it
	delete sys;

	cout << "System tests completed successfully." << endl;
}
