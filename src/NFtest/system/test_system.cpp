#include "test_system.hh"
#include "../../NFreactions/reactions/reaction.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

using namespace std;
using namespace NFcore;

void NFtest_system::run()
{
	cout << "Running System tests..." << endl;

	cout << "  Testing System Constructors..." << endl;

	// Test System(string name)
	System* sys1 = new System("TestSystem1");
	if (sys1->getName() != "TestSystem1") {
		throw std::runtime_error("System(string) did not initialize name correctly.");
	}
	if (sys1->isUsingComplex() != false) {
		throw std::runtime_error("System(string) did not initialize useComplex to false.");
	}
	if (sys1->getGlobalMoleculeLimit() != 100000) {
		throw std::runtime_error("System(string) did not initialize globalMoleculeLimit to 100000.");
	}
	delete sys1;

	// Test System(string name, bool useComplex)
	System* sys2 = new System("TestSystem2", true);
	if (sys2->getName() != "TestSystem2") {
		throw std::runtime_error("System(string, bool) did not initialize name correctly.");
	}
	if (sys2->isUsingComplex() != true) {
		throw std::runtime_error("System(string, bool) did not initialize useComplex correctly.");
	}
	if (sys2->getGlobalMoleculeLimit() != 100000) {
		throw std::runtime_error("System(string, bool) did not initialize globalMoleculeLimit to 100000.");
	}
	delete sys2;

	// Test System(string name, bool useComplex, int globalMoleculeLimit)
	System* sys3 = new System("TestSystem3", false, 500);
	if (sys3->getName() != "TestSystem3") {
		throw std::runtime_error("System(string, bool, int) did not initialize name correctly.");
	}
	if (sys3->isUsingComplex() != false) {
		throw std::runtime_error("System(string, bool, int) did not initialize useComplex correctly.");
	}
	if (sys3->getGlobalMoleculeLimit() != 500) {
		throw std::runtime_error("System(string, bool, int) did not initialize globalMoleculeLimit correctly.");
	}
	delete sys3;

	// Test edge cases for globalMoleculeLimit
	System* sys4 = new System("TestSystem4", true, 0);
	if (sys4->getGlobalMoleculeLimit() != 0) {
		throw std::runtime_error("System(string, bool, int) did not initialize globalMoleculeLimit to 0.");
	}
	delete sys4;

	System* sys5 = new System("TestSystem5", false, -100);
	if (sys5->getGlobalMoleculeLimit() != -100) {
		throw std::runtime_error("System(string, bool, int) did not initialize globalMoleculeLimit to negative correctly.");
	}
	delete sys5;

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

	cout << "  Testing System::getReactionByName..." << endl;

	// Create an empty TransformationSet for dummy reactions
	vector<TemplateMolecule*> emptyTemplates;
	TransformationSet* ts1 = new TransformationSet(emptyTemplates);
	ts1->finalize();

	TransformationSet* ts2 = new TransformationSet(emptyTemplates);
	ts2->finalize();

	ReactionClass* rxn1 = new BasicRxnClass("Rxn1", 1.0, "", ts1, sys);
	ReactionClass* rxn2 = new BasicRxnClass("Rxn2", 2.0, "", ts2, sys);

	sys->addReaction(rxn1);
	sys->addReaction(rxn2);

	// Verify we can find reactions by exact name
	if (sys->getReactionByName("Rxn1") != rxn1) {
		throw std::runtime_error("System::getReactionByName did not return correct reaction for 'Rxn1'.");
	}

	if (sys->getReactionByName("Rxn2") != rxn2) {
		throw std::runtime_error("System::getReactionByName did not return correct reaction for 'Rxn2'.");
	}

	// Verify we get NULL for nonexistent reactions
	if (sys->getReactionByName("Rxn3") != NULL) {
		throw std::runtime_error("System::getReactionByName did not return NULL for nonexistent reaction 'Rxn3'.");
	}

	if (sys->getReactionByName("") != NULL) {
		throw std::runtime_error("System::getReactionByName did not return NULL for empty string name.");
	}

	cout << "  System::getReactionByName tests passed!" << endl;

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
	// It also cleans up reactions added to it.
	delete sys;

	cout << "System tests completed successfully." << endl;
}
