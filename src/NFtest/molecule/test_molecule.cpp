#include "test_molecule.hh"
#include "../../NFfunction/NFfunction.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <sstream>

using namespace std;
using namespace NFcore;

void NFtest_molecule::run()
{
	cout << "Running NFcore::Molecule tests..." << endl;

	cout << "  Testing Molecule::setLocalFunctionValue..." << endl;

	// Need a system to create a molecule type
	System* s = new System("test");

	// Create a component
	vector<string> compNames;
	compNames.push_back("c");

	vector<string> defaultStates;
	defaultStates.push_back("s");

	vector<vector<string>> allowedStates;
	vector<string> compAllowedStates;
	compAllowedStates.push_back("s");
	compAllowedStates.push_back("p");
	allowedStates.push_back(compAllowedStates);

	// Create a dummy MoleculeType
	MoleculeType* mt = new MoleculeType("testMT", compNames, defaultStates, allowedStates, s);

	// In order for Molecule::setUpLocalFunctionList() to allocate array,
	// parentMoleculeType->getNumOfTypeIFunctions() must be > 0.
	// So we need to add a dummy local function of type I to MoleculeType.

	vector<string> args;
	args.push_back("x"); // Avoid global function check
	vector<string> varRefNames;
	vector<string> varObservableNames;
	vector<Observable *> varObservables;
	vector<int> varRefScope;
	vector<string> paramNames;

	LocalFunction* lf1 = new LocalFunction(s, "testFunc1", "1.0", "1.0", args, varRefNames, varObservableNames, varObservables, varRefScope, paramNames);
	LocalFunction* lf2 = new LocalFunction(s, "testFunc2", "2.0", "2.0", args, varRefNames, varObservableNames, varObservables, varRefScope, paramNames);

	mt->addLocalFunc_TypeI(lf1);
	mt->addLocalFunc_TypeI(lf2);

	// Create a molecule
	Molecule* m = new Molecule(mt, 0, NULL);

	m->setUpLocalFunctionList();

	// setLocalFunctionValue tests
	m->setLocalFunctionValue(5.5, 0);
	m->setLocalFunctionValue(10.2, 1);

	// getLocalFunctionValue to verify
	if (m->getLocalFunctionValue(0) != 5.5) {
		throw runtime_error("getLocalFunctionValue did not return expected value for index 0");
	}
	if (m->getLocalFunctionValue(1) != 10.2) {
		throw runtime_error("getLocalFunctionValue did not return expected value for index 1");
	}

	cout << "  Molecule::setLocalFunctionValue tests passed!" << endl;

	cout << "  Testing Molecule::printDetails..." << endl;

	ostringstream oss;
	m->printDetails(oss);
	string output = oss.str();

	if (output.find("++ Molecule instance of type: testMT") == string::npos) {
		throw runtime_error("printDetails did not print the correct molecule type name");
	}
	if (output.find("testFunc1(x)=5.5") == string::npos) {
		throw runtime_error("printDetails did not print the correct local function value 1");
	}
	if (output.find("testFunc2(x)=10.2") == string::npos) {
		throw runtime_error("printDetails did not print the correct local function value 2");
	}

	Molecule* m2 = new Molecule(mt, 0, NULL);
	m2->setUpLocalFunctionList();

	Molecule::bind(m, 0, m2, 0);

	ostringstream oss2;
	m->printDetails(oss2);
	string output2 = oss2.str();

	if (output2.find("c=s") == string::npos) {
		throw runtime_error("printDetails did not print the bonded component state");
	}
	if (output2.find("bond=testMT_") == string::npos) {
		throw runtime_error("printDetails did not print the bonded molecule type name correctly");
	}

	cout << "  Molecule::printDetails tests passed!" << endl;

	cout << "NFcore::Molecule tests completed successfully." << endl;

	// System destructor will free local functions, molecule types, and molecules instantiated.
	delete s;
}
