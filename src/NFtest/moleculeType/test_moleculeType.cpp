#include "test_moleculeType.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

using namespace std;
using namespace NFcore;

void NFtest_moleculeType::run()
{
	cout << "Running NFcore::MoleculeType tests..." << endl;

	cout << "  Testing MoleculeType::getCompIndexFromName..." << endl;

	// Need a system to create a molecule type
	System* s = new System("testSys");

	// Create a component
	vector<string> compNames;
	compNames.push_back("siteA");
	compNames.push_back("siteB");

	vector<string> defaultStates;
	defaultStates.push_back("u");
	defaultStates.push_back("unbound");

	vector<vector<string>> allowedStates;
	vector<string> compAllowedStatesA;
	compAllowedStatesA.push_back("u");
	compAllowedStatesA.push_back("p");
	allowedStates.push_back(compAllowedStatesA);

	vector<string> compAllowedStatesB;
	compAllowedStatesB.push_back("unbound");
	compAllowedStatesB.push_back("bound");
	allowedStates.push_back(compAllowedStatesB);

	// Create a MoleculeType
	MoleculeType* mt = new MoleculeType("testMT", compNames, defaultStates, allowedStates, s);

	// Happy path tests
	int indexA = mt->getCompIndexFromName("siteA");
	if (indexA != 0) {
		throw runtime_error("getCompIndexFromName did not return 0 for 'siteA', returned: " + to_string(indexA));
	}

	int indexB = mt->getCompIndexFromName("siteB");
	if (indexB != 1) {
		throw runtime_error("getCompIndexFromName did not return 1 for 'siteB', returned: " + to_string(indexB));
	}

	// Exception path test
	bool threw = false;
	try {
		mt->getCompIndexFromName("invalidSite");
	} catch (const std::runtime_error& e) {
		threw = true;
	}

	if (!threw) {
		throw runtime_error("getCompIndexFromName did not throw runtime_error for invalid site name 'invalidSite'");
	}

	cout << "  MoleculeType::getCompIndexFromName tests passed!" << endl;
	cout << "NFcore::MoleculeType tests completed successfully." << endl;

	// System destructor will free molecule types instantiated.
	delete s;
}
