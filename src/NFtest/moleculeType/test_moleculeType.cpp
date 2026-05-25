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

	vector<bool> isIntegerComponent;
	isIntegerComponent.push_back(false);
	isIntegerComponent.push_back(true);

	// Create a MoleculeType
	MoleculeType* mt = new MoleculeType("testMT", compNames, defaultStates, allowedStates, isIntegerComponent, s);

	// isIntegerComponent tests
	cout << "  Testing MoleculeType::isIntegerComponent..." << endl;
	if (mt->isIntegerComponent("siteA") != false) {
		throw runtime_error("isIntegerComponent did not return false for 'siteA'");
	}
	if (mt->isIntegerComponent("siteB") != true) {
		throw runtime_error("isIntegerComponent did not return true for 'siteB'");
	}
	if (mt->isIntegerComponent(0) != false) {
		throw runtime_error("isIntegerComponent did not return false for index 0");
	}
	if (mt->isIntegerComponent(1) != true) {
		throw runtime_error("isIntegerComponent did not return true for index 1");
	}
	bool threwIntComp = false;
	try {
		mt->isIntegerComponent("invalidSite");
	} catch (const std::runtime_error& e) {
		threwIntComp = true;
	}
	if (!threwIntComp) {
		throw runtime_error("isIntegerComponent did not throw runtime_error for invalid site name 'invalidSite'");
	}
	threwIntComp = false;
	try {
		mt->isIntegerComponent(99);
	} catch (const std::runtime_error& e) {
		threwIntComp = true;
	}
	if (!threwIntComp) {
		throw runtime_error("isIntegerComponent did not throw runtime_error for invalid site index 99");
	}

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


	cout << "  Testing MoleculeType::isEquivalentComponent..." << endl;

	vector<vector<string>> eqComps;
	vector<string> eqCompGroup;
	eqCompGroup.push_back("siteA");
	eqCompGroup.push_back("siteB");
	eqComps.push_back(eqCompGroup);

	mt->addEquivalentComponents(eqComps);

	if (!mt->isEquivalentComponent("siteA")) {
		throw runtime_error("isEquivalentComponent(string) did not return true for 'siteA'");
	}

	if (!mt->isEquivalentComponent("siteB")) {
		throw runtime_error("isEquivalentComponent(string) did not return true for 'siteB'");
	}

	bool threw_equiv = false;
	try {
		mt->isEquivalentComponent("invalidSite");
	} catch (const std::runtime_error& e) {
		threw_equiv = true;
	}
	if (!threw_equiv) {
		throw runtime_error("isEquivalentComponent(string) did not throw for 'invalidSite'");
	}

	if (!mt->isEquivalentComponent(mt->getCompIndexFromName("siteA"))) {
		throw runtime_error("isEquivalentComponent(int) did not return true for index of 'siteA'");
	}

	if (!mt->isEquivalentComponent(mt->getCompIndexFromName("siteB"))) {
		throw runtime_error("isEquivalentComponent(int) did not return true for index of 'siteB'");
	}

	if (mt->isEquivalentComponent(-1)) {
		throw runtime_error("isEquivalentComponent(int) returned true for index -1");
	}

	if (mt->isEquivalentComponent(999)) {
		throw runtime_error("isEquivalentComponent(int) returned true for index 999");
	}

	cout << "  MoleculeType::isEquivalentComponent tests passed!" << endl;

	cout << "NFcore::MoleculeType tests completed successfully." << endl;

	// System destructor will free molecule types instantiated.
	delete s;
}
