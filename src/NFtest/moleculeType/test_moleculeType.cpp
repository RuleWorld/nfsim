#include "test_moleculeType.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <sstream>

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

	cout << "  Testing MoleculeType::printDetails..." << endl;

	// Create an integer component
	vector<string> compNames2;
	compNames2.push_back("siteInt");
	compNames2.push_back("siteA");

	vector<string> defaultStates2;
	defaultStates2.push_back("0");
	defaultStates2.push_back("u");

	vector<vector<string>> allowedStates2;
	vector<string> compAllowedStatesInt;
	for (int i=0; i<=5; i++) compAllowedStatesInt.push_back(to_string(i));
	allowedStates2.push_back(compAllowedStatesInt);

	vector<string> compAllowedStatesA2;
	compAllowedStatesA2.push_back("u");
	compAllowedStatesA2.push_back("p");
	allowedStates2.push_back(compAllowedStatesA2);

	vector<bool> isIntegerComponent2;
	isIntegerComponent2.push_back(true);
	isIntegerComponent2.push_back(false);

	MoleculeType* mt2 = new MoleculeType("testMT_print", compNames2, defaultStates2, allowedStates2, isIntegerComponent2, s);

	// Redirect cout
	stringstream buffer;
	streambuf* old = cout.rdbuf(buffer.rdbuf());

	mt2->printDetails();

	// Restore cout
	cout.rdbuf(old);

	string output = buffer.str();

	if (output.find("Molecule Type: testMT_print type ID: " + to_string(mt2->getTypeID())) == string::npos) {
		throw runtime_error("printDetails did not print correct Molecule Type and ID. Output:\n" + output);
	}
	if (output.find("siteInt~integer[0-5]") == string::npos) {
		throw runtime_error("printDetails did not print correct integer component details. Output:\n" + output);
	}
	if (output.find("siteA~u~p") == string::npos) {
		throw runtime_error("printDetails did not print correct standard component details. Output:\n" + output);
	}

	cout << "  MoleculeType::printDetails tests passed!" << endl;

	cout << "  Testing MoleculeType::addEquivalentComponents..." << endl;

	vector<string> compNames3;
	compNames3.push_back("site1");
	compNames3.push_back("site2");
	compNames3.push_back("site3");
	compNames3.push_back("otherSite");

	vector<string> defaultStates3;
	defaultStates3.push_back("u");
	defaultStates3.push_back("u");
	defaultStates3.push_back("u");
	defaultStates3.push_back("u");

	vector<vector<string>> allowedStates3;
	vector<string> compAllowedStates3;
	compAllowedStates3.push_back("u");
	compAllowedStates3.push_back("p");
	allowedStates3.push_back(compAllowedStates3);
	allowedStates3.push_back(compAllowedStates3);
	allowedStates3.push_back(compAllowedStates3);
	allowedStates3.push_back(compAllowedStates3);

	MoleculeType* mt3 = new MoleculeType("testMT_eq", compNames3, defaultStates3, allowedStates3, s);

	vector<vector<string>> identicalComponents;
	vector<string> eqGroup;
	eqGroup.push_back("site1");
	eqGroup.push_back("site2");
	eqGroup.push_back("site3");
	identicalComponents.push_back(eqGroup);

	mt3->addEquivalentComponents(identicalComponents);

	if (mt3->getNumOfEquivalencyClasses() != 1) {
		throw runtime_error("addEquivalentComponents did not set the correct number of equivalency classes. Expected 1, got " + to_string(mt3->getNumOfEquivalencyClasses()));
	}

	if (mt3->getEquivalencyClassCompNames()[0] != "site") {
		throw runtime_error("addEquivalentComponents did not set the correct generic component name. Expected 'site', got '" + mt3->getEquivalencyClassCompNames()[0] + "'");
	}

	if (mt3->getEquivalencyClassNumber("site") != 0) {
		throw runtime_error("getEquivalencyClassNumber('site') returned " + to_string(mt3->getEquivalencyClassNumber("site")) + " instead of 0");
	}

	if (mt3->getEquivalenceClassNumber(0) != 0 ||
	    mt3->getEquivalenceClassNumber(1) != 0 ||
	    mt3->getEquivalenceClassNumber(2) != 0) {
		throw runtime_error("getEquivalenceClassNumber did not map site1, site2, site3 to class 0 properly.");
	}

	if (mt3->getEquivalenceClassNumber(3) != -1) {
		throw runtime_error("getEquivalenceClassNumber did not map otherSite to class -1 properly, got " + to_string(mt3->getEquivalenceClassNumber(3)));
	}

	int* components;
	int n_components;
	mt3->getEquivalencyClass(components, n_components, "site");

	if (n_components != 3) {
		throw runtime_error("getEquivalencyClass returned " + to_string(n_components) + " components for 'site' instead of 3");
	}

	if (components[0] != 0 || components[1] != 1 || components[2] != 2) {
		throw runtime_error("getEquivalencyClass did not return the correct component indices for 'site'");
	}

	cout << "  MoleculeType::addEquivalentComponents tests passed!" << endl;

	cout << "NFcore::MoleculeType tests completed successfully." << endl;

	// System destructor will free molecule types instantiated.
	delete s;
}
