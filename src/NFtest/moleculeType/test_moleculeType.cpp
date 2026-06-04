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


	cout << "NFcore::MoleculeType tests completed successfully." << endl;

	// System destructor will free molecule types instantiated.
	delete s;
}
