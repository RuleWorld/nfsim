#include "test_templateMolecule.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <sstream>

using namespace std;
using namespace NFcore;

void NFtest_templateMolecule::run()
{
	cout << "Running NFcore::TemplateMolecule tests..." << endl;

	cout << "  Testing TemplateMolecule::printDetails..." << endl;

	// Set up a basic system and molecule type to test against
	System* s = new System("test");

	vector<string> compNames;
	compNames.push_back("c");

	vector<string> defaultStates;
	defaultStates.push_back("s");

	vector<vector<string>> allowedStates;
	vector<string> compAllowedStates;
	compAllowedStates.push_back("s");
	compAllowedStates.push_back("p");
	allowedStates.push_back(compAllowedStates);

	MoleculeType* mt = new MoleculeType("testMT", compNames, defaultStates, allowedStates, s);

	// Create a TemplateMolecule
	TemplateMolecule* tm = new TemplateMolecule(mt);

	// Add an empty site constraint
	tm->addEmptyComponent("c");

	// Redirect cout to capture string stream
	ostringstream oss;
	tm->printDetails(oss);

	string output = oss.str();

	// Check for expected output substrings
	if (output.find("TemplateMolecule of type:   testMT") == string::npos) {
		throw runtime_error("printDetails did not output the expected MoleculeType name");
	}

	if (output.find("Connected-to:                       none") == string::npos) {
		throw runtime_error("printDetails did not output expected 'Connected-to' state");
	}

	if (output.find("Empty Binding Site Constraints:      c(index=0)") == string::npos) {
		throw runtime_error("printDetails did not output expected empty component constraint");
	}

	if (output.find("Occupied Binding Site Constraints:   none") == string::npos) {
		throw runtime_error("printDetails did not output expected occupied component state");
	}

	cout << "  TemplateMolecule::printDetails tests passed!" << endl;
	cout << "NFcore::TemplateMolecule tests completed successfully." << endl;

	delete s; // s destructor cascade-deletes tm
}
