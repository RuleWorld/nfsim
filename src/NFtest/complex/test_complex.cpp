#include "test_complex.hh"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <string>

using namespace std;
using namespace NFcore;

void NFtest_complex::run()
{
	cout << "Running NFcore::Complex tests..." << endl;

	cout << "  Testing Complex::printDetails..." << endl;

	// Redirect cout to capture output
	streambuf* oldCoutStreamBuf = cout.rdbuf();
	ostringstream strCout;
	cout.rdbuf(strCout.rdbuf());

	// Create a System, MoleculeType, and Molecule
	System* s = new System("test");

	vector<string> compNames;
	compNames.push_back("c");

	vector<string> defaultStates;
	defaultStates.push_back("s");

	vector<vector<string>> allowedStates;
	vector<string> compAllowedStates;
	compAllowedStates.push_back("s");
	allowedStates.push_back(compAllowedStates);

	MoleculeType* mt = new MoleculeType("testMT", compNames, defaultStates, allowedStates, s);

	Molecule* m = new Molecule(mt, 0, NULL);

	// Create a complex and test printDetails
	Complex* c = new Complex(s, 123, m);

	c->printDetails();

	// Restore cout
	cout.rdbuf(oldCoutStreamBuf);

	string output = strCout.str();
	string expected = "   -Complex 123: (1) - testMT__u" + to_string(m->getUniqueID()) + "\n";

	if (output != expected) {
		throw runtime_error("printDetails output did not match expected output.\nExpected: '" + expected + "'\nGot: '" + output + "'");
	}

	cout << "  Complex::printDetails tests passed!" << endl;
	cout << "NFcore::Complex tests completed successfully." << endl;

	delete s;
}
