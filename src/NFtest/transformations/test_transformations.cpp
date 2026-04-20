#include "test_transformations.hh"
#include "../../NFreactions/transformations/transformation.hh"
#include "../../NFcore/NFcore.hh"
#include "../../NFreactions/mappings/mapping.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

using namespace std;
using namespace NFcore;

void NFtest_transformations::run()
{
	cout << "Running transformations tests..." << endl;

	cout << "  Testing StateChangeTransform..." << endl;

	// Create a minimal system
	System *sys = new System("TestSystem");

	// Create a molecule type with 1 component that has states
	string mtName = "TestMol";
	vector<string> compNames;
	compNames.push_back("site");

	vector<string> defaultStates;
	defaultStates.push_back("S0");

	vector<vector<string> > compStateNames;
	vector<string> sNames;
	sNames.push_back("S0");
	sNames.push_back("S1");
	sNames.push_back("S2");
	compStateNames.push_back(sNames);

	MoleculeType *mt = new MoleculeType(mtName, compNames, defaultStates, compStateNames, sys);

	// Create a molecule of this type
	Molecule *mol = mt->genDefaultMolecule();

	// Initially, the state of component 0 should be 0
	if (mol->getComponentState(0) != 0) {
		throw runtime_error("Molecule default state is not 0.");
	}

	// Create a mapping to this molecule
	// type is TransformationFactory::STATE_CHANGE = 1, index is 0
	Mapping *m = new Mapping(1, 0);
	m->setMolecule(mol);

	// Create the StateChangeTransform
	// We want to change the state of component 0 to 2
	int cIndex = 0;
	int newValue = 2;
	StateChangeTransform *sct = new StateChangeTransform(cIndex, newValue);

	// Apply the transform
	sct->apply(m, NULL);

	// Check if the state changed
	if (mol->getComponentState(0) != 2) {
		throw runtime_error("StateChangeTransform failed to change the state.");
	}

	// Test the logstr signature
	string logstr = "initial"; // Initialize so that !logstr.empty() is true
	StateChangeTransform *sct2 = new StateChangeTransform(cIndex, 1);
	sct2->apply(m, NULL, logstr);

	if (mol->getComponentState(0) != 1) {
		throw runtime_error("StateChangeTransform failed to change the state with logstr.");
	}

	if (logstr.find("StateChange") == string::npos) {
		throw runtime_error("StateChangeTransform did not properly populate logstr.");
	}

	// Cleanup
	m->clear();
	delete sct;
	delete sct2;
	delete m;
	// mol is cleaned up when system/molecule type are deleted
	delete sys;

	cout << "  StateChangeTransform tests passed!" << endl;
	cout << "Transformations tests completed successfully." << endl;
}
