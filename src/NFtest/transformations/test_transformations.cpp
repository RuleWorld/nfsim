#include "test_transformations.hh"
#include "../../NFreactions/transformations/moleculeCreator.hh"
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

	System *sys = new System("TestSystem");
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
	MoleculeType *mtState = new MoleculeType(mtName, compNames, defaultStates, compStateNames, sys);
	Molecule *mol = mtState->genDefaultMolecule();
	if (mol->getComponentState(0) != 0) {
		throw runtime_error("Molecule default state is not 0.");
	}
	Mapping *m = new Mapping(1, 0);
	m->setMolecule(mol);
	StateChangeTransform *sct = new StateChangeTransform(0, 2);
	sct->apply(m, NULL);
	if (mol->getComponentState(0) != 2) {
		throw runtime_error("StateChangeTransform failed to change the state.");
	}
	string logstr = "initial";
	StateChangeTransform *sct2 = new StateChangeTransform(0, 1);
	sct2->apply(m, NULL, logstr);
	if (mol->getComponentState(0) != 1) {
		throw runtime_error("StateChangeTransform failed to change the state with logstr.");
	}
	if (logstr.find("StateChange") == string::npos) {
		throw runtime_error("StateChangeTransform did not properly populate logstr.");
	}
	m->clear();
	delete sct;
	delete sct2;
	delete m;
	delete sys;
	cout << "  StateChangeTransform tests passed!" << endl;

	cout << "  Testing TransformationSet..." << endl;

	System *s = new System("TestSystem");

	vector<string> xComps;
	xComps.push_back("a");
	xComps.push_back("p");
	vector<string> xStates;
	xStates.push_back("No State");
	xStates.push_back("U");
	vector<vector<string> > xAllowedStates;
	vector<string> aStates;
	xAllowedStates.push_back(aStates);
	vector<string> pStates;
	pStates.push_back("U");
	pStates.push_back("P");
	xAllowedStates.push_back(pStates);
	MoleculeType *molX = new MoleculeType("X", xComps, xStates, xAllowedStates, s);
	s->addMoleculeType(molX);

	vector<string> yComps;
	yComps.push_back("a");
	vector<string> yStates;
	yStates.push_back("No State");
	vector<vector<string> > yAllowedStates;
	vector<string> aStatesY;
	yAllowedStates.push_back(aStatesY);
	MoleculeType *molY = new MoleculeType("Y", yComps, yStates, yAllowedStates, s);
	s->addMoleculeType(molY);

	TemplateMolecule *tx = new TemplateMolecule(molX);
	TemplateMolecule *ty = new TemplateMolecule(molY);

	vector<TemplateMolecule*> reactants;
	reactants.push_back(tx);
	reactants.push_back(ty);

	TransformationSet *ts = new TransformationSet(reactants);

	if (ts->getNreactants() != 2) {
		throw runtime_error("TransformationSet getNreactants failed");
	}

	if (ts->getNmappingSets() != 2) {
		throw runtime_error("TransformationSet getNmappingSets failed");
	}

	ts->addStateChangeTransform(tx, "p", "P");
	if (ts->getNumOfTransformations(0) != 1) {
	    throw runtime_error("TransformationSet getNumOfTransformations failed for state change");
	}

	ts->addBindingTransform(tx, "a", ty, "a");
	if (ts->getNumOfTransformations(0) != 2) {
	    throw runtime_error("TransformationSet getNumOfTransformations failed for binding (first reactant)");
	}
	if (ts->getNumOfTransformations(1) != 1) {
	    throw runtime_error("TransformationSet getNumOfTransformations failed for binding (second reactant)");
	}

	if (ts->hasSymBindingTransform()) {
	    throw runtime_error("TransformationSet hasSymBindingTransform failed: X and Y are not symmetric");
	}

	TemplateMolecule *newT = new TemplateMolecule(molX);
	vector<pair<int, int> > stateValues;
	MoleculeCreator *mc = new MoleculeCreator(newT, molX, stateValues, NULL);
	ts->addAddMolecule(mc);
	if (ts->getNumOfAddMoleculeTransforms() != 1) {
	    throw runtime_error("TransformationSet getNumOfAddMoleculeTransforms failed");
	}

	ts->finalize();
	if (!ts->isFinalized()) {
		throw runtime_error("TransformationSet finalize failed");
	}

	ts->setComplexBookkeeping(true);
	if (!ts->getComplexBookkeeping()) {
	    throw runtime_error("TransformationSet getComplexBookkeeping failed");
	}

	delete ts;

    TemplateMolecule *tx3 = new TemplateMolecule(molX);
	TemplateMolecule *ty3 = new TemplateMolecule(molX);
	vector<TemplateMolecule*> reactants3;
	reactants3.push_back(tx3);
	reactants3.push_back(ty3);
	TransformationSet *ts3 = new TransformationSet(reactants3);
	ts3->addBindingTransform(tx3, "a", ty3, "a");
	if (!ts3->hasSymBindingTransform()) {
	    throw runtime_error("TransformationSet hasSymBindingTransform failed: should be symmetric");
	}
	ts3->finalize();
	delete ts3;

    TemplateMolecule *tx2 = new TemplateMolecule(molX);
    TemplateMolecule *ty2 = new TemplateMolecule(molX);
    tx2->addBond("a", ty2, "a");
	ty2->addBond("a", tx2, "a");
	vector<TemplateMolecule*> reactants2;
	reactants2.push_back(tx2);
	TransformationSet *ts2 = new TransformationSet(reactants2);
	ts2->addUnbindingTransform(tx2, "a", ty2, "a");
	if (!ts2->hasSymUnbindingTransform()) {
	    throw runtime_error("TransformationSet hasSymUnbindingTransform failed");
	}
	ts2->finalize();
	delete ts2;

	cout << "  TransformationSet basic tests passed!" << endl;
	cout << "Transformations tests completed successfully." << endl;
}
