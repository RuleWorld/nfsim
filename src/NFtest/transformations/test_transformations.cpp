#include "test_transformations.hh"
#include "../../NFreactions/transformations/moleculeCreator.hh"
#include "../../NFreactions/transformations/transformation.hh"
#include "../../NFreactions/transformations/speciesCreator.hh"
#include "../../NFcore/NFcore.hh"
#include "../../NFreactions/mappings/mapping.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <sstream>

using namespace std;
using namespace NFcore;

class TestTransformationSet : public TransformationSet {
public:
	TestTransformationSet(vector<TemplateMolecule*> &reactants) : TransformationSet(reactants) {}
	bool testCanReach(Molecule *m1, Molecule *m2, int excludeComp) {
		return canReachExcludingBond(m1, m2, excludeComp);
	}
};


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

	{
		TemplateMolecule *tMissing = new TemplateMolecule(molX);
		std::ostringstream localCerr;
		std::streambuf* oldCerr = std::cerr.rdbuf(localCerr.rdbuf());
		bool result = ts->addStateChangeTransform(tMissing, "p", "P");
		std::cerr.rdbuf(oldCerr);

		if (result) {
			throw runtime_error("TransformationSet addStateChangeTransform should have failed for missing template");
		}
		if (localCerr.str().find("Couldn't find the template you gave me") == string::npos) {
			throw runtime_error("TransformationSet addStateChangeTransform did not output expected error message");
		}
		delete tMissing;
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






	// --- Testing canReachExcludingBond ---
	cout << "  Testing canReachExcludingBond..." << endl;
	vector<string> ringComps;
	ringComps.push_back("s1");
	ringComps.push_back("s2");
	ringComps.push_back("s3");
	vector<string> ringStates;
	ringStates.push_back("No State");
	ringStates.push_back("No State");
	ringStates.push_back("No State");
	vector<vector<string> > ringAllowedStates(3); // 3 components, empty lists means no states
	vector<string> noStates;
	ringAllowedStates[0] = noStates;
	ringAllowedStates[1] = noStates;
	ringAllowedStates[2] = noStates;
	MoleculeType *molRing = new MoleculeType("Ring", ringComps, ringStates, ringAllowedStates, s);
	s->addMoleculeType(molRing);

	Molecule *m1 = molRing->genDefaultMolecule();
	Molecule *m2 = molRing->genDefaultMolecule();
	Molecule *m3 = molRing->genDefaultMolecule();
	Molecule *m4 = molRing->genDefaultMolecule();

	// Topology:
	// m1(s1) - m2(s1)
	// m2(s2) - m3(s1)
	// m3(s2) - m4(s1)

	Molecule::bind(m1, 0, m2, 0);
	Molecule::bind(m2, 1, m3, 0);
	Molecule::bind(m3, 1, m4, 0);

	TemplateMolecule *tm1 = new TemplateMolecule(molRing);
	vector<TemplateMolecule*> ringReactants;
	ringReactants.push_back(tm1);
	TestTransformationSet *testTS = new TestTransformationSet(ringReactants);

	// Test on the line m1 - m2. Exclude bond at m1's s1 (index 0).
	if (testTS->testCanReach(m1, m2, 0) != false) {
		throw runtime_error("canReachExcludingBond failed on line topology (should be false)");
	}

	// Close the ring: m4(s2) - m1(s2)
	Molecule::bind(m4, 1, m1, 1);

	// Now m1, m2, m3, m4 are in a ring. Test excluding bond at m1's s1 (index 0).
	if (testTS->testCanReach(m1, m2, 0) != true) {
		throw runtime_error("canReachExcludingBond failed on ring topology (should be true)");
	}

	// Add another branch to test BFS robustness
	Molecule *m5 = molRing->genDefaultMolecule();
	Molecule::bind(m3, 2, m5, 0); // m3(s3) - m5(s1)

	if (testTS->testCanReach(m1, m2, 0) != true) {
		throw runtime_error("canReachExcludingBond failed on ring topology with branch (should be true)");
	}

	delete testTS;
	delete tm1;

	cout << "  canReachExcludingBond tests passed!" << endl;



	cout << "  TransformationSet basic tests passed!" << endl;

	cout << "  Testing SpeciesCreator..." << endl;

	vector<MoleculeType*> productMoleculeTypes;
    productMoleculeTypes.push_back(molX);
    productMoleculeTypes.push_back(molY);

    vector<vector<int> > stateInformation(3);
    stateInformation[0].push_back(0); // ndStateMolecule: molX
    stateInformation[1].push_back(1); // ndStateIndex: "p"
    stateInformation[2].push_back(1); // ndStateValue: "P"

    vector<vector<int> > bindingSiteInformation(4);
    bindingSiteInformation[0].push_back(0); // bMolecule1: molX
    bindingSiteInformation[1].push_back(0); // bSite1: "a"
    bindingSiteInformation[2].push_back(1); // bMolecule2: molY
    bindingSiteInformation[3].push_back(0); // bSite2: "a"

    SpeciesCreator *sc = new SpeciesCreator(productMoleculeTypes, stateInformation, bindingSiteInformation);

    s->prepareForSimulation();

    int initialX = molX->getMoleculeCount();
    int initialY = molY->getMoleculeCount();

    if (initialX != 0 || initialY != 0) {
        throw runtime_error("Initial counts should be 0");
    }

    sc->create();

    if (molX->getMoleculeCount() != 1) throw runtime_error("Expected 1 X molecule");
    if (molY->getMoleculeCount() != 1) throw runtime_error("Expected 1 Y molecule");

    string logstr_species = "initial";
    sc->create(logstr_species);

    if (molX->getMoleculeCount() != 2) throw runtime_error("Expected 2 X molecules");
    if (molY->getMoleculeCount() != 2) throw runtime_error("Expected 2 Y molecules");

    if (logstr_species.find("AddSpecies") == string::npos) {
        throw runtime_error("Log string should contain AddSpecies");
    }
    if (logstr_species.find("StateChange") == string::npos) throw runtime_error("Log string should contain StateChange");
    if (logstr_species.find("AddBond") == string::npos) throw runtime_error("Log string should contain AddBond");

    delete sc;
    // Do not delete s here to avoid segmentation fault; it appears molecules hold references and System destructor cascades deletes.

	cout << "  SpeciesCreator tests passed!" << endl;

    // Test addExcludeReactant
    TemplateMolecule *filterPattern = new TemplateMolecule(molX);
    map<string, TemplateMolecule*> parsedTemplates;

    // We can reuse ts3 or create a new one. Let's create a new one.
    TemplateMolecule *tx4 = new TemplateMolecule(molX);
    vector<TemplateMolecule*> reactants4;
    reactants4.push_back(tx4);
    TransformationSet *ts4 = new TransformationSet(reactants4);

    // The filter expects it to match to return false.
    ts4->addExcludeReactant(0, filterPattern, parsedTemplates);

    Molecule *molX_test = molX->genDefaultMolecule();
    bool checkFilter = ts4->checkReactantFilters(0, molX_test);
    if (checkFilter) {
        throw runtime_error("checkReactantFilters failed to exclude molecule matching pattern");
    }

    // Test that a filter on a different reactant index doesn't exclude it
    bool checkFilterDiffIndex = ts4->checkReactantFilters(1, molX_test);
    if (!checkFilterDiffIndex) {
        throw runtime_error("checkReactantFilters excluded molecule when index didn't match");
    }

    delete ts4;
    cout << "  TransformationSet::addExcludeReactant tests passed!" << endl;

	cout << "Transformations tests completed successfully." << endl;
}
