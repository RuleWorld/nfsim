#include "test_system.hh"
#include "../../NFcore/compartment.hh"
#include "../../NFreactions/reactions/reaction.hh"
#include "../../NFfunction/NFfunction.hh"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <string>

using namespace std;
using namespace NFcore;

void NFtest_system::run()
{
	cout << "Running System tests..." << endl;

	// A compartment size explicitly linked to a parameter must follow updates.
	{
		System linkedSizeSystem("LinkedSizeSystem");
		Compartment *compartment = new Compartment("c", 3, 2.0);
		linkedSizeSystem.addCompartment(compartment);
		linkedSizeSystem.addParameter("Vc", 2.0);
		linkedSizeSystem.addParameter("kp", 100.0);
		vector<TemplateMolecule *> noReactants;
		TransformationSet *zeroOrderSet = new TransformationSet(noReactants);
		zeroOrderSet->finalize();
		BasicRxnClass *synthesis = new BasicRxnClass("synthesis", 100.0, "kp", zeroOrderSet, &linkedSizeSystem);
		linkedSizeSystem.addReaction(synthesis);
		synthesis->applyZeroOrderVolumeScaling(compartment);
		linkedSizeSystem.prepareForSimulation();
		if (synthesis->getBaseRate() != 200.0)
			throw runtime_error("Zero-order rate was not scaled by its initial compartment size.");
		if (!linkedSizeSystem.setCompartmentSizeParameter("c", "Vc"))
			throw runtime_error("Failed to link compartment size to parameter.");
		linkedSizeSystem.setParameter("Vc", 4.0);
		linkedSizeSystem.updateSystemWithNewParameters();
		if (compartment->getSize() != 4.0)
			throw runtime_error("Compartment size did not follow its linked parameter.");
		if (synthesis->getBaseRate() != 400.0)
			throw runtime_error("Zero-order propensity rate did not follow its updated compartment size.");
		linkedSizeSystem.setParameter("Vc", 2.0);
		linkedSizeSystem.updateSystemWithNewParameters();
		if (compartment->getSize() != 2.0)
			throw runtime_error("Compartment size did not recover after a repeated parameter update.");
		if (synthesis->getBaseRate() != 200.0)
			throw runtime_error("Repeated zero-order size updates compounded volume scaling.");
		linkedSizeSystem.setParameter("kp", 250.0);
		linkedSizeSystem.updateSystemWithNewParameters();
		if (synthesis->getBaseRate() != 500.0)
			throw runtime_error("Zero-order rate parameter update was not preserved through volume scaling.");
	}
	{
		System parameterOnlySystem("ParameterOnlySystem");
		parameterOnlySystem.addParameter("kp", 3.0);
		vector<TemplateMolecule *> noReactants;
		TransformationSet *zeroOrderSet = new TransformationSet(noReactants);
		zeroOrderSet->finalize();
		BasicRxnClass *synthesis = new BasicRxnClass("synthesis", 3.0, "kp", zeroOrderSet, &parameterOnlySystem);
		parameterOnlySystem.addReaction(synthesis);
		parameterOnlySystem.prepareForSimulation();
		parameterOnlySystem.setParameter("kp", 7.0);
		parameterOnlySystem.updateSystemWithNewParameters();
		if (synthesis->getBaseRate() != 7.0)
			throw runtime_error("Zero-order rate parameter update without compartments was lost.");
	}
	{
		System bimolecularSystem("BimolecularVolumeSystem");
		Compartment *compartment = new Compartment("c", 3, 2.0);
		bimolecularSystem.addCompartment(compartment);
		vector<TemplateMolecule *> noReactants;
		TransformationSet *reactionSet = new TransformationSet(noReactants);
		reactionSet->finalize();
		BasicRxnClass *reaction = new BasicRxnClass("association", 8.0, "", reactionSet, &bimolecularSystem);
		reaction->setReactantVolumeCompartments(vector<Compartment *>(1, compartment));
		if (reaction->getBaseRate() != 4.0)
			throw runtime_error("Reactant-bearing count rate did not divide by the compartment size.");
		bimolecularSystem.addParameter("Vc", 4.0);
		bimolecularSystem.setCompartmentSizeParameter("c", "Vc");
		compartment->setSize(4.0);
		reaction->refreshReactantVolumeFactor();
		if (reaction->getBaseRate() != 2.0)
			throw runtime_error("Reactant-bearing count rate did not follow the updated compartment size.");
		delete reaction;
	}
	{
		System unitVolumeSystem("UnitVolumeSystem");
		unitVolumeSystem.addCompartment(new Compartment("unit", 3, 1.0));
		if (!unitVolumeSystem.hasStaticUnitSizedCompartments())
			throw runtime_error("Static unit-size compartments were not recognized as invariant volume conversion.");
		unitVolumeSystem.addParameter("unit", 1.0);
		if (!unitVolumeSystem.setCompartmentSizeParameter("unit", "unit") ||
			unitVolumeSystem.hasStaticUnitSizedCompartments())
			throw runtime_error("Parameter-linked unit-size compartments were incorrectly treated as static.");
	}

	cout << "  Testing System Constructors..." << endl;

	// Test System(string name)
	System* sys1 = new System("TestSystem1");
	if (sys1->getName() != "TestSystem1") {
		throw std::runtime_error("System(string) did not initialize name correctly.");
	}
	if (sys1->isUsingComplex() != false) {
		throw std::runtime_error("System(string) did not initialize useComplex to false.");
	}
	if (sys1->getGlobalMoleculeLimit() != 100000) {
		throw std::runtime_error("System(string) did not initialize globalMoleculeLimit to 100000.");
	}
	delete sys1;

	// Test System(string name, bool useComplex)
	System* sys2 = new System("TestSystem2", true);
	if (sys2->getName() != "TestSystem2") {
		throw std::runtime_error("System(string, bool) did not initialize name correctly.");
	}
	if (sys2->isUsingComplex() != true) {
		throw std::runtime_error("System(string, bool) did not initialize useComplex correctly.");
	}
	if (sys2->getGlobalMoleculeLimit() != 100000) {
		throw std::runtime_error("System(string, bool) did not initialize globalMoleculeLimit to 100000.");
	}
	delete sys2;

	// Test System(string name, bool useComplex, int globalMoleculeLimit)
	System* sys3 = new System("TestSystem3", false, 500);
	if (sys3->getName() != "TestSystem3") {
		throw std::runtime_error("System(string, bool, int) did not initialize name correctly.");
	}
	if (sys3->isUsingComplex() != false) {
		throw std::runtime_error("System(string, bool, int) did not initialize useComplex correctly.");
	}
	if (sys3->getGlobalMoleculeLimit() != 500) {
		throw std::runtime_error("System(string, bool, int) did not initialize globalMoleculeLimit correctly.");
	}
	delete sys3;

	// Test edge cases for globalMoleculeLimit
	System* sys4 = new System("TestSystem4", true, 0);
	if (sys4->getGlobalMoleculeLimit() != 0) {
		throw std::runtime_error("System(string, bool, int) did not initialize globalMoleculeLimit to 0.");
	}
	delete sys4;

	System* sys5 = new System("TestSystem5", false, -100);
	if (sys5->getGlobalMoleculeLimit() != -100) {
		throw std::runtime_error("System(string, bool, int) did not initialize globalMoleculeLimit to negative correctly.");
	}
	delete sys5;

	cout << "  Testing System name validation..." << endl;
	const vector<string> invalidNames = {
		"../escape", "/absolute", "C:\\escape", "C:escape",
		"nested/name", "nested\\name", "contains..dots"
	};
	for (const string &invalidName : invalidNames) {
		for (int constructor = 0; constructor < 3; ++constructor) {
			bool rejected = false;
			try {
				if (constructor == 0) {
					System invalidSystem(invalidName);
				} else if (constructor == 1) {
					System invalidSystem(invalidName, true);
				} else {
					System invalidSystem(invalidName, true, 100);
				}
			} catch (const runtime_error &) {
				rejected = true;
			}
			if (!rejected) {
				throw runtime_error("System accepted unsafe name: " + invalidName);
			}
		}
	}

	cout << "  Testing System::addMoleculeType..." << endl;

	// Instantiate a System
	System* sys = new System("TestSystem");

	// Initial checks
	if (sys->getNumOfMoleculeTypes() != 0) {
		throw std::runtime_error("System did not initialize with 0 MoleculeTypes.");
	}

	// Create first MoleculeType
	vector<string> compNameA;
	compNameA.push_back("x");
	MoleculeType* mtA = new MoleculeType("MolA", compNameA, sys);

	// Verify first MoleculeType addition
	if (mtA->getTypeID() != 0) {
		throw std::runtime_error("First MoleculeType type ID is not 0.");
	}
	if (sys->getNumOfMoleculeTypes() != 1) {
		throw std::runtime_error("System number of MoleculeTypes is not 1 after adding first MoleculeType.");
	}
	if (sys->getMoleculeType(0) != mtA) {
		throw std::runtime_error("System::getMoleculeType(0) did not return the correct MoleculeType pointer.");
	}

	// Create second MoleculeType
	vector<string> compNameB;
	compNameB.push_back("y");
	compNameB.push_back("z");
	MoleculeType* mtB = new MoleculeType("MolB", compNameB, sys);

	// Verify second MoleculeType addition
	if (mtB->getTypeID() != 1) {
		throw std::runtime_error("Second MoleculeType type ID is not 1.");
	}
	if (sys->getNumOfMoleculeTypes() != 2) {
		throw std::runtime_error("System number of MoleculeTypes is not 2 after adding second MoleculeType.");
	}
	if (sys->getMoleculeType(1) != mtB) {
		throw std::runtime_error("System::getMoleculeType(1) did not return the correct MoleculeType pointer.");
	}

	// Verify that the first MoleculeType is still accessible correctly
	if (sys->getMoleculeType(0) != mtA) {
		throw std::runtime_error("System::getMoleculeType(0) did not return the correct MoleculeType pointer after adding second MoleculeType.");
	}

	cout << "  System::addMoleculeType tests passed!" << endl;

	cout << "  Testing System::getReactionByName..." << endl;

	// Create an empty TransformationSet for dummy reactions
	vector<TemplateMolecule*> emptyTemplates;
	TransformationSet* ts1 = new TransformationSet(emptyTemplates);
	ts1->finalize();

	TransformationSet* ts2 = new TransformationSet(emptyTemplates);
	ts2->finalize();

	ReactionClass* rxn1 = new BasicRxnClass("Rxn1", 1.0, "", ts1, sys);
	ReactionClass* rxn2 = new BasicRxnClass("Rxn2", 2.0, "", ts2, sys);

	sys->addReaction(rxn1);
	sys->addReaction(rxn2);

	// Verify we can find reactions by exact name
	if (sys->getReactionByName("Rxn1") != rxn1) {
		throw std::runtime_error("System::getReactionByName did not return correct reaction for 'Rxn1'.");
	}

	if (sys->getReactionByName("Rxn2") != rxn2) {
		throw std::runtime_error("System::getReactionByName did not return correct reaction for 'Rxn2'.");
	}

	// Verify we get NULL for nonexistent reactions
	if (sys->getReactionByName("Rxn3") != NULL) {
		throw std::runtime_error("System::getReactionByName did not return NULL for nonexistent reaction 'Rxn3'.");
	}

	if (sys->getReactionByName("") != NULL) {
		throw std::runtime_error("System::getReactionByName did not return NULL for empty string name.");
	}

	cout << "  System::getReactionByName tests passed!" << endl;

	cout << "  Testing System::getGlobalFunctionByName..." << endl;

	vector<string> emptyVec;
	GlobalFunction* gf1 = new GlobalFunction("Func1", "1.0", emptyVec, emptyVec, emptyVec, sys);
	GlobalFunction* gf2 = new GlobalFunction("Func2", "2.0", emptyVec, emptyVec, emptyVec, sys);

	sys->addGlobalFunction(gf1);
	sys->addGlobalFunction(gf2);

	if (sys->getGlobalFunctionByName("Func1") != gf1) {
		throw std::runtime_error("System::getGlobalFunctionByName did not return correct function for 'Func1'.");
	}

	if (sys->getGlobalFunctionByName("Func2") != gf2) {
		throw std::runtime_error("System::getGlobalFunctionByName did not return correct function for 'Func2'.");
	}

	if (sys->getGlobalFunctionByName("Func3") != NULL) {
		throw std::runtime_error("System::getGlobalFunctionByName did not return NULL for nonexistent function 'Func3'.");
	}

	if (sys->getGlobalFunctionByName("") != NULL) {
		throw std::runtime_error("System::getGlobalFunctionByName did not return NULL for empty string name.");
	}

	cout << "  System::getGlobalFunctionByName tests passed!" << endl;

	cout << "  Testing System::getMoleculeTypeByName..." << endl;

	// Test happy path
	if (sys->getMoleculeTypeByName("MolA") != mtA) {
		throw std::runtime_error("System::getMoleculeTypeByName did not return correct MoleculeType for 'MolA'.");
	}
	if (sys->getMoleculeTypeByName("MolB") != mtB) {
		throw std::runtime_error("System::getMoleculeTypeByName did not return correct MoleculeType for 'MolB'.");
	}

	// Test error path
	bool caughtError = false;
	try {
		sys->getMoleculeTypeByName("NonExistentMol");
	} catch (const std::runtime_error& e) {
		caughtError = true;
	}

	if (!caughtError) {
		throw std::runtime_error("System::getMoleculeTypeByName did not throw std::runtime_error for non-existent molecule type.");
	}

	cout << "  System::getMoleculeTypeByName tests passed!" << endl;

	cout << "  Testing System::getLocalFunctionByName..." << endl;

	// Create dummy arguments for LocalFunction
	vector<string> args;
	args.push_back("arg1");
	vector<string> varRefNames;
	vector<string> varObservableNames;
	vector<Observable *> varObservables;
	vector<int> varRefScope;
	vector<string> paramNames;

	LocalFunction* lf1 = new LocalFunction(sys, "func1", "1.0", "1.0", args, varRefNames, varObservableNames, varObservables, varRefScope, paramNames);
	LocalFunction* lf2 = new LocalFunction(sys, "func2", "2.0", "2.0", args, varRefNames, varObservableNames, varObservables, varRefScope, paramNames);

	sys->addLocalFunction(lf1);
	sys->addLocalFunction(lf2);

	if (sys->getLocalFunctionByName("func1") != lf1) {
		throw std::runtime_error("System::getLocalFunctionByName did not return correct LocalFunction for 'func1'.");
	}
	if (sys->getLocalFunctionByName("func2") != lf2) {
		throw std::runtime_error("System::getLocalFunctionByName did not return correct LocalFunction for 'func2'.");
	}
	if (sys->getLocalFunctionByName("NonExistentFunc") != NULL) {
		throw std::runtime_error("System::getLocalFunctionByName did not return NULL for non-existent local function.");
	}
	if (sys->getLocalFunctionByName("") != NULL) {
		throw std::runtime_error("System::getLocalFunctionByName did not return NULL for empty string name.");
	}

	cout << "  System::getLocalFunctionByName tests passed!" << endl;

	cout << "  Testing System::getObservableByName..." << endl;

	// Create and add an observable
	vector<TemplateMolecule*> obsTemplates;
	TemplateMolecule* tmObs = new TemplateMolecule(mtA);
	obsTemplates.push_back(tmObs);
	Observable* obs1 = new MoleculesObservable("Obs1", obsTemplates);
	sys->addObservableForOutput(obs1);

	if (sys->getObservableByName("Obs1") != obs1) {
		throw std::runtime_error("System::getObservableByName did not return correct observable for 'Obs1'.");
	}

	// Capture cout/cerr for edge cases (avoid console spam during tests)
	std::streambuf* oldCout = std::cout.rdbuf();
	std::streambuf* oldCerr = std::cerr.rdbuf();
	std::ostringstream localCout;
	std::ostringstream localCerr;
	std::cout.rdbuf(localCout.rdbuf());
	std::cerr.rdbuf(localCerr.rdbuf());

	Observable* notFoundObs = sys->getObservableByName("NonExistentObs");

	std::cout.rdbuf(oldCout);
	std::cerr.rdbuf(oldCerr);

	if (notFoundObs != NULL) {
		throw std::runtime_error("System::getObservableByName did not return NULL for nonexistent observable.");
	}

	string expectedWarning = "!!Warning, the system could not identify the observable: NonExistentObs.\n";
	if (localCerr.str().find(expectedWarning) == string::npos) {
		throw std::runtime_error("System::getObservableByName did not print correct warning for nonexistent observable.");
	}

	cout << "  System::getObservableByName tests passed!" << endl;

	vector<string> funcs, argNames;
	CompositeFunction* cf1 = new CompositeFunction(sys, "Func1", "x + y", funcs, argNames, paramNames);
	CompositeFunction* cf2 = new CompositeFunction(sys, "Func2", "a * b", funcs, argNames, paramNames);
	sys->addCompositeFunction(cf1);
	sys->addCompositeFunction(cf2);

	if (sys->getCompositeFunctionByName("Func1") != cf1) {
		throw std::runtime_error("System::getCompositeFunctionByName did not return correct function for 'Func1'.");
	}
	if (sys->getCompositeFunctionByName("Func2") != cf2) {
		throw std::runtime_error("System::getCompositeFunctionByName did not return correct function for 'Func2'.");
	}
	if (sys->getCompositeFunctionByName("Func3") != NULL) {
		throw std::runtime_error("System::getCompositeFunctionByName did not return NULL for nonexistent function 'Func3'.");
	}
	if (sys->getCompositeFunctionByName("") != NULL) {
		throw std::runtime_error("System::getCompositeFunctionByName did not return NULL for empty string name.");
	}

	cout << "  System::getCompositeFunctionByName tests passed!" << endl;

	// Cleanup
	// Note: System destructor deletes MoleculeTypes added to it
	// It also cleans up reactions added to it.
	delete sys;

	cout << "System tests completed successfully." << endl;
}
