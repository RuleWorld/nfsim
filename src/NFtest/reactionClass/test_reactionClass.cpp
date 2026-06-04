#include "test_reactionClass.hh"
#include "../../NFreactions/reactions/reaction.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

using namespace std;
using namespace NFcore;

void NFtest_reactionClass::run()
{
	cout << "Running ReactionClass tests..." << endl;

	cout << "  Testing ReactionClass::fire..." << endl;

	// Instantiate a System
	System* sys = new System("TestSystem");

	// Create a dummy TransformationSet
	vector<TemplateMolecule*> emptyTemplates;
	TransformationSet* ts1 = new TransformationSet(emptyTemplates);
	ts1->finalize();

	// Create a ReactionClass
	ReactionClass* rxn1 = new BasicRxnClass("Rxn1", 1.0, "", ts1, sys);

	// Test fire(double random_A_number, bool track)
	int initialFireCounter = rxn1->getFireCounter();

	// Test without tracking
	rxn1->fire(0.5, false);
	if (rxn1->getFireCounter() != initialFireCounter + 1) {
		throw std::runtime_error("ReactionClass::fire without tracking did not increment fireCounter.");
	}

	// Test with tracking
	string log = rxn1->fire(0.5, true);
	if (rxn1->getFireCounter() != initialFireCounter + 2) {
		throw std::runtime_error("ReactionClass::fire with tracking did not increment fireCounter.");
	}

	// Because n_reactants is 0, transformationSet->checkMolecularity(mappingSet) might fail if it needs reactants,
	// or it might pass. We just need to check that fire() executes and updates the fireCounter.

	cout << "  ReactionClass::fire tests passed!" << endl;

	delete sys;

	cout << "ReactionClass tests completed successfully." << endl;
}
