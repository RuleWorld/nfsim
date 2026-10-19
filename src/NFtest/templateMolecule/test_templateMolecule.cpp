#include "test_templateMolecule.hh"
#include "../../NFinput/NFinput.hh"
#include <memory>
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

    // RuleWorld/nfsim#107: disjoint sites must not connect basic rules.
    {
    vector<string> sites{"a", "b"};
    vector<string> defaults(2, "none");
    vector<vector<string>> states(2);
    MoleculeType* polymer = new MoleculeType("polymer", sites, defaults, states, s);
    // MoleculeType owns and deletes registered templates with the System.
    // Keep them alive until then; stack instances would leave dangling entries.
    TemplateMolecule* left = new TemplateMolecule(polymer);
    TemplateMolecule* right = new TemplateMolecule(polymer);
    TemplateMolecule* same = new TemplateMolecule(polymer);
    TemplateMolecule* bare = new TemplateMolecule(polymer);
    left->addEmptyComponent("a");
    right->addEmptyComponent("b");
    same->addEmptyComponent("a");
    if (!left->isTemplateCompatible(right) || left->sharesSpecifiedComponentWith(right))
        throw runtime_error("disjoint sites were not pruned (nfsim#107)");
    if (!left->sharesSpecifiedComponentWith(same) || !left->sharesSpecifiedComponentWith(bare))
        throw runtime_error("overlap or bare template was pruned (nfsim#107)");
    cout << "Component overlap regression passed." << endl;
    }

    // Full XML connectivity matrix from the original BNG3 regression.
    int traversal = ReactionClass::NO_LIMIT;
    unique_ptr<System> connected(NFinput::initializeFromXML(
        "validate/regressions/issue107.xml", false, 100000, false,
        traversal, true, true));
    if (!connected) throw runtime_error("could not load issue107 fixture");
    connected->setUniversalTraversalLimit(traversal);
    connected->prepareForSimulation();
    const auto reactions = connected->getAllReactions();
    if (reactions.size() != 5) throw runtime_error("expected five issue107 rules");
    int edges = 0;
    for (auto* fired : reactions)
        for (auto* other : reactions)
            if (connected->areReactionsConnected(fired->getRxnId(), other->getRxnId())) ++edges;
    if (edges != 10) throw runtime_error("expected 10 issue107 edges, got " + to_string(edges));
    cout << "Sparse connectivity regression passed (10/25 edges)." << endl;

    // A compartment-changing rule can enable a target constrained to another
    // component in the destination compartment; retain that matrix edge.
    unique_ptr<System> moved(NFinput::initializeFromXML(
        "validate/regressions/issue107_move.xml", false, 100000, false,
        traversal, true, true));
    if (!moved) throw runtime_error("could not load issue107 move fixture");
    moved->setUniversalTraversalLimit(traversal);
    moved->prepareForSimulation();
    const auto moveReactions = moved->getAllReactions();
    if (moveReactions.size() != 2 ||
        !moved->areReactionsConnected(moveReactions[0]->getRxnId(),
                                      moveReactions[1]->getRxnId()))
        throw runtime_error("MOVE-to-compartment connectivity edge was pruned (nfsim#107)");
    cout << "Compartment-changing connectivity regression passed." << endl;

	delete s; // s destructor cascade-deletes tm
}
