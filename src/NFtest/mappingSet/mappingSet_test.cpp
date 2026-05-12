#include "mappingSet_test.hh"
#include "../../NFreactions/transformations/transformation.hh"
#include "../../NFreactions/mappings/mappingSet.hh"
#include <iostream>
#include <vector>

using namespace std;
using namespace NFcore;

void NFtest_mappingSet::run()
{
    cout << "Running MappingSet tests..." << endl;

    int failCount = 0;

    // Create an empty transformation set
    vector<Transformation*> transformations;

    // Add simple empty transformations to be able to create MappingSet
    transformations.push_back(TransformationFactory::genEmptyTransform());
    transformations.push_back(TransformationFactory::genEmptyTransform());

    // Create mapping set with ID 1
    MappingSet *ms = new MappingSet(1, transformations);

    // Verify initial state
    if (ms->getNumOfMappings() != 2) {
        cerr << "Failed MappingSet initialization, expected 2 mappings, got " << ms->getNumOfMappings() << endl;
        failCount++;
    }

    // Set clonedMappingSet manually for testing the clear function
    MappingSet *msClone = new MappingSet(2, transformations);
    MappingSet::clone(ms, msClone);
    if (ms->getClonedMapping() != msClone->getId()) {
        cerr << "Failed MappingSet clone, expected " << msClone->getId() << ", got " << ms->getClonedMapping() << endl;
        failCount++;
    }

    // Create a mock system and molecule to set on mapping
    System* sys = new System("testSys");
    vector<string> compNames;
    vector<string> defaultStates;
    vector<vector<string>> possibleStates;
    MoleculeType* mt = new MoleculeType("TestMol", compNames, defaultStates, possibleStates, sys);

    // We need to pass mt, listId=0, and compartment=NULL (or appropriate compartment)
    Molecule* mol = new Molecule(mt, 0, NULL);

    // set molecule
    ms->set(0, mol);
    ms->set(1, mol);

    if (ms->get(0)->getMolecule() != mol) {
        cerr << "Failed mapping set to molecule" << endl;
        failCount++;
    }

    // Call clear
    ms->clear();

    // Verify it doesn't crash and clears mappings
    for (unsigned int i = 0; i < ms->getNumOfMappings(); i++) {
        Mapping* m = ms->get(i);
        if (m->getMolecule() != NULL) {
            cerr << "Failed MappingSet::clear, expected NULL molecule, got a molecule pointer" << endl;
            failCount++;
        }
    }

    if (ms->getClonedMapping() != MappingSet::NO_CLONE) {
        cerr << "Failed MappingSet::clear, expected clonedMappingSet to be NO_CLONE, got " << ms->getClonedMapping() << endl;
        failCount++;
    }

    delete mol;
    delete sys; // Deletes molType too
    delete ms;
    delete msClone;

    // Clean up transformations
    for (auto t : transformations) {
        delete t;
    }

    if (failCount == 0) {
        cout << "All MappingSet tests passed successfully!" << endl;
    } else {
        cout << "MappingSet tests failed with " << failCount << " errors." << endl;
        exit(1);
    }
}
