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
    Molecule* mol1 = new Molecule(mt, 0, NULL);
    Molecule* mol2 = new Molecule(mt, 1, NULL);
    Molecule* mol3 = new Molecule(mt, 2, NULL);

    // set molecule
    ms->set(0, mol1);
    ms->set(1, mol2);

    if (ms->get(0)->getMolecule() != mol1) {
        cerr << "Failed mapping set to molecule" << endl;
        failCount++;
    }

    // Test MappingSet::checkForCollisions
    MappingSet *ms2 = new MappingSet(3, transformations);

    // Setup for NO collision test
    ms->set(0, mol1);
    ms->set(1, mol2);
    ms2->set(0, mol3);
    ms2->set(1, mol3);

    if (MappingSet::checkForCollisions(ms, ms2)) {
        cerr << "Failed MappingSet::checkForCollisions: incorrectly detected a collision when there was none" << endl;
        failCount++;
    }

    // Setup for collision test (both have mol2)
    ms2->set(0, mol3);
    ms2->set(1, mol2);

    if (!MappingSet::checkForCollisions(ms, ms2)) {
        cerr << "Failed MappingSet::checkForCollisions: failed to detect a collision" << endl;
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

    delete mol1;
    delete mol2;
    delete mol3;
    delete sys; // Deletes molType too
    delete ms;
    delete ms2;
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
