#include "systemSnapshot.hh"
#include "NFcore.hh"

using namespace std;
using namespace NFcore;

SystemSnapshot::SystemSnapshot() : valid(false) {}
SystemSnapshot::~SystemSnapshot() {}

void SystemSnapshot::capture(System *s) {
    complexes.clear();

    auto captureComplexMembers = [this](const list<Molecule*> &members) {
        if (members.empty()) {
            return;
        }

        ComplexSnapshot cs;
        cs.count = 1;

        // Map molecule pointers to indices within this connected component.
        map<Molecule*, int> molIndex;
        int idx = 0;
        for (auto molIter = members.begin(); molIter != members.end(); ++molIter, ++idx) {
            molIndex[*molIter] = idx;
        }

        // Snapshot each molecule.
        for (auto molIter = members.begin(); molIter != members.end(); ++molIter) {
            Molecule *mol = *molIter;
            MoleculeSnapshot ms;
            ms.moleculeTypeName = mol->getMoleculeTypeName();
            ms.compartmentId = mol->getCompartmentId();

            MoleculeType *mt = mol->getMoleculeType();
            int nComp = mt->getNumOfComponents();

            for (int ci = 0; ci < nComp; ci++) {
                ms.componentStates.push_back(mol->getComponentState(ci));

                if (mol->isBindingSiteBonded(ci)) {
                    Molecule *partner = mol->getBondedMolecule(ci);
                    int partnerIdx = molIndex.count(partner) ? molIndex[partner] : -1;
                    int partnerSite = mol->getBondedMoleculeBindingSiteIndex(ci);
                    ms.bondPartners.push_back(partnerIdx);
                    ms.bondPartnerSites.push_back(partnerSite);
                } else {
                    ms.bondPartners.push_back(-1);
                    ms.bondPartnerSites.push_back(-1);
                }
            }
            cs.molecules.push_back(ms);
        }

        complexes.push_back(cs);
    };

    // Traverse every live molecule directly so the snapshot works whether or not
    // the complex list is being maintained.
    set<Molecule*> visited;
    for (unsigned int mtIndex = 0; mtIndex < s->getNumOfMoleculeTypes(); ++mtIndex) {
        MoleculeType *molType = s->getMoleculeType(mtIndex);
        for (int molIndex = 0; molIndex < molType->getMoleculeCount(); ++molIndex) {
            Molecule *mol = molType->getMolecule(molIndex);
            if (!mol->isAlive()) continue;
            if (!visited.insert(mol).second) continue;

            list<Molecule*> members;
            mol->traverseBondedNeighborhood(members, ReactionClass::NO_LIMIT);
            for (auto memberIter = members.begin(); memberIter != members.end(); ++memberIter) {
                visited.insert(*memberIter);
            }
            captureComplexMembers(members);
        }
    }

    valid = true;
}

void SystemSnapshot::restore(System *s) {
    if (!valid) {
        cerr << "Error: no saved concentrations to restore." << endl;
        return;
    }

    // 1. Destroy all existing molecules and their bookkeeping
    s->destroyAllMolecules();

    // 2. Recreate from snapshot
    for (const auto &cs : complexes) {
        for (int copy = 0; copy < cs.count; copy++) {
            // Create molecules
            vector<Molecule*> newMols;
            for (const auto &ms : cs.molecules) {
                MoleculeType *mt = s->getMoleculeTypeByName(ms.moleculeTypeName);
                Compartment *comp = ms.compartmentId.empty() ?
                    nullptr : s->getCompartment(ms.compartmentId);
                Molecule *mol = mt->genDefaultMolecule(comp);

                // Set component states
                for (int ci = 0; ci < (int)ms.componentStates.size(); ci++) {
                    mol->setComponentState(ci, ms.componentStates[ci]);
                }
                newMols.push_back(mol);
            }

            // Form bonds
            for (int mi = 0; mi < (int)cs.molecules.size(); mi++) {
                for (int ci = 0; ci < (int)cs.molecules[mi].bondPartners.size(); ci++) {
                    int partnerIdx = cs.molecules[mi].bondPartners[ci];
                    int partnerSite = cs.molecules[mi].bondPartnerSites[ci];
                    if (partnerIdx > mi) {  // Only bind once per pair
                        Molecule::bind(newMols[mi], ci,
                                      newMols[partnerIdx], partnerSite);
                    }
                }
            }

        }
    }

    // 3. Register all recreated molecules with their applicable reactions
    for (unsigned int mt = 0; mt < s->getNumOfMoleculeTypes(); mt++) {
        MoleculeType *molType = s->getMoleculeType(mt);
        for (int m = 0; m < molType->getMoleculeCount(); m++) {
            molType->getMolecule(m)->updateRxnMembership();
        }
    }

    // 4. Rebuild selector, observables, and propensities
    s->prepareForSimulation();
}