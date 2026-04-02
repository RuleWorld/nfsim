#include "systemSnapshot.hh"
#include "NFcore.hh"

using namespace std;
using namespace NFcore;

SystemSnapshot::SystemSnapshot() : valid(false) {}
SystemSnapshot::~SystemSnapshot() {}

void SystemSnapshot::capture(System *s) {
    complexes.clear();

    // Iterate all complexes in the system
    s->getAllComplexes().resetComplexIter();
    Complex *c;
    while ((c = s->getAllComplexes().nextComplex()) != nullptr) {
        if (!c->isAlive() || c->getComplexSize() == 0) continue;

        ComplexSnapshot cs;
        cs.count = 1;

        // Map molecule pointers to indices within this complex
        map<Molecule*, int> molIndex;
        int idx = 0;
        for (auto molIter = c->complexMembers.begin();
             molIter != c->complexMembers.end(); ++molIter, ++idx) {
            molIndex[*molIter] = idx;
        }

        // Snapshot each molecule
        for (auto molIter = c->complexMembers.begin();
             molIter != c->complexMembers.end(); ++molIter) {
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
    }

    valid = true;
}

void SystemSnapshot::restore(System *s) {
    if (!valid) {
        cerr << "Error: no saved concentrations to restore." << endl;
        return;
    }

    // 1. Destroy all existing molecules
    // This requires a System method to clear all molecules
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

            // Add the created molecules to the system without full update since we bulk update below
            for (Molecule *mol : newMols) {
                MoleculeType *mt_local = mol->getMoleculeType();
                mt_local->addMoleculeToRunningSystemButDontUpdate(mol);
            }
        }
    }

    // 3. Recalculate all observables
    s->recalculateAllObservables();

    // 4. Update all reaction propensities
    s->updateAllReactionPropensities();
}
