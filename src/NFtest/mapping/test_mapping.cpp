#include "test_mapping.hh"
#include "../../NFreactions/mappings/mapping.hh"
#include "../../NFreactions/mappings/mappingSet.hh"
#include "../../NFreactions/transformations/transformation.hh"
#include "../../NFcore/NFcore.hh"
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;
using namespace NFcore;

void test_Mapping_clone() {
    cout << "  Testing Mapping::clone..." << endl;

    System *s = new System("test_sys");
    vector<string> compNames;
    compNames.push_back("site");
    MoleculeType *mt = new MoleculeType("A", compNames, s);
    Molecule *m1 = new Molecule(mt, 0, nullptr);
    Molecule *m2 = new Molecule(mt, 1, nullptr);

    Mapping *originalMapping = new Mapping(TransformationFactory::STATE_CHANGE, 0);
    originalMapping->setMolecule(m1);

    Mapping *clonedMapping = new Mapping(TransformationFactory::STATE_CHANGE, 0);
    clonedMapping->setMolecule(m2);

    Mapping::clone(originalMapping, clonedMapping);

    if (clonedMapping->getMolecule() != m1) {
        throw runtime_error("Mapping::clone failed to clone molecule.");
    }
    if (clonedMapping->getType() != TransformationFactory::STATE_CHANGE) {
        throw runtime_error("Mapping::clone failed to preserve type.");
    }
    if (clonedMapping->getIndex() != 0) {
        throw runtime_error("Mapping::clone failed to preserve index.");
    }

    Molecule *m3 = new Molecule(mt, 2, nullptr);
    originalMapping->setMolecule(m3);

    if (clonedMapping->getMolecule() != m1) {
        throw runtime_error("Mapping::clone original changed clone unexpectedly.");
    }

    delete originalMapping;
    delete clonedMapping;

    delete s; // This destroys mt and all molecules associated with mt

    cout << "  Mapping::clone tests passed!" << endl;
}

void test_MappingSet_clone() {
    cout << "  Testing MappingSet::clone..." << endl;

    System *s = new System("test_sys2");
    vector<string> compNames;
    compNames.push_back("site");
    MoleculeType *mt = new MoleculeType("A", compNames, s);
    Molecule *m1 = new Molecule(mt, 0, nullptr);
    Molecule *m2 = new Molecule(mt, 1, nullptr);

    vector<Transformation *> transformations;
    transformations.push_back(TransformationFactory::genStateChangeTransform(0, 1));

    MappingSet *originalSet = new MappingSet(1, transformations);
    originalSet->set(0, m1);

    MappingSet *clonedSet = new MappingSet(2, transformations);
    clonedSet->set(0, m2);

    MappingSet::clone(originalSet, clonedSet);

    if (clonedSet->get(0)->getMolecule() != m1) {
        throw runtime_error("MappingSet::clone failed to clone molecule.");
    }

    if (originalSet->getClonedMapping() != 2) {
        throw runtime_error("MappingSet::clone failed to set clone ID.");
    }

    // Change original, check if clone is unaffected pointer-wise
    Molecule *m3 = new Molecule(mt, 2, nullptr);
    originalSet->set(0, m3);

    if (clonedSet->get(0)->getMolecule() != m1) {
         throw runtime_error("MappingSet::clone original modification affected clone.");
    }

    delete originalSet;
    delete clonedSet;

    for (Transformation *t : transformations) delete t;
    delete s;

    cout << "  MappingSet::clone tests passed!" << endl;
}

void NFtest_mapping::run()
{
	cout << "Running mapping tests..." << endl;

    test_Mapping_clone();
    test_MappingSet_clone();

	cout << "Mapping tests completed successfully." << endl;
}
