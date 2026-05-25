#include "test_observable.hh"
#include "../../NFcore/NFcore.hh"
#include "../../NFfunction/NFfunction.hh"
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace NFcore;

// Create a concrete implementation of Observable for testing
class DummyObservable : public Observable {
public:
    DummyObservable(string name) : Observable(name) {}
    ~DummyObservable() {}
    Observable* clone() { return new DummyObservable(obsName + "_clone"); }
    int isObservable(Molecule *m) const { return 0; }
    int isObservable(Complex *c) const { return 0; }
};

void NFtest_observable::run()
{
	cout << "Running Observable tests..." << endl;

	// Create a dummy System
	System* sys = new System("TestSystem");

	// Create an Observable
	DummyObservable* obs = new DummyObservable("TestObs");

	vector<string> varRefNames;
	vector<string> varRefTypes;
	vector<string> paramNames;

	// Create a dummy GlobalFunction
	GlobalFunction* gf = new GlobalFunction("TestFunc", "1", varRefNames, varRefTypes, paramNames, sys);

	// Call the untested method
	obs->addReferenceToGlobalFunction(gf);

    // Verify
    if (gf->getCtrType() != "Observable") {
        throw std::runtime_error("addReferenceToGlobalFunction failed to set counter type correctly on GlobalFunction.");
    }

    // Test that the pointer matches the observable's count by changing the count and checking the function
    obs->straightAdd(); // Increases count by 1
    if (gf->getCounterValue() != 1.0) {
        throw std::runtime_error("addReferenceToGlobalFunction failed to set pointer correctly.");
    }

	delete gf;
	delete obs;
	delete sys;
	cout << "Observable tests completed successfully." << endl;
}
