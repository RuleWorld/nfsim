#include "test_compositeFunction.hh"
#include "../../NFfunction/NFfunction.hh"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

using namespace std;
using namespace NFcore;

void NFtest_compositeFunction::run() {
    cout << "Running CompositeFunction tests..." << endl;

    System *s = new System("test_sys");

    vector<string> functions;
    vector<string> argNames;
    vector<string> paramNames;

    bool anyFailed = false;

    // Test case: Valid reactant reference
    {
        bool threw = false;
        try {
            CompositeFunction *cf = new CompositeFunction(s, "f1", "reactant_1()", functions, argNames, paramNames);
            cf->finalizeInitialization(s);
            // Note: intentionally leaking cf here on success to avoid complex cascade deletion bugs in mock Systems
        } catch (std::runtime_error& e) {
            threw = true;
        } catch (...) {
            threw = true;
        }
        if (threw) {
            cout << "  FAILED: Valid reactant reference threw an exception." << endl;
            anyFailed = true;
        }
    }

    // Redirect cerr to suppress expected error output
    streambuf *old_cerr = cerr.rdbuf();
    stringstream ss;
    cerr.rdbuf(ss.rdbuf());

    // Test case: Invalid reactant reference (no number)
    {
        bool threw = false;
        try {
            CompositeFunction *cf = new CompositeFunction(s, "f2", "reactant_X()", functions, argNames, paramNames);
            cf->finalizeInitialization(s);
            // Note: intentionally leaking cf here on success to avoid complex cascade deletion bugs in mock Systems
        } catch (std::runtime_error& e) {
            threw = true;
        } catch (...) {
            threw = true;
        }
        if (!threw) {
            cerr.rdbuf(old_cerr);
            cout << "  FAILED: Invalid reactant reference did not throw an exception." << endl;
            anyFailed = true;
            cerr.rdbuf(ss.rdbuf());
        }
    }

    // Test case: Invalid reactant reference (two digit number)
    {
        bool threw = false;
        try {
            CompositeFunction *cf = new CompositeFunction(s, "f3", "reactant_10()", functions, argNames, paramNames);
            cf->finalizeInitialization(s);
            // Note: intentionally leaking cf here on success to avoid complex cascade deletion bugs in mock Systems
        } catch (std::runtime_error& e) {
            threw = true;
        } catch (...) {
            threw = true;
        }
        if (!threw) {
            cerr.rdbuf(old_cerr);
            cout << "  FAILED: Two digit reactant reference did not throw an exception." << endl;
            anyFailed = true;
            cerr.rdbuf(ss.rdbuf());
        }
    }

    // Restore cerr
    cerr.rdbuf(old_cerr);

    // Note: intentionally omitting delete s; to avoid cascade teardown segfaults with partial state

    if (anyFailed) {
        throw std::runtime_error("CompositeFunction tests failed!");
    }

    cout << "  CompositeFunction tests passed!" << endl;
}
