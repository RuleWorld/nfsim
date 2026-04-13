#include "util_test.hh"
#include <iostream>
#include <climits>
#include <string>

using namespace std;
using namespace NFutil;

void NFtest_util::run()
{
    cout << "Running NFutil tests" << endl;

    int failCount = 0;

    // Test 0
    if (toString(0) != "0") {
        cerr << "Failed toString(0), expected '0', got '" << toString(0) << "'" << endl;
        failCount++;
    }

    // Test positive integers
    if (toString(42) != "42") {
        cerr << "Failed toString(42), expected '42', got '" << toString(42) << "'" << endl;
        failCount++;
    }
    if (toString(123456789) != "123456789") {
        cerr << "Failed toString(123456789)" << endl;
        failCount++;
    }

    // Test negative integers
    if (toString(-42) != "-42") {
        cerr << "Failed toString(-42), expected '-42', got '" << toString(-42) << "'" << endl;
        failCount++;
    }
    if (toString(-123456789) != "-123456789") {
        cerr << "Failed toString(-123456789)" << endl;
        failCount++;
    }

    // Test INT_MAX and INT_MIN
    string maxStr = std::to_string(INT_MAX);
    if (toString(INT_MAX) != maxStr) {
        cerr << "Failed toString(INT_MAX), expected '" << maxStr << "', got '" << toString(INT_MAX) << "'" << endl;
        failCount++;
    }

    string minStr = std::to_string(INT_MIN);
    if (toString(INT_MIN) != minStr) {
        cerr << "Failed toString(INT_MIN), expected '" << minStr << "', got '" << toString(INT_MIN) << "'" << endl;
        failCount++;
    }

    if (failCount == 0) {
        cout << "All NFutil tests passed successfully!" << endl;
    } else {
        cout << "NFutil tests failed with " << failCount << " errors." << endl;
        exit(1);
    }
}
