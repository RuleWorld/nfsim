#include "test_util.hh"
#include <iostream>
#include <stdexcept>
#include <climits>
#include <string>
#include "../../NFcore/NFcore.hh"

using namespace std;
using namespace NFutil;
using namespace NFcore;

void test_trim() {
    cout << "  Testing trim..." << endl;

    // Test empty string
    string s1 = "";
    trim(s1);
    if (s1 != "") throw runtime_error("trim failed on empty string");

    // Test purely whitespace string
    string s2 = "   \t  ";
    trim(s2);
    if (s2 != "") throw runtime_error("trim failed on purely whitespace string, got '" + s2 + "'");

    // Test leading whitespace
    string s3 = "  \t Hello";
    trim(s3);
    if (s3 != "Hello") throw runtime_error("trim failed on leading whitespace, got '" + s3 + "'");

    // Test trailing whitespace
    string s4 = "World \t  ";
    trim(s4);
    if (s4 != "World") throw runtime_error("trim failed on trailing whitespace, got '" + s4 + "'");

    // Test leading and trailing whitespace
    string s5 = " \t NFsim \t ";
    trim(s5);
    if (s5 != "NFsim") throw runtime_error("trim failed on leading and trailing whitespace, got '" + s5 + "'");

    // Test internal whitespace (should not be trimmed)
    string s6 = "  Hello \t World  ";
    trim(s6);
    if (s6 != "Hello \t World") throw runtime_error("trim failed on internal whitespace, got '" + s6 + "'");

    // Test no whitespace
    string s7 = "Testing";
    trim(s7);
    if (s7 != "Testing") throw runtime_error("trim failed on no whitespace, got '" + s7 + "'");

    cout << "  trim tests passed!" << endl;
}

void test_toString() {
    cout << "  Testing toString..." << endl;

    // Test 0
    if (toString(0) != "0") {
        throw runtime_error("Failed toString(0), expected '0', got '" + toString(0) + "'");
    }

    // Test positive integers
    if (toString(42) != "42") {
        throw runtime_error("Failed toString(42), expected '42', got '" + toString(42) + "'");
    }
    if (toString(123456789) != "123456789") {
        throw runtime_error("Failed toString(123456789), got '" + toString(123456789) + "'");
    }

    // Test negative integers
    if (toString(-42) != "-42") {
        throw runtime_error("Failed toString(-42), expected '-42', got '" + toString(-42) + "'");
    }
    if (toString(-123456789) != "-123456789") {
        throw runtime_error("Failed toString(-123456789), got '" + toString(-123456789) + "'");
    }

    // Test INT_MAX and INT_MIN
    string maxStr = std::to_string(INT_MAX);
    if (toString(INT_MAX) != maxStr) {
        throw runtime_error("Failed toString(INT_MAX), expected '" + maxStr + "', got '" + toString(INT_MAX) + "'");
    }

    string minStr = std::to_string(INT_MIN);
    if (toString(INT_MIN) != minStr) {
        throw runtime_error("Failed toString(INT_MIN), expected '" + minStr + "', got '" + toString(INT_MIN) + "'");
    }

    cout << "  toString tests passed!" << endl;
}

void NFtest_util::run()
{
	cout << "Running NFutil and Core tests..." << endl;

	const int NUM_ITERATIONS = 100000;

	cout << "  Testing RANDOM_INT..." << endl;

	// Test case 1: range [0, 10)
	unsigned long min1 = 0;
	unsigned long max1 = 10;
	bool hit_min1 = false;
	bool hit_max_minus_1 = false;
	for (int i = 0; i < NUM_ITERATIONS; ++i) {
		int result = NFutil::RANDOM_INT(min1, max1);
		if (result < (int)min1 || result >= (int)max1) {
			throw std::runtime_error("RANDOM_INT generated a number out of range: " + std::to_string(result));
		}
		if (result == (int)min1) hit_min1 = true;
		if (result == (int)max1 - 1) hit_max_minus_1 = true;
	}
	if (!hit_min1) throw std::runtime_error("RANDOM_INT did not hit the minimum value");
	if (!hit_max_minus_1) throw std::runtime_error("RANDOM_INT did not hit the maximum - 1 value");

	// Test case 2: range [5, 20)
	unsigned long min2 = 5;
	unsigned long max2 = 20;
	bool hit_min2 = false;
	bool hit_max_minus_2 = false;
	for (int i = 0; i < NUM_ITERATIONS; ++i) {
		int result = NFutil::RANDOM_INT(min2, max2);
		if (result < (int)min2 || result >= (int)max2) {
			throw std::runtime_error("RANDOM_INT generated a number out of range: " + std::to_string(result));
		}
		if (result == (int)min2) hit_min2 = true;
		if (result == (int)max2 - 1) hit_max_minus_2 = true;
	}
	if (!hit_min2) throw std::runtime_error("RANDOM_INT did not hit the minimum value");
	if (!hit_max_minus_2) throw std::runtime_error("RANDOM_INT did not hit the maximum - 1 value");

	// Test case 3: very tight range [5, 6)
	unsigned long min3 = 5;
	unsigned long max3 = 6;
	for (int i = 0; i < NUM_ITERATIONS; ++i) {
		int result = NFutil::RANDOM_INT(min3, max3);
		if (result != 5) {
			throw std::runtime_error("RANDOM_INT generated a number out of range for tight bound: " + std::to_string(result));
		}
	}

	cout << "  RANDOM_INT tests passed!" << endl;

    cout << "  Testing Complex::mergeWithList..." << endl;
    System *s = new System("TestSys");
    vector<string> compNames;
    vector<string> defaultStates;
    vector<vector<string>> possibleStates;
    vector<bool> isb;

    MoleculeType *mt = new MoleculeType("A", compNames, defaultStates, possibleStates, isb, false, s);
    Molecule *m1 = new Molecule(mt, 0, NULL);
    Molecule *m2 = new Molecule(mt, 0, NULL);

    Complex *c1 = new Complex(s, 0, m1);
    Complex *c2 = new Complex(s, 1, m2);

    // self-merge
    c1->mergeWithList(c1);
    if (c1->getComplexSize() != 1) {
        throw runtime_error("Self merge changed complex size unexpectedly.");
    }

    // merge two different complexes
    c1->mergeWithList(c2);
    if (c1->getComplexSize() != 2) {
        throw runtime_error("Merge with other list failed, size is " + to_string(c1->getComplexSize()));
    }

    cout << "  Complex mergeWithList tests passed!" << endl;

    delete s;

    test_trim();
    test_toString();

	cout << "NFutil and Core tests completed successfully." << endl;
}
