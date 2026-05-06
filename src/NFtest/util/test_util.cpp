#include "test_util.hh"
#include "../../NFcore/moleculeLists/moleculeList.hh"
#include "../../NFcore/NFcore.hh"
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <climits>
#include <string>
#include <vector>

using namespace std;
using namespace NFutil;
using namespace NFcore;

void testTemplateMoleculeConstraint() {
    cout << "  Testing TemplateMolecule constraints..." << endl;

    System *s = new System("test");

    vector <string> compName;
    vector <string> defaultCompState;
    vector < vector <string> > possibleCompStates;

    compName.push_back("p");
    defaultCompState.push_back("0");
    vector <string> pStates;
    pStates.push_back("0");
    pStates.push_back("1");
    possibleCompStates.push_back(pStates);

    MoleculeType *molX = new MoleculeType("X", compName, defaultCompState, possibleCompStates, s);

    TemplateMolecule *xTemp = new TemplateMolecule(molX);
    xTemp->addComponentConstraint("p", "0");

    // Adding the same constraint should be fine (idempotent)
    xTemp->addComponentConstraint("p", "0");

    bool threw_exception = false;
    try {
        xTemp->addComponentConstraint("p", "1");
    } catch (const std::runtime_error& e) {
        // Expected
        cout << "  Caught expected runtime_error: " << e.what() << endl;
        threw_exception = true;
    } catch (...) {
        throw std::logic_error("TemplateMolecule::addComponentConstraint threw unexpected exception type!");
    }

    if (!threw_exception) {
        throw std::logic_error("TemplateMolecule::addComponentConstraint should have thrown an exception for mutually exclusive constraints!");
    }

    // NOTE: delete s deletes the MoleculeType, which deletes all its templates, including xTemp
    // delete xTemp;
    delete s;
    cout << "  TemplateMolecule constraints tests passed!" << endl;
}

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

	cout << "  Testing RANDOM_GAUSSIAN..." << endl;
	double sum = 0;
	double sum_sq = 0;
	for (int i = 0; i < NUM_ITERATIONS; ++i) {
		double result = NFutil::RANDOM_GAUSSIAN();
		sum += result;
		sum_sq += result * result;
	}
	double mean = sum / NUM_ITERATIONS;
	double variance = (sum_sq / NUM_ITERATIONS) - (mean * mean);
	double stddev = std::sqrt(variance);

	if (std::abs(mean) > 0.05) {
		throw std::runtime_error("RANDOM_GAUSSIAN generated a mean far from 0: " + std::to_string(mean));
	}
	if (std::abs(stddev - 1.0) > 0.05) {
		throw std::runtime_error("RANDOM_GAUSSIAN generated a standard deviation far from 1: " + std::to_string(stddev));
	}
	cout << "  RANDOM_GAUSSIAN tests passed!" << endl;

	cout << "  Testing Complex::mergeWithList..." << endl;
	System *s = new System("TestSys");
	vector<string> compNames;
	vector<string> defaultStates;
	vector<vector<string> > possibleStates;
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

	cout << "Running NFcore::MoleculeList tests..." << endl;
	System *sList = new System("test_system");
	vector<string> compName;
	vector<string> defaultCompState;
	MoleculeType *mtList = new MoleculeType("TestMol", compName, defaultCompState, sList);
	MoleculeList *molList = new MoleculeList(mtList, 10, MoleculeList::NO_LIMIT);
	Molecule *m = NULL;
	int listId = 0;
	cout << "  Testing remove out of bounds gracefully..." << endl;
	molList->remove(listId, m);
	cout << "  Testing remove out of bounds passed (didn't crash)." << endl;
	delete molList;
	delete sList;
	cout << "MoleculeList tests completed successfully." << endl;

	testTemplateMoleculeConstraint();
	test_trim();
	test_toString();

	cout << "  Testing loadTimeSeries..." << endl;
	auto create_temp_file = [](const std::string& filename, const std::string& content) {
		std::ofstream ofs(filename);
		ofs << content;
		ofs.close();
	};

	// Test case 1: Valid increasing time series
	create_temp_file("test_ts_valid_inc.txt", "0.0 10.0\n1.0 20.0\n2.0 30.0\n");
	NFutil::TimeSeries ts1 = NFutil::loadTimeSeries("test_ts_valid_inc.txt", "testValidInc");
	if (ts1.time.size() != 3 || ts1.values.size() != 3 || ts1.time[0] != 0.0 || ts1.values[2] != 30.0) {
		throw std::runtime_error("loadTimeSeries failed to parse valid increasing series.");
	}
	remove("test_ts_valid_inc.txt");

	// Test case 2: Valid decreasing time series
	create_temp_file("test_ts_valid_dec.txt", "2.0 30.0\n1.0 20.0\n0.0 10.0\n");
	NFutil::TimeSeries ts2 = NFutil::loadTimeSeries("test_ts_valid_dec.txt", "testValidDec");
	if (ts2.time.size() != 3 || ts2.time[1] != 1.0) {
		throw std::runtime_error("loadTimeSeries failed to parse valid decreasing series.");
	}
	remove("test_ts_valid_dec.txt");

	// Test case 3: Missing/non-existent file
	bool threwMissing = false;
	try {
		NFutil::loadTimeSeries("test_ts_nonexistent.txt", "testMissing");
	} catch (const std::runtime_error& e) {
		threwMissing = true;
		if (std::string(e.what()).find("doesn't look like it exists") == std::string::npos) {
			throw std::runtime_error("loadTimeSeries threw incorrect error for missing file.");
		}
	}
	if (!threwMissing) {
		throw std::runtime_error("loadTimeSeries did not throw on missing file.");
	}

	// Test case 4: Empty file
	create_temp_file("test_ts_empty.txt", "");
	bool threwEmpty = false;
	try {
		NFutil::loadTimeSeries("test_ts_empty.txt", "testEmpty");
	} catch (const std::runtime_error& e) {
		threwEmpty = true;
		if (std::string(e.what()).find("empty or invalid format") == std::string::npos) {
			throw std::runtime_error("loadTimeSeries threw incorrect error for empty file.");
		}
	}
	if (!threwEmpty) {
		throw std::runtime_error("loadTimeSeries did not throw on empty file.");
	}
	remove("test_ts_empty.txt");

	// Test case 5: Duplicate time values
	create_temp_file("test_ts_dup.txt", "0.0 10.0\n1.0 20.0\n1.0 30.0\n");
	bool threwDup = false;
	try {
		NFutil::loadTimeSeries("test_ts_dup.txt", "testDup");
	} catch (const std::runtime_error& e) {
		threwDup = true;
		if (std::string(e.what()).find("Found duplicate time") == std::string::npos) {
			throw std::runtime_error("loadTimeSeries threw incorrect error for duplicate time values.");
		}
	}
	if (!threwDup) {
		throw std::runtime_error("loadTimeSeries did not throw on duplicate time values.");
	}
	remove("test_ts_dup.txt");

	// Test case 6: Non-monotonic time values (increasing then decreasing)
	create_temp_file("test_ts_nonmono.txt", "0.0 10.0\n1.0 20.0\n0.5 30.0\n");
	bool threwNonMono = false;
	try {
		NFutil::loadTimeSeries("test_ts_nonmono.txt", "testNonMono");
	} catch (const std::runtime_error& e) {
		threwNonMono = true;
		if (std::string(e.what()).find("strictly monotonic") == std::string::npos) {
			throw std::runtime_error("loadTimeSeries threw incorrect error for non-monotonic time values.");
		}
	}
	if (!threwNonMono) {
		throw std::runtime_error("loadTimeSeries did not throw on non-monotonic time values.");
	}
	remove("test_ts_nonmono.txt");

	// Test case 7: Invalid number formatting
	create_temp_file("test_ts_invalid_fmt.txt", "0.0 10.0\nabc 20.0\n");
	bool threwInvalidFmt = false;
	try {
		NFutil::loadTimeSeries("test_ts_invalid_fmt.txt", "testInvalidFmt");
	} catch (const std::runtime_error& e) {
		threwInvalidFmt = true;
		if (std::string(e.what()).find("invalid number format") == std::string::npos && std::string(e.what()).find("error in NFutil::convertToDouble") == std::string::npos) {
			throw std::runtime_error("loadTimeSeries threw incorrect error for invalid number formatting: " + std::string(e.what()));
		}
	}
	if (!threwInvalidFmt) {
		throw std::runtime_error("loadTimeSeries did not throw on invalid number formatting.");
	}
	remove("test_ts_invalid_fmt.txt");

	cout << "  loadTimeSeries tests passed!" << endl;

	cout << "NFutil and Core tests completed successfully." << endl;
}
