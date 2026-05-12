#include "test_util.hh"
#include <iostream>
#include <stdexcept>
#include "../../NFutil/MTrand/mtrand.h"

using namespace std;
using namespace NFutil;

void NFtest_util::run()
{
	cout << "Running NFutil tests..." << endl;

	const int NUM_ITERATIONS = 100000;

	cout << "  Testing MTRand_int32 initialization error paths..." << endl;

	// Test case: null pointer
	try {
		unsigned long* null_array = nullptr;
		MTRand_int32 rng(null_array, 4);
		throw std::runtime_error("MTRand_int32 did not throw on null array");
	} catch (const std::invalid_argument& e) {
		// Expected
	}

	// Test case: size <= 0
	try {
		unsigned long valid_array[4] = {0x123, 0x234, 0x345, 0x456};
		MTRand_int32 rng(valid_array, 0);
		throw std::runtime_error("MTRand_int32 did not throw on size == 0");
	} catch (const std::invalid_argument& e) {
		// Expected
	}

	try {
		unsigned long valid_array[4] = {0x123, 0x234, 0x345, 0x456};
		MTRand_int32 rng(valid_array, -1);
		throw std::runtime_error("MTRand_int32 did not throw on size == -1");
	} catch (const std::invalid_argument& e) {
		// Expected
	}

	// Test case: valid initialization
	unsigned long valid_array[4] = {0x123, 0x234, 0x345, 0x456};
	MTRand_int32 valid_rng(valid_array, 4); // Should not throw

	cout << "  MTRand_int32 initialization error paths passed!" << endl;

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
	cout << "NFutil tests completed successfully." << endl;
}
