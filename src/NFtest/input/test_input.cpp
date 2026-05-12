#include "test_input.hh"
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

// Declare the external functions we want to test in the global namespace
// as they are defined in walk.cpp without a namespace.
int getInput(int min, int max);
double getInput(double min);

void NFtest_input::run() {
    cout << "Running NFinput tests..." << endl;

    // Save old buffers
    streambuf* oldCin = cin.rdbuf();
    streambuf* oldCout = cout.rdbuf();

    try {
        // Test getInput(double min)
        // 1. Valid input
        {
            stringstream in("3.14\n");
            stringstream out;
            cin.rdbuf(in.rdbuf());
            cout.rdbuf(out.rdbuf());

            double result = getInput(0.0);
            if (result != 3.14) {
                throw runtime_error("Expected 3.14, got " + to_string(result));
            }
        }

        // 2. Invalid input (runtime_error on convertToDouble) followed by valid input
        {
            stringstream in("abc\n2.71\n");
            stringstream out;
            cin.rdbuf(in.rdbuf());
            cout.rdbuf(out.rdbuf());

            double result = getInput(0.0);
            if (result != 2.71) {
                throw runtime_error("Expected 2.71 after error, got " + to_string(result));
            }

            string output = out.str();
            if (output.find("not an option") == string::npos) {
                throw runtime_error("Expected error message in output, got: " + output);
            }
        }

        // 3. Below minimum followed by valid input
        {
            stringstream in("-1.0\n5.0\n");
            stringstream out;
            cin.rdbuf(in.rdbuf());
            cout.rdbuf(out.rdbuf());

            double result = getInput(0.0);
            if (result != 5.0) {
                throw runtime_error("Expected 5.0 after out of bounds, got " + to_string(result));
            }

            string output = out.str();
            if (output.find("not an option") == string::npos) {
                throw runtime_error("Expected error message in output, got: " + output);
            }
        }

        // Test getInput(int min, int max)
        // 1. Valid input
        {
            stringstream in("3\n");
            stringstream out;
            cin.rdbuf(in.rdbuf());
            cout.rdbuf(out.rdbuf());

            int result = getInput(0, 10);
            if (result != 3) {
                throw runtime_error("Expected 3, got " + to_string(result));
            }
        }

        // 2. Invalid input (runtime_error on convertToInt) followed by valid input
        {
            stringstream in("abc\n7\n");
            stringstream out;
            cin.rdbuf(in.rdbuf());
            cout.rdbuf(out.rdbuf());

            int result = getInput(0, 10);
            if (result != 7) {
                throw runtime_error("Expected 7 after error, got " + to_string(result));
            }

            string output = out.str();
            if (output.find("not an option") == string::npos) {
                throw runtime_error("Expected error message in output, got: " + output);
            }
        }

        // 3. Below minimum followed by valid input
        {
            stringstream in("-1\n5\n");
            stringstream out;
            cin.rdbuf(in.rdbuf());
            cout.rdbuf(out.rdbuf());

            int result = getInput(0, 10);
            if (result != 5) {
                throw runtime_error("Expected 5 after out of bounds, got " + to_string(result));
            }

            string output = out.str();
            if (output.find("not an option") == string::npos) {
                throw runtime_error("Expected error message in output, got: " + output);
            }
        }

        // 4. Above maximum followed by valid input
        {
            stringstream in("15\n8\n");
            stringstream out;
            cin.rdbuf(in.rdbuf());
            cout.rdbuf(out.rdbuf());

            int result = getInput(0, 10);
            if (result != 8) {
                throw runtime_error("Expected 8 after out of bounds, got " + to_string(result));
            }

            string output = out.str();
            if (output.find("not an option") == string::npos) {
                throw runtime_error("Expected error message in output, got: " + output);
            }
        }

    } catch (...) {
        // Restore buffers before throwing
        cin.rdbuf(oldCin);
        cout.rdbuf(oldCout);
        throw;
    }

    // Restore buffers
    cin.rdbuf(oldCin);
    cout.rdbuf(oldCout);

    // Explicitly write the passed message to the real cout
    cout << "NFinput tests passed!" << endl;
}
