#include "test_NFstream.hh"
#include "../../NFscheduler/NFstream.h"
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <string>

using namespace std;

namespace NFtest_scheduler {

void run_test_NFstream() {
    cout << "  - testing NFstream..." << endl;

    // Test constructors
    NFstream sstream_empty;

    // Default constructor is not necessarily string mode because of MPI, so we'll test explicitly
    // the boolean and string constructors which set it up differently (usually string mode)

    // Test boolean constructor
    NFstream sstream(false);
    sstream << "Hello string stream." << endl;
    sstream << "This is a test of string stream mode.";

    string expectedStr = "Hello string stream.\nThis is a test of string stream mode.";
    if (sstream.str() != expectedStr) {
        cout << "    FAILED: string stream output mismatch." << endl;
        cout << "    Expected: '" << expectedStr << "'" << endl;
        cout << "    Got: '" << sstream.str() << "'" << endl;
        exit(1);
    }

    // Test getStrName and string constructor
    string str_name = "test_string_name";
    NFstream named_sstream(str_name);
    if (named_sstream.getStrName() != "test_string_name") {
        cout << "    FAILED: getStrName returned '" << named_sstream.getStrName() << "' expected 'test_string_name'." << endl;
        exit(1);
    }

    // Test write and flush
    NFstream w_stream(false);
    w_stream.write("test", 4);
    w_stream.flush();
    if (w_stream.str() != "test") {
        cout << "    FAILED: write/flush string stream output mismatch." << endl;
        cout << "    Got: '" << w_stream.str() << "'" << endl;
        exit(1);
    }

    // Test setf and precision
    NFstream f_stream(false);
    f_stream.precision(4);
    f_stream.setf(ios_base::fixed, ios_base::floatfield);
    double test_val = 3.14159265;
    f_stream << test_val;
    if (f_stream.str() != "3.1416") {
        cout << "    FAILED: setf/precision string stream output mismatch." << endl;
        cout << "    Got: '" << f_stream.str() << "'" << endl;
        exit(1);
    }

    // Test is_open
    if (!sstream.is_open()) { // For string mode, is_open always returns true
        cout << "    FAILED: is_open returned false for string stream." << endl;
        exit(1);
    }

    // Test writing to file
    const char* filename = "test_nfstream_output.txt";
    NFstream fstream;
    fstream.open(filename);

    if (!fstream.is_open()) {
        cout << "    FAILED: is_open returned false for opened file stream." << endl;
        exit(1);
    }

    fstream << "Hello file stream." << endl;
    fstream << "This is a test of file stream mode.";

    // Test write and flush for file
    fstream.write("\ntest_write", 11);
    fstream.flush();

    fstream.close();

    // Read the file and check
    ifstream infile(filename);
    if (!infile.is_open()) {
        cout << "    FAILED: Could not open output file for checking." << endl;
        exit(1);
    }

    stringstream buffer;
    buffer << infile.rdbuf();
    string expectedFileStr = "Hello file stream.\nThis is a test of file stream mode.\ntest_write";
    if (buffer.str() != expectedFileStr) {
        cout << "    FAILED: file stream output mismatch." << endl;
        cout << "    Expected: '" << expectedFileStr << "'" << endl;
        cout << "    Got: '" << buffer.str() << "'" << endl;
        exit(1);
    }
    infile.close();

    // Clean up temporary file
    remove(filename);

    // Test filename constructor
    const char* filename2 = "test_nfstream_output2.txt";
    NFstream fstream2(filename2);
    fstream2 << "Hello filename stream.";
    fstream2.close();

    ifstream infile2(filename2);
    if (!infile2.is_open()) {
        cout << "    FAILED: Could not open output file 2 for checking." << endl;
        exit(1);
    }
    stringstream buffer2;
    buffer2 << infile2.rdbuf();
    if (buffer2.str() != "Hello filename stream.") {
        cout << "    FAILED: file stream 2 output mismatch." << endl;
        cout << "    Got: '" << buffer2.str() << "'" << endl;
        exit(1);
    }
    infile2.close();
    remove(filename2);

    cout << "  - run_test_NFstream passed!" << endl;
}

}
