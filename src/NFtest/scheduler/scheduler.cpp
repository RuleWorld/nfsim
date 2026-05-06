#include "scheduler.hh"
#include "../../NFscheduler/Scheduler.h"
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>

using namespace std;

namespace NFtest_scheduler {

void testConvertModelScanToJobs() {
    cout << "  - testing convertModelScanToJobs..." << endl;

    // Create a model with replicates
    model* m = new model();
    m->filename = "test_model.xml";
    m->processors = 2;
    m->replicates = 3;
    m->argument.push_back("-test");
    m->argval.push_back("val1");

    // Create a scan
    scan* s = new scan();
    s->parameter.push_back("param1");
    s->min.push_back(1.0);
    s->max.push_back(5.0);
    s->steps.push_back(3); // Should generate 1.0, 3.0, 5.0

    s->parameter.push_back("param2");
    s->min.push_back(10.0);
    s->max.push_back(20.0);
    s->steps.push_back(2); // Should generate 10.0, 20.0

    vector<job*> joblist;
    convertModelScanToJobs(m, s, joblist);

    // After function execution, m should be set to NULL
    if (m != NULL) {
        cout << "    FAILED: currentModel was not set to NULL" << endl;
        exit(1);
    }

    // Number of jobs = replicates * steps[0] * steps[1] = 3 * 3 * 2 = 18
    if (joblist.size() != 18) {
        cout << "    FAILED: expected 18 jobs, got " << joblist.size() << endl;
        exit(1);
    }

    // Verify properties of the first job
    job* firstJob = joblist[0];
    if (firstJob->filename != "test_model.xml") {
        cout << "    FAILED: filename not copied correctly" << endl;
        exit(1);
    }
    if (firstJob->processors != 2) {
        cout << "    FAILED: processors not copied correctly" << endl;
        exit(1);
    }
    if (firstJob->argument.size() != 1 || firstJob->argument[0] != "-test") {
        cout << "    FAILED: argument not copied correctly" << endl;
        exit(1);
    }
    if (firstJob->argval.size() != 1 || firstJob->argval[0] != "val1") {
        cout << "    FAILED: argval not copied correctly" << endl;
        exit(1);
    }
    if (firstJob->parameters.size() != 2) {
        cout << "    FAILED: parameters size is incorrect" << endl;
        exit(1);
    }
    if (firstJob->values.size() != 2) {
        cout << "    FAILED: values size is incorrect" << endl;
        exit(1);
    }

    // Clean up
    for (size_t i = 0; i < joblist.size(); i++) {
        delete joblist[i];
    }
    delete s;

    cout << "  - testConvertModelScanToJobs passed!" << endl;
}

void testNFstream() {
    cout << "  - testing NFstream..." << endl;

    // Test writing to string
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

    // Test writing to file
    const char* filename = "test_nfstream_output.txt";
    NFstream fstream;
    fstream.open(filename);
    fstream << "Hello file stream." << endl;
    fstream << "This is a test of file stream mode.";
    fstream.close();

    // Read the file and check
    ifstream infile(filename);
    if (!infile.is_open()) {
        cout << "    FAILED: Could not open output file for checking." << endl;
        exit(1);
    }

    stringstream buffer;
    buffer << infile.rdbuf();
    string expectedFileStr = "Hello file stream.\nThis is a test of file stream mode.";
    if (buffer.str() != expectedFileStr) {
        cout << "    FAILED: file stream output mismatch." << endl;
        cout << "    Expected: '" << expectedFileStr << "'" << endl;
        cout << "    Got: '" << buffer.str() << "'" << endl;
        exit(1);
    }
    infile.close();

    // Clean up temporary file
    remove(filename);

    cout << "  - testNFstream passed!" << endl;
}

void run() {
    cout << "Running NFtest_scheduler..." << endl;
    testConvertModelScanToJobs();
    testNFstream();
    cout << "NFtest_scheduler complete!" << endl;
}

}
