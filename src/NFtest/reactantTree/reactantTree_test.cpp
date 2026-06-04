#include "reactantTree_test.hh"
#include "../../NFreactions/reactantLists/reactantTree.hh"
#include "../../NFreactions/transformations/transformationSet.hh"
#include <iostream>
#include <vector>

#ifndef _WIN32
#include <unistd.h>
#include <sys/wait.h>
#endif

using namespace std;
using namespace NFcore;

void NFtest_reactantTree::run()
{
    cout << "Running ReactantTree tests..." << endl;

    int failCount = 0;

    // Test for exit(1) on removeMappingSet from empty ReactantTree
    cout << "  Testing removeMappingSet on empty tree (expecting exit(1))..." << endl;

#ifndef _WIN32
    // We will fork a process to test that it calls exit(1)
    pid_t pid = fork();
    if (pid == 0) {
        // In the child process
        // Redirect cerr so we don't spam the console if not needed, but here it's expected
        if (freopen("/dev/null", "w", stderr) == nullptr) {
            // Ignore if freopen fails
        }

        // Create an empty transformation set
        vector<TemplateMolecule*> tempMols;
        TransformationSet ts(tempMols);
        ts.finalize(); // Finalize to prevent "TransformationSet cannot generate blank mapping if it is not finalized!"

        // Create ReactantTree
        ReactantTree* tree = new ReactantTree(0, &ts, 10);

        // Call the method that should exit
        tree->removeMappingSet(123);

        // If we get here, the test failed (exit was not called)
        exit(0); // Return 0 to indicate failure of the test
    } else if (pid > 0) {
        // In the parent process
        int status;
        waitpid(pid, &status, 0);

        if (WIFEXITED(status)) {
            int exit_status = WEXITSTATUS(status);
            if (exit_status == 1) {
                cout << "    Success: Empty tree removeMappingSet exited with code 1." << endl;
            } else {
                cout << "    Failure: Child process exited with code " << exit_status << " instead of 1." << endl;
                failCount++;
            }
        } else {
            cout << "    Failure: Child process did not exit normally." << endl;
            failCount++;
        }
    } else {
        cerr << "Fork failed!" << endl;
        failCount++;
    }
#else
    cout << "    Skipping exit(1) test on Windows as fork() is not available." << endl;
#endif

    if (failCount == 0) {
        cout << "All ReactantTree tests passed successfully!" << endl;
    } else {
        cout << "ReactantTree tests failed with " << failCount << " errors." << endl;
        exit(1);
    }
}
