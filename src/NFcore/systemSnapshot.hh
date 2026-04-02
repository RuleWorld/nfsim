#ifndef SYSTEMSNAPSHOT_HH_
#define SYSTEMSNAPSHOT_HH_

#include <string>
#include <vector>
#include <map>

using namespace std;

namespace NFcore {
    class System;
    class Molecule;
    class MoleculeType;
    class Compartment;

    // Stores the state of a single molecule
    struct MoleculeSnapshot {
        string moleculeTypeName;
        vector<int> componentStates;     // State index for each component
        vector<int> bondPartners;        // -1 for unbound, index into complex's molecule list
        vector<int> bondPartnerSites;    // Which site the bond partner is on
        string compartmentId;            // "" if no compartment
    };

    // Stores the state of a single complex (connected species)
    struct ComplexSnapshot {
        vector<MoleculeSnapshot> molecules;
        int count;  // How many copies of this complex exist
    };

    // Stores the full system state
    class SystemSnapshot {
    public:
        SystemSnapshot();
        ~SystemSnapshot();

        void capture(System *s);
        void restore(System *s);
        bool isValid() const { return valid; }

    private:
        bool valid;
        vector<ComplexSnapshot> complexes;
        // Also store parameter values for completeness
        map<string, double> parameterValues;
    };
}

#endif
