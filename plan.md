1. **Refactor `LocalFunction` members in `src/NFfunction/NFfunction.hh`**:
   - Change `vector <MoleculeType *> typeI_mol;` to `MoleculeType ** typeI_mol;`.
   - Change `vector <int> typeI_localFunctionIndex;` to `int * typeI_localFunctionIndex;`.
   - Add `int n_typeImolecules;` to track the number of elements (with a max capacity bounded by total molecule types).
   - Change `vector <int> typeII_localFunctionIndex;` to `int * typeII_localFunctionIndex;`.

2. **Update `LocalFunction` initialization in `src/NFfunction/localFunction.cpp`**:
   - In the constructor, initialize `n_typeImolecules = 0`.
   - Allocate `typeI_mol` and `typeI_localFunctionIndex` as arrays with max capacity `s->getNumOfMoleculeTypes()`.
   - Update type II arrays initialization to use `new int[n_typeIImolecules]` instead of `push_back`.
   - Update the destructor to `delete [] typeI_mol`, `typeI_localFunctionIndex`, and `typeII_localFunctionIndex`.

3. **Update `LocalFunction` usages in `src/NFfunction/localFunction.cpp`**:
   - Replace loops `for(unsigned int ti=0; ti<typeI_mol.size(); ti++)` with `for(unsigned int ti=0; ti<n_typeImolecules; ti++)`.
   - Replace `.at(ti)` with array index `[ti]`.
   - Update `addTypeIMoleculeDependency` to assign values to the array and increment `n_typeImolecules`.

4. **Update `LocalFunctionException` in `src/NFcore/NFcore.hh`**:
   - Change `vector <MoleculeType *> * type1_Mol;` to `MoleculeType ** type1_Mol;` and add `int n_type1_Mol;`.
   - Update the getters and setters accordingly (`setType1_Mol(MoleculeType **m, int count)`).

5. **Update `DORRxnClass` usages in `src/NFreactions/reactions/DORreaction.cpp` and `.hh`**:
   - Change `pickLocalFunctionParameter(MappingSet* ms, int index, vector <MoleculeType *>* type1_Mol, int* reactantCounts)` to use `MoleculeType ** type1_Mol, int n_type1_Mol`.
   - Update `auto it: *(type1_Mol)` logic to loop through the array.

6. **Verify Build and Tests**:
   - Build with `cmake .. && make` in `build/` directory.
   - Run tests with `./build/NFsim -test util` and other test targets to make sure there are no regressions.
