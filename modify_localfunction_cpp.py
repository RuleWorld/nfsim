import re
with open("src/NFfunction/localFunction.cpp", "r") as f:
    content = f.read()

# Constructor updates
content = content.replace("this->isEverEvaluatedOnSpeciesScope=false;", "this->isEverEvaluatedOnSpeciesScope=false;\n\tthis->n_typeImolecules=0;\n\tthis->typeI_mol=new MoleculeType *[system->getNumOfMoleculeTypes()];\n\tthis->typeI_localFunctionIndex=new int[system->getNumOfMoleculeTypes()];")
content = content.replace("this->typeII_localFunctionIndex.push_back(index);", "this->typeII_localFunctionIndex[m]=index;")
content = content.replace("typeII_mol = new MoleculeType * [n_typeIImolecules];", "typeII_mol = new MoleculeType * [n_typeIImolecules];\n\ttypeII_localFunctionIndex = new int[n_typeIImolecules];")

# Destructor updates
content = content.replace("delete [] typeII_mol;", "delete [] typeII_mol;\n\tdelete [] typeII_localFunctionIndex;\n\tdelete [] typeI_mol;\n\tdelete [] typeI_localFunctionIndex;")

# Loop updates
content = content.replace("for(unsigned int ti=0; ti<typeI_mol.size(); ti++)", "for(int ti=0; ti<n_typeImolecules; ti++)")
content = content.replace("for ( unsigned int ti=0; ti<typeI_mol.size(); ti++)", "for(int ti=0; ti<n_typeImolecules; ti++)")
content = content.replace("typeI_mol.at(ti)", "typeI_mol[ti]")
content = content.replace("typeI_localFunctionIndex.at(ti)", "typeI_localFunctionIndex[ti]")
content = content.replace("typeI_mol.at(i)", "typeI_mol[i]")
content = content.replace("typeI_localFunctionIndex.at(i)", "typeI_localFunctionIndex[i]")

# addTypeIMoleculeDependency
content = content.replace("for(unsigned int i=0; i<this->typeI_mol.size(); i++)", "for(int i=0; i<this->n_typeImolecules; i++)")
content = content.replace("this->typeI_mol.push_back(mt);\n\tthis->typeI_localFunctionIndex.push_back(index);", "this->typeI_mol[this->n_typeImolecules]=mt;\n\tthis->typeI_localFunctionIndex[this->n_typeImolecules]=index;\n\tthis->n_typeImolecules++;")

# Exceptions
content = content.replace("lfe.setType1_Mol(&typeI_mol);", "lfe.setType1_Mol(typeI_mol, n_typeImolecules);")

with open("src/NFfunction/localFunction.cpp", "w") as f:
    f.write(content)
