import re
with open("src/NFcore/NFcore.hh", "r") as f:
    content = f.read()

content = content.replace("void setType1_Mol(vector <MoleculeType *>* type1_Mol){", "void setType1_Mol(MoleculeType ** type1_Mol, int n_type1_Mol){\n\t\t\tthis->n_type1_Mol = n_type1_Mol;")
content = content.replace("vector<MoleculeType*>* getType1_Mol() const{", "MoleculeType** getType1_Mol() const{")
content = content.replace("vector<MoleculeType*>* type1_Mol;", "MoleculeType** type1_Mol;\n\t\tint n_type1_Mol;")
content = content.replace("int getIndex() const{\n\t\t\treturn index;\n\t\t}", "int getIndex() const{\n\t\t\treturn index;\n\t\t}\n\n\t\tint get_n_type1_Mol() const{\n\t\t\treturn n_type1_Mol;\n\t\t}")

with open("src/NFcore/NFcore.hh", "w") as f:
    f.write(content)
