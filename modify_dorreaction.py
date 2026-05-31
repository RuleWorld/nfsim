import re
with open("src/NFreactions/reactions/DORreaction.cpp", "r") as f:
    content = f.read()

content = content.replace("double DORRxnClass::pickLocalFunctionParameter(MappingSet* ms, int index, vector <MoleculeType *>* type1_Mol, int* reactantCounts)", "double DORRxnClass::pickLocalFunctionParameter(MappingSet* ms, int index, MoleculeType ** type1_Mol, int n_type1_Mol, int* reactantCounts)")
content = content.replace("for (auto it: *(type1_Mol)){", "for(int type1_i=0; type1_i<n_type1_Mol; type1_i++){\n\t\t\t\tauto it = type1_Mol[type1_i];")
content = content.replace("return this->pickLocalFunctionParameter(ms, lfe.getIndex(), lfe.getType1_Mol(), reactantCounts);", "return this->pickLocalFunctionParameter(ms, lfe.getIndex(), lfe.getType1_Mol(), lfe.get_n_type1_Mol(), reactantCounts);")

with open("src/NFreactions/reactions/DORreaction.cpp", "w") as f:
    f.write(content)

with open("src/NFreactions/reactions/reaction.hh", "r") as f:
    content = f.read()

content = content.replace("virtual double pickLocalFunctionParameter(MappingSet *ms, int, vector <MoleculeType *>*, int*);", "virtual double pickLocalFunctionParameter(MappingSet *ms, int, MoleculeType **, int, int*);")

with open("src/NFreactions/reactions/reaction.hh", "w") as f:
    f.write(content)
