/*
 * NFinput_energy.cpp
 *
 * Integration between the energy rule expansion engine and NFsim's
 * XML input pipeline. Contains:
 *   - parseEnergyPatterns(): reads <ListOfEnergyPatterns> from XML
 *   - createExpandedBindingReactions(): builds BasicRxnClass instances from
 *     expanded energy rules
 *
 * These are called from NFinput.cpp at the appropriate points in
 * the parsing pipeline.
 */

#include "../NFcore/NFcore.hh"
#include "../NFcore/energyPattern.hh"
#include "../NFreactions/reactions/reaction.hh"
#include "NFinput.hh"
#include "TinyXML/tinyxml.h"

#include <map>
#include <string>
#include <vector>
#include <iostream>

using namespace std;
using namespace NFcore;


namespace NFinput {

/*
 * Parse <ListOfEnergyPatterns> from the XML model element.
 */
bool parseEnergyPatterns(
    TiXmlElement *pModel,
    System *s,
    map<string, double> &parameter,
    bool verbose)
{
    TiXmlElement *pList = pModel->FirstChildElement("ListOfEnergyPatterns");
    if (!pList) return true;

    if (verbose) cout << "\n\tReading list of Energy Patterns..." << endl;

    double phi = 0.5;
    double RT = 2.478;

    if (parameter.find("phi") != parameter.end())
        phi = parameter.find("phi")->second;
    if (parameter.find("RT") != parameter.end())
        RT = parameter.find("RT")->second;

    EnergyFunction *ef = new EnergyFunction(phi, RT);

    TiXmlElement *pEP;
    for (pEP = pList->FirstChildElement("EnergyPattern"); pEP != 0; pEP = pEP->NextSiblingElement("EnergyPattern"))
    {
        if (!pEP->Attribute("id") || !pEP->Attribute("expression")) {
            cerr << "Error: EnergyPattern missing 'id' or 'expression'." << endl;
            delete ef;
            return false;
        }

        string epId = pEP->Attribute("id");
        string epExpr = pEP->Attribute("expression");

        double energyVal = 0.0;
        if (parameter.find(epExpr) != parameter.end()) {
            energyVal = parameter.find(epExpr)->second;
        } else {
            try { energyVal = NFutil::convertToDouble(epExpr); }
            catch (...) {
                cerr << "Error: cannot resolve energy '" << epExpr << "' for pattern " << epId << endl;
                delete ef;
                return false;
            }
        }

        EnergyPatternInfo epInfo;
        epInfo.id = epId;
        epInfo.energyValue = energyVal;

        TiXmlElement *pPattern = pEP->FirstChildElement("Pattern");
        TiXmlElement *pListOfMols = pPattern ? pPattern->FirstChildElement("ListOfMolecules") : pEP->FirstChildElement("ListOfMolecules");
        if (!pListOfMols) {
            cerr << "Error: EnergyPattern " << epId << " has no ListOfMolecules." << endl;
            delete ef;
            return false;
        }

        map<string, pair<int,int>> compIdMap;
        TiXmlElement *pMol;
        for (pMol = pListOfMols->FirstChildElement("Molecule"); pMol != 0; pMol = pMol->NextSiblingElement("Molecule"))
        {
            EpMolecule mol;
            mol.xmlId = pMol->Attribute("id") ? pMol->Attribute("id") : "";
            mol.typeName = pMol->Attribute("name") ? pMol->Attribute("name") : "";
            int molIdx = (int)epInfo.molecules.size();

            TiXmlElement *pListOfComps = pMol->FirstChildElement("ListOfComponents");
            if (pListOfComps) {
                TiXmlElement *pComp;
                for (pComp = pListOfComps->FirstChildElement("Component"); pComp != 0; pComp = pComp->NextSiblingElement("Component"))
                {
                    EpMolecule::CompInfo ci;
                    ci.name = pComp->Attribute("name") ? pComp->Attribute("name") : "";
                    string compId = pComp->Attribute("id") ? pComp->Attribute("id") : "";
                    string numBonds = pComp->Attribute("numberOfBonds") ? pComp->Attribute("numberOfBonds") : "0";
                    ci.isBound = (numBonds != "0" && numBonds != "");
                    ci.bondPartnerId = "";

                    if (pComp->Attribute("state")) ci.stateConstraint = pComp->Attribute("state");

                    int compIdx = (int)mol.components.size();
                    compIdMap[compId] = make_pair(molIdx, compIdx);
                    mol.components.push_back(ci);
                }
            }
            epInfo.molecules.push_back(mol);
        }

        TiXmlElement *pListOfBonds = pPattern ? pPattern->FirstChildElement("ListOfBonds") : pEP->FirstChildElement("ListOfBonds");
        if (pListOfBonds) {
            TiXmlElement *pBond;
            for (pBond = pListOfBonds->FirstChildElement("Bond"); pBond != 0; pBond = pBond->NextSiblingElement("Bond"))
            {
                string site1 = pBond->Attribute("site1") ? pBond->Attribute("site1") : "";
                string site2 = pBond->Attribute("site2") ? pBond->Attribute("site2") : "";
                if (compIdMap.count(site1) && compIdMap.count(site2)) {
                    EnergyPatternInfo::Bond bond;
                    bond.mol1 = compIdMap[site1].first; bond.comp1 = compIdMap[site1].second;
                    bond.mol2 = compIdMap[site2].first; bond.comp2 = compIdMap[site2].second;
                    epInfo.bonds.push_back(bond);
                    epInfo.molecules[bond.mol1].components[bond.comp1].bondPartnerId = epInfo.molecules[bond.mol2].xmlId;
                    epInfo.molecules[bond.mol2].components[bond.comp2].bondPartnerId = epInfo.molecules[bond.mol1].xmlId;
                }
            }
        }
        ef->addEnergyPattern(epInfo);
    }

    if (verbose)
        cout << "\n\tParsed " << ef->getNumPatterns() << " energy pattern(s) with RT=" << RT << endl;
    s->setEnergyFunction(ef);
    return true;
}

/*
 * Create expanded BasicRxnClass instances from an energy binding rule.
 */
bool createExpandedBindingReactions(
    const string &rxnName,
    double phi_val,
    double Ea0,
    MoleculeType *molType1, const string &bindSite1,
    MoleculeType *molType2, const string &bindSite2,
    System *s,
    map<string, double> &parameter,
    map<string, int> &allowedStates,
    bool blockSameComplexBinding,
    bool verbose,
    int &reaction_count)
{
    EnergyFunction *ef = s->getEnergyFunction();
    if (!ef) return false;

    string mt1Name = molType1->getName();
    string mt2Name = molType2->getName();

    // Run the expansion algorithm
    vector<ExpandedRuleInfo> expanded = ef->expandBindingRule(
        rxnName, Ea0, phi_val, mt1Name, bindSite1, mt2Name, bindSite2);

    for (const auto &rule : expanded) {
        TemplateMolecule *t1, *t2;

        if (rule.isForward) {
            t1 = new TemplateMolecule(molType1);
            t1->addEmptyComponent(bindSite1);
            t2 = new TemplateMolecule(molType2);
            t2->addEmptyComponent(bindSite2);
        } else {
            t1 = new TemplateMolecule(molType1);
            t2 = new TemplateMolecule(molType2);
            TemplateMolecule::bind(t1, bindSite1, "", t2, bindSite2, "");
        }

        for (const auto &cc : rule.constraints) {
            TemplateMolecule *target = (cc.reactantIdx == 0) ? t1 : t2;
            if (cc.mustBeBound) target->addBoundComponent(cc.compName);
            else target->addEmptyComponent(cc.compName);
        }

        vector<TemplateMolecule *> templates;
        templates.push_back(t1);
        if (rule.isForward) templates.push_back(t2);

        TransformationSet *ts = new TransformationSet(templates);
        if (rule.isForward) ts->addBindingTransform(t1, bindSite1, t2, bindSite2);
        else ts->addUnbindingTransform(t1, bindSite1, t2, bindSite2);

        // Wire complex bookkeeping for blockSameComplexBinding flag
        ts->setComplexBookkeeping(blockSameComplexBinding);

        ts->finalize();

        // Rate is already correctly calculated inside expandBindingRule using phi
        double rate = rule.rate;
        BasicRxnClass *r = new BasicRxnClass(rule.name, rate, "", ts, s);
        
        // Register with system!
        s->addReaction(r);
        reaction_count++;

        if (verbose) {
            cout << "\t  Created " << (rule.isForward ? "forward" : "reverse")
                 << " rule: " << rule.name << "  rate=" << rate << endl;
        }
    }
    return true;
}

} // namespace NFinput
