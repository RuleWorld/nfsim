/*
 * NFinput_energy.cpp
 *
 * Integration between the energy rule expansion engine and NFsim's
 * XML input pipeline. Contains:
 *   - parseEnergyPatterns(): reads <ListOfEnergyPatterns> from XML
 *   - createExpandedBindingReactions(): builds compact mapping-local or
 *     materialized reaction classes from energy rules
 *
 * These are called from NFinput.cpp at the appropriate points in
 * the parsing pipeline.
 */

#include "../NFcore/NFcore.hh"
#include "../NFcore/energyPattern.hh"
#include "../NFfunction/NFfunction.hh"
#include "../NFreactions/reactions/reaction.hh"
#include "NFinput.hh"
#include "TinyXML/tinyxml.h"

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;
using namespace NFcore;


namespace NFinput {

static double evaluateEnergyExpression(const string &expression,
                                       map<string, double> &parameters)
{
    mu::Parser parser;
    for (auto &parameter : parameters)
        parser.DefineVar(parameter.first, &parameter.second);
    double value;
    try {
        parser.SetExpr(expression);
        value = parser.Eval();
    } catch (const mu::ParserError &error) {
        throw std::runtime_error(error.GetMsg());
    }
    if (!std::isfinite(value))
        throw std::runtime_error("energy expression is not finite");
    return value;
}

static void parseEnergyMolecules(TiXmlElement *pListOfMols, EnergyPatternInfo &epInfo, map<string, pair<int,int>> &compIdMap)
{
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
}

static void parseEnergyBonds(TiXmlElement *pListOfBonds, EnergyPatternInfo &epInfo, map<string, pair<int,int>> &compIdMap)
{
    if (!pListOfBonds) return;
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
    // BioNetGen stores energies in RT units and applies Arrhenius exponents
    // without a second RT divisor. A model parameter named RT is an ordinary
    // parameter used by expressions such as Gf/RT; it does not change this
    // convention.
    const double RT = 1.0;

    if (parameter.find("phi") != parameter.end())
        phi = parameter.find("phi")->second;
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
        try { energyVal = evaluateEnergyExpression(epExpr, parameter); }
        catch (const std::exception &error) {
            cerr << "Error: cannot evaluate energy expression '" << epExpr
                 << "' for pattern " << epId << ": " << error.what() << endl;
            delete ef;
            return false;
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
        parseEnergyMolecules(pListOfMols, epInfo, compIdMap);

        TiXmlElement *pListOfBonds = pPattern ? pPattern->FirstChildElement("ListOfBonds") : pEP->FirstChildElement("ListOfBonds");
        parseEnergyBonds(pListOfBonds, epInfo, compIdMap);

        ef->addEnergyPattern(epInfo);
    }

    if (verbose)
        cout << "\n\tParsed " << ef->getNumPatterns() << " energy pattern(s) with RT=" << RT << endl;
    s->setEnergyFunction(ef);
    return true;
}

/*
 * Create compact or materialized reaction classes from an energy binding
 * rule.
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
    int &reaction_count,
    bool includeForward,
    bool includeReverse)
{
    EnergyFunction *ef = s->getEnergyFunction();
    if (!ef) return false;

    string mt1Name = molType1->getName();
    string mt2Name = molType2->getName();

    /* A mapping-local context can be evaluated from the selected reaction
     * mapping. Keep the legacy materialized expansion for contexts that span
     * both reactants or cannot be represented by one weighted molecule. */
    EnergyBindingContext compactContext;
    bool useCompact = ef->getBindingContext(
        mt1Name, bindSite1, mt2Name, bindSite2, compactContext);
    if (mt1Name == mt2Name) useCompact = false;
    if (useCompact) {
        int contextReactant = -1;
        for (const auto &condition : compactContext.conditions) {
            if (contextReactant < 0) contextReactant = condition.reactantIdx;
            if (condition.reactantIdx != contextReactant) {
                useCompact = false;
                break;
            }
        }

        /* A term may require several occupancy predicates simultaneously.
         * EnergyRxnClass evaluates the full mask against the selected
         * weighted molecule, so conjunctions remain compact rather than
         * expanding into one rule per boolean context combination. */
        for (const auto &term : compactContext.conditionalTerms) {
            if (term.conditionMask == 0) {
                useCompact = false;
                break;
            }
        }

        /* The reverse binding transformation is represented by the first
         * connected template molecule.  Restrict the initial compact path to
         * contexts on reactant 0 so both directions use the same direct
         * molecule mapping without changing any other Arrhenius behavior. */
        if (contextReactant != 0) useCompact = false;
    }

    if (useCompact) {
        const double activationRate = std::exp(-Ea0 / ef->getRT());
        if (!std::isfinite(activationRate) || activationRate < 0.0) {
            cerr << "Error: Arrhenius activation rate for " << rxnName
                 << " is not finite and non-negative." << endl;
            return false;
        }
        for (int direction = 0; direction < 2; direction++) {
            const bool isForward = (direction == 0);
            if ((isForward && !includeForward) ||
                (!isForward && !includeReverse)) continue;
            TemplateMolecule *t1 = new TemplateMolecule(molType1);
            TemplateMolecule *t2 = new TemplateMolecule(molType2);

            if (isForward) {
                t1->addEmptyComponent(bindSite1);
                t2->addEmptyComponent(bindSite2);
            } else {
                TemplateMolecule::bind(t1, bindSite1, "", t2, bindSite2, "");
            }

            vector<TemplateMolecule *> templates;
            templates.push_back(t1);
            if (isForward) templates.push_back(t2);

            TransformationSet *ts = new TransformationSet(templates);
            if (isForward) ts->addBindingTransform(t1, bindSite1, t2, bindSite2);
            else ts->addUnbindingTransform(t1, bindSite1, t2, bindSite2);
            ts->setComplexBookkeeping(blockSameComplexBinding);
            ts->finalize();

            /* Pull out the activation-only term as DOR's base rate.  The
             * EnergyRxnClass supplies the context factor for each mapping. */
            string directionName = rxnName + (isForward ? "_fwd" : "_rev");
            EnergyRxnClass *r = new EnergyRxnClass(
                directionName, activationRate, "", ts, 0, compactContext,
                phi_val, ef->getRT(), isForward, s);

            s->addReaction(r);
            reaction_count++;

            if (verbose) {
                cout << "\t  Created compact "
                     << (isForward ? "forward" : "reverse")
                     << " rule: " << directionName
                     << " (mapping-local Arrhenius context)" << endl;
            }
        }
        return true;
    }

    // Run the legacy expansion algorithm for non-factorized contexts.
    vector<ExpandedRuleInfo> expanded = ef->expandBindingRule(
        rxnName, Ea0, phi_val, mt1Name, bindSite1, mt2Name, bindSite2);

    for (const auto &rule : expanded) {
        if ((rule.isForward && !includeForward) ||
            (!rule.isForward && !includeReverse)) continue;
        if (!std::isfinite(rule.rate) || rule.rate < 0.0) {
            cerr << "Error: Arrhenius rate for " << rule.name
                 << " is not finite and non-negative." << endl;
            return false;
        }
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

bool createExpandedStateChangeReactions(
    const string &rxnName,
    double phi_val,
    double Ea0,
    MoleculeType *molType,
    const string &component,
    const string &stateFrom,
    const string &stateTo,
    System *s,
    bool blockSameComplexBinding,
    bool verbose,
    int &reaction_count,
    bool includeForward,
    bool includeReverse)
{
    EnergyFunction *ef = s->getEnergyFunction();
    if (!ef || !molType || stateFrom.empty() || stateTo.empty()) return false;

    vector<ExpandedRuleInfo> expanded = ef->expandStateChangeRule(
        rxnName, Ea0, phi_val, molType->getName(), component, stateFrom, stateTo);
    for (const auto &rule : expanded) {
        if ((rule.isForward && !includeForward) ||
            (!rule.isForward && !includeReverse)) continue;

        const string &fromState = rule.isForward ? stateFrom : stateTo;
        const string &toState = rule.isForward ? stateTo : stateFrom;
        TemplateMolecule *templateMolecule = new TemplateMolecule(molType);
        TransformationSet *transformations = nullptr;
        try {
            templateMolecule->addComponentConstraint(component, fromState);
            for (const auto &constraint : rule.constraints) {
                if (constraint.reactantIdx != 0 ||
                    molType->isEquivalentComponent(constraint.compName)) {
                    return false;
                }
                if (constraint.mustBeBound)
                    templateMolecule->addBoundComponent(constraint.compName);
                else
                    templateMolecule->addEmptyComponent(constraint.compName);
            }

            vector<TemplateMolecule *> templates(1, templateMolecule);
            transformations = new TransformationSet(templates);
            if (!transformations->addStateChangeTransform(
                    templateMolecule, component, toState)) {
                delete transformations;
                return false;
            }
            transformations->setComplexBookkeeping(blockSameComplexBinding);
            transformations->finalize();
        } catch (const std::exception &) {
            delete transformations;
            return false;
        }

        if (!std::isfinite(rule.rate) || rule.rate < 0.0) {
            delete transformations;
            return false;
        }
        BasicRxnClass *reaction = new BasicRxnClass(
            rule.name, rule.rate, "", transformations, s);
        s->addReaction(reaction);
        ++reaction_count;
        if (verbose)
            cout << "\t  Created " << (rule.isForward ? "forward" : "reverse")
                 << " state-change rule: " << rule.name
                 << "  rate=" << rule.rate << endl;
    }
    return !expanded.empty();
}

} // namespace NFinput
