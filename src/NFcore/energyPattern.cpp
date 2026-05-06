/*
 * energyPattern.cpp
 *
 * Implementation of the Sekar energy rule expansion algorithm.
 *
 * The expansion works by:
 * 1. Identifying energy patterns that overlap with the reaction center
 *    (Corollary 3.3-43: only these contribute to ΔG)
 * 2. Extracting "context conditions" — components whose bond/state
 *    determines whether a pattern matches
 * 3. Enumerating all combinations of context conditions
 * 4. For each combination, computing ΔG and the Arrhenius rate
 * 5. Producing a conventional rule for each combination
 */

#include "energyPattern.hh"
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace NFcore;
using namespace std;

EnergyFunction::EnergyFunction(double phi, double RT) : phi(phi), RT(RT) {}
EnergyFunction::~EnergyFunction() {}

void EnergyFunction::addEnergyPattern(const EnergyPatternInfo &ep) {
    patterns.push_back(ep);
}

/*
 * Find energy patterns relevant to a binding rule.
 *
 * A pattern is "relevant" if it contains a bond between molType1.site1
 * and molType2.site2 — i.e., the bond being formed/broken by the rule.
 * Only such patterns have different match counts in reactant vs product,
 * so only they contribute to ΔG (Sekar Corollary 3.3-43).
 */
vector<int> EnergyFunction::findRelevantPatternsForBinding(
    const string &molType1, const string &site1,
    const string &molType2, const string &site2
) const {
    vector<int> relevant;

    for (int i = 0; i < (int)patterns.size(); i++) {
        const EnergyPatternInfo &ep = patterns[i];

        // Check if this pattern contains a bond between molType1.site1 and molType2.site2
        for (const auto &bond : ep.bonds) {
            const EpMolecule &m1 = ep.molecules[bond.mol1];
            const EpMolecule &m2 = ep.molecules[bond.mol2];
            const string &c1 = m1.components[bond.comp1].name;
            const string &c2 = m2.components[bond.comp2].name;

            bool match_forward = (m1.typeName == molType1 && c1 == site1 &&
                                  m2.typeName == molType2 && c2 == site2);
            bool match_reverse = (m1.typeName == molType2 && c1 == site2 &&
                                  m2.typeName == molType1 && c2 == site1);

            if (match_forward || match_reverse) {
                relevant.push_back(i);
                break;
            }
        }
    }

    return relevant;
}

/*
 * Find energy patterns relevant to a state-change rule.
 *
 * A pattern is relevant if it constrains the state of molType.comp,
 * because a state change on that component changes the match count.
 */
vector<int> EnergyFunction::findRelevantPatternsForStateChange(
    const string &molType, const string &comp
) const {
    vector<int> relevant;

    for (int i = 0; i < (int)patterns.size(); i++) {
        const EnergyPatternInfo &ep = patterns[i];
        for (const auto &mol : ep.molecules) {
            if (mol.typeName == molType) {
                for (const auto &c : mol.components) {
                    if (c.name == comp && !c.stateConstraint.empty()) {
                        relevant.push_back(i);
                        goto next_pattern;
                    }
                }
            }
        }
        next_pattern:;
    }

    return relevant;
}

/*
 * Extract context conditions from relevant energy patterns.
 *
 * For each relevant pattern, identify components beyond the reaction
 * center that the pattern constrains. These become context conditions
 * that must be resolved during expansion.
 *
 * Example: Pattern S(A!1,B!2).A(s!1).B(s!2) with reaction center A-s bond
 *          → Context condition: S.B must be bound to B.s
 */
vector<ContextCondition> EnergyFunction::extractContextConditions(
    const vector<int> &relevantPatternIndices,
    const string &molType1, const string &site1,
    const string &molType2, const string &site2
) const {
    // Collect unique context conditions across all relevant patterns.
    // Key: (reactantMolType, compName) → condition info
    map<pair<string,string>, ContextCondition> condMap;

    for (int pi : relevantPatternIndices) {
        const EnergyPatternInfo &ep = patterns[pi];

        // For each molecule in the pattern that matches a reactant type,
        // check for components beyond the reaction center
        for (int mi = 0; mi < (int)ep.molecules.size(); mi++) {
            const EpMolecule &mol = ep.molecules[mi];

            int reactantIdx = -1;
            if (mol.typeName == molType1) reactantIdx = 0;
            else if (mol.typeName == molType2) reactantIdx = 1;
            else continue;  // molecule from outside the reactants — skip for now

            for (const auto &comp : mol.components) {
                // Skip the reaction center component
                if (reactantIdx == 0 && comp.name == site1) continue;
                if (reactantIdx == 1 && comp.name == site2) continue;

                // This component is a context condition
                if (comp.isBound) {
                    auto key = make_pair(mol.typeName, comp.name);
                    if (condMap.find(key) == condMap.end()) {
                        ContextCondition cc;
                        cc.molType = mol.typeName;
                        cc.reactantIdx = reactantIdx;
                        cc.compName = comp.name;

                        // Find the bond partner
                        for (const auto &bond : ep.bonds) {
                            if (bond.mol1 == mi && bond.comp1 == (int)(&comp - &mol.components[0])) {
                                cc.partnerType = ep.molecules[bond.mol2].typeName;
                                cc.partnerComp = ep.molecules[bond.mol2].components[bond.comp2].name;
                                break;
                            }
                            if (bond.mol2 == mi && bond.comp2 == (int)(&comp - &mol.components[0])) {
                                cc.partnerType = ep.molecules[bond.mol1].typeName;
                                cc.partnerComp = ep.molecules[bond.mol1].components[bond.comp1].name;
                                break;
                            }
                        }

                        condMap[key] = cc;
                    }
                    condMap[key].gatedPatternIndices.push_back(pi);
                }
            }
        }
    }

    vector<ContextCondition> result;
    for (auto &kv : condMap) {
        // Deduplicate gated pattern indices
        sort(kv.second.gatedPatternIndices.begin(), kv.second.gatedPatternIndices.end());
        kv.second.gatedPatternIndices.erase(
            unique(kv.second.gatedPatternIndices.begin(), kv.second.gatedPatternIndices.end()),
            kv.second.gatedPatternIndices.end());
        result.push_back(kv.second);
    }
    return result;
}

/*
 * Expand a binding energy rule into conventional rules.
 *
 * Algorithm (Sekar §3.4):
 * 1. Find energy patterns containing the reaction center bond
 * 2. Separate into "always-matching" (no extra context) and "conditional"
 * 3. Extract context conditions from conditional patterns
 * 4. Enumerate all 2^n combinations of n boolean context conditions
 * 5. For each combination, compute ΔG and rates
 */
vector<ExpandedRuleInfo> EnergyFunction::expandBindingRule(
    const string &rxnName,
    double Ea0,
    double phi,
    const string &molType1, const string &bindSite1,
    const string &molType2, const string &bindSite2
) const {
    vector<ExpandedRuleInfo> expanded;

    // Step 1: Find relevant patterns
    vector<int> relevant = findRelevantPatternsForBinding(
        molType1, bindSite1, molType2, bindSite2);

    if (relevant.empty()) {
        // No energy patterns overlap with the reaction center.
        // ΔG = 0 for all contexts → single rule with rate = exp(-Ea0/RT)
        cerr << "Warning: Arrhenius rule " << rxnName
             << " has no overlapping energy patterns. ΔG=0." << endl;

        ExpandedRuleInfo fwd;
        fwd.name = rxnName + "_fwd";
        fwd.deltaG = 0.0;
        fwd.rate = computeForwardRate(Ea0, 0.0, phi);
        fwd.isForward = true;
        expanded.push_back(fwd);

        ExpandedRuleInfo rev;
        rev.name = rxnName + "_rev";
        rev.deltaG = 0.0;
        rev.rate = computeReverseRate(Ea0, 0.0, phi);
        rev.isForward = false;
        expanded.push_back(rev);
        return expanded;
    }

    // Step 2: Classify relevant patterns into "always" and "conditional"
    // "Always" patterns: contain ONLY the reaction center molecules+sites,
    //                    no extra context → always match when the bond exists
    // "Conditional" patterns: have additional context beyond the center

    vector<int> alwaysPatterns;    // indices into 'relevant'
    vector<int> conditionalPatterns;

    for (int ri = 0; ri < (int)relevant.size(); ri++) {
        int pi = relevant[ri];
        const EnergyPatternInfo &ep = patterns[pi];

        bool hasExtraContext = false;
        for (const auto &mol : ep.molecules) {
            bool isReactantType = (mol.typeName == molType1 || mol.typeName == molType2);
            if (!isReactantType) {
                // Pattern involves a third molecule type → conditional
                hasExtraContext = true;
                break;
            }
            for (const auto &comp : mol.components) {
                // Skip the reaction center components
                if (mol.typeName == molType1 && comp.name == bindSite1) continue;
                if (mol.typeName == molType2 && comp.name == bindSite2) continue;
                // Any other component with a constraint → extra context
                if (comp.isBound || !comp.stateConstraint.empty()) {
                    hasExtraContext = true;
                    break;
                }
            }
            if (hasExtraContext) break;
        }

        if (hasExtraContext)
            conditionalPatterns.push_back(pi);
        else
            alwaysPatterns.push_back(pi);
    }

    // Base ΔG from always-matching patterns (forward direction: bond forms → +G_e)
    double baseG = 0.0;
    for (int pi : alwaysPatterns) {
        baseG += patterns[pi].energyValue;
    }

    // Step 3: Extract context conditions from conditional patterns
    vector<ContextCondition> conditions = extractContextConditions(
        conditionalPatterns, molType1, bindSite1, molType2, bindSite2);

    if (conditions.empty()) {
        // All relevant patterns are "always" — single forward + reverse rule
        ExpandedRuleInfo fwd;
        fwd.name = rxnName + "_fwd";
        fwd.deltaG = baseG;
        fwd.rate = computeForwardRate(Ea0, baseG, phi);
        fwd.isForward = true;
        expanded.push_back(fwd);

        ExpandedRuleInfo rev;
        rev.name = rxnName + "_rev";
        rev.deltaG = baseG;
        rev.rate = computeReverseRate(Ea0, baseG, phi);
        rev.isForward = false;
        expanded.push_back(rev);

        cout << "  Expanded " << rxnName << " → 1 forward + 1 reverse rule"
             << "  (ΔG=" << baseG << ", k_fwd=" << fwd.rate
             << ", k_rev=" << rev.rate << ")" << endl;
        return expanded;
    }

    // Step 4: Enumerate all 2^n combinations of context conditions
    int nCond = (int)conditions.size();
    int nCombinations = 1 << nCond;  // 2^n

    cout << "  Expanding " << rxnName << " with " << nCond
         << " context condition(s), " << nCombinations
         << " variant(s) per direction:" << endl;

    for (int combo = 0; combo < nCombinations; combo++) {
        // Determine which conditional patterns are active in this combination
        set<int> activePatterns;  // indices into patterns[]
        for (int pi : alwaysPatterns) activePatterns.insert(pi);

        // Build context constraints for this combination
        vector<ExpandedRuleInfo::ContextConstraint> constraints;
        for (int ci = 0; ci < nCond; ci++) {
            bool conditionMet = (combo >> ci) & 1;

            ExpandedRuleInfo::ContextConstraint cc;
            cc.reactantIdx = conditions[ci].reactantIdx;
            cc.compName = conditions[ci].compName;
            cc.mustBeBound = conditionMet;
            constraints.push_back(cc);

            if (conditionMet) {
                for (int pi : conditions[ci].gatedPatternIndices) {
                    activePatterns.insert(pi);
                }
            }
        }

        // Step 5: Compute ΔG for this combination
        double deltaG = 0.0;
        for (int pi : activePatterns) {
            deltaG += patterns[pi].energyValue;
        }

        // Create forward rule
        {
            ExpandedRuleInfo rule;
            stringstream ss;
            ss << rxnName << "_fwd_v" << combo;
            rule.name = ss.str();
            rule.deltaG = deltaG;
            rule.rate = computeForwardRate(Ea0, deltaG, phi);
            rule.isForward = true;
            rule.constraints = constraints;
            expanded.push_back(rule);
        }

        // Create reverse rule
        {
            ExpandedRuleInfo rule;
            stringstream ss;
            ss << rxnName << "_rev_v" << combo;
            rule.name = ss.str();
            rule.deltaG = deltaG;
            rule.rate = computeReverseRate(Ea0, deltaG, phi);
            rule.isForward = false;
            // Reverse rule has same context constraints
            rule.constraints = constraints;
            expanded.push_back(rule);
        }

        cout << "    v" << combo << ": ΔG=" << deltaG
             << "  k_fwd=" << computeForwardRate(Ea0, deltaG, phi)
             << "  k_rev=" << computeReverseRate(Ea0, deltaG, phi)
             << "  context=[";
        for (int ci = 0; ci < nCond; ci++) {
            if (ci > 0) cout << ", ";
            cout << conditions[ci].molType << "." << conditions[ci].compName
                 << "=" << (((combo >> ci) & 1) ? "bound" : "free");
        }
        cout << "]" << endl;
    }

    return expanded;
}

/*
 * Expand a state-change energy rule.
 * Same algorithm but for unimolecular reactions.
 */
vector<ExpandedRuleInfo> EnergyFunction::expandStateChangeRule(
    const string &rxnName,
    double Ea0,
    double phi,
    const string &molType, const string &comp,
    const string &stateFrom, const string &stateTo
) const {
    vector<ExpandedRuleInfo> expanded;

    vector<int> relevant = findRelevantPatternsForStateChange(molType, comp);

    // For state change, a pattern is relevant if it constrains the state
    // of the changing component. We need to determine for each relevant
    // pattern whether it matches the "from" state or "to" state.

    // Patterns matching the "to" state contribute +G_e to ΔG
    // Patterns matching the "from" state contribute -G_e to ΔG
    // (because they match in reactant but not product)

    // Classify relevant patterns into "always" and "conditional"
    vector<int> alwaysPatterns;
    vector<int> conditionalPatterns;

    for (int ri = 0; ri < (int)relevant.size(); ri++) {
        int pi = relevant[ri];
        const EnergyPatternInfo &ep = patterns[pi];

        bool hasExtraContext = false;
        for (const auto &mol : ep.molecules) {
            if (mol.typeName != molType) {
                // Pattern involves another molecule type -> conditional
                hasExtraContext = true;
                break;
            }
            for (const auto &c : mol.components) {
                // Skip the reaction center component (the one changing state)
                if (c.name == comp) continue;
                // Any other component with a constraint -> extra context
                if (c.isBound || !c.stateConstraint.empty()) {
                    hasExtraContext = true;
                    break;
                }
            }
            if (hasExtraContext) break;
        }

        if (hasExtraContext)
            conditionalPatterns.push_back(pi);
        else
            alwaysPatterns.push_back(pi);
    }

    // Base ΔG from always-matching patterns
    double baseG = 0.0;
    for (int pi : alwaysPatterns) {
        const EnergyPatternInfo &ep = patterns[pi];
        for (const auto &mol : ep.molecules) {
            if (mol.typeName == molType) {
                for (const auto &c : mol.components) {
                    if (c.name == comp) {
                        if (c.stateConstraint == stateTo) {
                            baseG += ep.energyValue;
                        } else if (c.stateConstraint == stateFrom) {
                            baseG -= ep.energyValue;
                        }
                    }
                }
            }
        }
    }

    // Extract context conditions
    // Unimolecular: only molType is involved, use extractContextConditions passing it as molType1
    vector<ContextCondition> conditions = extractContextConditions(
        conditionalPatterns, molType, comp, "", "");

    if (conditions.empty()) {
        ExpandedRuleInfo fwd;
        fwd.name = rxnName + "_fwd";
        fwd.deltaG = baseG;
        fwd.rate = computeForwardRate(Ea0, baseG, phi);
        fwd.isForward = true;
        expanded.push_back(fwd);

        ExpandedRuleInfo rev;
        rev.name = rxnName + "_rev";
        rev.deltaG = baseG;
        rev.rate = computeReverseRate(Ea0, baseG, phi);
        rev.isForward = false;
        expanded.push_back(rev);

        cout << "  Expanded state-change " << rxnName << " → 1 forward + 1 reverse rule"
             << "  (ΔG=" << baseG << ", k_fwd=" << fwd.rate
             << ", k_rev=" << rev.rate << ")" << endl;
        return expanded;
    }

    // Enumerate all 2^n combinations of context conditions
    int nCond = (int)conditions.size();
    int nCombinations = 1 << nCond;  // 2^n

    cout << "  Expanding state-change " << rxnName << " with " << nCond
         << " context condition(s), " << nCombinations
         << " variant(s) per direction:" << endl;

    for (int combo = 0; combo < nCombinations; combo++) {
        set<int> activePatterns;
        for (int pi : alwaysPatterns) activePatterns.insert(pi);

        vector<ExpandedRuleInfo::ContextConstraint> constraints;
        for (int ci = 0; ci < nCond; ci++) {
            bool conditionMet = (combo >> ci) & 1;

            ExpandedRuleInfo::ContextConstraint cc;
            cc.reactantIdx = conditions[ci].reactantIdx;
            cc.compName = conditions[ci].compName;
            cc.mustBeBound = conditionMet;
            constraints.push_back(cc);

            if (conditionMet) {
                for (int pi : conditions[ci].gatedPatternIndices) {
                    activePatterns.insert(pi);
                }
            }
        }

        double deltaG = 0.0;
        for (int pi : activePatterns) {
            const EnergyPatternInfo &ep = patterns[pi];
            for (const auto &mol : ep.molecules) {
                if (mol.typeName == molType) {
                    for (const auto &c : mol.components) {
                        if (c.name == comp) {
                            if (c.stateConstraint == stateTo) {
                                deltaG += ep.energyValue;
                            } else if (c.stateConstraint == stateFrom) {
                                deltaG -= ep.energyValue;
                            }
                        }
                    }
                }
            }
        }

        // Create forward rule
        {
            ExpandedRuleInfo rule;
            stringstream ss;
            ss << rxnName << "_fwd_v" << combo;
            rule.name = ss.str();
            rule.deltaG = deltaG;
            rule.rate = computeForwardRate(Ea0, deltaG, phi);
            rule.isForward = true;
            rule.constraints = constraints;
            expanded.push_back(rule);
        }

        // Create reverse rule
        {
            ExpandedRuleInfo rule;
            stringstream ss;
            ss << rxnName << "_rev_v" << combo;
            rule.name = ss.str();
            rule.deltaG = deltaG;
            rule.rate = computeReverseRate(Ea0, deltaG, phi);
            rule.isForward = false;
            rule.constraints = constraints;
            expanded.push_back(rule);
        }

        cout << "    v" << combo << ": ΔG=" << deltaG
             << "  k_fwd=" << computeForwardRate(Ea0, deltaG, phi)
             << "  k_rev=" << computeReverseRate(Ea0, deltaG, phi)
             << "  context=[";
        for (int ci = 0; ci < nCond; ci++) {
            if (ci > 0) cout << ", ";
            cout << conditions[ci].molType << "." << conditions[ci].compName
                 << "=" << (((combo >> ci) & 1) ? "bound" : "free");
        }
        cout << "]" << endl;
    }

    return expanded;
}
