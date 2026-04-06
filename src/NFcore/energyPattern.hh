/*
 * energyPattern.hh
 *
 * Energy-based BioNetGen (eBNGL) support for NFsim.
 *
 * Implements the energy rule expansion algorithm from:
 *   Sekar, J. A. P. (2015). "Rule-based Modeling of Cell Signaling."
 *   PhD Dissertation, University of Pittsburgh, Chapter 3.
 *
 * The key insight (Corollary 3.3-43): only energy pattern embeddings
 * that overlap with the reaction center contribute to ΔG. This allows
 * each energy rule to be expanded into a finite set of conventional
 * rules with pre-computed rate constants, eliminating any runtime
 * energy computation or rejection sampling.
 *
 * Architecture:
 *   1. EnergyPatternInfo  — lightweight descriptor parsed from XML,
 *                           used for expansion analysis
 *   2. EnergyFunction     — holds all energy patterns + parameters,
 *                           drives the expansion
 *   3. ExpandedRule       — output of expansion: a conventional rule
 *                           with a pre-computed rate constant
 *
 * Authors: Achyudhan Kutuva, James R. Faeder
 */

#ifndef ENERGYPATTERN_HH_
#define ENERGYPATTERN_HH_

#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>

namespace NFcore {

    class System;
    class TemplateMolecule;
    class MoleculeType;
    class TransformationSet;

    /*
     * Lightweight description of one molecule within an energy pattern.
     * Parsed directly from XML — no TemplateMolecule needed.
     */
    struct EpMolecule {
        std::string typeName;     // e.g. "S"
        std::string xmlId;        // molecule id within the pattern XML

        // Components mentioned in this molecule within the pattern.
        // Key: component name (e.g. "A"), Value: bond partner xmlId or "" if unbound requirement
        // A component listed here with a bond means the pattern REQUIRES that bond.
        struct CompInfo {
            std::string name;           // component name, e.g. "A"
            std::string bondPartnerId;  // xmlId of partner comp, or "" if free/wildcard
            bool        isBound;        // true if pattern requires this comp bound
            std::string stateConstraint;// internal state constraint, or "" if none
        };
        std::vector<CompInfo> components;
    };

    /*
     * Lightweight description of an energy pattern parsed from XML.
     * Contains enough info to analyze overlaps with reaction centers
     * and determine context conditions for rule expansion.
     */
    struct EnergyPatternInfo {
        std::string id;           // pattern id from XML
        double      energyValue;  // Gibbs free energy of formation
        std::vector<EpMolecule> molecules;

        // Bond pairs within this pattern: (mol_idx_1, comp_idx_1, mol_idx_2, comp_idx_2)
        struct Bond {
            int mol1, comp1, mol2, comp2;
        };
        std::vector<Bond> bonds;
    };

    /*
     * Describes a "context condition" — an additional constraint on a
     * reactant component that must be resolved during expansion.
     *
     * For example, if an energy pattern requires S(A!1,B!2).A(s!1).B(s!2),
     * and the reaction center is the A-s bond, then the context condition
     * is: "S.B must be bound to B.s".
     */
    struct ContextCondition {
        std::string molType;      // molecule type name (e.g. "S")
        int         reactantIdx;  // which reactant (0 or 1) this is on
        std::string compName;     // component name (e.g. "B")
        // The partner type and site, if the condition requires a specific bond
        std::string partnerType;  // e.g. "B"
        std::string partnerComp;  // e.g. "s"
        // Which energy pattern(s) this condition gates
        std::vector<int> gatedPatternIndices;
    };

    /*
     * Output of the expansion: one conventional rule with a fixed rate.
     */
    struct ExpandedRuleInfo {
        std::string name;
        double      rate;
        double      deltaG;
        bool        isForward;  // true = binding, false = unbinding (for reversible rules)

        // Context constraints to add to reactant templates.
        // Key: (reactantIdx, compName), Value: true=must be bound, false=must be empty
        struct ContextConstraint {
            int         reactantIdx;
            std::string compName;
            bool        mustBeBound;
        };
        std::vector<ContextConstraint> constraints;
    };

    /*
     * EnergyFunction: holds all parsed energy patterns and parameters,
     * and implements the Sekar expansion algorithm.
     */
    class EnergyFunction {
    public:
        EnergyFunction(double phi, double RT);
        ~EnergyFunction();

        double getPhi() const { return phi; }
        double getRT()  const { return RT; }

        void addEnergyPattern(const EnergyPatternInfo &ep);
        int  getNumPatterns() const { return (int)patterns.size(); }
        const EnergyPatternInfo& getPattern(int i) const { return patterns[i]; }

        /*
         * Core expansion algorithm (Sekar §3.4).
         *
         * Given a binding energy rule:
         *   reactantMolTypes[0](bindSites[0]) + reactantMolTypes[1](bindSites[1])
         *     <-> reactantMolTypes[0](bindSites[0]!1).reactantMolTypes[1](bindSites[1]!1)
         *
         * Returns the set of expanded rules (forward + reverse) with
         * pre-computed Arrhenius rate constants.
         *
         * @param rxnName       Base name for the rule
         * @param Ea0           Activation energy parameter
         * @param molType1      Name of first reactant molecule type
         * @param bindSite1     Binding site on first reactant
         * @param molType2      Name of second reactant molecule type
         * @param bindSite2     Binding site on second reactant
         * @return              Vector of expanded rules
         */
        std::vector<ExpandedRuleInfo> expandBindingRule(
            const std::string &rxnName,
            double Ea0,
            double phi,
            const std::string &molType1, const std::string &bindSite1,
            const std::string &molType2, const std::string &bindSite2
        ) const;

        /*
         * Expansion for unimolecular state-change rules:
         *   molType(comp~stateFrom) <-> molType(comp~stateTo)
         */
        std::vector<ExpandedRuleInfo> expandStateChangeRule(
            const std::string &rxnName,
            double Ea0,
            double phi,
            const std::string &molType, const std::string &comp,
            const std::string &stateFrom, const std::string &stateTo
        ) const;

    private:
        double phi;   // Default distribution parameter (typically 0.5)
        double RT;    // Gas constant × Temperature

        std::vector<EnergyPatternInfo> patterns;

        /*
         * Find energy patterns that overlap with a binding reaction center.
         * Returns indices into the patterns vector.
         */
        std::vector<int> findRelevantPatternsForBinding(
            const std::string &molType1, const std::string &bindSite1,
            const std::string &molType2, const std::string &bindSite2
        ) const;

        /*
         * Find energy patterns that overlap with a state-change reaction center.
         */
        std::vector<int> findRelevantPatternsForStateChange(
            const std::string &molType, const std::string &comp
        ) const;

        /*
         * Extract context conditions from a set of relevant energy patterns.
         * These are the additional component states that need to be resolved
         * to determine which patterns contribute to ΔG.
         */
        std::vector<ContextCondition> extractContextConditions(
            const std::vector<int> &relevantPatternIndices,
            const std::string &molType1, const std::string &site1,
            const std::string &molType2, const std::string &site2
        ) const;

        /*
         * Compute Arrhenius rate constant.
         *   k_fwd = exp(-(Ea0 + phi * deltaG) / RT)
         *   k_rev = exp(-(Ea0 + (phi - 1) * deltaG) / RT)
         */
        double computeForwardRate(double Ea0, double deltaG, double phi) const {
            return std::exp(-(Ea0 + phi * deltaG) / RT);
        }
        double computeReverseRate(double Ea0, double deltaG, double phi) const {
            return std::exp(-(Ea0 + (phi - 1.0) * deltaG) / RT);
        }
    };

} // namespace NFcore

#endif /* ENERGYPATTERN_HH_ */
