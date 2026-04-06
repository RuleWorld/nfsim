/*
 * NFinput_energy.hh
 *
 * Declarations for energy pattern parsing and Arrhenius rule expansion.
 * Called from NFinput.cpp during model initialization.
 */

#ifndef NFINPUT_ENERGY_HH_
#define NFINPUT_ENERGY_HH_

#include "../NFcore/NFcore.hh"
#include "../NFcore/energyPattern.hh"
#include "TinyXML/tinyxml.h"

#include <map>
#include <string>

namespace NFinput {

    /*
     * Parse <ListOfEnergyPatterns> from the model XML.
     * Creates an EnergyFunction on the System if patterns exist.
     * Call this AFTER initParameters() and BEFORE initReactionRules().
     */
    bool parseEnergyPatterns(
        TiXmlElement *pModel,
        NFcore::System *s,
        std::map<std::string, double> &parameter,
        bool verbose);

    /*
     * Create expanded BasicRxnClass instances for a binding energy rule.
     * Call this from initReactionRules() when rateLawType=="Arrhenius".
     */
    bool createExpandedBindingReactions(
        const std::string &rxnName,
        double phi_val,
        double Ea0,
        NFcore::MoleculeType *molType1, const std::string &bindSite1,
        NFcore::MoleculeType *molType2, const std::string &bindSite2,
        NFcore::System *s,
        std::map<std::string, double> &parameter,
        std::map<std::string, int> &allowedStates,
        bool blockSameComplexBinding,
        bool verbose,
        int &reaction_count);

} // namespace NFinput

#endif /* NFINPUT_ENERGY_HH_ */
