


#include "NFcore.hh"
#include <unordered_set>


using namespace std;
using namespace NFcore;



ReactionClass::ReactionClass(string name, double baseRate, string baseRateParameterName, TransformationSet *transformationSet, System *s)
{
	this->system=s;
	this->tagged = false;
	this->useRuleMonkey = false;
	this->useConnectivity = false;
	this->directProductMolecules = 0;

	totalRateFlag=false;
	isDimerStyle=false;
	//Setup the basic properties of this reactionClass
	this->name = name;
	this->baseRate = baseRate;
	this->baseRateParameterName=baseRateParameterName;
	this->fireCounter = 0;
	this->a = 0;
	this->volumeConversionFactor = 1.0;
	this->traversalLimit = ReactionClass::NO_LIMIT;
	this->transformationSet = transformationSet;
	this->productComponentsTruncated = false;


	//Set up the template molecules from the transformationSet
	this->n_reactants   = transformationSet->getNreactants();
	this->n_mappingsets = transformationSet->getNmappingSets();
    this->reactantTemplates = new TemplateMolecule *[n_reactants];
	vector <TemplateMolecule*> tmList;
	vector <int> hasMapGenerator;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		//The main reactant should be the one that is getting modified...
		//In other words, we select the reactant that has at least one map generator, and
		//to minimize mistakes, with the least sym sites...
		TemplateMolecule *curTemplate = transformationSet->getTemplateMolecule(r);
		TemplateMolecule::traverse(curTemplate,tmList,TemplateMolecule::FIND_ALL);

		//First, single out all the templates that have at least one map generator
		for(unsigned int i=0, n=tmList.size(); i<n; i++) {

			if(tmList.at(i)->getN_mapGenerators()>0) {
				hasMapGenerator.push_back(i);
			}
		}

		//Find the one with the least sym comp bonds...
		int minSymSites = 999999;
		for(unsigned int k=0, n=hasMapGenerator.size(); k<n; k++) {
			if(tmList.at(hasMapGenerator.at(k))->getN_symCompBonds()<minSymSites) {
				curTemplate = tmList.at(hasMapGenerator.at(k));
				minSymSites = curTemplate->getN_symCompBonds();
			}
		}

		reactantTemplates[r] = curTemplate;
		tmList.clear(); hasMapGenerator.clear();
	}
	mappingSet = new MappingSet *[n_mappingsets];



	/* create blank mappingSets for the added molecules. These will be used
	 * to hold mappings to added molecules, which is useful for rules that create
	 *  molecules and then perform other transformations.  --Justin, 1Mar2011
	 */
	for ( unsigned int r = n_reactants; r < n_mappingsets; ++r )
	{
		mappingSet[r] = transformationSet->generateBlankMappingSet(r,0);
	}



	//Here, if we identify that there are disjoint sets in this pattern, from
	//the connected-to syntax, then we have to flag the ones that we actually
	//have to traverse down...
	for(unsigned int r=0; r<n_reactants; r++)
	{
		tmList.clear();

		// Get the connected set of molecules
		TemplateMolecule *curTemplate = reactantTemplates[r];
		TemplateMolecule::traverse(curTemplate,tmList,TemplateMolecule::FIND_ALL);

		//Label the unique sets, and only continue if we have more than one set
		vector <vector <TemplateMolecule *> > sets;
		vector <int> uniqueSetId;
		int setCount = TemplateMolecule::getNumDisjointSets(tmList,sets,uniqueSetId);
		if(setCount<=1) continue;


		//count the number of map generators (rxn centers) in each set
		vector <int> numMapGenerators;
		for(int s=0; s<setCount; s++) { numMapGenerators.push_back(0); }

		int curTemplateSetId = -1;
		for(unsigned int t=0, n=tmList.size(); t<n; t++) {
			if(tmList.at(t)==curTemplate) {
				curTemplateSetId = uniqueSetId.at(t);
			}
			int n_maps = numMapGenerators.at(uniqueSetId.at(t));
			numMapGenerators.at(uniqueSetId.at(t)) = n_maps+tmList.at(t)->getN_mapGenerators();
		}



		//Lets rearrange the connected-to elements so that the one head is listed as
		//connected to all other molecules.  This will better suit our needs.

		// first, clear out the old connections
		for(unsigned int i=0, n=tmList.size(); i<n; i++) {
			tmList.at(i)->clearConnectedTo();
		}

		// add back the connections, but always through the head template
		int rxnCenterSets = 1;
		int curSet=0;
		for(unsigned int i=0, n=uniqueSetId.size(); i<n; i++) {
			if(uniqueSetId.at(i)==curTemplateSetId) {
				if(curSet==curTemplateSetId) curSet++;
				continue;
			}
			if(uniqueSetId.at(i)==curSet) {
				bool otherHasRxnCenter = false;
				if(numMapGenerators.at(curSet)>0) {
					otherHasRxnCenter=true;
					rxnCenterSets++;
				}
				TemplateMolecule *otherTemplate = tmList.at(i);
				int ctIndex1=curTemplate->getN_connectedTo();
				int ctIndex2=otherTemplate->getN_connectedTo();
				curTemplate->addConnectedTo(otherTemplate,ctIndex2,otherHasRxnCenter);
				otherTemplate->addConnectedTo(curTemplate,ctIndex1);
				curSet++;
			}
		}

		if(rxnCenterSets>2) {
			cout.flush();
			cerr<<"\n\n   Error in Reaction Rule: "<<name<<endl;
			cerr<<"   You created a reaction with a pattern that includes the connected-to\n";
			cerr<<"   syntax (ie: A().B()).  You included 3 or more disjoint sets of molecules\n";
			cerr<<"   where there are more than 2 sets with rxn centers.  This may work ok, \n";
			cerr<<"   but you really shouldn't ever do something this crazy, so I'm just going\n";
			cerr<<"   to stop you now.  Goodbye.\n"<<endl;
			exit(1);
		}




		//Finally, clear out the data structures.
		for(unsigned int i=0, n=sets.size(); i<n; i++) sets.at(i).clear();
		sets.clear(); uniqueSetId.clear();
		numMapGenerators.clear();
	}


	//Check here to see if we have molecule types that are the same across different reactants
	//Because if so, we will give a warning
	if(n_reactants>2) cerr<<"Warning!! You created a reaction ("<< name <<") that has more than 2 reactants.  This has not been extensively tested!"<<endl;

	if(n_reactants==2)
	{	//If the reactants are of the same type, then we have to make a few special considerations
		if(reactantTemplates[0]->getMoleculeType()==reactantTemplates[1]->getMoleculeType())
		{
			cout<<endl;
			cout<<"Warning! You have a binding rxn (" << name << ") that allows a moleculeType to bind another of the same type."<<endl;
			cout<<"Make sure that is correct, because this can potentially make long polymers or large aggregates."<<endl;
			cout<<endl;
		}
	}

	// Correct the rate for reaction center symmetry.  Note that these all assign
	// through 'this->': the constructor argument 'baseRate' shadows the member,
	// and the member was copied from it above, so scaling the argument here would
	// discard the correction.  Only rate laws that route through setBaseRate()
	// (which applies the factor itself) used to recover it; every other rate law
	// -- global function, local function (DOR), function product, MM -- is
	// constructed with baseRate=1 and never calls setBaseRate, so a symmetric
	// rule fired at 1/symmetryFactor times its intended rate (2x for a homodimer).
	if ( this->transformationSet->usingSymmetryFactor() )
	{	// new general method for handling reaction center symmetry
		this->baseRate *= this->transformationSet->getSymmetryFactor();
	}
	else
	{	// old method for handling symmetric binding and unbinding
		if(n_reactants==2)
		{
			//If the binding is symmetric
			if(transformationSet->hasSymBindingTransform()) {
				cout<<endl;
				cout<<"Warning! You have an binding rxn (" << name << ") that is symmetric."<<endl;
				cout<<"Make sure that is correct."<<endl;

				cout<<endl;
				this->baseRate = this->baseRate*0.5;  //We have to correct the rate to get the proper factor
				isDimerStyle=true;
			}
		}
		if(n_reactants==1)
		{
			if(transformationSet->hasSymUnbindingTransform())
			{
				cout<<endl;
				cout<<"Warning! You have an unbinding rxn (" << name << ") that is symmetric."<<endl;
				cout<<"Make sure that is correct."<<endl;
				cout<<endl;
				this->baseRate = this->baseRate*0.5;  //We have to correct the rate to get the proper factor
				isDimerStyle=true;
			}
		}
	}


	onTheFlyObservables=true;


	// check for population type reactants
	isPopulationType = new bool[n_reactants];
	matchOncePerReactant = new bool[n_reactants];
	contextCountsPerComplex = new bool[n_reactants];
	for( unsigned int i=0; i < n_reactants; ++i )
	{
		isPopulationType[i] = reactantTemplates[i]->getMoleculeType()->isPopulationType();
		matchOncePerReactant[i] = false;
		contextCountsPerComplex[i] = false;
	}

	// Count a reactant the rule never transforms once per complex, as BNG does.
	//
	// BNG gives a pure context pattern one reaction instance per matching complex
	// no matter how many molecules inside it match, because the reaction produced
	// does not depend on which embedding was chosen -- same reactants, same
	// products, same transformation.  Verified against BNG's generated network for
	// a homodimer, a single subunit of that homodimer, a heterodimer, and a
	// scaffold holding two *distinguishable* copies: all four get a bare rate
	// constant.  A pattern the rule does transform is different and is left alone
	// here: binding one of a homodimer's two sites really is two reactions, and
	// BNG emits the factor of two for it.
	//
	// Only meaningful while complexes are tracked; see contextCountsPerComplex.
	if( system->isUsingComplex() )
	{
		for( unsigned int i=0; i < n_reactants; ++i )
		{
			if( isPopulationType[i] ) continue;
			contextCountsPerComplex[i] = transformationSet->isPureContextReactant(i);
		}
	}


	// calculate discrete count corrections for symmetric population reactants
	//  e.g. number of reactant pairs = A*(A-1)/2.  Note that the factor of two
	//  is part of the symmetry factor above.
	identicalPopCountCorrection = new int[n_reactants];
	for ( int i=0; i < (int)n_reactants; ++i )
	{
		identicalPopCountCorrection[i] = 0;
		if ( isPopulationType[i] )
		{
			for ( int j=i-1; j >= 0; --j )
			{
				if ( reactantTemplates[i]->getMoleculeType() == reactantTemplates[j]->getMoleculeType() )
				{
					identicalPopCountCorrection[i] = identicalPopCountCorrection[j] + 1;
					break;
				}
			}
		}
	}
}


void ReactionClass::appendConnectedRxn(ReactionClass * rxn) {
	this->connectedReactions.push_back(rxn);
}

bool ReactionClass::isReactionConnected(ReactionClass * rxn) {
	// First check if any of the operations share MoleculeType and components with
	// one of the reactant templates of rxn.
	if (this->transformationSet->checkConnection(rxn)) return true;

	// Full membership refresh revisits every explicit reactant template in the
	// fired rule, not only templates that carry direct transformations.
	for (unsigned int i=0, n=allReactantTemplates.size(); i<n; i++) {
		if (rxn->isTemplateCompatible(allReactantTemplates[i])) return true;
	}

	// Product templates can also create new compatible mappings, but avoid
	// broadening pure-synthesis rules where this over-connects add-only paths.
	if (n_reactants > 0) {
		for (unsigned int i=0, n=allProductTemplates.size(); i<n; i++) {
			if (rxn->isTemplateCompatible(allProductTemplates[i])) return true;
		}
	}
	return false;
}

ReactionClass::~ReactionClass()
{
	delete [] reactantTemplates;
	delete transformationSet;
	for ( unsigned int r = n_reactants; r < n_mappingsets; ++r )
	{
		delete mappingSet[r];
	}

	delete [] mappingSet;
	delete [] isPopulationType;
	delete [] matchOncePerReactant;
	delete [] contextCountsPerComplex;
	delete [] identicalPopCountCorrection;
	delete directProductMolecules;
	connectedReactions.clear();
}

/** Fill the reactant and product templates for inferring reaction connectivity matrix
 * @author Arvind Rasi Subramaniam
 */
void ReactionClass::setAllReactantAndProductTemplates(map <string,TemplateMolecule *> reactants,
					map <string,TemplateMolecule *> products) {
	map <string, TemplateMolecule *>::iterator it;
	// Fill the reactant template pattern
	for (it = reactants.begin(); it != reactants.end(); ++it)
		this->allReactantTemplates.push_back(it->second);
	// Fill the product template pattern
	for (it = products.begin(); it != products.end(); ++it)
		this->allProductTemplates.push_back(it->second);
}


void ReactionClass::setBaseRate(double newBaseRate,string newBaseRateName) {
	if ( this->transformationSet->usingSymmetryFactor() )
	{	this->baseRate = this->transformationSet->getSymmetryFactor() * newBaseRate;   }
	else if (isDimerStyle)
	{	this->baseRate = 0.5 * newBaseRate;   }
	else
	{	this->baseRate = newBaseRate;   }

	this->baseRateParameterName = newBaseRateName;
	update_a();
};


void ReactionClass::resetBaseRateFromSystemParamter() {

	if(!this->baseRateParameterName.empty()) {
		if ( transformationSet->usingSymmetryFactor() ) {
			this->baseRate = transformationSet->getSymmetryFactor() * system->getParameter(this->baseRateParameterName);
		}
		else if (isDimerStyle) {
			this->baseRate = 0.5 * system->getParameter(this->baseRateParameterName);
		}
		else {
			this->baseRate=system->getParameter(this->baseRateParameterName);
		}
		this->update_a();
	}

}


/** For use in MoleculeTye::updateRxnMembership
 * @author Arvind Rasi Subramaniam
 */
MoleculeType *ReactionClass::getMoleculeTypeOfReactantTemplate(int pos) const {
	// return reactantTemplates.at(pos)->getMoleculeType();
	return reactantTemplates[pos]->getMoleculeType();
}

bool ReactionClass::isDirectProductMolecule(Molecule *molecule,
		bool compactDirectProducts) const
{
	if (compactDirectProducts)
		return std::find(directProductMoleculeList.begin(),
				directProductMoleculeList.end(), molecule) !=
			directProductMoleculeList.end();
	return directProductMolecules != 0 &&
		directProductMolecules->find(molecule) != directProductMolecules->end();
}


void ReactionClass::printDetails() const {
	cout << name << "  (id=" << this->rxnId << ", baseRate=" << baseRate
			<< ",  a=" << a << ", fired=" << fireCounter << " times )" << endl;
	for (unsigned int r = 0; r < n_reactants; r++) {
		cout << "      -|" << this->getReactantCount(r) << " mappings|\t";
		cout << this->reactantTemplates[r]->getPatternString() << "\n";
	}
	if (n_reactants == 0)
		cout
				<< "      >No Reactants: so this rule either creates new species or does nothing."
				<< endl;
	cout << "\n";
}

void ReactionClass::fire(double random_A_number) {
	this->fire(random_A_number, false);
}

// AS2023 - Alternative call signature to tell fire call when we are tracking 
// each firing for the rxnlog argument
string ReactionClass::fire(double random_A_number, bool track) {
	fireCounter++;
	struct ProfileScope {
		System *system;
		int rxnId;
		const string *rxnName;
		clock_t start;
		int nullEventsBefore;
		bool enabled;

		ProfileScope(System *s, int id, const string &name)
			: system(s), rxnId(id), rxnName(&name), start(0),
			  nullEventsBefore(System::NULL_EVENT_COUNTER), enabled(false) {
			enabled = system != 0 && system->isProfilingEnabled();
			if (enabled) {
				start = clock();
				system->beginProfileReactionFire(rxnId, *rxnName);
			}
		}

		~ProfileScope() {
			if (enabled) {
				system->recordProfileReactionFire(
					this->rxnId, *rxnName, clock() - start,
					System::NULL_EVENT_COUNTER > nullEventsBefore);
			}
		}
	} profileScope(system, rxnId, name);

	productComponentSizes.clear();
	productComponentsTruncated = false;

	// First randomly pick the reactants to fire by selecting the MappingSets
	this->pickMappingSets(random_A_number);

	// Defensive check: a picked MappingSet can occasionally contain an unmapped entry
	// (null molecule) in edge cases involving internal bond reconnection/symmetry.
	// Treat this as a null event and skip firing to avoid dereferencing null mappings.
	// This has to come before checkMolecularity, which reads the complex of each
	// reactant's first mapped molecule when complex bookkeeping is on.
	for (unsigned int k=0; k<n_reactants; k++) {
		Mapping *picked = mappingSet[k]->get(0);
		if (picked == 0 || picked->getMolecule() == 0) {
			++(System::NULL_EVENT_COUNTER);
			return string("");
		}
	}

	// Check reactants for correct molecularity:
	if ( ! transformationSet->checkMolecularity(mappingSet) ) {
		// wrong molecularity!  this is a NULL event
		++(System::NULL_EVENT_COUNTER);
		// AS2023 - we need to return a string now that this can return 
		// an event log if track is true
		return string("");
	}

	/* Compact EnergyPattern binding rules can detect an occupied endpoint
	 * before the generic product/membership pipeline.  This preserves the
	 * transformation's null-event semantics while avoiding all work that
	 * cannot change the state. */
	if (!this->checkPreFireConditions(mappingSet)) {
		++(System::NULL_EVENT_COUNTER);
		return string("");
	}

	// Check reactant-side filters (include_reactants / exclude_reactants).
	for (unsigned int k=0; k<n_reactants; k++) {
		Molecule *mol = mappingSet[k]->get(0)->getMolecule();
		if (!transformationSet->checkReactantFilters(k, mol)) {
			++(System::NULL_EVENT_COUNTER);
			return string("");
		}
	}


	// Generate the set of possible products that we need to update
	// (excluding new molecules, we'll get those later --Justin).  A compact
	// energy rule can use only its explicitly mapped endpoints when it has
	// already proven that no observable, Type-II function, or indirect
	// membership dependency needs the rest of the bonded complex.
	bool directProductsPrepared = false;
	if (this->canUseDirectProductList()) {
		ProfileTime directProductStart = system->isProfileReactionActive()
			? profileNow() : ProfileTime();
		directProductMoleculeList.clear();
		for (unsigned int msIndex = 0; msIndex < n_mappingsets; ++msIndex) {
			MappingSet *ms = mappingSet[msIndex];
			if (ms == 0) continue;
			for (unsigned int mapIndex = 0;
					mapIndex < ms->getNumOfMappings(); ++mapIndex) {
				Mapping *mapping = ms->get(mapIndex);
				if (mapping == 0 || mapping->getMolecule() == 0) continue;
				Molecule *molecule = mapping->getMolecule();
				if (std::find(directProductMoleculeList.begin(),
						directProductMoleculeList.end(), molecule) ==
						directProductMoleculeList.end())
				{
					directProductMoleculeList.push_back(molecule);
					products.push_back(molecule);
				}
			}
		}
		directProductsPrepared = true;
		if (system->isProfileReactionActive())
			system->recordProfileProductPreparation(
					profileElapsedSeconds(directProductStart),
					static_cast<unsigned long long>(products.size()));
	} else {
		this->transformationSet->getListOfProducts(
				mappingSet, products, traversalLimit, &productComponentSizes,
				&productComponentsTruncated);
	}

	// Check product-side filters (include_products / exclude_products).
	// If the resulting products don't pass the filter, treat this as a null event.
	if (!transformationSet->checkProductFilters(products)) {
		products.clear();
		++(System::NULL_EVENT_COUNTER);
		return string("");
	}

	// Loop through the products (excluding added molecules) and remove from observables
	if (this->onTheFlyObservables) {
		bool profileObservables = system->isProfileReactionActive();
		ProfileTime profileObservablesStart = profileObservables
			? profileNow() : ProfileTime();
		unsigned long long profileObservableMolecules = profileObservables
			? static_cast<unsigned long long>(products.size()) : 0;
		std::unordered_set<int> updatedComplexIds;

		// molecule observables..
		for ( molIter = products.begin(); molIter != products.end(); molIter++ )
			(*molIter)->removeFromObservables();

		// species observables..
		if(system->getNumOfSpeciesObs()>0) {
			// we can find reactant complexes by following mappingSets to target molecules
			int matches = 0;
			Complex * c;
			for ( unsigned int k=0; k<transformationSet->getNreactants(); k++) {
				// get complexID and check if we've already updated that complex
				int complexId = mappingSet[k]->get(0)->getMolecule()->getComplexID();
				if ( updatedComplexIds.insert(complexId).second ) {
					// complex has not been updated, so do it now.
					c = mappingSet[k]->get(0)->getMolecule()->getComplex();
					for(int i=0; i<system->getNumOfSpeciesObs(); i++) {
						matches = system->getSpeciesObs(i)->isObservable(c);
						if (matches > 0) {
							system->getSpeciesObs(i)->subtract(matches);
						}
					}
				}
			}

			// grab added molecules that are represented as populations and remove from observables
			for ( int k=0; k<transformationSet->getNumOfAddMoleculeTransforms(); k++)
			{
				Molecule * addmol = transformationSet->getPopulationPointer((unsigned int)k);
				if ( addmol == NULL ) continue;

				// get complexID and check if we've already updated that complex
				int complexId = addmol->getComplexID();
				if ( updatedComplexIds.insert(complexId).second ) {
					// complex has not been updated, so do it now.
					c = addmol->getComplex();
					for (int i=0; i < system->getNumOfSpeciesObs(); i++) {
						matches = system->getSpeciesObs(i)->isObservable(c);
						if (matches > 0) {
							system->getSpeciesObs(i)->subtract(matches);
						}
					}
				}
			}
		}
		if (profileObservables)
			system->recordProfileObservableRemoval(
					profileElapsedSeconds(profileObservablesStart),
					profileObservableMolecules);
	}

	// Through the MappingSet, transform all the molecules as neccessary
	//  This will also create new molecules, as required.  As a side effect,
	//  deleted molecules will be removed from observables.
	// AS2023 - if tracking is turned on, transform needs a string to build up
	string logstr;
	if (this->system->getReactionTrackingStatus()) {
		logstr = this->transformationSet->transform(this->mappingSet, true);
		
	} else {
		logstr = this->transformationSet->transform(this->mappingSet);
	}

	// Add newly created molecules to the list of products
	this->transformationSet->getListOfAddedMolecules(mappingSet,products,traversalLimit);

	// Track molecules that were explicitly mapped by this firing when either
	// connectivity-aware membership or a compact energy reaction may use the
	// endpoint-local membership filter. Products added through bonded-neighborhood
	// traversal remain conservative and use the full updater.
	bool hasIndirectProducts = false;
	bool trackDirectProducts = useConnectivity || this->usesIncrementalMembership();
	if (trackDirectProducts) {
		indirectMembershipDecisions.clear();
		if (!directProductsPrepared) {
			if (directProductMolecules == 0)
				directProductMolecules = new unordered_set<Molecule *>();
			directProductMolecules->clear();
			for (unsigned int msIndex=0; msIndex<n_mappingsets; msIndex++) {
				MappingSet *ms = mappingSet[msIndex];
				if (ms==0) continue;
				for (unsigned int mapIndex=0; mapIndex<ms->getNumOfMappings(); mapIndex++) {
					Mapping *mapping = ms->get(mapIndex);
					if (mapping==0) continue;
					Molecule *directMol = mapping->getMolecule();
					if (directMol!=0) directProductMolecules->insert(directMol);
				}
			}
		}
		if (this->usesIncrementalMembership()) {
			for (molIter = products.begin(); molIter != products.end(); ++molIter) {
				Molecule *mol = *molIter;
				if (!mol->isAlive() ||
						isDirectProductMolecule(mol, directProductsPrepared))
					continue;
				MoleculeType *mt = mol->getMoleculeType();
				if (indirectMembershipDecisions.find(mt) ==
						indirectMembershipDecisions.end()) {
					indirectMembershipDecisions.emplace(
							mt, mt->canSkipIndirectMembership(this));
				}
			}
		}
		if (useConnectivity) {
			for (molIter = products.begin(); molIter != products.end(); molIter++) {
				Molecule *mol = *molIter;
				if (!mol->isAlive()) continue;
				if (!isDirectProductMolecule(mol, directProductsPrepared)) {
					hasIndirectProducts = true;
					break;
				}
			}
		}
	}

	// One flag controls all per-fire profiler work below. Without -profile,
	// hot product and local-function loops never enter profiler wrappers.
	bool profileReaction = system->isProfileReactionActive();

	// if complex bookkeeping is on, find all product complexes
	// (this is useful for updating Species Observables and TypeII functions, so keep the info handy).
	// NOTE: this is a brute force approach: check complex of each molecule. there may be a more
	//  elegant way to do this, but it's tricky to get it right.
	if (system->isUsingComplex()) {
		std::unordered_set<Complex*> productComplexSet;
		Complex * complex;
		for ( molIter = products.begin(); molIter != products.end(); molIter++ ) {
			// skip dead molecules
			if ( ! (*molIter)->isAlive() ) continue;
			// get complexID and check if we've already updated that complex
			complex = (*molIter)->getComplex();
			if ( productComplexSet.insert(complex).second )
				productComplexes.push_back(complex);
		}
		if (profileReaction) {
			unsigned long long productComplexMolecules = 0;
			for (std::unordered_set<Complex*>::const_iterator it = productComplexSet.begin();
					it != productComplexSet.end(); ++it) {
				productComplexMolecules +=
						static_cast<unsigned long long>((*it)->getComplexSize());
			}
			system->recordProfileAffectedComplexes(
					static_cast<unsigned long long>(productComplexSet.size()),
					productComplexMolecules);
		}
	}


	// If we're handling observables on the fly, tell each molecule to add itself to observables.
	if (onTheFlyObservables) {
		bool profileObservables = profileReaction;
		ProfileTime profileObservablesStart = profileObservables
			? profileNow() : ProfileTime();
		unsigned long long profileObservableMolecules = profileObservables
			? static_cast<unsigned long long>(products.size()) : 0;

		// molecule observables..
		for ( molIter = products.begin(); molIter != products.end(); molIter++ ) {
			// skip dead molecules
			if ( ! (*molIter)->isAlive() ) continue;
			(*molIter)->addToObservables();
		}

		// species observables..
		if (system->getNumOfSpeciesObs()>0) {
			Complex * c;
			int matches;
			// we can assume that complex bookkeeping is enabled..
			for ( complexIter = productComplexes.begin(); complexIter != productComplexes.end(); ++complexIter ) {
				// update all species observables for this complex
				c = *complexIter;
				matches = 0;
				for ( int i=0; i < system->getNumOfSpeciesObs(); i++ ) {
					matches = system->getSpeciesObs(i)->isObservable(c);
					if (matches > 0) {
						system->getSpeciesObs(i)->add(matches);
					}
				}
			}

			// NOTE: we don't need to handle added population types separately since they are
			//  among the product molecules
		}
		if (profileObservables)
			system->recordProfileObservableAddition(
					profileElapsedSeconds(profileObservablesStart),
					profileObservableMolecules);
	}

	// Now update reaction membership, functions, and update any DOR Groups
	//  also, gather a list of typeII dependencies that will require updating
	typeII_products.clear();
	std::unordered_set<MoleculeType*> typeIIProductSet;
	bool profileMembership = profileReaction;
	ProfileTime profileMembershipStart = profileMembership
		? profileNow() : ProfileTime();
	bool deferMembershipPropensityUpdates =
		directProductsPrepared && this->usesIncrementalMembership() &&
		!useConnectivity;
	if (deferMembershipPropensityUpdates)
		this->system->beginDeferredMembershipPropensityUpdates();
	for ( molIter = products.begin(); molIter != products.end(); molIter++ ) {
		Molecule * mol = *molIter;
		MoleculeType * mt = mol->getMoleculeType();

		// If this moleculeType has typeII dependencies, add it to the list
		// (do this for alive and dead molecules, since molecule deletion may influence
		//    the value of a local function)
		if ( mt->getNumOfTypeIIFunctions() > 0 ) {
			if ( typeIIProductSet.insert(mt).second )
				typeII_products.push_back( mt );
		}

		//Update this molcule's reaction membership
		//  NOTE: as a side-effect, DORreactions that depend on molecule-scoped local functions
		//   (typeI relationship) will be updated as long as UTL is set appropriately.
		if ( mol->isAlive() ) {
			bool useConnectedUpdate =
					useConnectivity &&
					!hasIndirectProducts &&
					isDirectProductMolecule(mol, directProductsPrepared);
			bool directProduct = !trackDirectProducts ||
					isDirectProductMolecule(mol, directProductsPrepared);
			if (!directProduct) {
				auto decision = indirectMembershipDecisions.find(mt);
				if (decision != indirectMembershipDecisions.end() &&
						decision->second)
					continue;
			}
			mol->updateRxnMembership(this, useConnectedUpdate, directProduct);
		}
	}
	if (deferMembershipPropensityUpdates)
		this->system->endDeferredMembershipPropensityUpdates();
	if (profileMembership)
		system->recordProfileMembershipPhase(
			profileElapsedSeconds(profileMembershipStart));

	// update complex-scoped local functions for typeII dependencies
	// NOTE: as a side-effect, dependent DOR reactions (via typeI molecule dependencies) will be updated

	unsigned long long productComponentMoleculeCount = 0;
	bool productComponentSizesValid = !productComponentSizes.empty();
	for (vector<unsigned int>::iterator componentSize = productComponentSizes.begin();
			componentSize != productComponentSizes.end(); ++componentSize) {
		if (*componentSize == 0)
			productComponentSizesValid = false;
		productComponentMoleculeCount += *componentSize;
	}
	bool canReuseProductComponents =
			!transformationSet->hasTopologyChangingTransform() &&
			!productComponentsTruncated &&
			productComponentSizesValid &&
			productComponentMoleculeCount == products.size();

	// No Type-II local function depends on these products: avoid traversing
	// connected components solely to execute an empty update loop.
	if (!typeII_products.empty()) {
		// for each typeII product molecule, update all dependent local functions
		if (system->isUsingComplex()) {
			// this is the easy way: update all typeI molecules on each complex
			for ( typeII_iter = typeII_products.begin(); typeII_iter != typeII_products.end(); ++typeII_iter ) {
				MoleculeType * mt = *typeII_iter;
				for (int i=0; i < mt->getNumOfTypeIIFunctions(); i++) {
					for ( complexIter = productComplexes.begin(); complexIter != productComplexes.end(); ++complexIter )
						mt->getTypeIILocalFunction(i)->evaluateOn( *complexIter );
				}
			}
		}
		else {
			if (canReuseProductComponents) {
				// State-only rules leave the product components unchanged. Reuse the
				// component boundaries collected during product preparation.
				list<Molecule*>::iterator componentMolecule = products.begin();
				for (vector<unsigned int>::iterator componentSize = productComponentSizes.begin();
						componentSize != productComponentSizes.end(); ++componentSize) {
					list<Molecule*> connectedMols;
					Molecule *mol = *componentMolecule;
					for (unsigned int i = 0; i < *componentSize; ++i) {
						connectedMols.push_back(*componentMolecule);
						++componentMolecule;
						if (profileReaction)
							system->recordProfileLocalFunctionComponentCandidate(i != 0);
					}

					for ( typeII_iter = typeII_products.begin(); typeII_iter != typeII_products.end(); ++typeII_iter ) {
						MoleculeType * mt = *typeII_iter;
						for (int i=0; i<mt->getNumOfTypeIIFunctions(); i++)
							mt->getTypeIILocalFunction(i)->evaluateOn( mol, connectedMols );
					}
				}
			}
			else {
				// this is the hard way: find a representative molecule from each connected set
				//  and evaluate TypeII functions on that representative.
				std::unordered_set<Molecule*> allMols;
				Molecule * mol;
				for ( molIter = products.begin(); molIter != products.end(); molIter++ ) {
					mol = *molIter;
					bool isNewComponent = allMols.insert(mol).second;
					if (profileReaction)
						system->recordProfileLocalFunctionComponentCandidate(!isNewComponent);
					if ( isNewComponent ) {
						// remember everything connected to this molecule so we don't
						// evaluate this connected set multiple times.
						list<Molecule*> connectedMols;
						ProfileConnectivityScope profileConnectivityScope(
							system, PROFILE_CONNECTIVITY_LOCAL_FUNCTION);
						mol->traverseBondedNeighborhood( connectedMols, ReactionClass::NO_LIMIT );
						for ( list<Molecule*>::iterator cm = connectedMols.begin(); cm != connectedMols.end(); ++cm )
							allMols.insert(*cm);

						// evaluate typeII local functions on this connected set
						for ( typeII_iter = typeII_products.begin(); typeII_iter != typeII_products.end(); ++typeII_iter ) {
							MoleculeType * mt = *typeII_iter;
							for (int i=0; i<mt->getNumOfTypeIIFunctions(); i++)
								mt->getTypeIILocalFunction(i)->evaluateOn( mol, connectedMols );
						}
					}
				}
			}
		}
	}

	// update the last reaction firing time
	// this is written to molecule_type_list.tsv at the end of the simulation
	// @author: Arvind R. Subramaniam
	// @date: 13 Nov 2019
	this->system->setLastRxnTime(this->system->getCurrentTime());
	// output to a JSON if the reaction was tagged
	if (this->system->getReactionTrackingStatus()) {
		if (tagged && track) {
			string track_str = "";
			int level = 6; // indentation level
			this->system->current_cpu_time = ((double) (clock() - this->system->start) / (double) CLOCKS_PER_SEC);
			// we need the correct number of commas
			if (this->system->getGlobalEventCounter() != 1) {
				track_str += ",\n";
			} 
			// open firing and write info
			track_str += std::string(level,' ') + "{\n" +
				std::string(level+2,' ') + "\"props\": [";
			if (this->system->getRxnNumberTrack()) {
				track_str += string("\"") + to_string(rxnId) + "\",";
			} else {
				track_str += string("\"") + name + "\",";
			}
			track_str += to_string(this->system->getGlobalEventCounter()) +
					"," + to_string(this->system->getCurrentTime()) + "],\n";

			// write transformation log
			track_str += logstr;
					
			// close firing 
			track_str += std::string(level,' ') + "}";
			//Tidy up
			products.clear();
			productComplexes.clear();
			return track_str;
		}
	}
	//Tidy up
	products.clear();
	productComplexes.clear();
	// AS2023 - returning empty, if we are here logging was off
	return "";
}

void ReactionClass::identifyConnectedReactions() {
	ReactionClass * rxn;
	vector <ReactionClass *> allReactions;
	allReactions = system->getAllReactions();
	for (unsigned int r=0, n=allReactions.size(); r<n; r++) {
		rxn = allReactions.at(r);
		if (this->isReactionConnected(rxn)) this->appendConnectedRxn(rxn);
	}
}

bool ReactionClass::areMoleculeTypeAndComponentPresent(MoleculeType * mt, int cIndex) {
	TemplateMolecule * t2;
	for (unsigned int i=0, n=allReactantTemplates.size(); i<n; i++) {
		t2 = allReactantTemplates[i];
		if (t2->isMoleculeTypeAndComponentPresent(mt, cIndex)) return true;
	}

	return false;
}

bool ReactionClass::isTemplateCompatible(TemplateMolecule * t) {
	TemplateMolecule * t2;
	for (unsigned int i=0, n=allReactantTemplates.size(); i<n; i++) {
		t2 = allReactantTemplates[i];
		if (t->isTemplateCompatible(t2)) return true;
	}

	return false;
}
