



#include "reaction.hh"


using namespace std;
using namespace NFcore;




FunctionalRxnClass::FunctionalRxnClass(string name, GlobalFunction *gf, TransformationSet *transformationSet, System *s) :
	BasicRxnClass(name,1,"",transformationSet,s)
{
	this->reactionType = ReactionClass::OBS_DEPENDENT_RXN;
	this->cf=0;
	this->gf=gf;
	for(int vr=0; vr<gf->getNumOfVarRefs(); vr++) {
		if(gf->getVarRefType(vr)=="Observable") {
			Observable *obs = s->getObservableByName(gf->getVarRefName(vr));
			obs->addDependentRxn(this);
		} else {
			cerr<<"When creating a FunctionalRxnClass of name: "+name+" you provided a function that\n";
			cerr<<"depends on an observable type that I can't yet handle! (which is "+gf->getVarRefType(vr)+"\n";
			cerr<<"try using type: 'MoleculeObservable' for now.\n";
			cerr<<"quiting..."<<endl; exit(1);
		}
	}
}

FunctionalRxnClass::FunctionalRxnClass(string name, CompositeFunction *cf, TransformationSet *transformationSet, System *s) :
	BasicRxnClass(name,1, "", transformationSet,s)
{
	this->reactionType = ReactionClass::OBS_DEPENDENT_RXN;
	this->gf=0;
	this->cf=cf;
	this->cf->setGlobalObservableDependency(this,s);
}


FunctionalRxnClass::~FunctionalRxnClass() {};

double FunctionalRxnClass::update_a() {
	//cout<<"udpating a"<<endl;
	if(this->onTheFlyObservables==false) {
		cerr<<"Warning!!  You have on the fly observables turned off, but you are using functional\n";
		cerr<<"reactions which depend on observables.  Therefore, you cannot turn off onTheFlyObservables!\n";
		cerr<<"exiting now."<<endl;
		exit(1);
	}



	//	cout<<"here"<<endl;
	if(gf!=0) {
	//	cout<<"in here"<<endl;
		// AS-2021
		if(gf->fileFunc==true) {
			gf->fileUpdate();
		}
		// AS-2021
		a=FuncFactory::Eval(gf->p);
	} else if(cf!=0) {
		if (reactantCountBuffer.size() != n_reactants) {
			reactantCountBuffer.resize(n_reactants);
		}
		for(unsigned int r=0; r<n_reactants; r++) {
			reactantCountBuffer[r] = (int)getReactantCount(r);
		}
		int *reactantCountsPtr = (n_reactants > 0) ? &reactantCountBuffer[0] : 0;
		a=cf->evaluateOn(0,0, reactantCountsPtr, n_reactants);
	//	cout<<"and here"<<endl;
	} else {
		cout<<"Error!  Functional rxn is not properly initialized, but is being used!"<<endl;
		exit(1);
	}

	a *= this->volumeConversionFactor;

	// Scale by baseRate, exactly as BasicRxnClass/DORRxnClass/DOR2RxnClass do.
	// A FunctionalRxnClass is constructed with baseRate=1 and never routes
	// through setBaseRate(), so this factor is the reaction center symmetry
	// correction and nothing else.  Without it a symmetric rule with a
	// functional rate fires at 1/symmetryFactor times its intended rate.
	//
	// Not under TotalRate, though.  The symmetry factor corrects a counting
	// problem -- a reactant pattern with a non-trivial automorphism matches the
	// same reaction more than once -- and TotalRate states the whole propensity
	// of the rule outright, so there is no count to correct.  It has to be
	// skipped here rather than where the factor is folded into baseRate:
	// NFinput calls setTotalRateFlag() well after both this class' constructor
	// and setBaseRate(), so totalRateFlag is still false at both of those
	// points and a guard placed there would never fire.
	if (!this->totalRateFlag)
	{	a *= this->baseRate;   }

	if(a<0) {
		cout<<"Warning!!  The function you provided for functional rxn: '"<<name<<"' evaluates\n";
		cout<<"to a value less than zero!  You cannot have a negative propensity (here "<<a<<")!\n";
		cout<<"here is the offending function: \n";
		// A rate law that is a function call in an expression (k*f()) is a
		// CompositeFunction, and gf is null for it: printing gf crashed here.
		// The evaluation above has already refused a reaction with neither.
		if (gf != 0)
			gf->printDetails(this->system);
		else
			cf->printDetails(this->system);
		cout<<"\nhere is the offending reaction: \n";
		this->printDetails();
		cout<<"\n\nquitting."<<endl;
		exit(1);
	}


	// check here for the total rate flag - if this is set to true, then
	// use the rate exactly as given by the function, but if it is false,
	// then we have to multiply here by the reactant counts
	if(!this->totalRateFlag) {
		for(unsigned int i=0; i<n_reactants; i++)
			a*=(double)getCorrectedReactantCount(i);
	}
	else
	{
		// Check that we have at least one set of reactants!
		for(unsigned int i=0; i<n_reactants; i++) {
			if(getCorrectedReactantCount(i)==0) {
				a=0.0;
				break;
				//cout<<"Warning!  Function evaluates to positive rate for a reaction, but"<<endl;
				//cout<<"one of the reactant lists is empty!"<<endl;
				//this->printDetails();
				//cf->printDetails(reactantTemplates[0]->getMoleculeType()->getSystem());
				//exit(1);
			}
		}
	}
	
	return a;
}

void FunctionalRxnClass::printDetails() const {

	string trate = "off";
	if(this->totalRateFlag) trate = "on";

	if(gf!=0) {
		// AS-2021
		if (gf->fileFunc==true) {
			gf->fileUpdate();
		}
		// AS-2021
		cout<<"ReactionClass: " << name <<"  ( baseFunction="<<gf->getNiceName()<<"="<<FuncFactory::Eval(gf->p)<<",  a="<<a<<", fired="<<fireCounter<<" times, TotalRate="<<trate<<" )"<<endl;
	} else if(cf!=0) {
		int * reactantCounts = new int[this->n_reactants];
		for(unsigned int r=0; r<n_reactants; r++) {
			reactantCounts[r]=getReactantCount(r);
		}
		double value=cf->evaluateOn(0,0, reactantCounts, n_reactants);
		delete [] reactantCounts;
		cout<<"ReactionClass: " << name <<"  ( baseFunction="<<cf->getName()<<"="<<value<<",  a="<<a<<", fired="<<fireCounter<<" times, TotalRate="<<trate<<" )"<<endl;

	}

	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
		cout<<"	(count="<< this->getReactantCount(r) <<")."<<endl;
	}
	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}









MMRxnClass::MMRxnClass(string name, double kcat, double Km, TransformationSet *transformationSet,System *s,
		string kcatParameterName, string KmParameterName) :
	BasicRxnClass(name,1,"",transformationSet,s)
{
	this->Km = Km;
	this->kcat = kcat;
	this->sFree=0;
	this->kcatParameterName = kcatParameterName;
	this->KmParameterName = KmParameterName;
	if(n_reactants!=2) {
		cerr<<"You have tried to create a reaction with a Michaelis-Menten rate law (named: '"+name+"'\n')";
		cerr<<"but you don't have the correct number of reactants!  Michaelis-Menten reactions require\n";
		cerr<<"exactly 2 reactants.  A substrate (always given first) and an enzyme (always given second)\n";
		cerr<<"Read your tutorial next time... now I will quit."<<endl;
		exit(1);
	}
}
MMRxnClass::~MMRxnClass() {};

// kcat and Km are held apart from the base rate, which for this class carries
// only the symmetry factor (and any volume factor) and has no parameter name,
// so the base class refresh leaves them as they were read at load. Read them
// again from their parameters.
void MMRxnClass::resetBaseRateFromSystemParamter() {
	ReactionClass::resetBaseRateFromSystemParamter();
	if (!kcatParameterName.empty()) kcat = system->getParameter(kcatParameterName);
	if (!KmParameterName.empty()) Km = system->getParameter(KmParameterName);
	update_a();
}

double FunctionalRxnClass::exactRuleMonkey_a()
{
	return update_a();
}

double MMRxnClass::exactRuleMonkey_a()
{
	return update_a();
}

double MMRxnClass::update_a()
{
	// See FunctionalRxnClass::update_a() -- an MMRxnClass is likewise built with
	// baseRate=1 and never calls setBaseRate(), so baseRate carries the reaction
	// center symmetry correction.  BNG emits symmetry_factor with an MM rate law
	// whenever the substrate pattern has a non-trivial automorphism.
	//
	// It belongs on the substrate COUNT, inside the law -- not on the finished
	// propensity.  What the factor corrects is a match multiplicity:
	// getCorrectedReactantCount(0) counts pattern embeddings, and a substrate
	// pattern with a non-trivial automorphism matches each complex more than
	// once, so the law is being handed more substrate than exists.  Michaelis-
	// Menten is not linear in that count, so scaling the finished propensity
	// instead is exact only below saturation; scaling the count is exact
	// everywhere, and the two agree wherever the law is linear.
	//
	// Only the substrate needs it. The enzyme is pure context here -- an MM rule
	// does not transform it -- and BNG's MultScale counts reaction-center
	// automorphisms, so it emits symmetry_factor=1 for an enzyme-side symmetry
	// and this factor is always the substrate's. (NFsim does over-count a
	// symmetric context pattern, but that is a separate defect with no
	// symmetry_factor attached to it; see issue #87.)
	double S = (double)getCorrectedReactantCount(0) * this->baseRate;
	double E = (double)getCorrectedReactantCount(1);
	// Free substrate is the non-negative root of sFree^2 - b*sFree - Km*S = 0, b = S-Km-E.
	// The textbook form 0.5*(b + sqrt(b*b + 4*Km*S)) cancels catastrophically when b < 0,
	// and rounds to exactly zero once 4*Km*S drops below about 1e-16*b*b, which silently
	// stops the reaction when Km is small and the enzyme is in excess. The two roots
	// multiply to -Km*S, so the b < 0 case is written as a sum of like-signed quantities.
	// See https://github.com/RuleWorld/bionetgen/issues/323
	double b = S-Km-E;
	double q = sqrt(b*b + 4.0*Km*S);
	sFree = (b>=0.0) ? 0.5*(b+q) : 2.0*Km*S/(q-b);
	a=kcat*sFree*E/(Km+sFree);
	return a;
}

void MMRxnClass::printDetails() const {
	cout<<"ReactionClass: " << name <<"  ( Km="<<Km<<", kcat="<<kcat<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
		cout<<"	(count="<< this->getReactantCount(r) <<")."<<endl;
	}
	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}









/* --- pure context counting, shared by every reaction class ------------------
 *
 * BioNetGen gives a reactant pattern the rule does not transform one reaction
 * instance per matching COMPLEX, however many molecules in that complex match
 * it: every embedding produces the identical reaction, so there is only one.
 * NFsim enumerates matches per molecule, so these collapse the difference.
 * Both read the complex off mapping 0, which for such a reactant is the root
 * molecule placed by the placeholder transform TransformationSet::finalize()
 * adds. */

/* Scratch buffer for the counts below.  These run on the propensity hot path --
 * once per reactant per update_a(), which is millions of calls on a model with a
 * catalytic rule -- so they must not allocate.  A reused vector plus sort beats
 * building a std::set by more than an order of magnitude, and thread_local keeps
 * concurrent sessions from sharing it. */
static thread_local std::vector<int> s_complexScratch;

template <typename Container>
static int distinctComplexesOf(Container *rc, int size)
{
	s_complexScratch.clear();
	if ((int)s_complexScratch.capacity() < size) s_complexScratch.reserve(size);
	for (int i = 0; i < size; ++i) {
		MappingSet *ms = rc->getMappingSetByIndex(i);
		if (ms && ms->getNumOfMappings() > 0) {
			Mapping *mapping = ms->get(0);
			if (mapping && mapping->getMolecule()) {
				s_complexScratch.push_back(mapping->getMolecule()->getComplexID());
			}
		}
	}
	std::sort(s_complexScratch.begin(), s_complexScratch.end());
	return (int)(std::unique(s_complexScratch.begin(), s_complexScratch.end())
	             - s_complexScratch.begin());
}

int NFcore::countDistinctComplexes(ReactantList *rl)
{
	if (!rl->mayShareComplexes()) return rl->size();
	return distinctComplexesOf(rl, rl->size());
}

int NFcore::countDistinctComplexes(ReactantTree *tree)
{
	if (!tree->mayShareComplexes()) return tree->size();
	return distinctComplexesOf(tree, tree->size());
}

double NFcore::perComplexRateFactorSum(ReactantTree *tree)
{
	// One representative term per complex.  getRateFactor() is indexed by the
	// same flat array position as getMappingSetByIndex(), so this walks the live
	// entries and keeps the first seen for each complex.
	std::set<int> seen;
	double sum = 0.0;
	int size = tree->size();
	for (int i = 0; i < size; ++i) {
		MappingSet *ms = tree->getMappingSetByIndex(i);
		if (!ms || ms->getNumOfMappings() == 0) continue;
		Mapping *mapping = ms->get(0);
		if (!mapping || !mapping->getMolecule()) continue;
		if (seen.insert(mapping->getMolecule()->getComplexID()).second) {
			sum += tree->getRateFactor(i);
		}
	}
	return sum;
}

void NFcore::collectReactantRepresentatives(ReactantList *rl, bool perComplex,
                                            std::vector<MappingSet*> &out,
                                            std::vector<int> *flatIndices)
{
	out.clear();
	if (flatIndices) flatIndices->clear();
	int size = rl->size();

	// Nothing to collapse: either this reactant is transformed by the rule, or
	// every match sits in a single-molecule complex and one match already is one
	// complex.  Enumerate every live mapping set, exactly as before.
	if (!perComplex || !rl->mayShareComplexes()) {
		for (int i = 0; i < size; ++i) {
			MappingSet *ms = rl->getMappingSetByIndex(i);
			if (ms) {
				out.push_back(ms);
				if (flatIndices) flatIndices->push_back(i);
			}
		}
		return;
	}

	std::set<int> seen;
	for (int i = 0; i < size; ++i) {
		MappingSet *ms = rl->getMappingSetByIndex(i);
		if (!ms || ms->getNumOfMappings() == 0) continue;
		Mapping *mapping = ms->get(0);
		if (!mapping || !mapping->getMolecule()) continue;
		if (seen.insert(mapping->getMolecule()->getComplexID()).second) {
			out.push_back(ms);
			if (flatIndices) flatIndices->push_back(i);
		}
	}
}


BasicRxnClass::BasicRxnClass(string name, double baseRate, string baseRateName, TransformationSet *transformationSet, System *s) :
	ReactionClass(name,baseRate,baseRateName,transformationSet,s)
{
	this->reactionType = BASIC_RXN;  //set as normal reaction here, but deriving reaction classes can change this
	reactantLists = new ReactantList *[n_reactants];
	//Set up the reactantLists
	for(unsigned int r=0; r<n_reactants; r++)
		reactantLists[r]=(new ReactantList(r,transformationSet,25,this->system));
	
	this->connectivityFlag = s->getConnectivityFlag();
	
	msPairBuffer = new MappingSet*[2];
}


BasicRxnClass::~BasicRxnClass()
{
    //cout<<"  -------------------------------\n  ----------------------------\n";
	//cout<<"Reaction: "<<name<<endl;
	//this->reactantLists[0]->printDetails();

	//this->reactantLists[0]->removeMappingSet(0);
	//this->reactantLists[0]->removeMappingSet(3);
	//this->reactantLists[0]->removeMappingSet(6);
	//this->reactantLists[0]->removeMappingSet(9);

	//cout<<endl<<endl<<endl;
	//this->reactantLists[0]->printDetails();




	if(DEBUG) cout<<"Destroying rxn: "<<name<<endl;

	for(unsigned int r=0; r<n_reactants; r++)
	{
		//delete reactantTemplates[r]; DO NOT DELETE HERE (MoleculeType has responsibility of
		//deleting all template molecules of its type now.
		delete reactantLists[r];
	}
	delete [] reactantLists;
	
	delete [] msPairBuffer;
}

void BasicRxnClass::init()
{
	for(unsigned int r=0; r<n_reactants; r++)
	{
		reactantTemplates[r]->getMoleculeType()->addReactionClass(this,r);
	}
}


void BasicRxnClass::prepareForSimulation()
{

}


int BasicRxnClass::checkForEquality(Molecule *m, MappingSet* ms, int rxnIndex, ReactantList* reactantList){
	/*
	Check if mapping set clashes with any of the mapping sets already in reactantList
	*/
	const MappingIdSet& tempSet = m->getRxnListMappingSet(rxnIndex);
	for(MappingIdSet::const_iterator it= tempSet.begin();it!= tempSet.end(); ++it){
		MappingSet* ms2 = reactantList->getMappingSet(*it);
		if(MappingSet::checkForEquality(ms,ms2)){
			return *it;
		}
	}
	return -1;


}

/**
 * Updates a molecule's reaction membership and a reaction's reactant list
 *
 * If a molecule matches the TemplateMolecule of a reaction, then it gets added.
 * Note that NFsim has a single arbitrary TemplateMolecule for each reactant
 * complex for a reaction. So it assumes that you will always check the molecule
 * that matches the designated TemplateMolecule of a reaction for updates. So if
 * I modify the identification of reactant molecules in
 * TransformationSet::getListOfProducts(), I need to make sure that every molecule
 * that is part of a reaction whose propensity might change is included in the
 * listOfProducts.
 *
 * This is the gateway function to TemplateMolecule::compare().
 * @author Arvind Rasi Subramaniam
 * @param m - molecule that is being compared against the TemplateMolecule of a reaction
 * @param reactantPos - the reactant number among the reaction's reacant complexes
 * (note that every connected complex gets only one reactant number.)
 * @return true if there are no errors.
 */
bool BasicRxnClass::tryToAdd(Molecule *m, unsigned int reactantPos)
{
	if (system != 0 && system->isProfilingEnabled())
		system->recordProfileMatchCandidate();

	//First a bit of error checking, that you should skip unless we are debugging...
	//	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL)
	//	{
	//		cout<<"Error adding molecule to reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
	//		exit(1);
	//	}

	//Get the specified reactantList
	ReactantList *rl = reactantLists[reactantPos];

	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
	//cout<<" got mappingSetId: " << m->getRxnListMappingId(rxnIndex)<<" size: " <<rl->size()<<endl;
	//cout<< " testing whether to add molecule ";
	//m->printDetails();
	//cout<<" ... as a mormal reaction "<<this->name<<endl;


	//If this reaction has multiple instances, we always remove them all!
	// then we remap because other mappings may have changed.  Yes, this may
	// be more ineffecient, but it is the fast implementation
	if(rl->getHasClonedMappings()) {
		while(m->getRxnListMappingId(rxnIndex)>=0) {
			rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->deleteRxnListMappingId(rxnIndex,m->getRxnListMappingId(rxnIndex));
			//m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}
	}

	//Here we get the standard update...
	MappingIdSet deleteMs = m->getRxnListMappingSet(rxnIndex);

	// Cheap note for countDistinctComplexes(): while every matched molecule is a
	// singleton complex, no complex can hold two matches and the distinct-complex
	// count is just size().
	if (contextCountsPerComplex[reactantPos]) {
		rl->noteMappedComplexSize(m->getComplex()->getComplexSize());
	}

	//Try to map it!
	MappingSet *ms = rl->pushNextAvailableMappingSet();
	symmetricMappingSet.clear();
	comparisonResult = reactantTemplates[reactantPos]->compare(m,rl,ms,false,&symmetricMappingSet);
	if(!comparisonResult) {
		//cout << "no mapping in normal reaction, remove"<<endl;
		//we must remove, if we did not match.  This will also remove
		//everything that was cloned off of the mapping set
		rl->removeMappingSet(ms->getId());
		//JJT: removes any symmetric mapping sets that might have been added since we are not using them
		for(vector<MappingSet *>::iterator it=symmetricMappingSet.begin();it!=symmetricMappingSet.end();++it){
			rl->removeMappingSet((*it)->getId());
		}
	} else {
		//cout << "should be in normal reaction, confirm push"<<endl;
		//ms->printDetails();
		
		if (symmetricMappingSet.size() > 0){
            rl->removeMappingSet(ms->getId());
			for(vector<MappingSet *>::iterator it=symmetricMappingSet.begin();it!=symmetricMappingSet.end();++it){
					int mapIndex = checkForEquality(m,*it,rxnIndex,rl);
					if(mapIndex >= 0){
						deleteMs.erase(mapIndex);
						rl->removeMappingSet((*it)->getId());
					}
					else{
						m->setRxnListMappingId(rxnIndex,(*it)->getId());
					}
            }
		}
		else{
			int mapIndex = checkForEquality(m,ms,rxnIndex,rl);
			if(mapIndex >= 0){
				deleteMs.erase(mapIndex);
				rl->removeMappingSet(ms->getId());
			}
			else{
				m->setRxnListMappingId(rxnIndex,ms->getId());
			}
		}
		
	}

	for (MappingIdSet::iterator it = deleteMs.begin(); it != deleteMs.end(); ++it) {
		rl->removeMappingSet(*it);
		m->deleteRxnListMappingId(rxnIndex, *it);
	}

	return true;
}



void BasicRxnClass::remove(Molecule *m, unsigned int reactantPos)
{

	//First a bit of error checking...
	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL)
	{
		cout<<"Error removing molecule from a reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
		exit(1);
	}


	//Get the specified reactantList
	ReactantList *rl = reactantLists[reactantPos];

	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
	bool isInRxn = (m->getRxnListMappingId(rxnIndex)>=0);


	if(isInRxn)
	{
		rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
		m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
	}
}




void BasicRxnClass::notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex)
{
	cerr<<"You are trying to notify a Basic Reaction of a rate Factor Change!!! You should only use this"<<endl;
	cerr<<"function for DORrxnClass rules!  For this offense, I must abort now."<<endl;
	exit(1);
}

double BasicRxnClass::exactRuleMonkey_a()
{
	if(this->totalRateFlag) {
		double exact_a = baseRate;
		for(unsigned int i=0; i<n_reactants; i++) {
			if(getCorrectedReactantCount(i)==0) exact_a = 0.0;
		}
		return exact_a;
	}

	double validCombinations = 0.0;
	if (n_reactants == 0) {
		validCombinations = 1.0;
	} else if (n_reactants == 1) {
		validCombinations = getCorrectedReactantCount(0);
	} else if (n_reactants == 2) {
		// Exact calculation: subtract null events.  Enumerate one representative
		// per complex for a reactant the rule does not transform, so both the
		// total and the invalid pairs subtracted from it are counted on the same
		// footing as getCorrectedReactantCount(); for every other reactant this
		// is still every live mapping set, and the arithmetic is unchanged.
		static thread_local std::vector<MappingSet*> reps0, reps1;
		collectReactantRepresentatives(reactantLists[0], contextCountsPerComplex[0], reps0);
		collectReactantRepresentatives(reactantLists[1], contextCountsPerComplex[1], reps1);

		// Use raw counts here because invalid self-pairs are removed explicitly below.
		double totalCombinations = (double)reps0.size() * (double)reps1.size();
		double invalidCombinations = 0;

		for (size_t i = 0; i < reps0.size(); ++i) {
			msPairBuffer[0] = reps0[i];
			for (size_t j = 0; j < reps1.size(); ++j) {
				msPairBuffer[1] = reps1[j];

				// check for collision
				if (!transformationSet->checkMolecularity(msPairBuffer)) {
					invalidCombinations++;
				}
			}
		}
		validCombinations = totalCombinations - invalidCombinations;
		if (validCombinations < 0) validCombinations = 0;
	} else {
		// fallback to standard approximation
		validCombinations = 1.0;
		for(unsigned int i=0; i<n_reactants; i++) {
			validCombinations *= getCorrectedReactantCount(i);
		}
	}

	return validCombinations * baseRate;
}

double BasicRxnClass::update_a()
{
	if (useRuleMonkey) {
		a = exactRuleMonkey_a();
		return a;
	}
	// Use the total rate law convention (macroscopic rate)
	if(this->totalRateFlag) {
		a=baseRate;
		for(unsigned int i=0; i<n_reactants; i++)
			if(getCorrectedReactantCount(i)==0) a=0.0;

	// Use the standard microscopic rate
	} else {
		a = 1.0;
		for(unsigned int i=0; i<n_reactants; i++) {
			a*=getCorrectedReactantCount(i);
		}
		a*=baseRate;
	}
	return a;
}


int BasicRxnClass::getReactantCount(unsigned int reactantIndex) const
{
	return isPopulationType[reactantIndex] ?
			   reactantLists[reactantIndex]->getPopulation()
			 : reactantLists[reactantIndex]->size();
}


int BasicRxnClass::getCorrectedReactantCount(unsigned int reactantIndex) const
{
	/*cerr << "  getCorrectedReactantCount rindex: " << reactantIndex << "  isPop? " << isPopulationType[reactantIndex] << endl;
	if ( isPopulationType[reactantIndex] )
	{
		cerr << "  corr:  " << identicalPopCountCorrection[reactantIndex];
		cerr << "  pop:   " << reactantLists[reactantIndex]->getPopulation() << endl;
		cerr << "  final: " << std::max( reactantLists[reactantIndex]->getPopulation()
	             - identicalPopCountCorrection[reactantIndex], 0 ) << endl;
	}
	else
	{
		cerr << "  count: " << reactantLists[reactantIndex]->size() << endl;
	}
	*/

	// MatchOnce is the user's explicit request; contextCountsPerComplex is the
	// same counting applied automatically to a reactant the rule never transforms,
	// which is how BNG counts one.
	if ((matchOncePerReactant[reactantIndex] || contextCountsPerComplex[reactantIndex])
	    && !isPopulationType[reactantIndex]) {
		return countDistinctComplexes(reactantLists[reactantIndex]);
	}

	return isPopulationType[reactantIndex] ?
			   std::max( reactantLists[reactantIndex]->getPopulation()
			             - identicalPopCountCorrection[reactantIndex], 0 )
			 : reactantLists[reactantIndex]->size();
}

void BasicRxnClass::printFullDetails() const
{
	cout<<"BasicRxnClass: "<<name<<endl;
	for(unsigned int i=0; i<n_reactants; i++)
		reactantLists[i]->printDetails();
}


void BasicRxnClass::pickRuleMonkeyMappingSets(double random_A_number) const
{
	if (n_reactants != 2 || totalRateFlag) {
		for(unsigned int i=0; i<n_reactants; i++) {
			if ( isPopulationType[i] ) {
				reactantLists[i]->pickRandomFromPopulation(mappingSet[i], system->getRNG());
			} else {
				reactantLists[i]->pickRandom(mappingSet[i], system->getRNG());
			}
		}
		return;
	}

	// For molecularity=2, we have to find a valid pair (no null events)
	int size0 = getReactantCount(0);
	int size1 = getReactantCount(1);
	
	validPairsBuffer.clear();
	for (int i = 0; i < size0; ++i) {
		msPairBuffer[0] = reactantLists[0]->getMappingSet(i);
		for (int j = 0; j < size1; ++j) {
			msPairBuffer[1] = reactantLists[1]->getMappingSet(j);
			
			if (transformationSet->checkMolecularity(msPairBuffer)) {
				validPairsBuffer.push_back(make_pair(i, j));
			}
		}
	}
	
	if (validPairsBuffer.empty()) {
		// Safety fallback: this should be unreachable when exactRuleMonkey_a() > 0.
		// If reached, preserve legacy behavior by allowing the standard chooser path.
		for(unsigned int i=0; i<n_reactants; i++) {
			if ( isPopulationType[i] ) {
				reactantLists[i]->pickRandomFromPopulation(mappingSet[i], system->getRNG());
			} else {
				reactantLists[i]->pickRandom(mappingSet[i], system->getRNG());
			}
		}
		return;
	}
	
	// Select a valid pair
	int selectedIndex = system->getRNG().random_int(0, validPairsBuffer.size());
	int i = validPairsBuffer[selectedIndex].first;
	int j = validPairsBuffer[selectedIndex].second;
	
	mappingSet[0] = reactantLists[0]->getMappingSet(i);
	mappingSet[1] = reactantLists[1]->getMappingSet(j);
}


void BasicRxnClass::pickMappingSets(double random_A_number) const
{
	if (useRuleMonkey) {
		pickRuleMonkeyMappingSets(random_A_number);
		return;
	}
	//Note here that we completely ignore the argument.  The argument is only
	//used for DOR reactions because we need that number to select the reactant to fire

	//Select a reactant from each list
	for(unsigned int i=0; i<n_reactants; i++)
	{
		if ( isPopulationType[i] ) {
			reactantLists[i]->pickRandomFromPopulation(mappingSet[i], system->getRNG());
		} else {
			reactantLists[i]->pickRandom(mappingSet[i], system->getRNG());
		}
	}
}
