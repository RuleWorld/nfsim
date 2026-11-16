#include <iostream>
#include "NFcore.hh"

#if defined(_MSC_VER)
#include <intrin.h>
#endif


using namespace std;
using namespace NFcore;

namespace {

unsigned int compactMembershipBitCount(std::uint64_t value)
{
	unsigned int count = 0;
	while (value != 0) {
		value &= (value - 1);
		++count;
	}
	return count;
}

void setCompactMembershipCandidateBit(
		vector<std::uint64_t> &bits, unsigned int reactionIndex)
{
	const std::size_t wordIndex = reactionIndex >> 6;
	if (bits.size() <= wordIndex)
		bits.resize(wordIndex + 1, 0);
	bits[wordIndex] |=
			(std::uint64_t(1) << (reactionIndex & 63));
}

unsigned int compactMembershipTrailingZeroCount(std::uint64_t value)
{
#if defined(_MSC_VER)
	unsigned long bit = 0;
	_BitScanForward64(&bit, value);
	return static_cast<unsigned int>(bit);
#else
	return static_cast<unsigned int>(__builtin_ctzll(value));
#endif
}

bool compactMembershipBitIsSet(
		const vector<std::uint64_t> &bits, unsigned int reactionIndex)
{
	const std::size_t wordIndex = reactionIndex >> 6;
	return wordIndex < bits.size() &&
			(bits[wordIndex] &
				(std::uint64_t(1) << (reactionIndex & 63))) != 0;
}

}




MoleculeType::MoleculeType(
	string name,
	vector <string> &compName,
	System *s)
 : population_type( false ), isFixed_(false), fixedCount_(0), fixedCompartment_(nullptr)
{
	vector <string> defaultCompState;
	vector < vector <string> > possibleCompStates;
	vector <bool> isIntegerComponent;
	unsigned int compNameSize = compName.size();
	for(unsigned int i=0; i<compNameSize; i++) {
		vector <string> v;
		possibleCompStates.push_back(v);
		defaultCompState.push_back("NO_STATE");
		isIntegerComponent.push_back(false);
	}
	init(name, compName, defaultCompState, possibleCompStates, isIntegerComponent, s);

}

MoleculeType::MoleculeType(
	string name,
	vector <string> &compName,
	vector <string> &defaultCompState,
	System *s)
 : population_type( false ), isFixed_(false), fixedCount_(0), fixedCompartment_(nullptr)
{
	vector < vector <string> > possibleCompStates;
	vector <bool> isIntegerComponent;
	unsigned int compNameSize = compName.size();
	for(unsigned int i=0; i<compNameSize; i++) {
		vector <string> v;
		possibleCompStates.push_back(v);
		isIntegerComponent.push_back(false);
	}
	init(name, compName, defaultCompState, possibleCompStates, isIntegerComponent, s);
}







MoleculeType::MoleculeType(
		string name,
		vector <string> &compName,
		vector <string> &defaultCompState,
		vector < vector<string> > &possibleCompStates,
		System *system)
 : population_type( false ), isFixed_(false), fixedCount_(0), fixedCompartment_(nullptr)
{
	vector <bool> isIntegerComponent;
	unsigned int compNameSize = compName.size();
	for(unsigned int i=0; i<compNameSize; i++) {
		isIntegerComponent.push_back(false);
	}
	init(name, compName, defaultCompState, possibleCompStates, isIntegerComponent, system);
}

MoleculeType::MoleculeType(
		string name,
		vector <string> &compName,
		vector <string> &defaultCompState,
		vector < vector<string> > &possibleCompStates,
		vector <bool> isIntegerComponent,
		System *system)
 : population_type( false ), isFixed_(false), fixedCount_(0), fixedCompartment_(nullptr)
{
	init(name, compName, defaultCompState, possibleCompStates, isIntegerComponent, system);
}


MoleculeType::MoleculeType(
		string name,
		vector <string> &compName,
		vector <string> &defaultCompState,
		vector < vector<string> > &possibleCompStates,
		vector <bool> isIntegerComponent,
		bool pop_type,
		System *system)
 : population_type( pop_type ), isFixed_(false), fixedCount_(0), fixedCompartment_(nullptr)
{
	init(name, compName, defaultCompState, possibleCompStates, isIntegerComponent, system);
}


const string& MoleculeType::getComponentName(int cIndex) const {
	if(cIndex>=this->numOfComponents) {
		cerr<<"in MoleculeType:getComponentName, can't get that component!"<<endl;
		cerr<<"looking for cIndex: "<<cIndex<<endl;
		cerr<<"I only have: "<<numOfComponents<<" components."<<endl;
		exit(2);
	}
	return this->compName[cIndex];
}

void MoleculeType::init(
	string name,
	vector <string> &compName,
	vector <string> &defaultCompState,
	vector < vector<string> > &possibleCompStates,
	vector <bool> isIntegerComponent,
	System *system)
{
	//Basics...
	this->name=name;
	this->numOfComponents=compName.size();

	//First, some quick error checks
	if((int)defaultCompState.size()!=numOfComponents || (int)possibleCompStates.size()!=numOfComponents ||
			(int)isIntegerComponent.size()!=numOfComponents) {
		cout<<"Error creating MoleculeType: '"<<name<<"': The length of the input vectors\n";
		cout<<"do not match, so I can't initialize this object.\n";
		cout<<"quitting now."<<endl; exit(1);
	}

	//Now we can get on with initializing the MoleculeType information
	this->compName=new string [numOfComponents];
	this->defaultCompState = new int [numOfComponents];
	this->isIntegerCompState = new bool [numOfComponents];

	int nostate = Molecule::NOSTATE;
	for(int c=0; c<numOfComponents; c++) {
		this->compName[c]=compName.at(c);
		this->compNameMap[compName.at(c)] = c;
		this->isIntegerCompState[c]=isIntegerComponent.at(c);

		bool foundDefaultState=false;
		vector <string> p;
		for(unsigned int i=0; i<possibleCompStates.at(c).size(); i++) {
			p.push_back(possibleCompStates.at(c).at(i));
			if(possibleCompStates.at(c).at(i) == defaultCompState.at(c)) {
				this->defaultCompState[c]=i; foundDefaultState=true;
			}
		}
		if(!foundDefaultState) this->defaultCompState[c]=Molecule::NOSTATE;
		this->possibleCompStates.push_back(p);
	}


	reactionMappingIndices.clear();
	reactionMappingCount = 0;

	//Register myself with the system, and get an ID number
	this->system = system;
	this->type_id = this->system->addMoleculeType(this);
	compactPartnerPools.assign(numOfComponents,
			static_cast<CompactPartnerPool *>(0));
	compactEnergyCenterCandidateBits.assign(
			numOfComponents, vector<std::uint64_t>());
	compactEnergyContextCandidateBits.assign(
			numOfComponents, vector<std::uint64_t>());
	compactPartnerCandidateBits.assign(
			numOfComponents, vector<std::uint64_t>());
	compactPartnerReactionIndices.assign(
			numOfComponents, vector<unsigned int>());
	compactEnergyContextMinimumRequiredBits.assign(numOfComponents, 0);
	nonCompactMembershipCandidateBits.clear();
	hasCompactEnergyMembershipIndex = false;


	mList = new MoleculeList(this,2,system->getGlobalMoleculeLimit());
	n_eqComp = 0;
	eqCompOriginalName = nullptr;
	eqCompSizes = nullptr;
	eqCompName = nullptr;
	eqCompIndex = nullptr;
	indexToEqClass = nullptr;
}






MoleculeType::~MoleculeType()
{
	if(DEBUG) cout << "Destroying MoleculeType " << name << endl;

	//Delete freestore component information
	delete [] compName;
	delete [] defaultCompState;
	delete [] isIntegerCompState;

	//Delete details about equivalent components
	delete [] eqCompSizes;
	for(int i=0; i<n_eqComp; i++) {
		delete [] eqCompName[i];
		delete [] eqCompIndex[i];
	}
	delete [] eqCompName;
	delete [] eqCompIndex;
	delete [] eqCompOriginalName;
	if (indexToEqClass) {
		delete [] indexToEqClass;
	}
	for (vector<CompactPartnerPool *>::iterator it = compactPartnerPools.begin();
			it != compactPartnerPools.end(); ++it) {
		delete *it;
	}




	//Delete all template molecules of this type that exist
	TemplateMolecule *t;
	while(allTemplates.size()>0)
	{
		t = allTemplates.back();
		allTemplates.pop_back();
		delete t;
	}




	delete mList;
}

void MoleculeType::addEquivalentComponents(vector <vector <string> > &identicalComponents)
{
	this->n_eqComp = identicalComponents.size();
	eqCompOriginalName = new string [n_eqComp];
	eqCompName=new string * [n_eqComp];
	eqCompIndex=new int *[n_eqComp];
	eqCompSizes=new int [n_eqComp];

	if (indexToEqClass == nullptr) {
		indexToEqClass = new int[numOfComponents];
		for(int c=0; c<numOfComponents; c++) {
			indexToEqClass[c] = -1;
		}
	}

	for(int i=0; i<n_eqComp; i++) {
		eqCompSizes[i]=identicalComponents.at(i).size();
		eqCompName[i] = new string [eqCompSizes[i]];
		eqCompIndex[i] = new int [eqCompSizes[i]];
		for(int k=0; k<eqCompSizes[i]; k++) {
			if(k==0) {
				string tempString = identicalComponents.at(i).at(k);
				eqCompOriginalName[i] = tempString.substr(0,tempString.size()-1);
			}
			eqCompName[i][k] = identicalComponents.at(i).at(k);
			eqCompIndex[i][k] = getCompIndexFromName(eqCompName[i][k]);

			// Map the component index to the equivalency class index
			if(eqCompIndex[i][k] >= 0 && eqCompIndex[i][k] < numOfComponents) {
				indexToEqClass[eqCompIndex[i][k]] = i;
			}
		}
	}
}


bool MoleculeType::isIntegerComponent(const string& cName) const {
	auto it = compNameMap.find(cName);
	if (it != compNameMap.end()) {
		return this->isIntegerCompState[it->second];
	}
	cerr<<"!!! error !!! cannot find site name "<< cName << " in MoleculeType: "<<name;
	cerr<<"in function isIntegerComponent(string cName).  "<<endl;
	this->printDetails();
	exit(1);
}
bool MoleculeType::isIntegerComponent(int cIndex) const {
	if(cIndex>=0 && cIndex<numOfComponents) {
		return this->isIntegerCompState[cIndex];
	} else {
		cerr<<"!!! error !!! "<< cIndex << " is not a valid component index in MoleculeType: "<<name;
		cerr<<"in function isIntegerComponent(int cIndex).  "<<endl;
		this->printDetails();
		exit(1);
	}
}


bool MoleculeType::isEquivalentComponent(const string& cName) const {
	for(int i=0; i<n_eqComp; i++) {
		if(eqCompOriginalName[i].compare(cName)==0) {
			return true;
		}
	}
	return false;
}
bool MoleculeType::isEquivalentComponent(int cIndex) const {
	if (indexToEqClass && cIndex >= 0 && cIndex < numOfComponents) {
		return indexToEqClass[cIndex] != -1;
	}
	return false;
}

int MoleculeType::getEquivalenceClassNumber(int cIndex) const {
	if (indexToEqClass && cIndex >= 0 && cIndex < numOfComponents) {
		return indexToEqClass[cIndex];
	}
	return -1;
}

string MoleculeType::getEquivalenceClassComponentNameFromComponentIndex(int cIndex) const {
	if (indexToEqClass && cIndex >= 0 && cIndex < numOfComponents) {
		int eqClassIndex = indexToEqClass[cIndex];
		if (eqClassIndex != -1) {
			return eqCompOriginalName[eqClassIndex];
		}
	}
	cerr<<"Could not find equivalency class component string for component number: "<<cIndex<<"!!!"<<endl;
	exit(1);
}

void MoleculeType::getEquivalencyClass(int *&components, int &n_components, const string& cName) const {
	for(int i=0; i<n_eqComp; i++) {
		if(eqCompOriginalName[i].compare(cName)==0) {
			components = eqCompIndex[i];
			n_components=eqCompSizes[i];
			return;
		}
	}
}
int MoleculeType::getEquivalencyClassNumber(const string& cName) const {
	for(int i=0; i<n_eqComp; i++) {
		if(eqCompOriginalName[i].compare(cName)==0) {
				return i;
		}
	}
	cerr<<"Could not find equivalency class number for component named: "<<cName<<"!!!"<<endl;
	exit(1);
}




string MoleculeType::getComponentStateName(int cIndex, int cValue) {
	if(cValue==Molecule::NOSTATE) return "NO_STATE";
	if( cIndex>=(int)possibleCompStates.size() || cIndex<0 ){
		cerr<<"Component index out of range (moltype="<<getName()<<" cIndex="<<cIndex<<")!!!"<<endl;
		exit(1);
	}
	if( cValue>=(int)possibleCompStates.at(cIndex).size() || cValue<0 ){
		cerr<<"State index out of range (moltype="<<getName()<<" cIndex="<<cIndex<<" cValue="<<cValue<<")!!!"<<endl;
		exit(1);
	}
	return possibleCompStates.at(cIndex).at(cValue);
}




Molecule *MoleculeType::genDefaultMolecule(Compartment *c)
{
	Molecule *m;
	mList->create(m);
	m->setAlive(true);
	m->setCompartment(c);
	//cout<<"adding molecule: "<<m->getMoleculeTypeName()<<"_"<<m->getUniqueID()<<endl;

	return m;
}


void MoleculeType::addMoleculeToRunningSystem(Molecule *&mol)
{
	//cout<<"adding molecule: "<<mol->getMoleculeTypeName()<<"_"<<mol->getUniqueID()<<endl;
	//First prepare the molecule for simulation
	mol->setUpLocalFunctionList();
	mol->prepareForSimulation();
	mol->setAlive(true);

	mol->addToObservables();
	this->updateRxnMembership(mol);
}


void MoleculeType::addMoleculeToRunningSystemButDontUpdate(Molecule *&mol)
{
	//First prepare the molecule for simulation

	mol->setUpLocalFunctionList();
	mol->prepareForSimulation();
	mol->setAlive(true);

	//We assume observables and reaction membership will be updated later
	// (this is now the case for reaction firing)
}


void MoleculeType::removeMoleculeFromRunningSystem(Molecule *&m)
{
	if (system->isUsingComplex())
		// Need to manually unset canonical flag since we're not calling a Complex method
		m->getComplex()->unsetCanonical();

	//Remove this guy from the list, the observables list, and from all rxns
	mList->remove(m->getMolListId(), m);
	removeFromObservables(m);
	removeFromRxns(m);


	//We also have to remove all bonds
	for(int c=0; c<getNumOfComponents(); c++) {
		if(m->isBindingSiteBonded(c)) {
			Molecule::unbind(m,c);
		}
	}

	m->setAlive(false);

}


void MoleculeType::removeAllMolecules()
{
	// Iterate through all molecules and remove them. Loop backwards because
	// remove() removes by swapping with the last element.
	//
	// Do NOT delete the Molecule objects here: mList is a fixed-capacity pool
	// that owns every Molecule across [0, capacity) and recycles them on
	// genDefaultMolecule()/create(). remove() only unbinds, drops the molecule
	// from observables/reactions, marks it dead, and decrements the live count
	// — the object stays in the pool for reuse. Deleting it leaves a dangling
	// pointer in the pool (use-after-free on the next create()) and a double
	// free in ~MoleculeList(). This is the crash behind resetConcentrations().
	for (int m = mList->size() - 1; m >= 0; m--) {
		Molecule *mol = mList->at(m);
		removeMoleculeFromRunningSystem(mol);
	}
}

void MoleculeType::removeMoleculeFromRunningSystemButDontUpdate(Molecule *&m)
{
	//Remove this guy from the list, the observables list, and from all rxns
	mList->remove(m->getMolListId(), m);
	//removeFromObservables(m);
	//removeFromRxns(m);

	//We also have to remove all bonds
	for(int c=0; c<getNumOfComponents(); c++) {
		if(m->isBindingSiteBonded(c)) {
			Molecule::unbind(m,c);
		}
	}

	m->setAlive(false);
}


Molecule * MoleculeType::getMolecule(int ID_molecule) const {
	return mList->at(ID_molecule);
}
int MoleculeType::getMoleculeCount() const {
	return mList->size();
}

CompactPartnerPool *MoleculeType::getOrCreateCompactPartnerPool(
		int componentIndex)
{
	if (componentIndex < 0 || componentIndex >= numOfComponents)
		return 0;
	if (compactPartnerPools[componentIndex] == 0)
		compactPartnerPools[componentIndex] = new CompactPartnerPool();
	return compactPartnerPools[componentIndex];
}


void MoleculeType::addTemplateMolecule(TemplateMolecule *t)
{
	if(t->getMoleculeType()==this)
		allTemplates.push_back(t);
	else
		cout<<"!!!!Error: trying to add molecule of type " << t->getMoleculeTypeName() << " to MoleculeType " << name << endl;
}

string MoleculeType::getMolObsName(int obsIndex) const {
	return molObs.at(obsIndex)->getName();
}

int MoleculeType::getMolObsCount(int obsIndex) const {
	return molObs.at(obsIndex)->getCount();
}



int MoleculeType::getCompIndexFromName(const string& cName) const
{
	auto it = compNameMap.find(cName);
	if(it != compNameMap.end()) return it->second;

	cerr<<"!!! warning !!! cannot find site name "<< cName << " in MoleculeType: "<<name<<endl;
	this->printDetails();
	throw std::runtime_error("Cannot find site name " + cName + " in MoleculeType: " + name);
}

int MoleculeType::getStateValueFromName(int cIndex, string stateName) const
{
	for(unsigned int s=0; s<possibleCompStates.at(cIndex).size(); s++) {
		if(possibleCompStates.at(cIndex).at(s)==stateName) {
			return s;
		}
	}
	cerr<<"Error!  '"<<stateName<<" is not a recognized possible state for '"<<compName[cIndex]<<"' in MoleculeType: '"<<name<<"'"<<endl;
	cerr<<"For that, I'm quitting!";
	printDetails();
	exit(1);
}




void MoleculeType::addReactionClass(ReactionClass * r, int rPosition)
{
	/* The decision vectors follow the reaction list, so invalidate them if a
	 * reaction is registered after the cache has been populated. */
	directMembershipDecisionCache.clear();
	directMembershipDecisionCacheSafe.clear();
	unsigned int reactionIndex = static_cast<unsigned int>(reactions.size());
	this->reactions.push_back(r);
	this->reactionPositions.push_back(rPosition);
	int partnerComponent = -1;
	bool compactPartnerRegistration =
		r->getCompactPartnerPoolInfo(rPosition, partnerComponent) &&
		partnerComponent >= 0 && partnerComponent < numOfComponents;
	reactionMappingIndices.push_back(compactPartnerRegistration
			? -1 : reactionMappingCount++);

	int reactionCenterComponent = -1;
	std::uint64_t contextComponentMask = 0;
	unsigned int minimumContextComponents = 0;
	bool indexed = r->getCompactMembershipIndexInfo(
			rPosition, reactionCenterComponent, contextComponentMask,
			minimumContextComponents) &&
			reactionCenterComponent >= 0 &&
			reactionCenterComponent < numOfComponents;
	if (indexed) {
		hasCompactEnergyMembershipIndex = true;
		setCompactMembershipCandidateBit(
				compactEnergyCenterCandidateBits[reactionCenterComponent],
				reactionIndex);
		for (int componentIndex = 0;
				componentIndex < numOfComponents && componentIndex < 64;
				++componentIndex) {
			if ((contextComponentMask &
					(std::uint64_t(1) << componentIndex)) == 0)
				continue;
			setCompactMembershipCandidateBit(
					compactEnergyContextCandidateBits[componentIndex],
					reactionIndex);
			unsigned int &minimum =
					compactEnergyContextMinimumRequiredBits[componentIndex];
			if (minimum == 0 || minimumContextComponents < minimum)
				minimum = minimumContextComponents;
		}
	} else {
		setCompactMembershipCandidateBit(
				nonCompactMembershipCandidateBits, reactionIndex);
	}
	if (compactPartnerRegistration) {
		CompactPartnerPool *partnerPool = r->getCompactPartnerPool();
		if (partnerPool != 0)
			partnerPool->registerReaction(
					r, r->supportsCompactPartnerPoolUpdate());
		setCompactMembershipCandidateBit(
				compactPartnerCandidateBits[partnerComponent], reactionIndex);
		compactPartnerReactionIndices[partnerComponent].push_back(reactionIndex);
	}

	//We also have to check to make sure that if the reaction is a DOR reaction,
	//we remember it so we can updated it
	if(r->getRxnType()==ReactionClass::DOR_RXN) {
		if( r->getDORreactantPosition()==rPosition) {
			indexOfDORrxns.push_back(reactions.size()-1);
		}
	}
	else if(r->getRxnType()==ReactionClass::DOR2_RXN) {
		if( r->getDORreactantPosition()==rPosition) {
			indexOfDORrxns.push_back(reactions.size()-1);
		}
		else if( r->getDORreactantPosition2()==rPosition) {
			indexOfDORrxns.push_back(reactions.size()-1);
		}
	}
}

bool MoleculeType::canSkipIndirectMembership(
		ReactionClass *firedReaction) const
{
	for (vector<ReactionClass *>::const_iterator it = reactions.begin();
			it != reactions.end(); ++it) {
		if (!(*it)->canSkipIndirectMembership(firedReaction))
			return false;
	}
	return true;
}



void MoleculeType::populateWithDefaultMolecules(int moleculeCount)
{
	if(DEBUG) cout<< " Populating "<< this->name << " with " << moleculeCount << " molecule(s)";
	if(DEBUG) cout<< " for a total of " << mList->size()+moleculeCount << " molecule(s)."<<endl;
	//mInstances.reserve(mInstances.size()+moleculeCount);
	for(int m=0; m<moleculeCount; m++)
	{
		if(DEBUG) cout<<" ("<<m+1<<") ";

		//Create the molecule (which knows how many components to make)
		this->genDefaultMolecule();
		//new Molecule(this);

		//Add the molecule to the list of molecules so we save it (does this automatically now!!!! )
		//mInstances.push_back(mol);
	}
}



void MoleculeType::setUpLocalFunctionListForMolecules()
{
	Molecule *mol;
	for(int m=0; m<mList->size(); m++ )
	{
	  	mol = mList->at(m);
	  	mol->setUpLocalFunctionList();
	}
}

void MoleculeType::prepareForSimulation()
{
	//cout<<"Preparing: "<<name<<endl;
	//Check each reaction and add this molecule as a reactant if we have to
	int r=0;
	for(rxnIter = reactions.begin(), r=0; rxnIter != reactions.end(); rxnIter++, r++ )
	{
		system->registerRxnIndex((*rxnIter)->getRxnId(), reactionPositions.at(r),r);
  	}


	//Our iterators that we will use to loop through every molecule
	Molecule *mol;
  	for( int m=0; m<mList->size(); m++ )
  	{
  		//First prepare the molecule for simulation
  		mol = mList->at(m);
  		mol->prepareForSimulation();

  		//Check each observable and see if this molecule should be counted
  		this->addToObservables(mol);

  		//Check each reaction and add this molecule as a reactant if we have to
		for(rxnIter = reactions.begin(), r=0; rxnIter != reactions.end(); rxnIter++, r++ )
		{
			if ((*rxnIter)->usesIncrementalMembership())
				(*rxnIter)->tryToAddWithIndex(
						mol, reactionPositions.at(r), r);
			else
				(*rxnIter)->tryToAdd(mol, reactionPositions.at(r));
  		}
	}
}

void MoleculeType::updateRxnMembership(Molecule * m,
		ReactionClass * firedReaction, bool directProduct)
{
	const vector<unsigned char> *cachedDecisions = 0;
	IncrementalMembershipChange membershipChange;
	bool hasMembershipChange = directProduct && firedReaction != 0 &&
		firedReaction->usesIncrementalMembership() &&
		firedReaction->getIncrementalMembershipChange(membershipChange);
	bool refineMembershipChange = hasMembershipChange &&
		m->getMoleculeType() == membershipChange.moleculeType1;
	if (refineMembershipChange && membershipChange.componentIndex1 >= 0 &&
			membershipChange.componentIndex1 < 64) {
		/* When the weighted molecule is full immediately before or after the
		 * event, every accepted context dependency on the changed site crosses
		 * its predicate.  The endpoint cache is therefore already exact and the
		 * per-reaction mask test would only add overhead. */
		std::uint64_t changedBit = std::uint64_t(1) <<
				membershipChange.componentIndex1;
		int componentCount = m->getMoleculeType()->getNumOfComponents();
		if (componentCount <= 64) {
			std::uint64_t fullMask = componentCount == 64
					? ~std::uint64_t(0)
					: ((std::uint64_t(1) << componentCount) - 1);
			std::uint64_t newMask = m->getBoundComponentMask();
			std::uint64_t oldMask = membershipChange.isBoundAfter1
					? (newMask & ~changedBit) : (newMask | changedBit);
			if (newMask == fullMask || oldMask == fullMask)
				refineMembershipChange = false;
		}
	}
	bool useCompactMembershipIndex = hasMembershipChange &&
			m->getMoleculeType() == membershipChange.moleculeType1 &&
			membershipChange.componentIndex1 >= 0 &&
			membershipChange.componentIndex1 < 64 &&
			static_cast<unsigned int>(membershipChange.componentIndex1) <
				compactEnergyCenterCandidateBits.size() &&
			static_cast<unsigned int>(membershipChange.componentIndex1) <
				compactEnergyContextCandidateBits.size() &&
			hasCompactEnergyMembershipIndex;
	bool useCompactPartnerPoolIndex = hasMembershipChange && directProduct &&
			membershipChange.moleculeType1 != membershipChange.moleculeType2 &&
			m->getMoleculeType() == membershipChange.moleculeType2 &&
			membershipChange.componentIndex2 >= 0 &&
			static_cast<unsigned int>(membershipChange.componentIndex2) <
				compactPartnerCandidateBits.size() &&
			static_cast<unsigned int>(membershipChange.componentIndex2) <
				compactPartnerReactionIndices.size() &&
			!compactPartnerReactionIndices[membershipChange.componentIndex2].empty();
	/* The compact weighted-side candidate index contains exactly the endpoint
	 * dependencies for simple EnergyPattern rules.  When every registered
	 * reaction on this type uses that index, the generic per-fire decision cache
	 * would only repeat the same dependency test.  Keep the fallback for mixed
	 * molecule types and same-type binding, where the partner endpoint can also
	 * be represented by this molecule type. */
	bool compactMembershipDecisionsComplete = useCompactMembershipIndex &&
		membershipChange.moleculeType1 != membershipChange.moleculeType2 &&
		nonCompactMembershipCandidateBits.empty();
	bool compactPartnerPoolChanged = false;
	bool compactPartnerPoolBatchScheduled = false;
	CompactPartnerPool *compactPartnerPool = 0;
	int oldCompactPartnerPoolSize = 0;
	const vector<std::uint64_t> *partnerCandidates = 0;
	if (useCompactPartnerPoolIndex) {
		int partnerComponent = membershipChange.componentIndex2;
		partnerCandidates = &compactPartnerCandidateBits[partnerComponent];
		const vector<unsigned int> &partnerReactions =
				compactPartnerReactionIndices[partnerComponent];
		unsigned int firstPartnerReaction = partnerReactions.front();
		compactPartnerPool = reactions[firstPartnerReaction]->
				getCompactPartnerPool();
		if (compactPartnerPool != 0)
			oldCompactPartnerPoolSize = compactPartnerPool->size();
		compactPartnerPoolChanged = reactions[firstPartnerReaction]->
				refreshCompactPartnerPool(
						m, reactionPositions[firstPartnerReaction]);
		if (compactPartnerPoolChanged && compactPartnerPool != 0 &&
					this->system->isDeferringMembershipPropensityUpdates() &&
					compactPartnerPool->supportsBatchUpdate()) {
			this->system->deferCompactPartnerPoolUpdate(
					compactPartnerPool, oldCompactPartnerPoolSize,
					compactPartnerPool->size());
			compactPartnerPoolBatchScheduled = true;
		}
	}
	bool allReactionsUseCompactPartnerPool =
			useCompactPartnerPoolIndex &&
			compactPartnerReactionIndices[membershipChange.componentIndex2].size() ==
				reactions.size();
	if (allReactionsUseCompactPartnerPool) {
		if (compactPartnerPoolChanged) {
			bool defer = this->system->isDeferringMembershipPropensityUpdates();
			if (defer && compactPartnerPoolBatchScheduled)
				return;
			const vector<unsigned int> &partnerReactions =
					compactPartnerReactionIndices[membershipChange.componentIndex2];
			for (vector<unsigned int>::const_iterator it =
					partnerReactions.begin(); it != partnerReactions.end(); ++it) {
				ReactionClass *rxn = reactions[*it];
				if (defer) {
					this->system->deferMembershipPropensityUpdate(rxn);
				} else {
					double oldA = rxn->get_a();
					double newA = rxn->update_a();
					this->system->update_A_tot(rxn, oldA, newA);
				}
			}
		}
		return;
	}
	if (!compactMembershipDecisionsComplete && directProduct && firedReaction != 0 &&
		firedReaction->usesIncrementalMembership()) {
		unordered_map<ReactionClass *, bool>::iterator safe =
				directMembershipDecisionCacheSafe.find(firedReaction);
		if (safe == directMembershipDecisionCacheSafe.end()) {
			bool typeInvariant = true;
			for (unsigned int r = 0; r < reactions.size(); ++r) {
				if (!reactions[r]->membershipDecisionIsTypeInvariant()) {
					typeInvariant = false;
					break;
				}
			}
			safe = directMembershipDecisionCacheSafe.emplace(
					firedReaction, typeInvariant).first;
		}
		if (safe->second) {
			unordered_map<ReactionClass *,
					DirectMembershipDecisionCacheEntry>::iterator cached =
					directMembershipDecisionCache.find(firedReaction);
			if (cached == directMembershipDecisionCache.end()) {
				cached = directMembershipDecisionCache.emplace(
						firedReaction, DirectMembershipDecisionCacheEntry()).first;
				DirectMembershipDecisionCacheEntry &entry = cached->second;
				entry.decisions.reserve(reactions.size());
				for (unsigned int r = 0; r < reactions.size(); ++r) {
					entry.decisions.push_back(
							reactions[r]->shouldUpdateMembership(
									m, firedReaction, true));
				}
				std::size_t affectedReactionCount = 0;
				for (unsigned int r = 0; r < entry.decisions.size(); ++r) {
					if (entry.decisions[r]) ++affectedReactionCount;
				}
				if (affectedReactionCount * 2 < entry.decisions.size()) {
					entry.reactionIndices.reserve(affectedReactionCount);
					for (unsigned int r = 0; r < entry.decisions.size(); ++r) {
						if (entry.decisions[r])
							entry.reactionIndices.push_back(r);
					}
					vector<unsigned char>().swap(entry.decisions);
					entry.useReactionIndices = true;
				}
			}
			const DirectMembershipDecisionCacheEntry &entry = cached->second;
			if (entry.useReactionIndices && !useCompactMembershipIndex) {
				for (vector<unsigned int>::const_iterator it =
						entry.reactionIndices.begin();
						it != entry.reactionIndices.end(); ++it) {
					unsigned int r = *it;
					ReactionClass *rxn = reactions.at(r);
					if (useCompactPartnerPoolIndex &&
							compactMembershipBitIsSet(*partnerCandidates, r)) {
						if (!compactPartnerPoolChanged) continue;
						bool defer = this->system->isDeferringMembershipPropensityUpdates();
						if (defer) {
							if (compactPartnerPoolBatchScheduled)
								continue;
							this->system->deferMembershipPropensityUpdate(rxn);
						} else {
							double oldA = rxn->get_a();
							double newA = rxn->update_a();
							this->system->update_A_tot(rxn, oldA, newA);
						}
						continue;
					}
					if (refineMembershipChange &&
							!rxn->shouldUpdateMembershipForChange(
									m, membershipChange))
						continue;
					bool useIndexedMembership =
						rxn->supportsDeferredMembershipUpdate();
					bool defer = this->system->isDeferringMembershipPropensityUpdates() &&
						useIndexedMembership;
					if (defer) {
						bool changed = useIndexedMembership
							? rxn->tryToAddAndReportChangeWithIndex(
									m, reactionPositions.at(r), r)
							: rxn->tryToAddAndReportChange(
									m, reactionPositions.at(r));
						if (changed)
							this->system->deferMembershipPropensityUpdate(rxn);
					} else {
						double oldA = rxn->get_a();
						if (useIndexedMembership)
							rxn->tryToAddWithIndex(
									m, reactionPositions.at(r), r);
						else
							rxn->tryToAdd(m, reactionPositions.at(r));
						double newA = rxn->update_a();
						this->system->update_A_tot(rxn, oldA, newA);
					}
				}
				return;
			}
			if (!entry.useReactionIndices)
				cachedDecisions = &entry.decisions;
		}
	}

	if (useCompactMembershipIndex) {
		int changedComponent = membershipChange.componentIndex1;
		std::uint64_t changedBit = std::uint64_t(1) << changedComponent;
		std::uint64_t newMask = m->getBoundComponentMask();
		std::uint64_t oldMask = membershipChange.isBoundAfter1
				? (newMask & ~changedBit) : (newMask | changedBit);
		unsigned int minimumContextBits =
				compactEnergyContextMinimumRequiredBits[changedComponent];
		bool includeContext = minimumContextBits == 0 ||
			compactMembershipBitCount(newMask | oldMask) >= minimumContextBits;
		const vector<std::uint64_t> &centerCandidates =
				compactEnergyCenterCandidateBits[changedComponent];
		const vector<std::uint64_t> &contextCandidates =
				compactEnergyContextCandidateBits[changedComponent];
		const std::size_t wordCount = std::max(
				nonCompactMembershipCandidateBits.size(),
				std::max(centerCandidates.size(),
						includeContext ? contextCandidates.size() : std::size_t(0)));

		for (std::size_t wordIndex = 0; wordIndex < wordCount; ++wordIndex) {
			std::uint64_t candidates = wordIndex <
					nonCompactMembershipCandidateBits.size()
				? nonCompactMembershipCandidateBits[wordIndex] : 0;
			std::uint64_t contextCandidateBits = 0;
			if (wordIndex < centerCandidates.size())
				candidates |= centerCandidates[wordIndex];
			if (includeContext && wordIndex < contextCandidates.size())
				contextCandidateBits = contextCandidates[wordIndex];
			candidates |= contextCandidateBits;
			while (candidates != 0) {
				unsigned int bit = compactMembershipTrailingZeroCount(candidates);
				unsigned int r = static_cast<unsigned int>((wordIndex << 6) + bit);
				candidates &= candidates - 1;
				ReactionClass *rxn = reactions[r];
				std::uint64_t reactionBit = std::uint64_t(1) << bit;
				bool contextCandidate =
						(contextCandidateBits & reactionBit) != 0;
				bool centerCandidate = wordIndex < centerCandidates.size() &&
						(centerCandidates[wordIndex] & reactionBit) != 0;
				/* A context-only change cannot create a weighted-side mapping:
				 * the reaction center occupancy is unchanged.  Avoid probing an
				 * inactive compact rule and retain the normal path for center and
				 * non-compact candidates. */
				if (contextCandidate && !centerCandidate &&
						m->getRxnListMappingId(r) < 0)
					continue;
				if (useCompactPartnerPoolIndex &&
						compactMembershipBitIsSet(*partnerCandidates, r)) {
					if (!compactPartnerPoolChanged) continue;
					bool defer = this->system->isDeferringMembershipPropensityUpdates();
					if (defer) {
						if (compactPartnerPoolBatchScheduled)
							continue;
						this->system->deferMembershipPropensityUpdate(rxn);
					} else {
						double oldA = rxn->get_a();
						double newA = rxn->update_a();
						this->system->update_A_tot(rxn, oldA, newA);
					}
					continue;
				}
				if (!compactMembershipDecisionsComplete) {
					if (cachedDecisions != 0) {
						if (!(*cachedDecisions)[r]) continue;
					} else if (!rxn->shouldUpdateMembership(
							m, firedReaction, directProduct))
						continue;
				}
				if (refineMembershipChange &&
						!rxn->shouldUpdateMembershipForChange(
								m, membershipChange))
					continue;
				bool useIndexedMembership =
					rxn->supportsDeferredMembershipUpdate();
				bool defer = this->system->isDeferringMembershipPropensityUpdates() &&
						useIndexedMembership;
				if (defer) {
					bool changed = useIndexedMembership
						? rxn->tryToAddAndReportChangeWithIndex(
								m, reactionPositions[r], r)
						: rxn->tryToAddAndReportChange(
								m, reactionPositions[r]);
					if (changed)
						this->system->deferMembershipPropensityUpdate(rxn);
				} else {
					double oldA = rxn->get_a();
					if (useIndexedMembership)
						rxn->tryToAddWithIndex(m, reactionPositions[r], r);
					else
						rxn->tryToAdd(m, reactionPositions[r]);
					double newA = rxn->update_a();
					this->system->update_A_tot(rxn, oldA, newA);
				}
			}
		}
		return;
	}

	for( unsigned int r=0; r<reactions.size(); r++ )
	{
		ReactionClass * rxn=reactions.at(r);
		if (useCompactPartnerPoolIndex &&
				compactMembershipBitIsSet(*partnerCandidates, r)) {
			if (!compactPartnerPoolChanged) continue;
			bool defer = this->system->isDeferringMembershipPropensityUpdates();
			if (defer) {
				if (compactPartnerPoolBatchScheduled)
					continue;
				this->system->deferMembershipPropensityUpdate(rxn);
			} else {
				double oldA = rxn->get_a();
				double newA = rxn->update_a();
				this->system->update_A_tot(rxn, oldA, newA);
			}
			continue;
		}
		if (!compactMembershipDecisionsComplete) {
			if (cachedDecisions != 0) {
				if (!(*cachedDecisions)[r]) continue;
			}
			else if (!rxn->shouldUpdateMembership(m, firedReaction, directProduct))
				continue;
		}
		if (refineMembershipChange &&
				!rxn->shouldUpdateMembershipForChange(m, membershipChange))
			continue;
		bool useIndexedMembership =
			rxn->supportsDeferredMembershipUpdate();
		bool defer = this->system->isDeferringMembershipPropensityUpdates() &&
			useIndexedMembership;
		if (defer) {
			bool changed = useIndexedMembership
				? rxn->tryToAddAndReportChangeWithIndex(
						m, reactionPositions.at(r), r)
				: rxn->tryToAddAndReportChange(
						m, reactionPositions.at(r));
			if (changed)
				this->system->deferMembershipPropensityUpdate(rxn);
		} else {
			double oldA = rxn->get_a();
			if (useIndexedMembership)
				rxn->tryToAddWithIndex(
						m, reactionPositions.at(r), r);
			else
				rxn->tryToAdd(m, reactionPositions.at(r));
			double newA = rxn->update_a();
			this->system->update_A_tot(rxn,oldA,newA);
		}
  	}

}

void MoleculeType::updateConnectedRxnMembership(Molecule * m,
		ReactionClass * firedReaction, bool directProduct)
{
	// Preserve the MoleculeType's native reaction order so the connectivity path
	// mutates reactant containers in the same sequence as a full membership
	// refresh, while still using the precomputed connectivity matrix.
	for (unsigned int r=0; r<reactions.size(); r++) {
		rxn = reactions.at(r);
		if (!this->system->areReactionsConnected(
				firedReaction->getRxnId(), rxn->getRxnId())) {
			continue;
		}
		if (!rxn->shouldUpdateMembership(m, firedReaction, directProduct))
			continue;
		int pos = reactionPositions.at(r);
		double oldA = rxn->get_a();
		// Only -trackconnected reports the propensity just before the update;
		// computing it costs a second update_a() per rule and molecule.
		bool trackConnected = this->system->getTrackConnected();
		double oldAwithTotal = trackConnected ? rxn->update_a() : oldA;
		if (rxn->usesIncrementalMembership())
			rxn->tryToAddWithIndex(m, pos, r);
		else
			rxn->tryToAdd(m, pos);
		double newA = rxn->update_a();
		this->system->update_A_tot(rxn,oldA,newA);
		// Used for debugging to see which reaction rates changed
		// upon updating molecule membership
		// Arvind Rasi Subramaniam Nov 21, 2018
		if (!trackConnected) continue;
		if (oldAwithTotal != newA) {
			this->system->getConnectedRxnFileStream() <<
			this->system->getGlobalEventCounter() << "\t" <<
			firedReaction->getName() << "\t" <<
					m->getMoleculeTypeName() << "\t" <<
					m->getUniqueID() << "\t" <<
					rxn->getName() << "\t" <<
					oldAwithTotal << "\t" << newA << endl;
		}
  	}
}


int MoleculeType::getRxnIndex(ReactionClass * rxn, int rxnPosition)
{
	return system->getRxnIndex(rxn->getRxnId(),rxnPosition);

	//The old way!!  (that is slow if we have many rxns of course!)
	int r=0;
	for(rxnIter = reactions.begin(); rxnIter != reactions.end(); rxnIter++, r++ )
	{
		if((*rxnIter)==rxn)
			if(reactionPositions.at(r) == rxnPosition)
				return r;
	}
	cerr<<"Could not find this rxn: " << rxn->getName() << " in molecule Type: "<<name<<endl;
	exit(1);
}




void MoleculeType::removeFromObservables(Molecule *m)
{
	//cout<<"removing from observables:"<<m->getMoleculeTypeName()<<"_"<<m->getUniqueID()<<endl;
	//m->printDetails();

	//Check each observable and see if this molecule was counted, and if so, remove
	int ind=0;
  	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++ )
  	{
  		//Only subtract if m happened to be an observable... this saves us a compare call
  		//int matches = (*molObsIter)->isObservable(m);

  		// How many times does this observable match the molecule?
  		int matches = m->isObs(ind);
  		// subtract matches from observable
  		(*molObsIter)->subtract(matches);
		// set IsObs to zero, so if remove is called twice, we don't get negative counts.
  		m->setIsObs(ind,0);

  		ind++;
	}
}

void MoleculeType::removeFromRxns(Molecule * m)
{
	int r=0;
	for(rxnIter = reactions.begin(); rxnIter != reactions.end(); rxnIter++, r++ )
	{
		double oldA = (*rxnIter)->get_a();
		(*rxnIter)->remove(m, reactionPositions.at(r));
		double newA = (*rxnIter)->update_a();
		this->system->update_A_tot((*rxnIter),oldA,newA);
  	}
}




//TypeI local function: this molecule type depends on the value of this
//evaluated function
int MoleculeType::addLocalFunc_TypeI(LocalFunction *lf) {
	locFuncs_typeI.push_back(lf);
	return locFuncs_typeI.size()-1;

}

//TypeII local function: this molecule type, when updated, changes the
//value of this function
int MoleculeType::addLocalFunc_TypeII(LocalFunction *lf) {
	locFuncs_typeII.push_back(lf);
	return locFuncs_typeII.size()-1;

}
















void MoleculeType::addAllToObservables()
{
	/////  WARNING:: when calling this function, be sure to clear all observables
	/////  first, because this function will not clear observables.

//	cout<<"+++++++++ "<<this->getName()<<endl;

	//Check each observable and see if this molecule should be counted
	Molecule *mol;  int o=0;  int matches=0;
  	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++)
  	{
  		//cout<<"comparing to obs: "<<(*molObsIter)->getName()<<endl;

  		for( int m=0; m<mList->size(); m++ )
  		{
  			mol = mList->at(m);
  			matches = (*molObsIter)->isObservable(mol);
  			(*molObsIter)->add(matches);
  			mol->setIsObs(o,matches);
  			//cout<<"matches:"<<matches<<endl;
  		}
  		o++;
	}

}




void MoleculeType::addToObservables(Molecule *m)
{
	//Check each observable and see if this molecule should be counted
	int o=0;
  	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++)
  	{
		//cout<<"Comparing(in add: ";
		//cout<<m->getUniqueID()<<")"<<endl;

		int matches = (*molObsIter)->isObservable(m);
		m->setIsObs(o,matches);

		(*molObsIter)->add(matches);
		o++;
	}
	
}


void MoleculeType::outputMolObsNames(NFstream &fout)
{
	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++ )
		fout<<"\t"<<(*molObsIter)->getName();
}

void MoleculeType::outputMolObsCounts(NFstream &fout)
{
	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++ )
		fout<<"\t"<<(*molObsIter)->getCount();
}

void MoleculeType::printMolObsNames()
{
	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++)
		cout<<"\t"<<(*molObsIter)->getName();
}

void MoleculeType::printMolObsCounts()
{
	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++ )
		cout<<"\t"<<(*molObsIter)->getCount();
}


void MoleculeType::printAllMolecules()
{
	for( int m=0; m<mList->size(); m++ ) {
		mList->at(m)->printDetails();
	}

}


void MoleculeType::printDetails() const
{
	cout<<"Molecule Type: "<< name << " type ID: " << type_id <<endl;

	cout<<"   -components ( ";
	for(int c=0; c<numOfComponents; c++) {

		cout<<compName[c];
		if(!isIntegerCompState[c]) {
			for(unsigned int s=0; s<possibleCompStates.at(c).size(); s++) {
				cout<<"~"<<possibleCompStates.at(c).at(s);
			}
		} else {
			cout<<"~integer[0-"<<possibleCompStates.at(c).at(possibleCompStates.at(c).size()-1)<<"]";
		}
		if(c<(numOfComponents-1)) cout<<", ";
	}
	cout<<" )"<<endl;

	//Output the local functions...
	cout<<"  Type I local functions include:";
	if(locFuncs_typeI.size()==0) cout<<"  none.";
	for(unsigned int ti=0; ti<locFuncs_typeI.size(); ti++) {
		cout<<"  "<<locFuncs_typeI.at(ti)->getNiceName();
	} cout<<endl;
	cout<<"  Type II local functions include:";
	if(locFuncs_typeII.size()==0) cout<<"  none.";
	for(unsigned int tii=0; tii<locFuncs_typeII.size(); tii++) {
		cout<<"  "<<locFuncs_typeII.at(tii)->getNiceName();
	} cout<<endl;


	cout<<"   -has "<< mList->size() <<" molecules."<<endl;
	cout<<"   -has "<< reactions.size() <<" reactions"<<endl;
//	cout<<"        of which "<< indexOfDORrxns.size() <<" are DOR rxns. "<<endl;
	cout<<"   -has "<< molObs.size() <<" molecules observables " <<endl;
}
