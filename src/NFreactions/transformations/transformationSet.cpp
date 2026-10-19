#include "transformationSet.hh"
#include "transformation.hh"
#include "../../NFcore/compartment.hh"
#include <queue>
#include <unordered_set>

using namespace NFcore;

namespace {
Compartment *selectReactionCompartment(MappingSet **mappingSets,
		unsigned int nReactants, bool &compatible)
{
	std::vector<Compartment *> volumes;
	std::vector<Compartment *> surfaces;
	bool missingCompartment = false;
	compatible = true;
	for (unsigned int r = 0; r < nReactants; ++r) {
		MappingSet *set = mappingSets ? mappingSets[r] : NULL;
		Mapping *mapping = set && set->getNumOfMappings() > 0 ? set->get(0) : NULL;
		Molecule *molecule = mapping ? mapping->getMolecule() : NULL;
		if (!molecule) { compatible = false; return NULL; }
		Compartment *compartment = molecule->getCompartment();
		if (!compartment) missingCompartment = true;
		else if (compartment->getSpatialDimensions() == 3) volumes.push_back(compartment);
		else if (compartment->getSpatialDimensions() == 2) surfaces.push_back(compartment);
		else { compatible = false; return NULL; }
	}
	if (missingCompartment && (!volumes.empty() || !surfaces.empty())) compatible = false;
	for (Compartment *compartment : volumes)
		if (compartment != volumes.front()) compatible = false;
	for (Compartment *compartment : surfaces)
		if (compartment != surfaces.front()) compatible = false;
	if (!volumes.empty() && !surfaces.empty() &&
			!volumes.front()->adjacent(surfaces.front())) compatible = false;
	if (!compatible) return NULL;
	// A surface-volume reaction is placed in the lower-dimensional surface.
	if (!surfaces.empty()) return surfaces.front();
	return volumes.empty() ? NULL : volumes.front();
}
}


list <Molecule *> TransformationSet::deleteList;
list <Molecule *> TransformationSet::updateAfterDeleteList;
list <Molecule *>::iterator TransformationSet::it;
void TransformationSet::initCommon()
{
	this->hasSymUnbinding = false;
	this->hasSymBinding = false;
	this->topologyChanging = false;
	this->hasProductSideCheck = false;
	this->productsCannotForm = false;

	// complex bookkeeping is off by default
	this->complex_bookkeeping = false;

	// for now, symmetry factor is off by default
	this->useSymmetryFactor = false;
	this->symmetryFactor = 1.0;

	// check collisions is off by default
	this->check_collisions = false;

	finalized = false;
}

TransformationSet::TransformationSet(vector <TemplateMolecule *> reactantTemplates) // @suppress("Class members should be properly initialized")
{
	this->initCommon();

	//cout<<"creating transformationSet..."<<endl;
	//Remember our reactants
	this->n_reactants = reactantTemplates.size();
	this->n_addmol  = 0;
	this->n_product_patterns = 1;

	this->reactants = new TemplateMolecule *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++)
		this->reactants[r] = reactantTemplates.at(r);

	this->addmol = new TemplateMolecule *[n_addmol];

	//Set up our transformation vectors
	this->transformations = new vector <Transformation *> [n_reactants];
}


TransformationSet::TransformationSet(vector <TemplateMolecule *> reactantTemplates,
		                             vector <TemplateMolecule *> addMoleculeTemplates )
{
	this->initCommon();

	//cout<<"creating transformationSet..."<<endl;
	//Remember our reactants
	this->n_reactants = reactantTemplates.size();
	this->n_product_patterns = 1;
	this->reactants = new TemplateMolecule *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++)
		this->reactants[r] = reactantTemplates.at(r);

	//Remember our add molecules
	this->n_addmol = addMoleculeTemplates.size();
	this->addmol = new TemplateMolecule *[n_addmol];
	for(unsigned int r=0; r<n_addmol; r++)
		this->addmol[r] = addMoleculeTemplates.at(r);

	//Set up our transformation vectors
	this->transformations = new vector <Transformation *> [ this->getNmappingSets() ];
}


TransformationSet::~TransformationSet()
{
	for(unsigned int r=0; r<getNmappingSets(); r++)  {
		Transformation *t;
		while(transformations[r].size()>0)
		{
			t = transformations[r].back();
			transformations[r].pop_back();
			delete t;
		}
	}

	Transformation *t;
	while(addMoleculeTransformations.size()>0)
	{
		t = addMoleculeTransformations.back();
		addMoleculeTransformations.pop_back();
		delete t;
	}

	while(addSpeciesTransformations.size()>0)
	{
		t = addSpeciesTransformations.back();
		addSpeciesTransformations.pop_back();
		delete t;
	}

	// Do NOT delete rf.parsedTemplates here. Every TemplateMolecule produced by
	// readPattern() for a reactant selector is registered with its MoleculeType
	// and freed when System deletes all MoleculeTypes (see system.cpp,
	// "Delete all MoleculeTypes (which deletes all molecules and templates)").
	// Deleting them here is a second free of the same objects, which double-frees
	// across the System teardown — intermittently on a single selector session,
	// deterministically once two selector-bearing sessions share a process
	// (SIGABRT/SIGSEGV). rf.pattern is one of these templates, so it must not be
	// deleted here either.
	reactantFilters.clear();
	productFilters.clear();

	delete [] transformations;
	delete [] reactants;
	delete [] addmol;

	this->n_reactants = 0;
	this->n_addmol = 0;
}


TemplateMolecule * TransformationSet::getTemplateMolecule( unsigned int reactantIndex ) const
{
	if ( reactantIndex < n_reactants )
	{
		return reactants[reactantIndex];
	}
	else if ( reactantIndex < getNmappingSets() )
	{
		return addmol[reactantIndex-n_reactants];
	}

	return 0;
}



bool TransformationSet::addStateChangeTransform(TemplateMolecule *t, string cName, int finalStateValue)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addStateChangeTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	// 2) Create a Transformation object to remember the information
	//cout<<"Adding state change transform to value: "<<finalStateValue<<endl;
	int cIndex = t->getMoleculeType()->getCompIndexFromName(cName);
	Transformation *transformation = TransformationFactory::genStateChangeTransform(cIndex, finalStateValue, t);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}

bool TransformationSet::addLocalFunctionReference(TemplateMolecule *t, string PointerName, int scope)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addStateChangeTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	Transformation *transformation = TransformationFactory::genLocalFunctionReference(PointerName,scope,t);
	transformations[reactantIndex].push_back(transformation);
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}


bool TransformationSet::addIncrementStateTransform(TemplateMolecule *t, string cName)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addIncrementStateTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	// 2) Create a Transformation object to remember the information
	int cIndex = t->getMoleculeType()->getCompIndexFromName(cName);
	Transformation *transformation = TransformationFactory::genIncrementStateTransform(cIndex, t);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}
bool TransformationSet::addDecrementStateTransform(TemplateMolecule *t, string cName)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addDecrementStateTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	// 2) Create a Transformation object to remember the information
	int cIndex = t->getMoleculeType()->getCompIndexFromName(cName);
	Transformation *transformation = TransformationFactory::genDecrementStateTransform(cIndex, t);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}



bool TransformationSet::addStateChangeTransform(TemplateMolecule *t, string cName, string finalStateValue)
{
	int cIndex = t->getMoleculeType()->getCompIndexFromName(cName);
	int fStateValue = t->getMoleculeType()->getStateValueFromName(cIndex,finalStateValue);
	return TransformationSet::addStateChangeTransform(t,cName, fStateValue);
}


bool TransformationSet::addBindingTransformImpl(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2, bool isNewMolecule)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	//Again, first find the reactants that the binding pertains to
	int reactantIndex1 = find(t1);
	int reactantIndex2 = find(t2);
	if(reactantIndex1==-1 || reactantIndex2==-1) {
		cerr<<"Couldn't find one of the templates you gave me!  In transformation set - addBindingTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	//Find the index of the respective binding sites
	unsigned int cIndex1 = t1->getMoleculeType()->getCompIndexFromName(bSiteName1);
	unsigned int cIndex2 = t2->getMoleculeType()->getCompIndexFromName(bSiteName2);

	//Check for symmetric binding
	bool isSymmetric = TemplateMolecule::checkSymmetry(t1,t2,bSiteName1,bSiteName2);
	if( isSymmetric )
		hasSymBinding = true;

	//Add transformation 1: Note that if both molecules involved with this bond are in the same reactant list, then
	//the mappingIndex will be size()+1.  But if they are on different reactant lists, then the mappingIndex will be exactly
	//equal to the size.
	Transformation *transformation1;
	unsigned int otherMappingIndex = transformations[reactantIndex2].size() + (reactantIndex1 == reactantIndex2 ? 1 : 0);

	if (isNewMolecule) {
		transformation1 = TransformationFactory::genNewMoleculeBindingTransform1(cIndex1, reactantIndex2, otherMappingIndex, t1);
	} else {
		transformation1 = TransformationFactory::genBindingTransform1(cIndex1, reactantIndex2, otherMappingIndex, t1);
	}

	Transformation *transformation2 = TransformationFactory::genBindingTransform2(cIndex2, t2);

	transformations[reactantIndex1].push_back(transformation1);
	MapGenerator *mg1 = new MapGenerator(transformations[reactantIndex1].size()-1);
	t1->addMapGenerator(mg1);

	transformations[reactantIndex2].push_back(transformation2);
	MapGenerator *mg2 = new MapGenerator(transformations[reactantIndex2].size()-1);
	t2->addMapGenerator(mg2);
	topologyChanging = true;

	return true;
}

bool TransformationSet::addBindingTransform(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2)
{
	return addBindingTransformImpl(t1, bSiteName1, t2, bSiteName2, false);
}

bool TransformationSet::addNewMoleculeBindingTransform(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2)
{
	return addBindingTransformImpl(t1, bSiteName1, t2, bSiteName2, true);
}



bool TransformationSet::addUnbindingTransform(TemplateMolecule *t, string bSiteName, TemplateMolecule *t2, string bSiteName2, bool separatesProducts, bool keepsProductsConnected)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }

	TemplateMolecule *tToTransform = 0;
	if(t==0 && t2==0) {
		cerr<<"Error in transformation set! when creating unbinding transform!"<<endl;
		cerr<<"Both molecules you gave me are null!\n";
		return false;
	} else if(t2==0) {
		tToTransform = t;
	} else if(t==0) {
		tToTransform = t2;
	} else {
		// they are both real, so randomly pick t1
		tToTransform=t;

		//Check for symmetric unbinding
		bool isSymmetric = TemplateMolecule::checkSymmetryAroundBond(t,t2,bSiteName,bSiteName2);
		if( isSymmetric ) {
			hasSymUnbinding = true;
		}

	}

	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(tToTransform);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set!"<<endl;
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you might get this error."<<endl;
		return false;
	}

	// 2) Create a Transformation object to remember the information
	unsigned int cIndex = tToTransform->getMoleculeType()->getCompIndexFromName(bSiteName);
	Transformation *transformation = TransformationFactory::genUnbindingTransform(cIndex, t);
	if (separatesProducts) {
		static_cast<UnbindingTransform *>( transformation )->setSeparatesProducts(true);
		hasProductSideCheck = true;
	}
	if (keepsProductsConnected) {
		static_cast<UnbindingTransform *>( transformation )->setKeepsProductsConnected(true);
		hasProductSideCheck = true;
	}

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 4) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	tToTransform->addMapGenerator(mg);

	// 4) Create an empty transformation for the binding partner so that
	// connectivity can be inferred. Arvind Rasi Subramaniam
	if (t2 != 0) {
		unsigned int cIndex2 = t2->getMoleculeType()->getCompIndexFromName(bSiteName2);
		Transformation *transformation2 = TransformationFactory::genUnbindingTransform2(cIndex2, t2);
		transformations[reactantIndex].push_back(transformation2);
		MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
		t2->addMapGenerator(mg);
	}
	topologyChanging = true;


	return true;
}





/*!
	Adds a delete rule to the given TemplateMolecule.
	@author Michael Sneddon
*/
bool TransformationSet::addDeleteMolecule(TemplateMolecule *t, int deletionType) {
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addDeleteMolecule!"<<endl;
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}
	Transformation *transformation = TransformationFactory::genRemoveMoleculeTransform(deletionType, t);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	topologyChanging = true;
	return true;
}


bool TransformationSet::addSurvivorOfDeletedMolecule(TemplateMolecule *survivor, int cluster, int productPattern)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	int reactantIndex = find(survivor);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addSurvivorOfDeletedMolecule!"<<endl;
		return false;
	}

	// one entry per molecule and cluster is enough
	for ( unsigned int i = 0; i < deletedMoleculeLinks.size(); ++i )
	{
		if ( deletedMoleculeLinks[i].survivor == survivor && deletedMoleculeLinks[i].cluster == cluster )
			return true;
	}

	// The rule may do nothing else to this molecule, in which case nothing maps
	// it. An empty transform gives it a mapping, so that the molecule it matched
	// can be found when the rule fires.
	transformations[reactantIndex].push_back( TransformationFactory::genEmptyTransform() );
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	survivor->addMapGenerator(mg);

	DeletedMoleculeLink link;
	link.survivor = survivor;
	link.reactantIndex = (unsigned int)reactantIndex;
	link.transformIndex = (unsigned int)(transformations[reactantIndex].size()-1);
	link.cluster = cluster;
	link.productPattern = productPattern;
	deletedMoleculeLinks.push_back(link);
	hasProductSideCheck = true;
	return true;
}


bool TransformationSet::addMoveTransform(TemplateMolecule *t, Compartment *c, bool moveConnected)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addMoveTransform!\n";
		return false;
	}

	std::vector<int> excludedMappingIndices;
	if (moveConnected) {
		std::vector<TemplateMolecule *> patternMolecules;
		TemplateMolecule::traverse(reactants[reactantIndex], patternMolecules,
				TemplateMolecule::FIND_ALL);
		for (TemplateMolecule *patternMolecule : patternMolecules) {
			if (patternMolecule == t) continue;
			transformations[reactantIndex].push_back(new EmptyTransform(-1, patternMolecule));
			int mappingIndex = static_cast<int>(transformations[reactantIndex].size() - 1);
			patternMolecule->addMapGenerator(new MapGenerator(mappingIndex));
			excludedMappingIndices.push_back(mappingIndex);
		}
	}
	Transformation *transformation = new MoveTransformation(
			c, t, moveConnected, reactantIndex, excludedMappingIndices);
	transformations[reactantIndex].push_back(transformation);

	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}


/*!
	Adds a decrement population rule to the given TemplateMolecule.
	@author Justin Hogg
*/
bool TransformationSet::addDecrementPopulation(TemplateMolecule *t)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addDecrementPopulation!"<<endl;
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}
	Transformation *transformation = TransformationFactory::genDecrementPopulationTransform(t);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}

bool TransformationSet::addAddSpecies( SpeciesCreator *sc )
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }

	// We don't need a polymorphic transform because AddTransforms are handled separately!
	//  But we do need to call some methods specific to AddMoleculeTransform.
	//  So we're modified TransformationFactory to return the specific object type  --JUstin
	AddSpeciesTransform * transformation = TransformationFactory::genAddSpeciesTransform( sc );

	// 3) Add the transformation object to the TransformationSet
	addSpeciesTransformations.push_back( transformation );

	// 3) No map generators needed for an add species!
	topologyChanging = true;
	return true;
}


bool TransformationSet::addAddMolecule( MoleculeCreator *mc )
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }

	// We don't need a polymorphic transform because AddTransforms are handled separately!
	//  But we do need to call some methods specific to AddMoleculeTransform.
	//  So we're modified TransformationFactory to return the specific object type  --JUstin
	AddMoleculeTransform * transformation = TransformationFactory::genAddMoleculeTransform( mc, mc->getTemplateMolecule() );

	// 3) Add the transformation object to the TransformationSet
	addMoleculeTransformations.push_back( transformation );

	// 3) No map generators needed for an add molecule!
	topologyChanging = true;
	return true;
}




bool TransformationSet::canReachExcludingBonds(Molecule *mol1, Molecule *mol2,
		const std::vector< std::pair<Molecule *, int> > &excludedBonds,
		const std::vector< std::pair<Molecule *, Molecule *> > *addedBonds)
{
	// Perform a BFS from mol2 to see if mol1 is reachable while refusing to
	// traverse any of the bonds in excludedBonds. Each excluded bond is recorded
	// as both of its (molecule, component-index) half-edges, so the bond is
	// skipped regardless of which endpoint the BFS reaches first.
	// Returns true if mol1 can still be reached, false otherwise.
	//
	// Excluding the full set of bonds that the firing rule deletes (rather than
	// one bond at a time) is what lets a multi-bond ring-opening dissociation be
	// recognized as genuinely separating its products: each individual bond may
	// leave the partners connected, while removing all of them does not.

	// Use a static queue for efficiency (reusing the BFS infrastructure)
	static queue <Molecule *> q;
	static list <Molecule *> visited;

	// Clear queues and visited list
	while(!q.empty()) q.pop();
	visited.clear();
	// Start BFS from mol2
	q.push(mol2);
	visited.push_back(mol2);
	mol2->setVisitedMolecule(true);

	bool found = false;

	while(!q.empty()) {
		Molecule *current = q.front();
		q.pop();

		// Check if we reached mol1
		if (current == mol1) {
			found = true;
			break;
		}

		// Explore neighbors through bonds
		int numComponents = current->getMoleculeType()->getNumOfComponents();
		for (int c = 0; c < numComponents; ++c) {
			// Skip any half-edge belonging to a bond the rule will delete.
			bool isExcluded = false;
			for (unsigned int e = 0; e < excludedBonds.size(); ++e) {
				if (excludedBonds[e].first == current && excludedBonds[e].second == c) {
					isExcluded = true;
					break;
				}
			}
			if (isExcluded) {
				continue;
			}

			// Check if this component is bonded
			if (current->isBindingSiteBonded(c)) {
				Molecule *neighbor = current->getBondedMolecule(c);

				// Check if we've already visited this neighbor
				if (!neighbor->getVisitedMolecule()) {
					// Haven't visited this neighbor yet
					neighbor->setVisitedMolecule(true);
					visited.push_back(neighbor);
					q.push(neighbor);
				}
			}
		}

		// Also cross any bond the rule is about to add at this molecule.
		if (addedBonds != 0) {
			for (unsigned int a = 0; a < addedBonds->size(); ++a) {
				Molecule *neighbor = 0;
				if ((*addedBonds)[a].first == current) neighbor = (*addedBonds)[a].second;
				else if ((*addedBonds)[a].second == current) neighbor = (*addedBonds)[a].first;
				if (neighbor != 0 && !neighbor->getVisitedMolecule()) {
					neighbor->setVisitedMolecule(true);
					visited.push_back(neighbor);
					q.push(neighbor);
				}
			}
		}
	}
	
	// Clear the hasVisitedMolecule values before returning
	for (list<Molecule *>::iterator it = visited.begin(); it != visited.end(); ++it) {
		(*it)->setVisitedMolecule(false);
	}

	// Return whether mol1 is reachable from mol2 when excluding the bond
	return found;
}


int TransformationSet::find(TemplateMolecule *t)
{
	if(finalized) { cerr<<"TransformationSet cannot search for a templateMolecule once it has been finalized!"<<endl; exit(1); }
	int findIndex = -1;
	for(unsigned int r=0; r<n_reactants; r++)  {
		if(this->reactants[r]->contains(t)) {
			if(findIndex==-1) {
				findIndex = r;
			}
			else {
				cerr<<"Found duplicate template molecule in two reaction lists!!  (in transformationSet)."<<endl;
				exit(1);
			}
		}
	}
	// also check add molecule templates
	for(unsigned int r=0; r<n_addmol; r++)  {
		if(this->addmol[r]->contains(t)) {
			if(findIndex==-1) {
				findIndex = r + n_reactants;
			}
			else {
				cerr<<"Found duplicate template molecule in two reaction lists!!  (in transformationSet)."<<endl;
				exit(1);
			}
		}
	}
	return findIndex;
}
// AS2023 - normal calls should have tracking off
string TransformationSet::transform(MappingSet **mappingSets)
{
	return this->transform(mappingSets, false);
}
// AS2023 - alternative call sig to store a log of the transform
string TransformationSet::transform(MappingSet **mappingSets, bool tracking)
{
	if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }
	System *profileSystem = 0;
	if (n_reactants > 0 && reactants[0] != 0 &&
			reactants[0]->getMoleculeType() != 0)
		profileSystem = reactants[0]->getMoleculeType()->getSystem();
	bool profile = profileSystem != 0 && profileSystem->isProfileReactionActive();
	ProfileConnectivityScope profileConnectivityScope(
		profileSystem, PROFILE_CONNECTIVITY_TRANSFORMATION);
	ProfileTime profileStart = profile ? profileNow() : ProfileTime();

	/*
	 * NOTE: Check for "null conditions" was moved to ReactionClass::fire. This allows rejection of a reaction 
	 * prior to removing molecules from observables. A general method TransformationSet::checkMolecularity has
	 * been implemented to check for incorrect molecularity or reaction center conflicts. --Justin
	 */

	// AS2023 - if we are tracking, initialize a log string
	string logstr;
	if (tracking) {
		logstr = "        \"ops\": [\n";
	} else {
		logstr = "";
	}

	// addMolecule transforms are applied before other transforms so the molecules exist
	//  for potential modification by other transforms.
	int size = addMoleculeTransformations.size();
	if(size>0) {
		for(int i=0; i<size; i++) {
			// AS2023 - since we are in the tracking call, track the application
			bool compatibleCompartment = true;
			Compartment *inheritedCompartment = selectReactionCompartment(
					mappingSets, n_reactants, compatibleCompartment);
			addMoleculeTransformations.at(i)->apply_and_map(
					mappingSets[n_reactants+i], logstr, inheritedCompartment);
		}
	}

	// apply addSpecies transforms, so we have all the molecules out there
	size = addSpeciesTransformations.size();
	if(size>0) {
		for(int i=0; i<size; i++) {
			// AS2023 - since we are in the tracking call, track the application
			addSpeciesTransformations.at(i)->apply(NULL, NULL, logstr);
		}
	}

	// loop over reactants and added molecules, apply transforms to each
	for(unsigned int r=0; r<getNmappingSets(); r++)
	{
		MappingSet *ms = mappingSets[r];
		for ( unsigned int t=0;  t<transformations[r].size();  t++ )
		{
			if( transformations[r].at(t)->getType()==(int)TransformationFactory::REMOVE )
			{	// handle deletions
				Mapping *mapObj = ms->get(t);
				if (!mapObj) {
					// This can happen when a MappingSet isn't properly populated due to
					// internal bond reconnection. Skip the deletion transform to avoid crashing.
					continue;
				}
				Molecule * mol = mapObj->getMolecule();
				if (!mol) {
					// Guard against null molecule mappings (see issue with internal bond reconnection).
					continue;
				}
				if ( transformations[r].at(t)->getRemovalType()==(int)TransformationFactory::COMPLETE_SPECIES_REMOVAL )
				{	// complex deletion: flag connected molecules for deletion
					// AS2023 - since we are in the tracking call, track the deletion events
					mol->traverseBondedNeighborhood(deleteList,ReactionClass::NO_LIMIT, logstr);
				}
				else
				{	// molecule deletion: flag this molecule for deletion
					// track deletions if tracking is on
					// this has to be done here
					if (tracking) {
						// AS2023 - since we are in the tracking call, track the deletion operations
						logstr += "          [\"Delete\"," + to_string(mol->getUniqueID()) + "],\n";
					}
					deleteList.push_back( mol );
				}
			}
			else
			{	// handle other transforms
				// AS2023 - since we are in the tracking call, track the operation
				transformations[r].at(t)->apply(ms->get(t), mappingSets, logstr);
			}
		}
	}


	//Each molecule that is on the delete list must be dealt with
	Molecule * mol;
	for( it = deleteList.begin(); it!=deleteList.end(); it++)
	{
		mol = *it;
		mol->getMoleculeType()->removeMoleculeFromRunningSystem(mol);
	}
	deleteList.clear();

	// finalize tracking
	if (tracking) {
		// need to delete the last comma if we actually added anything
		if (logstr.length() > 25) {  // 25 is roughly the length of "        \"ops\": [\n"
			logstr.erase(logstr.end()-2, logstr.end());
			logstr += "\n        ]\n";
		} else {
			logstr += "        ]\n";
		}
	}
	if (profile)
		profileSystem->recordProfileTransformation(profileElapsedSeconds(profileStart));
	return logstr;
}


bool TransformationSet::deletedBondsMatchProducts( MappingSet ** mappingSets )
{
	// First pass: gather every bond this rule's UNBINDING transforms will
	// delete. The connectivity test below must reflect the post-reaction
	// graph with ALL of them removed at once. Testing each bond in
	// isolation wrongly blocks multi-bond ring-opening dissociations
	// (e.g. a symmetric two-bond homodimer splitting into two monomers):
	// each single bond's removal leaves the partners connected through
	// the rule's other deleted bond, even though removing both genuinely
	// separates the products.
	//
	// Every deleted bond is excluded from the walk. What is then required of
	// its two molecules depends on how the rule writes them: with "+" they
	// must end up in different complexes, with "." they must still be in one.
	// Bonds the rule adds are walked as if already in place, so that a rule
	// which moves a bond from one pair of sites to another is judged on the
	// graph it leaves behind.
	std::vector< std::pair<Molecule *, int> > excludedBonds;          // half-edges
	std::vector< std::pair<Molecule *, Molecule *> > brokenPairs;     // must separate
	std::vector< std::pair<Molecule *, Molecule *> > connectedPairs;  // must stay connected
	std::vector< std::pair<Molecule *, Molecule *> > addedBonds;      // about to be created
	for ( unsigned int r = 0; r < n_reactants; ++r )
	{
		for ( unsigned int t = 0; t < getNumOfTransformations(r); ++t )
		{
			Transformation * transform = getTransformation(r, t);
			if ( transform->getType() == TransformationFactory::BINDING )
			{
				// Only bonds between molecules that already exist: a molecule the
				// rule creates is not there yet to be walked.
				BindingTransform * bindingTransform =
					static_cast<BindingTransform *>( transform );
				int otherReactant = bindingTransform->getOtherReactantIndex();
				if ( otherReactant < 0 || otherReactant >= (int)n_reactants ) continue;
				Mapping * mapping1 = mappingSets[r]->get(t);
				Mapping * mapping2 = mappingSets[otherReactant]->get(
						bindingTransform->getOtherMappingIndex() );
				if ( mapping1 == NULL || mapping2 == NULL ) continue;
				if ( mapping1->getMolecule() == NULL || mapping2->getMolecule() == NULL ) continue;
				addedBonds.push_back( std::make_pair(mapping1->getMolecule(), mapping2->getMolecule()) );
				continue;
			}
			if ( transform->getType() == TransformationFactory::REMOVE
					&& transform->getRemovalType() == (int)TransformationFactory::DELETE_MOLECULES )
			{
				// A molecule the rule deletes takes all of its bonds with it, so
				// none of them may be walked.
				Mapping * mapping = mappingSets[r]->get(t);
				if ( mapping == NULL || mapping->getMolecule() == NULL ) continue;
				Molecule * deletedMol = mapping->getMolecule();
				int numComponents = deletedMol->getMoleculeType()->getNumOfComponents();
				for ( int c = 0; c < numComponents; ++c )
				{
					if ( !deletedMol->isBindingSiteBonded(c) ) continue;
					Molecule * partner = deletedMol->getBondedMolecule(c);
					if ( partner == NULL ) continue;
					excludedBonds.push_back( std::make_pair(deletedMol, c) );
					excludedBonds.push_back( std::make_pair(partner,
							deletedMol->getBondedMoleculeBindingSiteIndex(c)) );
				}
				continue;
			}
			if ( transform->getType() != TransformationFactory::UNBINDING )
				continue;

			UnbindingTransform * unbindingTransform =
				static_cast<UnbindingTransform *>( transform );

			Mapping * mapping = mappingSets[r]->get(t);
			if (mapping == NULL) {
				continue;
			}

			// Get the molecule that will be unbound
			Molecule * mol = mapping->getMolecule();
			if (mol == NULL) {
				continue;
			}
			int componentIndex = unbindingTransform->getComponentIndex();

			// Get the molecule bonded at this site
			if ( !mol->isBindingSiteBonded(componentIndex) )
			{
				// Site is not bonded, so this won't actually unbind anything
				continue;
			}

			Molecule * bondedMol = mol->getBondedMolecule(componentIndex);
			if ( bondedMol == NULL )
			{
				// This shouldn't happen, but skip if it does
				continue;
			}

			int partnerIndex = mol->getBondedMoleculeBindingSiteIndex(componentIndex);
			excludedBonds.push_back( std::make_pair(mol, componentIndex) );
			excludedBonds.push_back( std::make_pair(bondedMol, partnerIndex) );
			if ( unbindingTransform->getSeparatesProducts() )
				brokenPairs.push_back( std::make_pair(mol, bondedMol) );
			else if ( unbindingTransform->getKeepsProductsConnected() )
				connectedPairs.push_back( std::make_pair(mol, bondedMol) );
		}
	}

	// Molecules left behind by a molecule the rule deletes. Two of them that
	// were held by the same deleted molecule, or by deleted molecules bonded
	// to one another, are judged like the two ends of a deleted bond: in
	// different product patterns they must separate, in the same one they
	// must stay connected.
	if ( !deletedMoleculeLinks.empty() )
	{
		std::vector< Molecule * > survivors( deletedMoleculeLinks.size(), (Molecule *)0 );
		for ( unsigned int i = 0; i < deletedMoleculeLinks.size(); ++i )
		{
			const DeletedMoleculeLink & link = deletedMoleculeLinks[i];
			Mapping * mapping = mappingSets[link.reactantIndex]->get(link.transformIndex);
			if ( mapping == NULL ) continue;
			survivors[i] = mapping->getMolecule();
		}
		for ( unsigned int i = 0; i < deletedMoleculeLinks.size(); ++i )
		{
			if ( survivors[i] == 0 ) continue;
			for ( unsigned int j = i + 1; j < deletedMoleculeLinks.size(); ++j )
			{
				if ( survivors[j] == 0 || survivors[j] == survivors[i] ) continue;
				if ( deletedMoleculeLinks[i].cluster != deletedMoleculeLinks[j].cluster ) continue;
				if ( deletedMoleculeLinks[i].productPattern != deletedMoleculeLinks[j].productPattern )
					brokenPairs.push_back( std::make_pair(survivors[i], survivors[j]) );
				else
					connectedPairs.push_back( std::make_pair(survivors[i], survivors[j]) );
			}
		}
	}

	// Second pass: with every deleted bond removed, each pair the rule
	// separates must end up in different connected components. If any pair
	// is still mutually reachable, the products do not actually separate,
	// which violates the product-side molecularity constraint -> refuse to fire.
	for ( unsigned int b = 0; b < brokenPairs.size(); ++b )
	{
		if ( canReachExcludingBonds( brokenPairs[b].first, brokenPairs[b].second, excludedBonds, &addedBonds ) )
		{
			return false;
		}
	}

	// And the other way round: each pair the rule keeps in one product
	// pattern must still be mutually reachable.
	for ( unsigned int b = 0; b < connectedPairs.size(); ++b )
	{
		if ( !canReachExcludingBonds( connectedPairs[b].first, connectedPairs[b].second, excludedBonds, &addedBonds ) )
		{
			return false;
		}
	}
	return true;
}


bool TransformationSet::checkMolecularity( MappingSet ** mappingSets )
{
	// A rule that writes molecules into one product pattern without ever
	// joining them cannot satisfy its own "." and so never fires.
	if ( complex_bookkeeping && productsCannotForm )
		return false;

	// Separate reactant species can interact only when they all have a
	// compartment, or all are unplaced. Compartmented reactions may span one
	// 3-D volume and one adjacent 2-D surface (BNG2 interactingSet semantics).
	if (n_reactants > 1) {
		bool compatibleCompartment = true;
		selectReactionCompartment(mappingSets, n_reactants, compatibleCompartment);
		if (!compatibleCompartment) return false;
	}
	if ( n_reactants < 2 )
	{	// unimolecular, so there's nothing to check on the reactant side.
		// HOWEVER: if the rule deletes a bond and writes "+" between the two
		// molecules it joined, we need to verify that they end up in separate
		// connected components. This enforces the BNGL semantics that
		// "A.B -> A + B" should only fire if breaking the bond actually
		// separates the molecules into different complexes, and that
		// "A.B -> A.B" with a bond deleted should only fire if it does not.
		// Issue #48: NFsim does not enforce product molecularity for unimolecular unbinding rules
		if ( n_reactants == 1 && complex_bookkeeping && hasProductSideCheck )
			return deletedBondsMatchProducts( mappingSets );
		return true;
	}
	else if ( complex_bookkeeping )
	{	// verify that each reactant pattern points to a unique complex
		if (n_reactants == 2 &&
				!reactants[0]->getMoleculeType()->isPopulationType() &&
				!reactants[1]->getMoleculeType()->isPopulationType()) {
			/* The general path below uses a hash set for arbitrary molecularity.
			 * For the overwhelmingly common bimolecular case, uniqueness is
			 * exactly one complex-ID comparison and does not need a container. */
			if (mappingSets[0]->getComplexID() == mappingSets[1]->getComplexID())
				return false;
		}
		else
		{
			complex_ids.clear();
			for ( unsigned int ir = 0;  ir < n_reactants;  ++ir )
			{
				// skip populations
				if ( reactants[ir]->getMoleculeType()->isPopulationType() ) continue;

				complex_id = mappingSets[ir]->getComplexID();
				if ( !complex_ids.insert( complex_id ).second )
				{   // two reactant patterns matched the same complex!
					return false;
				}
			}
		}

		// The same product-side rule as for one reactant pattern: a rule such as
		// "A.B + C -> A + B + C" may fire only if A and B really come apart.
		if ( hasProductSideCheck )
			return deletedBondsMatchProducts( mappingSets );
		return true;
	}
	else if ( check_collisions )
	{	// we won't do a proper check for molecularity, but we should ensure that mappingSets
		//  point to non-overlapping reaction centers.
		for ( collision_pair_iter = collision_pairs.begin(); collision_pair_iter != collision_pairs.end(); ++collision_pair_iter )
		{
			if ( MappingSet::checkForCollisions( mappingSets[ (*collision_pair_iter).first  ],
					                             mappingSets[ (*collision_pair_iter).second ]  ) )
			{	// reaction centers overlap!
				return false;
			}
		}
		return true;
	}
	else
	{   // do nothing
		return true;
	}
}


bool TransformationSet::getListOfProducts(
		MappingSet **mappingSets, list <Molecule *> &products, int traversalLimit,
		vector <unsigned int> *componentSizes, bool *componentsTruncated)
{
	System *profileSystem = 0;
	if (n_reactants > 0 && reactants[0] != 0 &&
			reactants[0]->getMoleculeType() != 0)
		profileSystem = reactants[0]->getMoleculeType()->getSystem();
	bool profile = profileSystem != 0 && profileSystem->isProfileReactionActive();
	ProfileConnectivityScope profileConnectivityScope(
		profileSystem, PROFILE_CONNECTIVITY_PRODUCT_PREPARATION);
	unsigned long long productsBefore = profile
		? static_cast<unsigned long long>(products.size()) : 0;
	ProfileTime profileStart = profile ? profileNow() : ProfileTime();
	std::unordered_set<Molecule*> product_set(products.begin(), products.end());
	list <Molecule *>::iterator molIter;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		// if we are deleting the entire complex, we don't have to track molecules in this complex
		if (mappingSets[r]->hasSpeciesDeletionTransform()) continue;

		//cout<<"Traversing:"<<endl;
		//mappingSets[r]->get(0)->getMolecule()->printDetails();
		//mappingSets[r]->get(0)->getMolecule()->traverseBondedNeighborhood(products,traversalLimit);

		/*
		 * I thought that making sure we don't go over the same molecule multiple
		 * times would make the code faster - but this is rarely used for most rxn
		 * systems, so it is commented out for now.  But actually, we do have to check
		 * because if have the same molecule in here twice, then it can mess up our
		 * observable lists...  --michael */
		else
		{
			// For each of the molecules that we possibly affect, traverse the neighborhood
			// Q: Is it sufficient to just look at the first mapping?
			// A: It should be if the traversal limit is set high enough, at least for
			// all standard reactions.  I'm wondering now, though, if it is enough in
			// all cases where you would use the connected-to syntax.  I think so, but
			// someone should test it.  --michael 9Mar2011
			Molecule * molecule = mappingSets[r]->get(0)->getMolecule();

			// is this molecule already on the product list?
			if ( product_set.find( molecule ) == product_set.end() )
			{	// Traverse neighbor and add molecules to list
				bool was_empty = products.empty();
				unsigned int productsBeforeTraversal = products.size();
				auto last = was_empty ? products.end() : std::prev(products.end());
				bool traversalTruncated = molecule->traverseBondedNeighborhood(
						products, traversalLimit);
				if (traversalTruncated && componentsTruncated != 0)
					*componentsTruncated = true;
				if (componentSizes != 0)
					componentSizes->push_back(
						static_cast<unsigned int>(products.size() - productsBeforeTraversal));
				// Sync only newly appended molecules into the set
				auto it = was_empty ? products.begin() : std::next(last);
				for (; it != products.end(); ++it) {
					product_set.insert(*it);
				}
				//molecule->traverseBondedNeighborhoodForUpdate(products,traversalLimit);
			}
		}
	}

	// Next, find added molecules that are treated as populations.
	//  Populations molecules have to be removed from observables, then incremented,
	//  and then added back to the observables (Add molecules treated as particles are handled later)
	vector <AddMoleculeTransform *>::iterator addmol_iter;
	for ( addmol_iter = addMoleculeTransformations.begin();
			addmol_iter != addMoleculeTransformations.end();  ++addmol_iter )
	{
		// get molecule creator
		AddMoleculeTransform * addmol = *addmol_iter;
		if ( !(addmol->isPopulationType()) ) continue;

		// Get the population molecule pointer
		Molecule * molecule = addmol->get_population_pointer();

		// is this molecule already on the product list?
		if ( product_set.find( molecule ) == product_set.end() )
		{	// Add molecule to list
			products.push_back( molecule );
			product_set.insert( molecule );
			if (componentSizes != 0)
				componentSizes->push_back(1);
		}
	}

	//cout<<"All together, we have: "<<products.size()<<endl;
	if (profile)
		profileSystem->recordProfileProductPreparation(
				profileElapsedSeconds(profileStart),
				static_cast<unsigned long long>(products.size()) - productsBefore);
	return true;
}


Molecule * TransformationSet::getPopulationPointer( unsigned int r ) const
{
	return addMoleculeTransformations.at(r)->isPopulationType()
			? addMoleculeTransformations.at(r)->get_population_pointer()
		    : NULL;
}

bool TransformationSet::getListOfAddedMolecules(MappingSet **mappingSets, list <Molecule *> &products, int traversalLimit)
{
	System *profileSystem = 0;
	if (n_reactants > 0 && reactants[0] != 0 &&
			reactants[0]->getMoleculeType() != 0)
		profileSystem = reactants[0]->getMoleculeType()->getSystem();
	bool profile = profileSystem != 0 && profileSystem->isProfileReactionActive();
	unsigned long long productsBefore = profile
		? static_cast<unsigned long long>(products.size()) : 0;
	ProfileTime profileStart = profile ? profileNow() : ProfileTime();
	if (addMoleculeTransformations.empty()) {
		if (profile)
			profileSystem->recordProfileProductCollection(0.0, 0);
		return true;
	}
	std::unordered_set<Molecule*> product_set(products.begin(), products.end());

	// Add new molecules (particle type) to the list of products
	list <Molecule *>::iterator molIter;
	for (unsigned int r=n_reactants; r<getNmappingSets(); r++)
	{
		//For each of the molecules that we possibly affect, traverse the neighborhood
		// NOTE: in this instance, it's okay to only look at the first mapping
		Molecule * molecule = mappingSets[r]->get(0)->getMolecule();

		// Skip populations
		if ( molecule->isPopulationType() ) continue;

		// Is the molecule already in the products list?  If not, add to list.
		if ( product_set.find( molecule ) == product_set.end() )
		{	// Add molecule to list.
			products.push_back( molecule );
			product_set.insert( molecule );
			// NOTE: we don't need to traverse neighbors. All new molecules will be put in this
			//  list separately and old molecules that bind to new molecules will be traversed elsewhere
		}
	}
	if (profile)
		profileSystem->recordProfileProductCollection(
				profileElapsedSeconds(profileStart),
				static_cast<unsigned long long>(products.size()) - productsBefore);

	return true;
}


MappingSet *TransformationSet::generateBlankMappingSet(unsigned int reactantIndex, unsigned int mappingSetId)
{
	if(!finalized) { cerr<<"TransformationSet cannot generate blank mapping if it is not finalized!"<<endl; exit(1); }
	if( reactantIndex>=getNmappingSets() ) {
		cerr<<"Gave me (a transformation Set) a reactant index that was too high!"<<endl;
		exit(1);
	}
	return new MappingSet(mappingSetId, transformations[reactantIndex]);
}


bool TransformationSet::isPureContextReactant(unsigned int reactantIndex) const
{
	if( reactantIndex >= n_reactants ) return false;
	return pureContextReactants.find(reactantIndex) != pureContextReactants.end();
}


void TransformationSet::finalize()
{
	//Be sure to add at least a blank transformation to every reactant if there is no transformation
	//specified so that we count the reactants even if we don't do anything to it.
	for(unsigned int r=0; r<getNmappingSets(); r++)  {

		// Record which reactants the rule does not transform, before the placeholder
		// below makes that indistinguishable.  EMPTY cannot be the test: it is also
		// what genBindingTransform2() puts on the second partner of a binding, which
		// very much is a reaction center.  A LOCAL_FUNCTION_REFERENCE, on the other
		// hand, only marks the molecule a local function reads, so a DOR reactant
		// can still be pure context.
		bool pureContext = true;
		for(unsigned int t=0; t<transformations[r].size(); t++) {
			if(transformations[r].at(t)->getType() != (int)TransformationFactory::LOCAL_FUNCTION_REFERENCE) {
				pureContext = false;
				break;
			}
		}
		if(pureContext) pureContextReactants.insert(r);

		if(transformations[r].size()==0) {
			transformations[r].push_back(TransformationFactory::genEmptyTransform());
			MapGenerator *mg = new MapGenerator(transformations[r].size()-1);
			getTemplateMolecule(r)->addMapGenerator(mg);
		}
	}

	// Determine if we need to do any reactant center overlap checks.
	// Currently we check a necessary (but not sufficient) condition for the
	//   possibility of reactant center overlap: are there a common molecule types in
	//   a pair of reactant templates. In the future, we could check a necessary and sufficent condition
	//   (e.g. pattern overlap) to avoid extra work.
	if ( (n_reactants>1)  &&  !complex_bookkeeping )
	{
		vector <TemplateMolecule *> tmList1;
		vector <TemplateMolecule *> tmList2;
		vector <TemplateMolecule *>::iterator tm_iter;

		vector <MoleculeType *> moltypes;
		vector <MoleculeType *>::iterator  found_iter;

		for ( unsigned int ir1 = 0;  ir1 < n_reactants;  ++ir1 )
		{
			// skip populations
			if ( reactants[ir1]->getMoleculeType()->isPopulationType() ) continue;

			tmList1.clear();
			moltypes.clear();

			// get all the template molecules in reactant pattern ir1
			TemplateMolecule::traverse(reactants[ir1], tmList1, false);

			// collect the molecule types included in the pattern
			for ( tm_iter = tmList1.begin(); tm_iter != tmList1.end();  ++tm_iter )
			{
				moltypes.push_back( (*tm_iter)->getMoleculeType() );
			}

			for ( unsigned  int ir2 = ir1 + 1;  ir2 < n_reactants;  ++ir2 )
			{
				// skip populations
				if ( reactants[ir2]->getMoleculeType()->isPopulationType() ) continue;

				tmList2.clear();

				// get all the template molecules in reactant pattern ir2
				TemplateMolecule::traverse(reactants[ir2], tmList2, false);

				// check if any moleculeTypes collide
				for ( tm_iter = tmList2.begin(); tm_iter != tmList2.end();  ++tm_iter )
				{
					found_iter = std::find( moltypes.begin(), moltypes.end(), (*tm_iter)->getMoleculeType() );
					if ( found_iter != moltypes.end() )
					{
						// we found a molecule type that is commont to patterns ir1 and ir2!
						check_collisions = true;
						collision_pairs.push_back( pair<int,int>(ir1,ir2) );
						break;
					}
				}
			}
		}

	}

	finalized = true;
}

bool TransformationSet::checkConnection(ReactionClass * rxn) {
	TemplateMolecule * t1;
	Transformation * transfn;
	for(unsigned int r=0; r<n_reactants; r++) {
		for (unsigned int i=0; i<transformations[r].size(); i++) {
			transfn = transformations[r].at(i);
			t1 = transfn->getTemplateMolecule();
			if (!t1) continue;
			// AS2023 - if this is not a removal, track connections, removal
			// doesn't give any reaction connections, so skip that
			if (transfn->getType()!=(int)TransformationFactory::REMOVE) {
				// A compartment-changing rule can enable a pattern in its
				// destination compartment even when the pattern constrains a
				// disjoint component. Preserve the broad fallback for MOVE.
				bool isCompatible = transfn->getType() == (int)TransformationFactory::MOVE
						? rxn->isTemplateCompatible(t1)
						: rxn->isTemplateCompatibleWithSharedComponent(t1);
				if (!isCompatible) continue;
				// Basic rules require a shared specified component; functional
				// rules retain the broader observable-dependent fallback.
				// Both checks passed for one op so return true
				return true;
			} else {
				return false;
			}
		}
	}
	// Do the same as above but now for the transformed product templates
	for(unsigned int r=0; r<n_reactants; r++) {
		for (unsigned int i=0; i<transformations[r].size(); i++) {
			transfn = transformations[r].at(i);
			// Use the transformed product here
			t1 = transfn->getTemplateMolecule();
			if (!t1) continue;
			t1 = t1->getMappedPartner();
			if (!t1) continue;
			bool isMove = transfn->getType() == (int)TransformationFactory::MOVE;
			bool isCompatible = isMove ? rxn->isTemplateCompatible(t1)
								: rxn->isTemplateCompatibleWithSharedComponent(t1);
			if (!isCompatible) continue;
			// Apply the same component gate to the mapped product.
			// Both checks passed for one op so return true
			return true;
		}
	}
	// Now consider the newly added molecules where components don't matter, but
	// compatibility
	for (unsigned int i=0; i<addMoleculeTransformations.size(); i++) {
		transfn = addMoleculeTransformations.at(i);
		// Use the transformed product here
		t1 = transfn->getTemplateMolecule();
		if (!t1) continue;
		// If the TemplateMolecule is 'incompatible' with any of the reactants
		// or products, then the reaction is not connected
		// Creating a molecule can add a matching target complex even when none
		// of its component sites overlap a target rule's specified sites.
		if (!rxn->isTemplateCompatible(t1)) continue;
		// Both checks passed for one op so return true
		return true;
	}
	// Both checks did not pass for any reactant or product template, so not connected
	return false;
}

void TransformationSet::addExcludeReactant(int reactantIndex, TemplateMolecule *pattern, const map<string, TemplateMolecule*>& parsedTemplates) {
	ReactantFilter rf;
	rf.reactantIndex = reactantIndex;
	rf.pattern = pattern;
	rf.isExclude = true;
	rf.parsedTemplates = parsedTemplates;
	reactantFilters.push_back(rf);
}

void TransformationSet::addIncludeReactant(int reactantIndex, TemplateMolecule *pattern, const map<string, TemplateMolecule*>& parsedTemplates) {
	ReactantFilter rf;
	rf.reactantIndex = reactantIndex;
	rf.pattern = pattern;
	rf.isExclude = false;
	rf.parsedTemplates = parsedTemplates;
	reactantFilters.push_back(rf);
}

bool TransformationSet::checkReactantFilters(int reactantIndex, Molecule *mol) const {
	if (reactantFilters.empty()) return true;

	for (const auto &rf : reactantFilters) {
		if (rf.reactantIndex != reactantIndex) continue;

		// Check if any molecule in the complex matches the filter pattern
		list<Molecule *> complexMembers;
		mol->traverseBondedNeighborhood(complexMembers, ReactionClass::NO_LIMIT);

		bool patternMatches = false;
		for (Molecule *cm : complexMembers) {
			if (rf.pattern->compare(cm)) {
				patternMatches = true;
				break;
			}
		}

		if (rf.isExclude && patternMatches) {
			return false;
		}
		if (!rf.isExclude && !patternMatches) {
			return false;
		}
	}
	return true;
}

void TransformationSet::addExcludeProduct(int productIndex, TemplateMolecule *pattern, const map<string, TemplateMolecule*>& parsedTemplates) {
	ProductFilter pf;
	pf.productIndex = productIndex;
	pf.pattern = pattern;
	pf.isExclude = true;
	pf.parsedTemplates = parsedTemplates;
	productFilters.push_back(pf);
}

void TransformationSet::addIncludeProduct(int productIndex, TemplateMolecule *pattern, const map<string, TemplateMolecule*>& parsedTemplates) {
	ProductFilter pf;
	pf.productIndex = productIndex;
	pf.pattern = pattern;
	pf.isExclude = false;
	pf.parsedTemplates = parsedTemplates;
	productFilters.push_back(pf);
}

bool TransformationSet::checkProductFilters(const list<Molecule *> &products) const {
	if (productFilters.empty()) return true;

	// Collect unique complexes from the product molecule list
	unordered_set<Complex*> productComplexes;
	for (Molecule *mol : products) {
		if (mol == 0 || !mol->isAlive()) continue;
		productComplexes.insert(mol->getComplex());
	}

	for (const auto &pf : productFilters) {
		bool anyComplexMatches = false;
		for (Complex *c : productComplexes) {
			bool patternMatches = false;
			for (Molecule *cm : c->complexMembers) {
				if (pf.pattern->compare(cm)) {
					patternMatches = true;
					break;
				}
			}
			if (patternMatches) {
				anyComplexMatches = true;
				break;
			}
		}

		if (pf.isExclude && anyComplexMatches) {
			return false;
		}
		if (!pf.isExclude && !anyComplexMatches) {
			return false;
		}
	}
	return true;
}
