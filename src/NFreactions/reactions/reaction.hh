#ifndef BASICREACTIONS_HH_
#define BASICREACTIONS_HH_

#include <cstdint>

#include "../NFreactions.hh"




using namespace std;

namespace NFcore
{

	/* Helpers for counting a pure context reactant the way BioNetGen does: one
	 * reaction instance per matching complex rather than one per matching
	 * molecule.  See ReactionClass::contextCountsPerComplex. */

	/*! Number of distinct complexes represented among a reactant list's matches. */
	int countDistinctComplexes(ReactantList *rl);

	/*! Same, for a DOR reactant's tree. */
	int countDistinctComplexes(ReactantTree *tree);

	/*! Sum of one representative rate factor per distinct complex in a DOR tree.
	    Exact: it sums the same terms BNG would count instances for. */
	double perComplexRateFactorSum(ReactantTree *tree);

	/*! Fills `out` with the mapping sets to enumerate for one reactant of a
	    RuleMonkey-exact pair count.  Ordinarily that is every live mapping set;
	    for a pure context reactant it is one representative per complex, so that
	    the enumerated pairs are counted on the same footing as
	    getCorrectedReactantCount().  `flatIndices`, when non-null, receives each
	    kept mapping set's flat array position, which a DOR caller needs to look
	    up the matching rate factor. */
	void collectReactantRepresentatives(ReactantList *rl, bool perComplex,
	                                    std::vector<MappingSet*> &out,
	                                    std::vector<int> *flatIndices = 0);

	class BasicRxnClass : public ReactionClass {
		public:
			BasicRxnClass(string name, double baseRate, string baseRateName, TransformationSet *transformationSet, System *s);
			virtual ~BasicRxnClass();

			virtual void init();
			virtual void prepareForSimulation();

			//JJT: checks if there's an existing mapping set in *m equal to *ms that maps to this reaction
			virtual int checkForEquality(Molecule *m, MappingSet* ms,int rxnIndex, ReactantList*);

			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
			virtual void remove(Molecule *m, unsigned int reactantPos);
			virtual double update_a();
			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
			virtual int getReactantCount(unsigned int reactantIndex) const;
			virtual int getCorrectedReactantCount(unsigned int reactantIndex) const;
			virtual ReactantContainer *getReactantContainer(unsigned int r) const { return reactantLists[r]; }

			virtual void printFullDetails() const;

				protected:
					virtual void pickRuleMonkeyMappingSets(double randNumber) const;
					virtual double exactRuleMonkey_a();

		protected:
			virtual void pickMappingSets(double randNumber) const;
			// AS-6/22
			bool connectivityFlag;

			ReactantList **reactantLists;


				// RuleMonkey buffers to avoid heap allocations
				MappingSet **msPairBuffer;
				mutable vector<pair<int, int> > validPairsBuffer;
	};


	class FunctionalRxnClass : public BasicRxnClass {

		public:
			FunctionalRxnClass(string name, GlobalFunction *gf, TransformationSet *transformationSet, System *s);
			FunctionalRxnClass(string name, CompositeFunction *cf, TransformationSet *transformationSet, System *s);

			virtual ~FunctionalRxnClass();

			virtual double update_a();
			virtual double exactRuleMonkey_a();
			virtual void pickRuleMonkeyMappingSets(double randNumber) const { BasicRxnClass::pickRuleMonkeyMappingSets(randNumber); }
			virtual void printDetails() const;

		protected:
			GlobalFunction *gf;
			CompositeFunction *cf;
			vector<int> reactantCountBuffer;
	};

	class MMRxnClass : public BasicRxnClass {

		public:
			MMRxnClass(string name, double kcat, double Km, TransformationSet *transformationSet, System *s);
			virtual ~MMRxnClass();

			virtual double update_a();
				virtual double exactRuleMonkey_a();
				virtual void pickRuleMonkeyMappingSets(double randNumber) const { BasicRxnClass::pickRuleMonkeyMappingSets(randNumber); }
			virtual void printDetails() const;

			// The substrate count enters the law through sFree, so a substrate
			// counted once per match cannot be thinned to once per complex.
			virtual bool canThinToOncePerComplex() const { return false; }

		protected:
			double Km;
			double kcat;
			double sFree;
	};


	class DORRxnClass : public ReactionClass {
		public:
			DORRxnClass(
					string name,
					double baseRate,
					string baseRateName,
					TransformationSet *transformationSet,
					CompositeFunction *function,
					vector <string> &lfArgumentPointerNameList,
					System *s);
			/* Internal constructor for reactions whose rate factor is supplied by
			 * another mapping-local evaluator rather than a BNGL local function. */
			DORRxnClass(
					string name,
					double baseRate,
					string baseRateName,
					TransformationSet *transformationSet,
					int dorReactantIndex,
					System *s,
					unsigned int reactantListInitialCapacity = 25,
					unsigned int reactantTreeInitialCapacity = 32,
					bool allocateReactantLists = true);
			virtual ~DORRxnClass();

			virtual void init();
			virtual void prepareForSimulation() {};
			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
			virtual void remove(Molecule *m, unsigned int reactantPos);
			virtual double update_a();

			virtual int getDORreactantPosition() const { return DORreactantIndex; };

			//JJT: checks if there's an existing mapping set in *m equal to *ms that maps to this reaction
			virtual int checkForCollision(Molecule *m, MappingSet* ms,int rxnIndex);

			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
			virtual int getReactantCount(unsigned int reactantIndex) const;
			virtual int getCorrectedReactantCount(unsigned int reactantIndex) const;
			virtual ReactantContainer *getReactantContainer(unsigned int r) const {
				return r == (unsigned)DORreactantIndex
					? static_cast<ReactantContainer *>(reactantTree)
					: static_cast<ReactantContainer *>(reactantLists[r]); }

			virtual void printDetails() const;
			virtual void printFullDetails() const {};

			void directAddForDebugging(Molecule *m);
			void printTreeForDebugging();

			static void test1(System *s);

				protected:
					virtual void pickRuleMonkeyMappingSets(double randNumber) const;
					virtual double exactRuleMonkey_a();

		protected:
			MappingSet **msPairBuffer;
			mutable vector<pair<int, int> > validPairsBuffer;

			virtual double evaluateLocalFunctions(MappingSet *ms);

			virtual void pickMappingSets(double randNumber) const;

			virtual double pickLocalFunctionParameter(MappingSet *ms, int, MoleculeType **, int, int*);

			mutable vector<double> validWeightsBuffer;

			ReactantList **reactantLists;
			ReactantTree *reactantTree;



			CompositeFunction *cf;

			//Parameters to keep track of local functions
			int DORreactantIndex;

			int n_argMolecules;
			int * argIndexIntoMappingSet;
			Molecule ** argMappedMolecule;
			int * argScope;

	};

	/*
	 * Compact Arrhenius energy reaction.  This uses the DOR mapping tree to
	 * select a reaction-center molecule with its context-dependent rate factor,
	 * but does not create one BasicRxnClass for every boolean context state.
	 * The input path only selects this class for factorized binding contexts
	 * whose conditional terms depend only on sites of the same reactant.
	 */
	class EnergyRxnClass : public DORRxnClass {
		public:
			EnergyRxnClass(
					string name,
					double baseRate,
					string baseRateName,
					TransformationSet *transformationSet,
					int dorReactantIndex,
					const EnergyBindingContext &context,
					double phi,
					double RT,
					bool isForward,
					System *s);
			virtual ~EnergyRxnClass();
			virtual bool usesIncrementalMembership() const { return simpleMembership; }
			virtual bool supportsSparseSelection() const {
				return simpleMembership;
			}
			virtual bool membershipDecisionIsTypeInvariant() const {
				return simpleMembership;
			}
			virtual bool getIncrementalMembershipChange(
					IncrementalMembershipChange &change) const;
			virtual bool getCompactMembershipIndexInfo(
					unsigned int reactantPos,
					int &reactionCenterComponent,
					std::uint64_t &contextComponentMask,
					unsigned int &minimumContextComponents) const;
			virtual bool getCompactPartnerPoolInfo(
					unsigned int reactantPos,
					int &partnerComponent) const;
			virtual bool refreshCompactPartnerPool(
					Molecule *m, unsigned int reactantPos);
			virtual bool supportsCompactPartnerPoolUpdate() const {
				return compactForwardPartnerPropensity;
			}
			virtual bool supportsCompactPartnerPoolScale() const {
				return compactForwardPartnerPropensity && !useRuleMonkey;
			}
			virtual CompactPartnerPool *getCompactPartnerPool() const {
				return partnerPool;
			}
			virtual double getCompactPartnerPoolCoefficient() const;
			virtual double update_a_for_compact_partner_pool(int poolSize);
			virtual bool shouldUpdateMembershipForChange(
					Molecule *m,
					const IncrementalMembershipChange &change) const;
			virtual bool supportsDeferredMembershipUpdate() const {
				return simpleMembership;
			}
			virtual double get_a() const;
			virtual bool tryToAddAndReportChange(
					Molecule *m, unsigned int reactantPos);
			virtual bool tryToAddWithIndex(
				Molecule *m, unsigned int reactantPos, int rxnIndex);
			virtual bool tryToAddAndReportChangeWithIndex(
				Molecule *m, unsigned int reactantPos, int rxnIndex);
			virtual void remove(Molecule *m, unsigned int reactantPos);
			virtual double update_a();
			virtual int getReactantCount(unsigned int reactantIndex) const;
			virtual int getCorrectedReactantCount(unsigned int reactantIndex) const;
			virtual bool canUseDirectProductList() const;
			virtual bool canSkipIndirectMembership(
					ReactionClass *firedReaction) const;
			virtual bool checkPreFireConditions(
					MappingSet **mappingSets) const;
			virtual void notifyRateFactorChange(
					Molecule *m, int reactantIndex, int rxnListIndex);
			virtual bool shouldUpdateMembership(Molecule *m,
					ReactionClass *firedReaction,
					bool directProduct) const;

		protected:
			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
			virtual double evaluateLocalFunctions(MappingSet *ms);
			virtual void pickRuleMonkeyMappingSets(double randNumber) const;
			virtual double exactRuleMonkey_a();

		private:
			vector<EnergyPatternTerm> conditionalTerms;
			vector<int> conditionComponentIndices;
			vector<std::uint64_t> conditionalComponentMasks;
			bool componentMaskFastPath;
			double baseEnergy;
			double phi;
			double RT;
			bool isForward;
			bool simpleMembership;
			bool compactFactorizedPropensity;
			bool compactForwardPartnerPropensity;
			bool compactReversePropensity;
			bool preFireBindingFastPath;
			int reactionCenterComponentIndex;
			int partnerComponentIndex;
			MoleculeType *partnerMoleculeType;
			CompactPartnerPool *partnerPool;
			MappingSet *compactPartnerMappingSet;
			std::uint64_t weightedDependencyMask;
			bool dependencyMaskValid;
			bool singleConditionalTermFastPath;
			double baseEnergyRateFactor;
			double conditionedEnergyRateFactor;
			bool multiConditionalTermFastPath;
			std::vector<double> conditionalRateFactors;
			double compactRateFactor;
			unsigned int minimumConditionalBits;
			mutable bool directProductListDecisionKnown;
			mutable bool directProductListSafe;

			bool tryToAddCompact(Molecule *m, unsigned int reactantPos,
					int rxnIndex = -1);
			void refreshCompactRateFactor();

			bool dependsOnEndpoint(MoleculeType *targetMoleculeType,
					MoleculeType *changedMoleculeType,
					int changedComponentIndex) const;

			virtual void pickMappingSets(double randNumber) const;
	};

	/* A reaction class with DOR calculations on two reactants.
	 * The rate function must be factored as  f(x,y) = g(x)*h(y)
	 * */
	class DOR2RxnClass : public ReactionClass {
		public:
			DOR2RxnClass(
					string name,
					double baseRate,
					string baseRateName,
					TransformationSet *transformationSet,
					CompositeFunction *function1,
					CompositeFunction *function2,
					vector <string> &lfArgumentPointerNameList1,
					vector <string> &lfArgumentPointerNameList2,
					System *s);
			virtual ~DOR2RxnClass();

			virtual void init();
			virtual void prepareForSimulation() {};
			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
			virtual void remove(Molecule *m, unsigned int reactantPos);
			virtual double update_a();

			virtual int getDORreactantPosition()  const { return DORreactantIndex1; };
			virtual int getDORreactantPosition2() const { return DORreactantIndex2; };

			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
			virtual int getReactantCount(unsigned int reactantIndex) const;
			virtual int getCorrectedReactantCount(unsigned int reactantIndex) const;
			virtual ReactantContainer *getReactantContainer(unsigned int r) const {
				if (r == (unsigned)DORreactantIndex1) return reactantTree1;
				if (r == (unsigned)DORreactantIndex2) return reactantTree2;
				return reactantLists[r]; }

			virtual void printDetails() const;
			virtual void printFullDetails() const {};

			void directAddForDebugging(Molecule *m);
			void printTreeForDebugging();

			static void test1(System *s);

			public:
				virtual void pickRuleMonkeyMappingSets(double randNumber) const;
				virtual double exactRuleMonkey_a();

		protected:
			MappingSet **msPairBuffer;
			mutable vector<pair<int, int> > validPairsBuffer;

			virtual double evaluateLocalFunctions1(MappingSet *ms);
			virtual double evaluateLocalFunctions2(MappingSet *ms);

			virtual void pickMappingSets(double randNumber) const;

			mutable vector<double> validWeightsBuffer;

			ReactantList **reactantLists;
			ReactantTree *reactantTree1;
			ReactantTree *reactantTree2;


			CompositeFunction *cf1;
			CompositeFunction *cf2;

			//Parameters to keep track of local functions
			int DORreactantIndex1;
			int DORreactantIndex2;

			int n_argMolecules1;
			int n_argMolecules2;
			int * argIndexIntoMappingSet1;
			int * argIndexIntoMappingSet2;
			Molecule ** argMappedMolecule1;
			Molecule ** argMappedMolecule2;
			int * argScope1;
			int * argScope2;

	};

}







#endif /*BASICREACTIONS_HH_*/
