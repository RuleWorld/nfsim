/*
 * reactionSelector.hh
 *
 *  Created on: Jul 23, 2009
 *      Author: msneddon
 */

#ifndef REACTIONSELECTOR_HH_
#define REACTIONSELECTOR_HH_



#include "../NFcore.hh"
//#include <vector>


using namespace std;


namespace NFcore
{
	//Forward Declarations
	class ReactionClass;
	class CompactPartnerPool;


	// Abstract Interface Class for the Reaction Selection Algorithm
	class ReactionSelector {

		public:

			//Initializations and basic functionality
			ReactionSelector() : sys_(nullptr) {};
			explicit ReactionSelector(System *sys) : sys_(sys) {};
			virtual ~ReactionSelector() {};

			virtual double refactorPropensities() = 0;


			virtual double update(ReactionClass *r,double oldA, double newA) = 0;
			virtual double updateBatch(vector<ReactionClass *> &rxns);
			virtual double updateCompactPartnerPoolBatch(
					const vector<ReactionClass *> &rxns,
					int oldPoolSize, int newPoolSize,
					unsigned long long deferredGeneration);
			virtual double getNextReactionClass(ReactionClass *&rc) = 0;
			virtual double getAtot() = 0;
			// The total to draw a waiting time from. A selector that keeps a
			// running total re-sums it here when rounding may have left it
			// inaccurate (issue #123).
			virtual double getCheckedAtot() { return getAtot(); }

		protected:
			System *sys_;  // Pointer to System for accessing per-instance RNG
	};


	class DirectSelector : public ReactionSelector {

		public:
			//Initializations and basic functionality
			DirectSelector(vector <ReactionClass *> &rxns, System *sys);
			virtual ~DirectSelector();

			virtual double refactorPropensities();


			virtual double update(ReactionClass *r,double oldA, double newA);
			virtual double updateBatch(vector<ReactionClass *> &rxns);
			virtual double updateCompactPartnerPoolBatch(
					const vector<ReactionClass *> &rxns,
					int oldPoolSize, int newPoolSize,
					unsigned long long deferredGeneration);
			virtual double getNextReactionClass(ReactionClass *&rc);
			virtual double getAtot();
			virtual double getCheckedAtot();


		protected:
			struct CompactPoolSelectionGroup {
				CompactPoolSelectionGroup() : pool(0), poolSize(0),
						totalCoefficient(0.0), blockCoefficients(), reactionIndices() {}
				CompactPartnerPool *pool;
				int poolSize;
				double totalCoefficient;
				vector<double> blockCoefficients;
				vector<int> reactionIndices;
			};

			int findCompactPoolGroup(CompactPartnerPool *pool) const;
			void applyCompactPoolGroupSize(int groupIndex, int newPoolSize);
			void synchronizeCompactPoolGroup(int groupIndex);
			void updateCompactPoolActiveBits(
					const CompactPoolSelectionGroup &group,
					bool active);
			double getSelectionPropensity(std::size_t reaction) const;
			double resumAtotFromCache();

			double Atot;
			// Largest value of Atot since it was last summed afresh.
			double peakAtot;
			// Running-sum updates to Atot since it was last summed afresh.
			unsigned int updatesSinceSum;
			int n_reactions;
			ReactionClass ** reactionClassList;
			bool sparseSelectionSafe;
			/* -1 until the first update; 1 when prepared reaction IDs match
			 * selector positions, 0 for direct-selector compatibility callers. */
			int reactionIndexMode;
			vector<std::uint64_t> activeReactionBits;
			unsigned int selectionBlockSize;
			vector<double> selectionBlockPropensities;
			/* Keep a contiguous actual-propensity mirror for non-group entries and
			 * diagnostics.  Grouped compact selection derives each value from its
			 * stored coefficient and shared pool size; dense legacy selectors retain
			 * the reaction-object read. */
			vector<double> reactionPropensities;
			/* For factorized compact forward rules, retain the weighted-side
			 * coefficient while the shared partner-pool size is selector state. */
			vector<int> compactPoolGroupByReaction;
			vector<double> compactPoolCoefficients;
			vector<CompactPoolSelectionGroup> compactPoolSelectionGroups;

		};


	class LogClassSelector : public ReactionSelector {

		public:
			//Initializations and basic functionality
			LogClassSelector(vector <ReactionClass *> &rxns, System *sys);
			virtual ~LogClassSelector();

			virtual double refactorPropensities();


			virtual double update(ReactionClass *r,double oldA, double newA);
			virtual double getNextReactionClass(ReactionClass *&rc);
			virtual double getAtot();


		protected:

			int calculateClass(double a);

			void place(ReactionClass *r,int logClass,double a);

			void setLogClassToActive(int logClass);
			void setLogClassToInactive(int logClass);


			// The absolute maximum log class
			int maxClassLimit;
			int minClassLimit;
			int totalLogClassCount;

			//the current range of log classes that we have to loop over
			int currentHighClass;
			int currentLowClass;


			// A 2d array of the logClasses
			ReactionClass *** logClassList;

			// A 1d array of the logClasses sizes
			int *logClassSize;
			int *logClassCapacity;


			//An array to keep track of the log classes with positive values
			int *activeLogClasses;
			bool *isLogClassActive;
			int n_activeLogClasses;



			// A 1d array of the propensity sum of the logClass
			double *logClassPropensity;



			int *mapRxnIdToLogClass;
			int *mapRxnIdToLogClassPosition;




			double Atot;
			int n_reactions;
			ReactionClass ** reactionClassList;

	};



}








#endif /* REACTIONSELECTOR_HH_ */
