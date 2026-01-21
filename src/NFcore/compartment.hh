/*!\file compartment.hh
    \brief Compartment class for cBNGL spatial models
*/

#ifndef COMPARTMENT_HH_
#define COMPARTMENT_HH_

#include <string>

using namespace std;

namespace NFcore
{
	//! Represents a spatial compartment in cBNGL models
	/*!
	    Compartments define spatial regions where molecules can exist.
	    Each compartment has a unique ID, spatial dimensions (2D or 3D),
	    and a size (area or volume).
	    
	    @author NFsim cBNGL Team
	*/
	class Compartment
	{
		public:
			//! Constructor
			/*!
			    @param id Unique compartment identifier (e.g., "c0", "cytoplasm")
			    @param spatialDimensions 2 for membrane/surface, 3 for volume
			    @param size Volume (3D) or area (2D) of the compartment
			    @param parent Optional parent compartment for nested structures
			*/
			Compartment(
				string id, 
				int spatialDimensions, 
				double size, 
				Compartment* parent = 0
			);
			
			//! Destructor
			~Compartment();
			
			//! Get the compartment ID
			string getId() const { return id; }
			
			//! Get spatial dimensions (2 or 3)
			int getSpatialDimensions() const { return spatialDimensions; }
			
			//! Get the size (volume or area)
			double getSize() const { return size; }
			
			//! Alias for getSize() - returnsibling volume or area
			double getVolume() const { return size; }
			
			//! Get parent compartment (for nested compartments)
			Compartment* getParent() const { return parent; }
			
			//! Set parent compartment
			void setParent(Compartment* p) { parent = p; }
			
			//! Check if this compartment is inside another
			/*!
			    For nested compartments, checks if this is a child of 'other' (directly or transitively)
			    @param other The potential parent compartment
			    @return true if this compartment is inside 'other'
			*/
			bool isInside(Compartment* other) const;
			
			//! Print compartment details for debugging
			void printDetails() const;
			
		private:
			string id;                    //!< Compartment identifier
			int spatialDimensions;        //!< 2 for membrane, 3 for volume
			double size;                  //!< Volume (3D) or area (2D)
			Compartment* parent;          //!< Parent compartment for nesting (optional)
	};
}

#endif /* COMPARTMENT_HH_ */
