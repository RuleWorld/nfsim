/*!\file compartment.cpp
    \brief Implementation of Compartment class
*/

#include "compartment.hh"
#include <iostream>

using namespace NFcore;
using namespace std;

Compartment::Compartment(
	string id, 
	int spatialDimensions, 
	double size, 
	Compartment* parent)
	: id(id),
	  spatialDimensions(spatialDimensions),
	  size(size),
	  parent(parent)
{
	// Validate inputs
	if (spatialDimensions != 2 && spatialDimensions != 3) {
		cerr << "Warning: Compartment '" << id << "' has unusual spatial dimensions: " 
		     << spatialDimensions << " (expected 2 or 3)" << endl;
	}
	
	if (size <= 0) {
		cerr << "Warning: Compartment '" << id << "' has non-positive size: " 
		     << size << endl;
	}
}

Compartment::~Compartment()
{
	// Parent is not owned by this compartment, so don't delete it
}

bool Compartment::isInside(Compartment* other) const
{
	if (!other) return false;
	if (this == other) return true;
	
	// Check parent chain
	Compartment* current = parent;
	while (current) {
		if (current == other) return true;
		current = current->getParent();
	}
	
	return false;
}

void Compartment::printDetails() const
{
	cout << "Compartment '" << id << "': ";
	cout << spatialDimensions << "D, ";
	cout << "size=" << size;
	if (parent) {
		cout << ", parent=" << parent->getId();
	}
	cout << endl;
}
