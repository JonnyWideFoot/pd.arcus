#ifndef __FFEXAMPLE_H
#define __FFEXAMPLE_H

#include "forcefields/forcefield.h" // Provides the base class

class PD_API Particle;






//-------------------------------------------------
//
/// \brief  Example forcefield & Template for new developments 
///
/// \details  This class simply exempifies the layout of a general forcefield component. 
///           THe functions are annotated to facilitate development of new forcefields.
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo 
///
/// \bug 
///
	class FF_Example: public Physics::ForcefieldBase
	{
	 public:
		
		/// Constructor does not need any arguments, but 
		/// you can add your own if you like.
		FF_Example();

		/// This is a function you must implement - it will be called
		/// automatically (by the ForcefieldContainer this forcefield component
		/// is added to) once before anyone tries to use this forcefield
		/// It should perform any essential setup functionality such as:
		///   - Reading in forcefield parameters
		///   - Setting up internal data strucutres (precalculated data,
		///       paramters, libraries, etc..)
		///   - Request a minimum cutoff from the neighbour list 
		/// Important: This function must be able to be run multiple times, i.e.
		///  you must clear any data structures and recreate them rather than
		///  adding to them.
		/// When this function finished
		virtual void setup();


		/// calcualtes energies AND forces
		virtual void calcForces(){}

		/// calcualtes energies only
		virtual void calcEnergies();   

		/// calculates and displays energies verbosely
		virtual void calcEnergiesVerbose(AtomicVerbosity level); 
	
	 protected:
		/// Put private data here such as internal stores of current step number, state,
		/// Library of structures, etc...
		/// 
		/// e.g.:

		/// Put private member function such as helpers and bits of algorithms
		/// etc.. here. You can use virtual functions to allow other people
		/// to modify your algorithm by overloading specific functional aspects 
		/// of your algorithm. Note that calling virtual functions does carry 
		/// a small overhead.

		/// prints a line of current energies/information/stepnumber etc..
		virtual void infoLine() const;       
		
		/// prints the headers for the above function
		virtual void infoLineHeader() const; 

	};

#endif
