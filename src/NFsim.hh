////////////////////////////////////////////////////////////////////////////////
//
//    NFsim: The Network Free Stochastic Simulator
//    A software platform for efficient simulation of biochemical reaction
//    systems with a large or infinite state space.
//
//    Copyright (C) 2016
//    Michael W. Sneddon, James R. Faeder, Thierry Emonet
//
//    Licensed under the MIT License. See LICENSE.txt for details.
//
//
//    For more information on NFsim, see http://emonet.biology.yale.edu/nfsim
//
////////////////////////////////////////////////////////////////////////////////



#ifndef NFSIM_HH_
#define NFSIM_HH_

//Include "mpi.h" in Scheduler.h first
#include "NFscheduler/Scheduler.h"

//Include the core files needed to run the simulation
#include "NFcore/NFcore.hh"
#include "NFutil/NFutil.hh"
#include "NFinput/NFinput.hh"
#include "NFreactions/NFreactions.hh"

//Include the specific tests
#include "NFfunction/NFfunction.hh"
//#include "NFtest/compare/compare.hh"
//#include  "NFtest/transformations/transformations.hh"
#include  "NFtest/simple_system/simple_system.hh"
#include  "NFtest/transcription/transcription.hh"
#include  "NFtest/tlbr/tlbr.hh"
#include  "NFtest/agentcell/agentcell.hh"
#include  "NFtest/util/util_test.hh"



//! Runs a given System with the specified arguments
/*!
  @author Michael Sneddon
*/
bool runFromArgs(System *s, map<string,string> argMap, bool verbose);


//! Initialize a system from command line flags
/*!
  @author Michael Sneddon
*/
System *initSystemFromFlags(map<string,string> argMap, bool verbose);








#endif /*NFSIM_HH_*/
