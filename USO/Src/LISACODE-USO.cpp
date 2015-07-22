/*
 *  LISACODE-USO.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 12/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-USO.h"

// *********************
// ***  Constructor  ***
// *********************

LCUSO::LCUSO()
{
	MT = new LCTools;
	MT->LocTools = true;
	initNULLBase(false);	
}


LCUSO::LCUSO(LCTools * MT_n)
{
	MT = MT_n;
	initNULLBase(false);
}



LCUSO::~LCUSO()
{
	initNULLBase(true);
}


void LCUSO::initNULLBase(bool CleanMem)
{
	if(CleanMem){
		if(tShift != NULL)
			delete tShift;
	}
	
	iSC = -1;
	tShift = NULL;
	
	t0tShift = 0.;
	tShift2Phase = 2.*M_PI*(LC::c_SI/1064.e-9);
	InterpType = LAG;
	InterpUtilValue = 7;
	
	N2Add = 1;
	
	IsFreqOutput = true;
	
	dtMes = 1.;
	dtPhy = 1.;
	tFirst = 0.;
	tLast = 30.;
	
}


// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************


void LCUSO::initNULL(bool CleanMem)
{
	
}


void LCUSO::config(ezxml_t noisexmlbloc)
{
	
}


void LCUSO::config(int iParam, double ParamVal)
{
	
}


void LCUSO::init()
{
	
}


void LCUSO::genertShift(int iStartBin)
{
	for(int i=iStartBin; i>=0; i--)
		tShift->addData(0.);
}


void LCUSO::DispInfo(char * BTab)
{
	
}


// ***********************
// ***  Local mehtods  ***
// ***********************


// ***************************************
// * Linking and initialization methods  *
// ***************************************


void LCUSO::initBase()
{
	double t0;
	int NDat;
	
	if(tShift != NULL)
		delete tShift;
	
	//! **** Creation of the serie
	t0 = tFirst - InterpUtilValue*dtPhy;
	NDat = MT->iceil((tLast-t0)/dtPhy);
	tShift = new LCSerie2(MT, t0, dtPhy, NDat);
	
	N2Add = MT->ifloor(dtMes/dtPhy);
}


// ********************
// *  Access methods  *
// ********************

void LCUSO::setTimeInfo(double dtMes_n, double dtPhy_n, double tFirst_n, double tLast_n)
{
	dtMes = dtMes_n;
	dtPhy = dtPhy_n;
	tFirst = tFirst_n;
	tLast = tLast_n;
}


void LCUSO::RefFrequency(double RefFreq)
{	
	//! *** Conversion between time shift and phase noise  : \f$ \Delta \phi = 2 \pi f \ \Delta t \f$
	tShift2Phase = 2.*M_PI*RefFreq;
}


// *********************
// *  Running methods  *
// *********************

void LCUSO::RunStep(double t)
{
	t0tShift = t;
	genertShift(N2Add-1);
}


double LCUSO::gT(double t)
{
	return(t + tShift->gData(t0tShift-t, InterpType, InterpUtilValue));
}


double LCUSO::gN(double tDelay)
{
	if(IsFreqOutput){
		//! *** Conversion between time shift and frequency noise : \f$  {\Delta \nu \over \nu } = {1 \over 2 \pi \nu } { d \Delta \phi / dt} = { d \Delta t / dt} \f$ 
		return( tShift->DerivBackOrder2Spe(tDelay, tShift->getRefStep(), InterpType, InterpUtilValue) );
	}else{
		//! *** Conversion between time shift and phase noise  : \f$ \Delta \phi = 2 \pi \nu \ \Delta t \f$
		return(tShift2Phase * tShift->gData(tDelay, InterpType, InterpUtilValue));
	}
}




// end of LISACODE-USO.cpp