/*
 *  LISACODE-LCTDIInt.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 20/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


#include "LISACODE-TDIInt.h"



// *********************
// ***  Constructor  ***
// *********************


LCTDIInt::LCTDIInt(LCTools * MT_n)
{
	MT = MT_n;
	initNULL(false);
}


LCTDIInt::~LCTDIInt()
{
	initNULL(true);
}



void LCTDIInt::initNULLBase(bool CleanMem)
{
	if(CleanMem){
		
		if(Sig!=NULL)
			MT->Free(Sig, NSig*sizeof(LCSerie2*));
		
		if(Dat!=NULL)
			delete Dat;
		
		if(Delays != NULL)
			MT->Free(Delays, NDelays*sizeof(LCSerie2*));
		
	}
	
	Sig = NULL;
	NSig = 0;
	Dat = NULL;
	Delays = NULL;
	NDelays = 0;
	
	InterpType = LAG;
	InterpUtilValue = 20;
	InterpType = LAG;
	InterpUtilValue = 6;
	
	strcpy(Name, "");
	iSC = -1;
	IndirectDir = -1;
	
	dtMes = 1.0;
	NData = 0;
	tShift = 0.;
	
}



// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCTDIInt::initNULL(bool CleanMem)
{
	initNULLBase(CleanMem);
}


void LCTDIInt::config(ezxml_t xmlbloc)
{
	
}


void LCTDIInt::LinkDetector(LCDetector * LISA)
{
	
}


void LCTDIInt::LinkDelays(LCSerie2** AllDelays, int AllNDelays)
{
	
}


void LCTDIInt::init()
{
	initBase();
}


int LCTDIInt::getNMaxDelay()
{
	return(0);
}


void LCTDIInt::RunStep(double t)
{
	
}


void LCTDIInt::DispInfo(char * BTab)
{
	if(MT->Disp()){
		DispInfoBase(BTab);
	}
}



// ***********************
// ***  Local mehtods  ***
// ***********************


// ***************************
// *  Configuration methods  *
// ***************************

// ********************
// *  Access methods  *
// ********************


int LCTDIInt::getNDatInterp()
{
	return(MT->iceil(InterpUtilValue));
}



void LCTDIInt::setNameEqPos(const char * NameEqPos)
{
	int ipos(strlen(NameEqPos)-1);
	
	//! *** Start from the end : detect if the direction of optical bench
	if(NameEqPos[ipos] == 's'){
		IndirectDir = 1;
		ipos--;
	}else{
		IndirectDir = 0;
	}
	
	//! *** Then : Spacecraft index
	char TmpiSC[2];
	TmpiSC[0] = NameEqPos[ipos];
	TmpiSC[1] = '\0';
	iSC = atoi(TmpiSC);
	if((iSC<1)||(iSC>3))
		throw std::invalid_argument("ERROR in LCNoise::setNamePos : Index of spacecraft should be in 1,2 or 3 !" );
	
	
	int NChar(MIN(strlen(NameEqPos), 8));
	for(int i=0; i<NChar; i++)
		Name[i] = NameEqPos[i];
	Name[NChar] = '\0';
	
}


void LCTDIInt::setTimeInfo(double dtMes_n, int NDataTDII_n, double tShiftTDII_n)
{
	dtMes = dtMes_n;
	NData = NDataTDII_n;
	tShift = tShiftTDII_n;
}


// ***************************************
// * Linking and initialization methods  *
// ***************************************

void LCTDIInt::initBase()
{
	//! *** Delete previous data
	if(Dat != NULL)
		delete Dat;
	
	//! *** Create data
	Dat = new LCSerie2(MT, 0., dtMes, NData);
	
}


// *********************
// *  Running methods  *
// *********************


// *******************
// *  Other methods  *
// *******************

double LCTDIInt::gN(int iN, double tDelay)
{
	if(Sig[iN] == NULL)
		return( 0. );
	else
		return( Sig[iN]->gData(tDelay, InterpType, InterpUtilValue) );
}


void LCTDIInt::DispInfoBase(char * BTab)
{
	if(MT->Disp()){
		Cout << BTab << "Intermediate TDI : " << Name << Endl; 
		Cout << BTab << "\t- position : spacecraft " << iSC << " and direction " << IndirectDir << Endl;
	}
	
}


// end of LISACODE-InterTDI.cpp
