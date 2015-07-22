/*
 *  LISACODE-Noise.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 05/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


#include "LISACODE-Noise.h"


// *********************
// ***  Constructor  ***
// *********************

/**
 *
 */
LCNoise::LCNoise()
{
	MT = new LCTools();
	MT->LocTools = true;
	initNULLBase(false);
}


LCNoise::LCNoise(LCTools * MT_n)
{
	MT = MT_n;
	initNULLBase(false);
}


LCNoise::~LCNoise()
{
	initNULLBase(true);
}



void LCNoise::initNULLBase(bool CleanMem)
{
	if(CleanMem){
		if(NoiseData != NULL)
			delete NoiseData;
		
	}
	
	NoiseData = NULL;
	InterpType = LAG;
	InterpUtilValue = 7;
	
	Name[0] = ' ';
	Name[1] = ' ';
	Name[2] = ' ';
	Name[3] = ' ';
	
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

void LCNoise::initNULL(bool CleanMem)
{
	
}


void LCNoise::config(ezxml_t noisexmlbloc)
{
	
}


void LCNoise::config(int iParam, double ParamVal)
{
	
}


void LCNoise::init()
{
	
}


void LCNoise::generNoise(int iStartBin)
{
	for(int i=iStartBin; i>=0; i--)
		NoiseData->addData(0.);
}


void LCNoise::DispInfo(char * BTab)
{

}	



// ***********************
// ***  Local mehtods  ***
// ***********************


// ***************************
// *  Configuration methods  *
// ***************************


// ***************************************
// * Linking and initialization methods  *
// ***************************************


void LCNoise::initBase()
{
	double t0;
	int NDat;
	
	if(NoiseData != NULL)
		delete NoiseData;
	
	//! **** Creation of the serie
	t0 = tFirst - InterpUtilValue*dtPhy;
	NDat = MT->iceil((tLast-t0)/dtPhy);
	NoiseData = new LCSerie2(MT, t0, dtPhy, NDat);
	
	N2Add = MT->ifloor(dtMes/dtPhy);
	
}




// ********************
// *  Access methods  *
// ********************


void LCNoise::setTimeInfo(double dtMes_n, double dtPhy_n, double tFirst_n, double tLast_n)
{
	dtMes = dtMes_n;
	dtPhy = dtPhy_n;
	tFirst = tFirst_n;
	tLast = tLast_n;
}


void LCNoise::setName(const char * NamePos)
{
	int ipos(strlen(NamePos)-1);
	int NChar(MIN(strlen(NamePos), 4));
	char TmpiSC[2];
	for(int i=0; i<NChar; i++)
		Name[i] = NamePos[i];
	
	//! *** Start from the end : detect if the direction of optical bench
	if(NamePos[ipos] == 's'){
		IndirectDir = 1;
		ipos--;
	}else{
		if(isdigit(NamePos[ipos])&&(!(isdigit(NamePos[ipos-1])))){
			IndirectDir = 0;
		}else{
			//! ** There are 2 digits at the end of the name , so it's not a LISA like detector (not 3 SCs)
			TmpiSC[0] = NamePos[ipos];
			TmpiSC[1] = '\0';
			IndirectDir = -1*atoi(TmpiSC);
			ipos--;
		}
	}
	
	//! *** Then : Spacecraft index
	TmpiSC[0] = NamePos[ipos];
	TmpiSC[1] = '\0';
	iSC = atoi(TmpiSC);
	//if((iSC<1)||(iSC>3))
	//	throw std::invalid_argument("ERROR in LCNoise::setNamePos : Index of spacecraft should be in 1,2 or 3 !" );
	
	
}

bool LCNoise::CmpName(const char * Name_cmp, int nchar)
{
	bool res(true);
	
	for(int i=0; i<nchar; i++)
		if(Name[i] != Name_cmp[i])
			res = false;
	return res;
}


void LCNoise::setInterpolation(INTERP InterpType_n, double InterpUtilValue_n)
{
	InterpType = InterpType_n;
	InterpUtilValue = InterpUtilValue_n;
}


// *********************
// *  Running methods  *
// *********************

void LCNoise::RunStep()
{
	generNoise(N2Add-1);
}


double LCNoise::gN(double tDelay)
{
	return( NoiseData->gData(tDelay, InterpType, InterpUtilValue) );
	
}


double LCNoise::gB(int iBin)
{
	return( NoiseData->getBinValue(iBin) );
	
}


// *******************
// *  Other methods  *
// *******************

void LCNoise::DispInfoBase(char * BTab)
{
	if(MT->Disp()){
		Cout << BTab << "Noise info : " << Name << " : " << Endl;
		Cout << BTab << "\t- position : spacecraft = " << iSC << "  , direction = ";
		if(IndirectDir)
			Cout << "indirect" << Endl;
		else 
			Cout << "direct" << Endl;
		Cout << BTab << "\t- time : dtMes = " << dtMes << " , dtPhy = " << dtPhy << " , tFisrt=" << tFirst << " , tLast=" << tLast << Endl;
	}
}


// end of LISACODE-Noise.cpp



