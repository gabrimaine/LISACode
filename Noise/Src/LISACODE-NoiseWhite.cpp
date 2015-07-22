/*
 *  LISACODE-NoiseWhite.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 09/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


#include "LISACODE-NoiseWhite.h"


// *********************
// ***  Constructor  ***
// *********************

/**
 *
 */
LCNoiseWhite::LCNoiseWhite()
: LCNoise()
{
	initNULL(false);
}


LCNoiseWhite::LCNoiseWhite(LCTools * MT_n)
: LCNoise(MT_n)
{
	initNULL(false);
}


LCNoiseWhite::~LCNoiseWhite()
{
	initNULL(true);
}



// ********************************************
// ***        Required methods              ***
// ********************************************

void LCNoiseWhite::initNULL(bool CleanMem)
{
	initNULLBase(CleanMem);
	
	Sigma = 1.;
	
	SqPSD = 1.;
}


void LCNoiseWhite::config(ezxml_t noisexmlbloc)
{
	ezxml_t param;
	char * SpectralType(NULL);
	
	for(param = ezxml_child(noisexmlbloc,"Param"); param; param = param->next){
		if(MT->wcmp(ezxml_attr(param,"Name"),"SpectralType")){
			MT->stripcopy((*param).txt, SpectralType);
			if(!MT->wcmp(SpectralType, "WhiteFrequency"))
				throw std::invalid_argument("ERROR in LCNoiseWhite::config : The xml bloc does not correspond to a spectral type of white noise.");
			
		}
		if(MT->wcmp(ezxml_attr(param,"Name"),"PowerSpectralDensity")){
			double PSD;
			PSD = atof((*param).txt);
			SqPSD = sqrt(PSD);
		}
	}
	if(SpectralType != NULL)
		MT->Free(SpectralType, (strlen(SpectralType)+1) * sizeof(char) );
	
	
}

void LCNoiseWhite::config(int iParam, double ParamVal)
{
	switch (iParam) {
		case 0:
			SqPSD = ParamVal;
			break;
		default:
			Cout << "ERROR in LCNoiseWhite::config : Unknow type of parameter : " << iParam << Endl;
			throw std::invalid_argument("ERROR in LCNoiseWhite::config : The type of parameters is unknown !");
			break;
	}
}


void LCNoiseWhite::init()
{
	initBase();
	
	Sigma = SqPSD/sqrt(2.0*dtMes);
	
	//! *** Fill the noise serie
	generNoise(NoiseData->getNmax()-1);
}


void LCNoiseWhite::generNoise(int iStartBin)
{
	double tmpN;
	for(int i=iStartBin; i>=0; i--){
		tmpN = Sigma * MT->RandGaussian(0., 1.0);
		NoiseData->addData(tmpN);
	}
}



void LCNoiseWhite::DispInfo(char * BTab)
{
	if(MT->Disp()){
		DispInfoBase(BTab);
		Cout << BTab << "\t- sqrt(PSD) = " << SqPSD << Endl;
	}
}	

// ***********************
// ***  Local mehtods  ***
// ***********************


void LCNoiseWhite::setSqPSD(double SqPSD_n)
{
    if(SqPSD_n<0.0)
        throw std::invalid_argument("LCNoiseWhite::setSqPSD : Power Spectral Density cannot be negative !");
	SqPSD = SqPSD_n;
}



// end of LISACODE-Noise.cpp



