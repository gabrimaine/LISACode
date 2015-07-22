/*
 *  LISACODE-USONoiseDeriv.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 12/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-USONoiseDeriv.h"

// *********************
// ***  Constructor  ***
// *********************

LCUSONoiseDeriv::LCUSONoiseDeriv()
: LCUSO()
{
	initNULL(false);
}


LCUSONoiseDeriv::LCUSONoiseDeriv(LCTools * MT_n)
: LCUSO(MT_n)
{
	initNULL(false);	
}


LCUSONoiseDeriv::~LCUSONoiseDeriv()
{
	initNULL(true);
}



// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************


void LCUSONoiseDeriv::initNULL(bool CleanMem)
{
	if(CleanMem){
		if(noise != NULL)
			delete noise;
		
	}
		
	tOffset = 0.;
	tSlope = 0.;
	noise = NULL;
}


void LCUSONoiseDeriv::config(ezxml_t usodata)
{
	ezxml_t param;
	char * USOType(NULL);
	
	//! **** Read the parameters
	
	for(param = ezxml_child(usodata,"Param"); param; param = param->next){
		//Cout << ezxml_attr(param,"Name") << Endl;
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"USOType")){
			MT->stripcopy((*param).txt, USOType);
			if(!MT->wcmp(USOType, "NoiseDeriv"))
				throw std::invalid_argument("ERROR in LCUSONoiseDeriv::config : The USO type of the xml bloc is not a 'NoiseDeriv' ! ");
		}
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"SCIndex")){
			iSC = atoi((*param).txt);
		}
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Offset")){
			tOffset = MT->gXMLTime(param);
		}
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Derivs")){
			tSlope = MT->gXMLTime(param);
		}
	}
	
	if(USOType!=NULL)
		MT->Free(USOType, (strlen(USOType)+1) * sizeof(char));
	
	
	//! **** Read the noise
	ezxml_t noisedata; 
	for (noisedata = ezxml_child(usodata, "XSIL"); noisedata; noisedata = noisedata->next) {
		ezxml_t param;
		char * SpectralType(NULL);
		for(param = ezxml_child(noisedata,"Param"); param; param = param->next){
			if(MT->wcmp(ezxml_attr(param,"Name"),"SpectralType")){
				MT->stripcopy((*param).txt, SpectralType);
				if(MT->wcmp(SpectralType, "White")){
					if(noise != NULL)
						delete noise;
					noise = new LCNoiseWhite(MT);
					noise->config(noisedata);
				}
				if((MT->wcmp(SpectralType, "Filter_f"))
				   ||(MT->wcmp(SpectralType, "Filter_1of"))
				   ||(MT->wcmp(SpectralType, "WhitePhase"))
				   ||(MT->wcmp(SpectralType, "BlueFrequency"))
				   ||(MT->wcmp(SpectralType, "WhiteFrequency"))
				   ||(MT->wcmp(SpectralType, "RedFrequency"))
				   ||(MT->wcmp(SpectralType, "PinkFrequency"))
				   ||(MT->wcmp(SpectralType, "PinkAcceleration"))
				   ||(MT->wcmp(SpectralType, "PinkPhase"))
				   ||(MT->wcmp(SpectralType, "Filter_1of_1of2"))){
					if(noise != NULL)
						delete noise;
					noise = new LCNoiseFilter(MT);
					noise->config(noisedata);
				}
			}
		}
		if(SpectralType != NULL)
			MT->Free(SpectralType, (strlen(SpectralType)+1) * sizeof(char));
	}
	
}


void LCUSONoiseDeriv::config(int iParam, double ParamVal)
{
	switch (iParam) {
		case 0:
			tOffset = ParamVal;
			break;
		case 1:
			tSlope = ParamVal;
			break;
		default:
			throw std::invalid_argument("ERROR in LCUSONoiseDeriv::config : The index of parameter is unknown (only 0 or 1 )! ");
			break;
	}
}


void LCUSONoiseDeriv::init()
{
	initBase();
	noise->setTimeInfo(dtMes, dtPhy, tFirst, tLast);
	noise->init();
	
	//! *** Fill the time shift serie
	genertShift(tShift->getNmax()-1);
	
}



void LCUSONoiseDeriv::genertShift(int iStartBin)
{
	double tmpT, tmpN;
	
	/*! *** Compute one nefile://localhost/Users/petiteau/Applications/src/LISACode/LISACode_2_0/LISACode/OptBenPhaMet/Include/LISACODE-OBPM.hw step of noise */
	noise->RunStep();
	
	/*! *** Compute the time shift values : \f$ \Delta t (i_B) = \Delta t_{offset} + \Delta t_{slope} * (t_0^{t_shift} - (i dt + x_0))  + Noise(i_B) \f$ */
	for(int i=iStartBin; i>=0; i--){
		tmpT = tSlope * ( t0tShift - tShift->getRef(i) );
		tmpN =  noise->gB(i);
		//Cout << tOffset << " + " << tmpT << " + " << tmpN << " = " <<  tOffset + tmpT + tmpN  << Endl;
		tShift->addData( tOffset + tmpT + tmpN/tShift2Phase );
	}
}


void LCUSONoiseDeriv::DispInfo(char * BTab)
{
	if(MT->Disp()){
		Cout << BTab << "USO with time evolution described by first time derivative and noise : spacecraft " << iSC << Endl;
		Cout << BTab << "\t- Offset (on time shift) = " << tOffset << " s" << Endl;
		Cout << BTab << "\t- First time derivative (on time shift) = " << tSlope << " s/s" << Endl;
		char BTab2[1024];
		strcpy(BTab2,BTab);
		strcat(BTab2, "\t");
		noise->DispInfo(BTab2);
	}
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


// ********************
// *  Access methods  *
// ********************


// *********************
// *  Running methods  *
// *********************





// end of LISACODE-USONoiseDeriv.cpp