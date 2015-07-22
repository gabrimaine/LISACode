/*
 *  LISACODE-GWFile.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 23/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-GWFile.h"


// *********************
// ***  Constructor  ***
// *********************

LCGWFile::LCGWFile()
: LCGW()
{
	initNULL(false);
}


LCGWFile::LCGWFile(LCTools * MT_n)
: LCGW(MT_n)
{
	initNULL(false);
}


LCGWFile::~LCGWFile()
{
	initNULL(true);
}




// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCGWFile::initNULL(bool CleanMem)
{
	initNULLBase(CleanMem);
	if(CleanMem){
		if(fDat != NULL)
			delete fDat;
	}
	
	fDat = NULL;
	InterpType = LAG;
	InterpUtilValue = 3;
}


// ***************************
// *  Configuration methods  *
// ***************************


void LCGWFile::config(ezxml_t xmlbloc)
{
	//! *** Delete previous allocation
	if(fDat != NULL)
		delete fDat;
	
	//! *** Read sky position and name
	configBase(xmlbloc);
	
	ezxml_t section;
	bool NoDatFound(true);
	for (section = ezxml_child(xmlbloc, "XSIL"); section; section = section->next) {
		if((NoDatFound)&&(MT->wcmp(ezxml_attr(section,"Type"),"TimeSeries"))){
			NoDatFound = false;
			fDat = new LCDataFileRead(MT);
			fDat->config(section);
		}
	}
}


void LCGWFile::config(int iParam, double ParamVal)
{
	
}


// ********************
// *  Access methods  *
// ********************


double LCGWFile::getParam(int iP)
{
	switch (iP) {
		case 0:
			return(Beta);
			break;
		case 1:
			return(Lambda);
			break;
		default:
			Cout << "WARNING : LCGWFile::getParam : The gravitaional wave is readed in a file so the parameters can not be modified (except the sky position : [0,1]) !";
			break;
	}
	return(0.);
}


void LCGWFile::setParam(int iP, double Param_n)
{
	switch (iP) {
		case 0:
			Beta = Param_n;
			break;
		case 1:
			Lambda = Param_n;
			break;
		default:
			Cout << "WARNING : LCGWFile::setParam : The gravitaional wave is readed in a file so the parameters can not be modified (except the sky position : [0,1]) !";
			break;
	}
	
}


void LCGWFile::getRange(int iP, double &Pmin, double &Pmax)
{
	switch (iP) {
		case 0:
			Pmin = -M_PI/2.0;
			Pmax = M_PI/2.0;
			break;
		case 1:
			Pmin = 0.0;
			Pmax = 2.0*M_PI;
			break;
		default:
			Cout << "WARNING : LCGWFile::getRange : The gravitaional wave is readed in a file so the range of parameters is unknown !";
			throw std::invalid_argument("ERROR in LCGWStochastic::getRange : Bad index of parameter.");
	}
}


double LCGWFile::getDelta(int iP)
{
	Cout << "WARNING : LCGWFile::getDelta : The gravitaional wave is readed in a file so the delta of parameters is unknown !";
	return(1.);
}


void LCGWFile::setSpecialParam(int iPS, double SpecParam_n)
{
	
}


double LCGWFile::getSpecialParam(int iPS)
{
	return(0.);
}


void LCGWFile::setTimeInfo(double t0_n, double dt_n, double TObs_n, double tAskMin_n, double tDOrbMax_n,  double tMaxDiff)
{
    t0Real = t0_n; 
	tAskMin = tAskMin_n;
	tDOrbMax = tDOrbMax_n;
}


void LCGWFile::RandParam(int iP)
{
	
}

// ***************************************
// * Linking and initialization methods  *
// ***************************************

int LCGWFile::init()
{
	initBase();
#ifdef _DEBUG_GW_    
	char fNCheck[512];
    int iRCheck(MT->ifloor(MT->RandUniform(0, 10000.)));
    sprintf(fNCheck,"CheckGWFile_%d.txt",iRCheck);
    std::cerr << "DEBUG:GW:LCGWFile : File = " << fNCheck << Endl;
    DEBUGfCheck = new std::ofstream(fNCheck);
#endif
	if(fDat==NULL)
		throw std::invalid_argument("ERROR in LCGWFile::init : There is no file to read !");
	fDat->init();
	return 0;
}


// *********************
// *  Running methods  *
// ********************* 


void LCGWFile::Computehpc(double t)
{
    double tR(t - t0Real + fDat->getx0());
    //Cout << t << "in ? [ " << fDat->getx0()+InterpUtilValue*fDat->getdx() << " = " << fDat->getx0() << " + " << InterpUtilValue*fDat->getdx() <<   "  :  "  << fDat->getxend()-InterpUtilValue*fDat->getdx() << " = " << fDat->getxend() << " - " << InterpUtilValue*fDat->getdx() << " ]" << Endl;  
    if ( (tR<fDat->getxend()-InterpUtilValue*fDat->getdx()) && (tR>fDat->getx0()+InterpUtilValue*fDat->getdx()) ) {
        hBpLast = fDat->gData(0, tR, InterpType, InterpUtilValue);
        hBcLast = fDat->gData(1, tR, InterpType, InterpUtilValue);
    }else{
        hBpLast = 0.;
        hBcLast = 0.;
    }
}




// *******************
// *  Other methods  *
// *******************

void LCGWFile::DispInfo(char * BTab)
{
	if(MT->Disp()){
		Cout << "Gravitaional wave read in file : " << Name << Endl;
		DispInfoBase(BTab);
		if(fDat != NULL)
			fDat->ControlDisplay();
	}
}

void LCGWFile::DispAllParam(std::ostream * out)
{
	for(int iP=0; iP<NParams; iP++)
		(*out) << " " << getParam(iP);
}


void LCGWFile::DispAllParamName(std::ostream * out)
{
	(*out) << "Beta Lambda Polarization";
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




// ***************************************
// * Linking and initialization methods  *
// ***************************************



// *********************
// *  Running methods  *
// *********************




// *******************
// *  Other methods  *
// *******************



// end of LISACODE-GWFile.cpp