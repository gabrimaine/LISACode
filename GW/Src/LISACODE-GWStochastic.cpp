/*
 *  LISACODE-GWStochastic.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 23/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-GWStochastic.h"


// *********************
// ***  Constructor  ***
// *********************

LCGWStochastic::LCGWStochastic()
: LCGW()
{
	initNULL(false);
}


LCGWStochastic::LCGWStochastic(LCTools * MT_n)
: LCGW(MT_n)
{
	initNULL(false);
}


LCGWStochastic::~LCGWStochastic()
{
	initNULL(true);
}




// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCGWStochastic::initNULL(bool CleanMem)
{
	initNULLBase(CleanMem);
	
	if(CleanMem){
		if(nhp!=NULL)
			delete nhp;
		
		if(nhc!=NULL)
			delete nhc;
	}
	
	nhp = NULL;
	nhc = NULL;
	
	Fact_hp = 1.;
	Fact_hc = 1.;
	
	dt = 1.0;
	tLast = 1000.;
	tLastRef = LC::DBLMINALLOW;
	
}


// ***************************
// *  Configuration methods  *
// ***************************


void LCGWStochastic::config(ezxml_t xmlbloc)
{
	
	//! *** Delete previous allocation
	if(nhp!=NULL)
		delete nhp;
	
	if(nhc!=NULL)
		delete nhc;
	
	//! *** Read sky position and name
	configBase(xmlbloc);
	
	//! *** Configuration of the two noises
	ezxml_t param;
	char * SpectralType(NULL);
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		if(MT->wcmp(ezxml_attr(param,"Name"),"SpectralType")){
			MT->stripcopy((*param).txt, SpectralType);
			if(MT->wcmp(SpectralType, "White")){
				nhp = new LCNoiseWhite(MT);
				nhp->config(xmlbloc);
				nhc = new LCNoiseWhite(MT);
				nhc->config(xmlbloc);
			}
			if((MT->wcmp(SpectralType, "Filter_f"))
			   ||(MT->wcmp(SpectralType, "Filter_1of"))
			   ||(MT->wcmp(SpectralType, "WhitePhase"))
			   ||(MT->wcmp(SpectralType, "BlueFrequency"))
			   ||(MT->wcmp(SpectralType, "WhiteFrequency"))
			   ||(MT->wcmp(SpectralType, "RedFrequency"))
			   ||(MT->wcmp(SpectralType, "PinkFrequency"))
			   ||(MT->wcmp(SpectralType, "PinkAcceleration"))
			   ||(MT->wcmp(SpectralType, "Filter_1of_1of2"))){
				nhp = new LCNoiseFilter(MT);
				nhp->config(xmlbloc);
				nhc = new LCNoiseFilter(MT);
				nhc->config(xmlbloc);
			}
		}
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"SpectralSlope")){
			
			//! *** Equivalent to white noise
			if(fabs(atof((*param).txt))<LC::PRECISION){
				nhp = new LCNoiseWhite(MT);
				nhp->config(xmlbloc);
				nhc = new LCNoiseWhite(MT);
				nhc->config(xmlbloc);
			}else{
				nhp = new LCNoiseFilter(MT);
				nhp->config(xmlbloc);
				nhc = new LCNoiseFilter(MT);
				nhc->config(xmlbloc);
				//throw std::invalid_argument("ERROR in LCGWStochastic::config : The noise with non null spectral slope are not yet defined.");
			}
		}
	}
	if(SpectralType != NULL)
		MT->Free(SpectralType, (strlen(SpectralType)+1) * sizeof(char));
}
	
	

void LCGWStochastic::config(int iParam, double ParamVal)
{
	
}


// ********************
// *  Access methods  *
// ********************


double LCGWStochastic::getParam(int iP)
{
	switch (iP) {
		case 0:
			return(Beta);
			break;
		case 1:
			return(Lambda);
			break;
		default:
			Cout << "WARNING : LCGWStochastic::getParam : For a stochastic gravitaional wave, the parameters can not be modified (except the sky position : [0,1]) !";
			break;
	}
	return(0.);
}


void LCGWStochastic::setParam(int iP, double Param_n)
{
	switch (iP) {
		case 0:
			Beta = Param_n;
			break;
		case 1:
			Lambda = Param_n;
			break;
		default:
			Cout << "WARNING : LCGWStochastic::setParam : The gravitaional wave is readed in a file so the parameters can not be modified (except the sky position : [0,1]) !";
			break;
	}
	
}


void LCGWStochastic::getRange(int iP, double &Pmin, double &Pmax)
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


double LCGWStochastic::getDelta(int iP)
{
	//! \todo TODO properly
	return(1.);
}


void LCGWStochastic::setSpecialParam(int iPS, double SpecParam_n)
{
	
}


double LCGWStochastic::getSpecialParam(int iPS)
{
	return(0.);
}


void LCGWStochastic::RandParam(int iP)
{
	
}


void LCGWStochastic::setTimeInfo(double t0_n, double dt_n, double TObs_n, double tAskMin_n, double tDOrbMax_n,  double tMaxDiff)
{
    t0Real = t0_n; 
	tAskMin = tAskMin_n;
	tDOrbMax = tDOrbMax_n;
	
	dt = dt_n;
	tLastRef = t0_n+tAskMin;
}


// ***************************************
// * Linking and initialization methods  *
// ***************************************

int LCGWStochastic::init()
{
	initBase();
#ifdef _DEBUG_GW_    
	char fNCheck[512];
    int iRCheck(MT->ifloor(MT->RandUniform(0, 10000.)));
    sprintf(fNCheck,"CheckGWStochastic_%d.txt",iRCheck);
    std::cerr << "DEBUG:GW:LCGWStochastic : File = " << fNCheck << Endl;
    DEBUGfCheck = new std::ofstream(fNCheck);
#endif
	
	if(nhp==NULL)
		throw std::invalid_argument("ERROR in LCGWStochastic::init : The noise associated to h+ has not been allocated !");
	
	if(nhc==NULL)
		throw std::invalid_argument("ERROR in LCGWStochastic::init : The noise associated to hx has not been allocated !");
	
	
	double tDeep( (tDOrbMax+2.0*MAX(nhp->getInterpUtilValue(),nhc->getInterpUtilValue())*dt) );
	
	if(MT->DispDet())
		Cout << "Deepness = " << tDeep << " s " << Endl;
	
	//! *** Initialization and stabilization of the two noise associated to \f$ h_+ \f$ \f$
	nhp->setTimeInfo(dt, dt, -2.0*nhp->getInterpUtilValue()*dt, tDeep); 
	nhp->init();
	
	//! *** Initialization and stabilization of the noise associated to \f$ h_{\times} \f$
	nhc->setTimeInfo(dt, dt, -nhc->getInterpUtilValue()*dt, tDeep);
	nhc->init();
	
	double t0(tLastRef);
	tLastRef -= tDeep;
	//for(int i=0; i<MT->iceil(tDeep/dt); i++)
	Computehpc(t0);
		
	//std::cerr << "LCGWStochastic::init : t0 = " << t0 << Endl;
	//DEBUGfCheckR.open("Check_hphc_R.txt");
	//DEBUGfCheckT.open("Check_hphc_T.txt");
	
	return 0;
}


// *********************
// *  Running methods  *
// *********************


void LCGWStochastic::Computehpc(double t)
{
	//! *** If needed add data
	while (tLastRef<t) {
		tLastRef += dt;
		nhp->RunStep();
		nhc->RunStep();
		//DEBUGfCheckR << tLastRef << " " << nhp->gN(0.) << " " << nhc->gN(0.) << Endl;
		//Cout << tLastRef << " " << nhp->gN(0.) << " " << nhc->gN(0.) << Endl;
	}
	
	//Cout << t << " " << tLastRef; MT->o->flush();
	hBpLast = nhp->gN(tLastRef-t);
	hBcLast = nhc->gN(tLastRef-t);
	//DEBUGfCheckT << t << " " << hBpLast << " " << hBcLast << Endl;
	//Cout << " " << hBpLast << " " << hBcLast << Endl;
}



// *******************
// *  Other methods  *
// *******************

void LCGWStochastic::DispInfo(char * BTab)
{
	if(MT->Disp()){
		char BTab2[1024];
		strcpy(BTab2,BTab);
		strcat(BTab2, "\t\t");
		Cout << "Stochastic gravitaional wave : " << Name << Endl;
		DispInfoBase(BTab);
		Cout << BTab << "\t- Noise modeling h+ : " << Endl;
		nhp->DispInfo(BTab2);
		Cout << BTab << "\t- Noise modeling hx : " << Endl;
		nhc->DispInfo(BTab2);
	}
	
}

void LCGWStochastic::DispAllParam(std::ostream * out)
{
	for(int iP=0; iP<NParams; iP++)
		(*out) << " " << getParam(iP);
}


void LCGWStochastic::DispAllParamName(std::ostream * out)
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



// end of LISACODE-GWStochastic.cpp