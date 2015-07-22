/*
 *  LISACODE-GWGalBin.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 06/05/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-GWGalBin.h"


// *********************
// ***  Constructor  ***
// *********************

LCGWGalBin::LCGWGalBin()
: LCGW()
{
	initNULL(false);
}


LCGWGalBin::LCGWGalBin(LCTools * MT_n)
: LCGW(MT_n)
{
	initNULL(false);
}


LCGWGalBin::~LCGWGalBin()
{
	initNULL(true);
}




// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCGWGalBin::initNULL(bool CleanMem)
{
	initNULLBase(CleanMem);
	
	NParams = 12;
	
	inc = 0.;
	Amp = 1.e-21;
	Freq = 1.e-3;
	DFreq = 0.;
	Phi0 = 0.;
	
	m1 = -1.;
	m2 = -1.;
	DL = -1.;
	ecc = 0.;
	DPeriod = 0.;
	
	DontChangeAmp = false;
	DontChangeDFreq = false;
	
	hp0 = 0.;
	hc0 = 0.;
	om = 1.e-3;
	Dom = 0.;
	
}


// ***************************
// *  Configuration methods  *
// ***************************


void LCGWGalBin::config(ezxml_t xmlbloc)
{
	//! *** Read sky position, polarization and name
	configBase(xmlbloc);
	
	ezxml_t param;
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Inclination"))
			setParam(3, MT->gXMLAngle(param));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Amplitude"))
			setParam(4, atof((*param).txt));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Frequency"))
			setParam(5, MT->gXMLFrequency(param));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Period")){
			double Period;
			Period = MT->gXMLTime(param);
			setParam(5, 2./Period);
		}
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"FrequencyDerivative"))
			setParam(6, 1.e16*atof((*param).txt));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"InitialPhase"))
			setParam(7, MT->gXMLAngle(param));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Mass1"))
			setParam(8, MT->gXMLAstroMass(param));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Mass2"))
			setParam(9, MT->gXMLAstroMass(param));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Distance"))
			setParam(10, MT->gXMLAstroDistance(param));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Eccentricity"))
			setParam(11, atof((*param).txt));
			
		if(MT->wcmp(ezxml_attr(param,"Name"),"PeriodDerivative"))
			setParam(12, atof((*param).txt));
	}
}


void LCGWGalBin::config(int iParam, double ParamVal)
{
	
}


// ********************
// *  Access methods  *
// ********************


double LCGWGalBin::getParam(int iP)
{
	switch (iP) {
		case 0:
			return(Beta);
			break;
		case 1:
			return(Lambda);
			break;
		case 2:
			return(Polarization);
			break;
		case 3:
			return(inc);
			break;
		case 4:
			return(Amp);
			break;
		case 5:
			return(Freq);
			break;
		case 6:
			return(DFreq);
			break;
		case 7:
			return(Phi0);
			break;
		case 8:
			return(m1);
			break;
		case 9:
			return(m2);
			break;
		case 10:
			return(DL);
			break;
		case 11:
			return(ecc);
			break;
		case 12:
			return(DPeriod);
			break;
		default:
			Cout << "ERROR in LCGWGalBin::getParam : The parameter " << iP << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWGalBin::getParam : The parameter is unknown !");
			break;
	}
	return(0.);
}


void LCGWGalBin::setParam(int iP, double Param_n)
{
	bool RecomputeAmp(false);
	bool RecomputeDFreq(false);
	switch (iP) {
		case 0:
			Beta = Param_n;
			break;
		case 1:
			Lambda = Param_n;
			break;
		case 2:
			Polarization = Param_n;
			break;
		case 3:
			inc = Param_n;
			break;
		case 4:
			Amp = Param_n;
			DontChangeAmp = true;
			break;
		case 5:
			Freq = Param_n;
			RecomputeDFreq = true;
			RecomputeAmp = true;
			break;
		case 6:
			DFreq = Param_n;
			DontChangeDFreq = true;
			break;
		case 7:
			Phi0 = Param_n;
			break;
		case 8:
			m1 = Param_n;
			RecomputeAmp = true;
			break;
		case 9:
			m2 = Param_n;
			RecomputeAmp = true;
			break;
		case 10:
			DL = Param_n;
			RecomputeAmp = true;
			break;
		case 11:
			ecc = Param_n;
			break;
		case 12:
			DPeriod = Param_n;
			RecomputeDFreq = true;
			break;
		default:
			Cout << "ERROR in LCGWGalBin::setParam : The parameter " << iP << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWGalBin::setParam : The parameter is unknown !");
			break;
	}
	
	//! *** If needed, recompute amplitude : \f[ A = { 2 (G M_{chirp})^{5/3} \over D c^4 } ( \pi f )^{2/3}  \f]
	if((!DontChangeAmp)&&(RecomputeAmp)){
		double Mtot, mu, eta, Mchirp; 
		Mtot = m1 + m2;
		mu = m1*m2/Mtot;
		eta = mu/Mtot;
		Mchirp = pow(mu, 0.6)*pow(Mtot, 0.4);
		Amp = 2. * pow(LC::G_SI * Mchirp * LC::MS_SI , 5./3.) * pow(M_PI*Freq, 2./3.) / ( DL * LC::kpc_m * pow(LC::c_SI,4.));
	}
	
	//! *** If needed, recompute first time derivative of frequency : \f[ {d f \over d t } = - f^2 {d P \over d t } \f]
	if((!DontChangeDFreq)&&(RecomputeDFreq)){
		DFreq = - 2.*DPeriod*(Freq*Freq);
		DFreq *= 1.e16;
	}
	
}


void LCGWGalBin::getRange(int iP, double &Pmin, double &Pmax)
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
		case 2:
			Pmin = 0.0;
			Pmax = 2.0*M_PI;
			break;
		case 3:
			Pmin = 0.0;
			Pmax = M_PI;
			break;
		case 4:
			Pmin = 1.e-25;
			Pmax = 1.e-10;
			break;
		case 5:
			Pmin = 1.e-5;
			Pmax = 1.;
			break;
		case 6:
			Pmin = -1e-13;
			Pmax = 1e-13;
			break;
		case 7:
			Pmin = 0.0;
			Pmax = 2.0*M_PI;
			break;
		case 8:
			Pmin = 0.001;
			Pmax = 10.;
			break;
		case 9:
			Pmin = 0.001;
			Pmax = 10.;
			break;
		case 10:
			Pmin = 0.01;
			Pmax = 100.;
			break;
		case 11:
			Pmin = 0.;
			Pmax = 1.;
			break;
		case 12:
			Pmin = -1e-10;
			Pmax = 1e-10;
		default:
			Cout << "ERROR in LCGWGalBin::setParam : The parameter " << iP << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWGalBin::setParam : The parameter is unknown !");
			break;
	}
}


double LCGWGalBin::getDelta(int iP)
{
	double tmpMin, tmpMax, DRes;
	getRange(iP, tmpMin, tmpMax);
	
	switch (iP){
		case 0:
			DRes = 2.e-6;
			break;
		case 1:
			DRes = 2.e-6;
			break;
		case 2:
			DRes = 1.e-5;
			break;
		case 3:
			DRes = 1.e-5;
			break;
		case 4:
			DRes = 1.e-26;
			break;
		case 6:
			DRes = 1.e-4;
			break;
		case 7:
			DRes = 1.e-4;
			break;
		default:
			DRes = (tmpMax-tmpMin)*5.0e-8;
	};
	
	//Cout << iP << " --> Delta = " << DRes<< Endl;
	return (DRes);
}


void LCGWGalBin::setSpecialParam(int iPS, double SpecParam_n)
{
	
}


double LCGWGalBin::getSpecialParam(int iPS)
{
	return(0.);
}


void LCGWGalBin::setTimeInfo(double t0_n, double dt_n, double TObs_n, double tAskMin_n, double tDOrbMax_n,  double tMaxDiff)
{	
    t0Real = t0_n; 
	tAskMin = tAskMin_n;
	tDOrbMax = tDOrbMax_n;
}


void LCGWGalBin::RandParam(int iP)
{
	double Pmin(0.0);
	double Pmax(1.0);
	
	getRange(iP, Pmin, Pmax);
	setParam(iP,MT->RandUniform(Pmin, Pmax));
	
	//! ** Specifical case : non-uniform
	if(iP == 0)
		setParam(0, M_PI/2.0 - acos(MT->RandUniform(-1.0, 1.0)));
	if(iP == 3)
		setParam(3, acos(MT->RandUniform(-1.0, 1.0)));
	
	if(iP == 4)
		DontChangeAmp = false;
	if(iP == 6)
		DontChangeDFreq = false;
}

// ***************************************
// * Linking and initialization methods  *
// ***************************************

int LCGWGalBin::init()
{
	initBase();
#ifdef _DEBUG_GW_    
	char fNCheck[512];
    int iRCheck(MT->ifloor(MT->RandUniform(0, 10000.)));
    sprintf(fNCheck,"CheckGWGalBin_%d.txt",iRCheck);
    std::cerr << "DEBUG:GW:LCGWGalBin : File = " << fNCheck << Endl;
    DEBUGfCheck = new std::ofstream(fNCheck);
#endif
	
	/*! *** Compute some fixed quantities :
	 * \f[ h{+,0} = A (1+ \cos^2 \iota)  \f]
	 * \f[ h{\times,0} = - 2 A (\cos \iota) \f]
	 * \f[ \omega = 2 \pi f \f]
	 * \f[ \dot{\omega} = \pi \dot{f} \f]
	 */
	
	
	hp0 = Amp * ( 1. + cos(inc)*cos(inc) ) ;
	hc0 = -2. * Amp * cos(inc) ;
	om = 2. * M_PI * Freq;
	Dom = M_PI * DFreq * 1.e-16;
	
	FreqMin = 1.e-6;
	FreqMax = 2.*Freq;
    
    return 0;
}


// *********************
// *  Running methods  *
// *********************


void LCGWGalBin::Computehpc(double t)
{
	/*! *** Compute the component of the waveform (in source frame) :
	 * \f[ h^S_{+}(t) = h{+,0} \cos [ \omega t + \dot{\omega} t^2 ] \f]
	 * \f[ h^S_{\times}(t) = h{\times,0} \sin [ \omega t + \dot{\omega} t^2 ] \f]
	 */
	hBpLast = hp0 * cos( om*t + Dom*t*t + Phi0);
	hBcLast = hc0 * sin( om*t + Dom*t*t + Phi0);
	
	
	/*
	if((t>2.1460e7)&&(t<2.1461e7)){
		if(DEBUGfCheck==NULL){
			DEBUGfCheck = new std::ofstream("CheckData.txt");
			DEBUGfCheck->precision(15);
		}
		(*DEBUGfCheck) << t << " " << hBpLast << " " << cos( om*t + Dom*t*t + Phi0) << " " << om*t + Dom*t*t << " " << (om*t + Dom*t*t + Phi0)  << Endl;
	}
	*/
}




// *******************
// *  Other methods  *
// *******************

void LCGWGalBin::DispInfo(char * BTab)
{
	if(MT->Disp()){
		DispInfoBase(BTab);
	}
	
	Cout << BTab << "\t- Inclination   = " << inc << " rad" << Endl;
	Cout << BTab << "\t- Amplitude     = " << Amp << Endl;
	Cout << BTab << "\t- Frequency     = " << Freq << " Hz" << Endl;
	Cout << BTab << "\t- dFrequency/dt = " << DFreq << " x 10^-16 Hz/s" << Endl;
	Cout << BTab << "\t- Initial phase = " << Phi0 << " rad" << Endl;
	
}


void LCGWGalBin::DispAllParam(std::ostream * out)
{
	for(int iP=0; iP<NParams; iP++)
		(*out) << " " << getParam(iP);
}


void LCGWGalBin::DispAllParamName(std::ostream * out)
{
	(*out) << " Bet Lam Pol Inc Amp Freq DFreq Phi0 m1 m2 DL ecc DPeriod";
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



// end of LISACODE-GWSpinBBH2.cpp