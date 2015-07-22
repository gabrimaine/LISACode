/*
 *  LISACODE-GWSpinBBHHHarm1.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 31/05/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-GWSpinBBHHHarm1.h"


// *********************
// ***  Constructor  ***
// *********************

LCGWSpinBBHHHarm1::LCGWSpinBBHHHarm1()
: LCGW()
{
	initNULL(false);
}


LCGWSpinBBHHHarm1::LCGWSpinBBHHHarm1(LCTools * MT_n)
: LCGW(MT_n)
{
	initNULL(false);
}


LCGWSpinBBHHHarm1::~LCGWSpinBBHHHarm1()
{
	
	initNULL(true);
}




// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCGWSpinBBHHHarm1::initNULL(bool CleanMem)
{
	initNULLBase(CleanMem);
	
	if(CleanMem){
		
		if (fOutCheckPrec != NULL){
			delete fOutCheckPrec;
			Cout << "Close LCCheckPrec.txt " << Endl;
		}
		
		if(coordI != NULL){
			for(int i=0; i<NMemItg; i++)
				MT->Free(coordI[i], 13*sizeof(double) );
			MT->Free(coordI, NMemItg*sizeof(double*) );
		}
		
	}
	
	NParams = 23;
	
	
	mass1 = 1.;
	mass2 = 1.;
	tc = 0.;
	DL = 0.;
	chi1 = 0.;
	chi2 = 0.;
	PolarAngleOfSpin1 = 0.;
	PolarAngleOfSpin2 = 0.;
	AzimuthalAngleOfSpin1 = 0.;
	AzimuthalAngleOfSpin2 = 0.;
	phic = 0.;
	Phi0 = 0.;
	InitialPolarAngleL = 0.;
	InitialAzimuthalAngleL = 0.;
	eta = 0.25;
	Mchirp = 1.;
	mum = 1.;
	Mtot = 1.;
	m1M2 = 1.;
	m2M2 = 1.;
	m1 = 1.;
	m2 = 1.;
	M = 1.;
	mu = 1.;
	dm = 0.;
	dist = 1.;
	
	x1 = 0.;
	x2 = 0.;
	iota = 0.;
	alpha = 0.;
	Lnx = 0.;
	Lny = 0.;
	Lnz = 0.;
	S1x = 0.;
	S1y = 0.;
	S1z = 0.;
	S2x = 0.;
	S2y = 0.;
	S2z = 0.;
	
	LS1 = 0.;
	LS2 = 0.;
	S1S2 = 0.;
    
    TypeExtraParamCompute = 2;
	
	
	taperQ = 3.0;
	tTaper = 0.;
	ApplyTaper = true;
    
    PolNR = 0.;
    ThdNR = 0.;
    PhdNR = 0.;
	
	
	t0Intg = 0.;
	t0IntgtoM = -1.1e30;
	omM0Intg = -1.1e30;
	t0Start = 0.;
	maxDur = 0.;
	back = false;
	dt = 1.;
	tIc = 0.;
	NMemItg = 0;
	coordI = NULL;
	tend = 0.;
	
	Amp = 0.;
	wk = 0.;
	v = 0.;
	theta = 0.;
	thetaS = 0.;
	phiS = 0.;
	cThS = 0.;
	sThS = 0.;
	psi = 0.;
	En_prev = 0.;
	om_prev = 0.;
	om0 =0.;
	om = 0.;
	phi =0.;
	Psi = 0.;
	nonspin = false;
	hp = 0.;
	hc = 0.;
	ex[0] = 0.; ex[1] = 0.; ex[2] = 0.;
	ey[0] = 0.; ey[1] = 0.; ey[2] = 0.;
	ez[0] = 0.; ez[1] = 0.; ez[2] = 0.;
	n[0] = 0.; n[1] = 0.; n[2] = 0.;
	LnBx = 0.;
	LnBy =0.;
	LnBz = 0.;
	PhLnB = 0.;
	ThLnB = 0.;
	PolB = 0.;
	c2PB = 0.;
	s2PB = 0.;
	h0PNp = 0.;
	h0PNc = 0.;
	h05PNp = 0.;
	h05PNc = 0.;
	h1PNp = 0.;
	h1PNc = 0.;
	beta = 0.;
	sigma = 0.;
	tau = 0.;
	tau38 = 0.;
	tau14 = 0.;
	tau58 = 0.;
	x = 0.;
	ci = 0.;
	si = 0.;
	cth = 0.;
	sth = 0.;
	ciby2 = 0.;
	siby2 = 0.;
	s2th = 0.;
	c2th = 0.;
	s3th = 0.;
	c3th = 0.;
	c4th = 0.;
	s4th = 0.;
	
	NRmrat = 0.;
	NReta = 0.;
	NRchi1 = 0.;
	NRchi2 = 0.;
	NRPhi0 = 0.;
	NRFreqMaxM = 0.; 
	NRomMInit = 0.; 
	NRtoMInitHyb = 0.;
	
	AmpL = 0.;
	AmpS1 = 0.;
	AmpS2 = 0.;
	S1B.setTools(MT);
	S1B.xyz(0.,0.,0.);
	S2B.setTools(MT);
	S2B.xyz(0.,0.,0.);
	LnB.setTools(MT);
	LnB.xyz(0.,0.,0.);
	S1N.setTools(MT);
	S1N.xyz(0.,0.,0.);
	S2N.setTools(MT);
	S2N.xyz(0.,0.,0.);
	LnN.setTools(MT);
	LnN.xyz(0.,0.,0.);
	JN.setTools(MT);
	JN.xyz(0.,0.,0.);
	
	if(ReadNRParam("NRdata.txt", NRmrat, NReta, NRchi1, NRchi2, NRPhi0, NRFreqMaxM, NRomMInit, NRtoMInitHyb, S1N, S2N, LnN))
		Cout << "WARNING : the file containing NR parameters ('NRdata.txt')cannot be openned ==> we cannot make conversion between (L, S1, S2) in SSB and (thetad, Phid, Polarization)";
	
	
	fOutCheckPrec = NULL;
	CheckData = false;
	
}


// ***************************
// *  Configuration methods  *
// ***************************


void LCGWSpinBBHHHarm1::config(ezxml_t xmlbloc)
{
	//! *** Read sky position, polarization and name
	configBase(xmlbloc);
	
	//! ** Reset declination to compute quantities related to sky position 
	setParam(0, Beta);
	
	
	ezxml_t param;
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Mass1"))
			setParam(2, MT->gXMLAstroMass(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Mass2"))
			setParam(3, MT->gXMLAstroMass(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"CoalescenceTime"))
			setParam(4, MT->gXMLTime(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Distance"))
			setParam(5, MT->gXMLAstroDistance(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Spin1"))
			setParam(6, atof(ezxml_txt(param)) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Spin2"))
			setParam(7, atof(ezxml_txt(param)) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"PolarAngleOfSpin1"))
			setParam(8, MT->gXMLAngle(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"PolarAngleOfSpin2"))
			setParam(9, MT->gXMLAngle(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"AzimuthalAngleOfSpin1"))
			setParam(10, MT->gXMLAngle(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"AzimuthalAngleOfSpin2"))
			setParam(11, MT->gXMLAngle(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"PhaseAtCoalescence"))
			setParam(18, MT->gXMLAngle(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"InitialPolarAngleL"))
			setParam(13, MT->gXMLAngle(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"InitialAzimuthalAngleL"))
			setParam(14, MT->gXMLAngle(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"taperQFactor"))
			taperQ = atof(ezxml_txt(param));
		
		//! *** Read polarization
		if(MT->wcmp(ezxml_attr(param,"Name"),"Polarization"))
			setParam(20, MT->gXMLAngle(param));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"TotalMass"))
			setParam(19, MT->gXMLAstroMass(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"PolarAngleTotalMomentuminSBB"))
			setParam(21, MT->gXMLAngle(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"AzimuthalAngleTotalMomentuminSBB"))
			setParam(22, MT->gXMLAngle(param) );
		
		
		//if(MT->wcmp(ezxml_attr(param,"Name"),"TaperApplied"))
		//	TaperApplied = atof(ezxml_txt(param));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"InitialPhase"))
			setParam(12, MT->gXMLAngle(param) );
		
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"InitialIntegrationTime"))
			setSpecialParam(0, MT->gXMLTime(param));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"InitialIntegrationTimeM"))
			setSpecialParam(1,  atof(ezxml_txt(param)));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"InitialIntegrationTimetoM"))
			setSpecialParam(2,  atof(ezxml_txt(param)));
		
		
	}
}


void LCGWSpinBBHHHarm1::config(int iParam, double ParamVal)
{
	
}


// ********************
// *  Access methods  *
// ********************


double LCGWSpinBBHHHarm1::getParam(int iP)
{
	switch (iP) {
		case 0:
			return (Beta);
			break;
		case 1:
			return (Lambda);
			break;
		case 2:
			return (mass1);
			break;
		case 3:
			return (mass2);
			break;
		case 4:
			return (tc);
			break;
		case 5:
			return (DL);
			break;
		case 6:
			return (chi1);
			break;
		case 7:
			return (chi2);
			break;
		case 8:
			return (PolarAngleOfSpin1);
			break;
		case 9:
			return (PolarAngleOfSpin2);
			break;
		case 10:
			return (AzimuthalAngleOfSpin1);
			break;
		case 11:
			return (AzimuthalAngleOfSpin2);
			break;
		case 12:
			return (Phi0);
			break;
		case 13:
			return (InitialPolarAngleL);
			break;
		case 14:
			return (InitialAzimuthalAngleL);
			break;
		case 15:
			return (Mchirp);
			break;
		case 16:
			return (mum);
			break;
		case 17:
			return (eta);
			break;
		case 18:
			return (phic);
			break;
		case 19:
			return (Mtot);
			break;
		case 20:
			return (PolNR);
			break;
		case 21:
			return (thetaJ);    // changed by Sofiane
			break;
		case 22:
			return (phiJ);   // changed by Sofiane
			break;
		default:
			Cout << "ERROR in LCGWSpinBBHHHarm1::getParam : The parameter " << iP << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWSpinBBHHHarm1::getParam : The parameter is unknown !");
			break;
	}
	return(0.);
}


void LCGWSpinBBHHHarm1::setParam(int iP, double Param_n)
{
	switch (iP) {
		case 0:
			Beta = Param_n;
			break;
		case 1:
			Lambda = Param_n;
			break;
		case 2:
			mass1 = Param_n;
			Mtot = mass1 + mass2;
			eta = mass1*mass2/(Mtot*Mtot);
			mum = mass1*mass2/Mtot;
			Mchirp = pow(mum, 0.6)*pow(Mtot, 0.4);
			break;
		case 3:
			mass2 = Param_n;
			Mtot = mass1 + mass2;
			eta = mass1*mass2/(Mtot*Mtot);
			mum = mass1*mass2/Mtot;
			Mchirp = pow(mum, 0.6)*pow(Mtot, 0.4);
			break;
		case 4:
			tc = Param_n;
			break;
		case 5:
			DL = Param_n;
			dist = DL*LC::kpc_s;
			break;
		case 6:
			chi1 = Param_n;
			x1 = chi1;
			break;
		case 7:
			chi2 = Param_n;
			x2 = chi2;
			break;
		case 8:
			PolarAngleOfSpin1 = Param_n;
			break;
		case 9:
			PolarAngleOfSpin2 = Param_n;
			break;
		case 10:
			AzimuthalAngleOfSpin1 = Param_n;
			break;
		case 11:
			AzimuthalAngleOfSpin2 = Param_n;
			break;
		case 12:
			Phi0 = Param_n;
			phic = -1.1e30;
			break;
		case 13:
			InitialPolarAngleL = Param_n;
			break;
		case 14:
			InitialAzimuthalAngleL = Param_n;
			break;
		case 15:
			Mchirp = Param_n;
			Mtot = Mchirp / pow(eta,0.6);
			mum = eta*Mtot;
			mass1 = Mtot*(1.0+sqrt(1.0-4.0*eta))/2.0;
			mass2 = Mtot*(1.0-sqrt(1.0-4.0*eta))/2.0;
			break;
		case 16:
			mum = Param_n;
			Mtot = pow(Mchirp,2.5) / pow(mum,1.5);
			eta = mum/Mtot;
			mass1 = Mtot*(1.0+sqrt(1.0-4.0*eta))/2.0;
			mass2 = Mtot*(1.0-sqrt(1.0-4.0*eta))/2.0;
			break;
		case 17:
			eta = Param_n;
			Mtot = Mchirp / pow(eta,0.6);
			mum = eta*Mtot;
			mass1 = Mtot*(1.0+sqrt(1.0-4.0*eta))/2.0;
			mass2 = Mtot*(1.0-sqrt(1.0-4.0*eta))/2.0;
			break;
		case 18:
			phic = Param_n;
			Phi0 = -1.1e30;
			break;
		case 19:
			Mtot = Param_n;
			Mchirp = Mtot * pow(eta,0.6);
			mum = eta*Mtot;
			mass1 = Mtot*(1.0+sqrt(1.0-4.0*eta))/2.0;
			mass2 = Mtot*(1.0-sqrt(1.0-4.0*eta))/2.0;
			break;
		case 20:
			PolNR = Param_n;
			break;
		case 21:
			thetaJ = Param_n;   // changed by Sofiane   ThdNR
			break;
		case 22:
			phiJ = Param_n;     // changed by Sofiane   PhdNR
			break;
		default:
			break;
	}
	
	if((iP==2)||(iP==3)||(iP==15)||(iP==16)||(iP==17)||(iP==19)){
		m1M2 = mass1*mass1/(Mtot*Mtot);
		m2M2 = mass2*mass2/(Mtot*Mtot);
		m1 = mass1 * LC::TSUN;
		m2 = mass2 * LC::TSUN; 
		M = Mtot * LC::TSUN;
		mu = m1*m2/M;
		dm = m1-m2;
	}
	
	if((iP==0)||(iP==1)){
		thetaS = M_PI/2 - Beta;
		phiS = Lambda;
		cThS = cos(thetaS);
		sThS = sin(thetaS);
	}
	
	if((iP!=4)&&(iP!=5)&&(iP!=12)&&(iP!=18))
		NeedExtraParamCompute = true;
	
	if((iP==8)||(iP==9)||(iP==10)||(iP==11)||(iP==13)||(iP==14))
		TypeExtraParamCompute = 2;
	
	if((iP==20)||(iP==21)||(iP==22))
		TypeExtraParamCompute = 1;
	
}


void LCGWSpinBBHHHarm1::getRange(int iP, double &Pmin, double &Pmax)
{
	double m1min(0.5e6),  m1max(5.0e6);
	double m1om2min(1.0), m1om2max(10.0);
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
			Pmin = m1min;
			Pmax = m1max;
			break;
		case 3:
			Pmin = m1min/m1om2max;
			Pmax = m1max/m1om2min;
			break;
		case 4:
			Pmin = 1.0e5;
			Pmax = 7.51e7;
			break;
		case 5:
			Pmin = 1.0e3;
			Pmax = 1.0e7;
			break;
		case 6:
			Pmin = 0.0;
			Pmax = 1.0;
			break;
		case 7:
			Pmin = 0.0;
			Pmax = 1.0;
			break;
		case 8:
			Pmin = 0.0;
			Pmax = M_PI;
			break;
		case 9:
			Pmin = 0.0;
			Pmax = M_PI;
			break;
		case 10:
			Pmin = 0.0;
			Pmax = 2.0*M_PI;
			break;
		case 11:
			Pmin = 0.0;
			Pmax = 2.0*M_PI;
			break;
		case 12:
			Pmin = 0.0;
			Pmax = 2.0*M_PI;
			break;
		case 13:
			Pmin = 0.0;
			Pmax = M_PI;
			break;
		case 14:
			Pmin = 0.0;
			Pmax = 2.0*M_PI;
			break;
		case 15:
			Pmin = m1min*pow(1.0+1.0/m1om2max,0.4)/pow(1+m1om2max,0.6);
			Pmax = m1max*pow(1.0+1.0/m1om2min,0.4)/pow(1+m1om2min,0.6);
			break;
		case 16:
			Pmin = m1min/(1.0+m1om2max);
			Pmax = m1max/(1.0+m1om2min);
			break;
		case 17:
			Pmin = 1.0/(2.0+m1om2max+1.0/m1om2max);
			Pmax = 0.25;
			break;
		case 18:
			Pmin = 0.0;
			Pmax = 2.0*M_PI;
			break;
		case 19:
			Pmin = m1min*(1.+1./m1om2max);
			Pmax = m1max*(1.+1./m1om2min);
			break;
		case 20:
			Pmin = 0.0;
			Pmax = 2.0*M_PI;
			break;
		case 21:
			Pmin = 0.0;
			Pmax = M_PI;
			break;
		case 22:
			Pmin = 0.0;
			Pmax = 2.0*M_PI;
			break;
		default:
			Cout << "ERROR in LCGWSpinBBHHHarm1::setParam : The parameter " << iP << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWSpinBBHHHarm1::setParam : The parameter is unknown !");
			break;
	}
}


double LCGWSpinBBHHHarm1::getDelta(int iP)
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
		case 4:
			DRes = 1.e-1;
			break;
		case 6:
			DRes = 1.e-5;
			break;
		case 7:
			DRes = 1.e-5;
			break;
		case 8:
			DRes = 1.e-4;
			break;
		case 9:
			DRes = 1.e-3;
			break;
		case 10:
			DRes = 1.e-4;
			break;
		case 11:
			DRes = 1.e-3;
			break;
		case 12:
			DRes = 1.e-3;
			break;
		case 13:
			DRes = 1.e-5;
			break;
		case 14:
			DRes = 1.e-5;
			break;
		case 15:
			DRes = 1.e-2;
			break;
		case 17:
			DRes = 1.e-6;
			break;
		case 19:
            DRes = 3.e-5*Mtot; // 3.e-1 changed by Sofiane
            break;
		case 20:
			DRes = 1.e-5;
			break;
		case 21:
			DRes = 2.e-4;  // changed by Sofiane
			break;
		case 22:
			DRes = 2.e-4;  // changed by Sofiane
			break;
		default:
			DRes = (tmpMax-tmpMin)*5.0e-8;
	};
	//Cout << iP << " --> Delta = " << DRes<< Endl;
	return (DRes);
}



double LCGWSpinBBHHHarm1::getSpecialParam(int iPS)
{
	switch (iPS){
		case 0:
			return(t0Intg);
			break;
		case 1:
			return(pow((tc-t0Intg)*(eta*256)/(5.*M),1./4.));
			break;
		default:
			Cout << "ERROR in LCGWSpinBBHHHarm1::getSpecialParam : The parameter " << iPS << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWSpinBBHHHarm1::getSpecialParam : The parameter is unknown !");
			break;
	}
	return(0.);
}


void LCGWSpinBBHHHarm1::setSpecialParam(int iPS, double SpecParam_n)
{
	switch (iPS){
		case 0:
			t0Intg = SpecParam_n;
			break;
		case 1:
			t0Intg = tc - pow(SpecParam_n,4.) * (5.*M)/(eta*256) ;
			break;
		default:
			Cout << "ERROR in LCGWSpinBBHHHarm1::setSpecialParam : The parameter " << iPS << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWSpinBBHHHarm1::setSpecialParam : The parameter is unknown !");
			break;
	}
}



void LCGWSpinBBHHHarm1::RandParam(int iP)
{
	double Pmin(0.0);
	double Pmax(1.0);
	
	getRange(iP, Pmin, Pmax);
	setParam(iP,MT->RandUniform(Pmin, Pmax));
	
	// ** Specifical case : non-uniform
	if(iP == 0)
		setParam(0, M_PI/2.0 - acos(MT->RandUniform(-1.0, 1.0)));
	if(iP == 8)
		setParam(8, acos(MT->RandUniform(-1.0, 1.0)));
	if(iP == 9)
		setParam(9, acos(MT->RandUniform(-1.0, 1.0)));
	if(iP == 13)
		setParam(13, acos(MT->RandUniform(-1.0, 1.0)));
}

void LCGWSpinBBHHHarm1::setTimeInfo(double t0_n, double dt_n, double TObs_n, double tAskMin_n, double tDOrbMax_n,  double tMaxDiff)
{
    t0Real = t0_n; 
	tAskMin = tAskMin_n;
	tDOrbMax = tDOrbMax_n;
	
	//! *** Free previous allocation of last steps of interation
	if(coordI != NULL){
		for(int i=0; i<NMemItg; i++)
			MT->Free(coordI[i], 13*sizeof(double) );
		MT->Free(coordI, NMemItg*sizeof(double*) );
	}
	
	//! *** Define the integration step and the number of steps to store
	dt = dt_n;
	NMemItg = MAX(2, MT->iceil( 2. * tMaxDiff / dt) ) + 1;
	
	/*! ** Set time for starting integration to be sure that we have enough data and be able to manage any kind of GW direction relative to detector position :
	 *	Actually the reference time for GW is defined at the solar system barycenter and the time of the simulation at the detector */
	t0Start = t0_n - tDOrbMax;
    t0Start = -4000;
	maxDur = TObs_n;
	
	//! *** Allocate memory for the last steps of interation
	coordI = (double **) MT->AllocMemory(NMemItg*sizeof(double*));
	for(int iM=0; iM<NMemItg; iM++){
		coordI[iM] = (double*) MT->AllocMemory(13*sizeof(double));
		for(int i=0; i<13; i++)
			coordI[iM][i] = 0.;
	}
	
}

// ***************************************
// * Linking and initialization methods  *
// ***************************************

int LCGWSpinBBHHHarm1::init()
{
	//DispInfo("> C > \t");
	
	initBase();
	
	double iota0(InitialPolarAngleL);
	
	if(MT->Disp()){
		Cout.precision(13);
		Cout << "   --> Initialisation of LCGWSpinBBHHHarm1 (compute precession for " << maxDur << " s starting at " << t0Start << " s ) ... " << Endl;
		//DispAllParam();
	}
	
	Amp = 2.*M*eta/dist;
	//Amp = pow(M,3./5.)*eta/dist;
	
	//ComputePositionConsThS();
	
	
	S1x = sin(PolarAngleOfSpin1)*cos(AzimuthalAngleOfSpin1);
	S1y = sin(PolarAngleOfSpin1)*sin(AzimuthalAngleOfSpin1);
	S1z = cos(PolarAngleOfSpin1);
	
	S2x = sin(PolarAngleOfSpin2)*cos(AzimuthalAngleOfSpin2);
	S2y = sin(PolarAngleOfSpin2)*sin(AzimuthalAngleOfSpin2);
	S2z = cos(PolarAngleOfSpin2);
	
	Lnx = sin(InitialPolarAngleL)*cos(InitialAzimuthalAngleL);
	Lny = sin(InitialPolarAngleL)*sin(InitialAzimuthalAngleL);
	Lnz = cos(InitialPolarAngleL);
	
	
	LS1 = Lnx*S1x + Lny*S1y + Lnz*S1z;
	LS2 = Lnx*S2x + Lny*S2y + Lnz*S2z;
	S1S2 = S1x*S2x + S1y*S2y + S1z*S2z;
	
	//! ** Estimation of initial frequency using 2PN expression and neglecThS precession 
	om = ComputeApprOrbFreq(t0Intg, tc); 
	om_prev = om;
	
	if((phic>-1.0e30)&&(Phi0<-1.0e30)){
		//! ** Estimation of initial phase using 2PN expression and neglecThS precession
		Phi0 = ComputeApprOrbPhase(t0Intg, tc); 
	}
	if(MT->Disp()){
		Cout << "\tIntitial frequency is " << om << "  " << om/M_PI << Endl;
		Cout << "\tIntitial phase is " << Phi0 << Endl;
	}
	
	
	//! *** Transform direction of vectors in the source frame: 
	
	//! ** Total angular momentum
	double Jx = eta*M*M*pow(M*om, -1./3)*Lnx + x1*m1*m1*S1x + x2*m2*m2*S2x;
	double Jy = eta*M*M*pow(M*om, -1./3)*Lny + x1*m1*m1*S1y + x2*m2*m2*S2y;
	double Jz = eta*M*M*pow(M*om, -1./3)*Lnz + x1*m1*m1*S1z + x2*m2*m2*S2z;
	double absJ = sqrt(Jx*Jx + Jy*Jy + Jz*Jz);
	
	/*
	 if(MT->Disp()){
	 Cout << "\tL in bar. ref. : " << eta*M*M*pow(M*om, -1./3)*Lnx << " " << eta*M*M*pow(M*om, -1./3)*Lny << " " << eta*M*M*pow(M*om, -1./3)*Lnz << Endl;
	 Cout << "\tJ in bar. ref. : " << Jx << " " << Jy << " " << Jz << Endl;
	 }*/
	
	if (!(absJ > 0.0))
		throw std::invalid_argument("ERROR in MFLModTSrcLWSpinBBH3PN::init : total angular momentum iz zero at t=0");
	double thetaJ = acos(Jz/absJ);
	double phiJ;
	if(fabs(thetaJ) <= 1.e-6 || fabs(thetaJ - M_PI) <= 1.e-6){
		phiJ = 0.0;
	}else{
		phiJ = atan2(Jy, Jx);
	}
	
	double stJ = sin(thetaJ);
	double ctJ = cos(thetaJ);
	//double down = cThS*stJ*cos(phiS-phiJ) - ctJ*sThS;
	//double up = stJ*sin(phiS - phiJ);
	//psi = atan2(-up, down);
	double up = cThS*stJ*cos(phiS-phiJ) - ctJ*sThS;
	double down = stJ*sin(phiS - phiJ);
	PolB = atan2(up, down);
	c2PB = cos(2.*PolB);
	s2PB = sin(2.*PolB);
	
	theta = acos(-ctJ*cThS - stJ*sThS*cos(phiS - phiJ));
	
	
	//! *** Compute components of the source frame basis in SSB
	
	n[0] = sThS*cos(phiS); 
	n[1] = sThS*sin(phiS); 
	n[2] = cThS;
	
	ez[0] = stJ*cos(phiJ);  
	ez[1] = stJ*sin(phiJ); 
	ez[2] = ctJ;
	
	double stheta = fabs(sin(theta));
	if (!(stheta != 0.0))
		throw std::invalid_argument("ERROR in MFLModTSrcLWSpinBBH3PN::init : source frame is undefined: source direction is colinear with J");
	
	
	ey[0] = (n[1]*ez[2] - n[2]*ez[1])/stheta;
	ey[1] = (n[2]*ez[0] - n[0]*ez[2])/stheta;
	ey[2] = (n[0]*ez[1] - n[1]*ez[0])/stheta;
	
	double nez = n[0]*ez[0] + n[1]*ez[1] + n[2]*ez[2];
	ex[0] = (ez[0]*nez - n[0])/stheta;
	ex[1] = (ez[1]*nez - n[1])/stheta;
	ex[2] = (ez[2]*nez - n[2])/stheta;
	
	
	//! *** Now we can compute the directions of Ln, S1, S2 in the source frame
	double iotaIn = acos(Lnx*ez[0] + Lny*ez[1] + Lnz*ez[2]);
	double LNey = Lnx*ey[0] + Lny*ey[1] + Lnz*ey[2]; 
	double LNex = Lnx*ex[0] + Lny*ex[1] + Lnz*ex[2];
	double alphaIn = atan2(LNey, LNex); 
	if(fabs(iotaIn) <= 1.e-6 || fabs(iotaIn - M_PI) <= 1.e-6){
		Cout << "Warning: Ln co(anti)-alligned with J, alpha will be put to zero" << Endl;
		alphaIn = 0.0;
	}
	double thetaS1In = acos(S1x*ez[0] + S1y*ez[1] + S1z*ez[2]);
	double S1ey = S1x*ey[0] + S1y*ey[1] + S1z*ey[2];
	double S1ex = S1x*ex[0] + S1y*ex[1] + S1z*ex[2];
	double phiS1In = atan2(S1ey, S1ex);
	if(fabs(thetaS1In) <= 1.e-6 || fabs(thetaS1In - M_PI) <= 1.e-6){
		Cout << "Warning: S1 co-alligned with J, phiS1 will be put to zero" << Endl;
		phiS1In = 0.0;
	}
	
	double thetaS2In = acos(S2x*ez[0] + S2y*ez[1] + S2z*ez[2]);
	double S2ey = S2x*ey[0] + S2y*ey[1] + S2z*ey[2];
	double S2ex = S2x*ex[0] + S2y*ex[1] + S2z*ex[2];
	double phiS2In = atan2(S2ey, S2ex);
	if(fabs(thetaS2In) <= 1.e-6 || fabs(thetaS2In - M_PI) <= 1.e-6){
		Cout << "Warning: S2 co-alligned with J, phiS2 will be put to zero" << Endl;
		phiS2In = 0.0;
	}
	
	iota = iotaIn;
	alpha = alphaIn;
	nonspin = false;
	if (chi1 == 0.0 && chi2 == 0.0){
        iota0 = 0.0; // \todo ASK TO STAS IF IT'S REALLY iota0 and not iota ?
        alpha = 0.0;
        nonspin = true;
	}
	
	Lnx = sin(iota)*cos(alpha);
	Lny = sin(iota)*sin(alpha);
	Lnz = cos(iota);
	
	S1z = cos(thetaS1In);
	S1x = sin(thetaS1In)*cos(phiS1In);
	S1y = sin(thetaS1In)*sin(phiS1In);
	S2z = cos(thetaS2In);
	S2x = sin(thetaS2In)*cos(phiS2In);
	S2y = sin(thetaS2In)*sin(phiS2In);
	
	/*
	 if(MT->Disp()){
	 Cout << "\tL in src. ref. : " << Lnx << " " << Lny << " " << Lnz << Endl;
	 Cout << "\tMatrix RotB2S : Xsrc = RotB2S Xbar : " << Endl;
	 Cout << "\t\t" << ex[0] << " " << ex[1] << " " << ex[2] << Endl ;
	 Cout << "\t\t" << ey[0] << " " << ey[1] << " " << ey[2] << Endl ;
	 Cout << "\t\t" << ez[0] << " " << ez[1] << " " << ez[2] << Endl ;
	 }*/
	
	/* can uncomment bellow and compare with original inner producThS
	 LS1 = Lnx*S1x + Lny*S1y + Lnz*S1z;
	 LS2 = Lnx*S2x + Lny*S2y + Lnz*S2z;
	 S1S2 = S1x*S2x + S1y*S2y + S1z*S2z; */
	
	
	ComputeInspiral();
	
	
	//! ** Variable depending on theta only
	cth = cos(theta);
	sth = sin(theta);
	c2th = cos(2.*theta);
	s2th = sin(2.*theta);
	c3th = cos(3.*theta);
	s3th = sin(3.*theta);
	c4th = cos(4.*theta);
	s4th = sin(4.*theta);
	cth2 = cth*cth;
	sth2 = sth*sth;
	cth3 = cth2*cth;
	cthsth2 = cth*sth*sth;
	sthcth3 = sth*cth3;
	
	
	ChpC1th = cos(theta);
	ChpS1th = sin(theta);
	ChpC1thE2 = ChpC1th * ChpC1th;        
	ChpS1thE2 = ChpS1th * ChpS1th;       
	ChpC1thE3 =  ChpC1thE2 * ChpC1th;
	ChpC1thE4 =  ChpC1thE3 * ChpC1th;
	ChpC2th = 2. * ChpC1thE2 - 1.;
	ChpC3th = 4. * ChpC1thE3 - 3. * ChpC1th;
	ChpC4th = 8 * ChpC1thE4 - 8. * ChpC1thE2 + 1.;
	ChpS2th = 2. * ChpS1th * ChpC1th;
	ChpS3th = ChpS1th * (4. * ChpC1thE2 - 1.);
	ChpS4th = ChpS1th * (8. * ChpC1thE3 - 4. * ChpC1th);
	Chp5S1thPS3th = 5. * ChpS1th + ChpS3th;
	Chp3PC2th = 3. + ChpC2th;
	Chp1P3C2th = 1. + 3. * ChpC2th;
	ChpS1thM3S3th = ChpS1th - 3. * ChpS3th;
	Chp5p4C2thP7C4th = 5. + 4. * ChpC2th + 7. * ChpC4th;
	ChpC1thXS1th = ChpC1th * ChpS1th;  
	
	
	//tTaper = tend - M_PI*taperQ/freq[ind_end-1];
	
	tTaper = tend - M_PI*taperQ/FreqMax;
	
	if(MT->Disp()){
		if(ApplyTaper)
			Cout << "\t Taper : applied" << Endl;
		else
			Cout << "\t Taper : not applied" << Endl;
	}
	
	//if(MT->Disp())
	//	Cout << " --> OK (Range of frequency = [" << FreqMin << "," << FreqMax << "] Hz)" << Endl;
    return 0 ;
}


// *********************
// *  Running methods  *
// *********************


void LCGWSpinBBHHHarm1::Computehpc(double t)
{
	double LinIO, LinIY;
	double rhs[13];
	double newstep = dt;
	double tr;
	int bin;
	
	/*! NOTE : that for the LW approximation we might want to 
	 *	compute the waveform at the un-even time intervals modified
	 *  by the Doppler modulation \vec{k}.\vec{R}
	 *  to do so we can use linear interpolation of smooth functions like 
	 *  spin and orbiatl evolution and then
	 * compute the waveform (which oscillatory)
	 */
	
	
	if (t <= tend){
		
		//! **** While required time #t is not reach by current time of integration then make forward integration
		while (tIc < t) {
			
			//! *** Store previous steps
			for(int iM=NMemItg-1; iM>0; iM--)
				for (int i=0; i<13; i++)
					coordI[iM][i] = coordI[iM-1][i];
			
			//! Compute current step
			Derivs(t, coordI[1], rhs, 13);
			newstep =  Integrator(coordI[1], rhs, coordI[0], tIc, dt, 13);
			tIc += dt;
		}
		
		
		//! *** Compute value at required time using linear interpolation
		//! Find the bin
		tr = (tIc-t)/dt ;
		bin = MT->ifloor(tr);
		if((bin<0)||(bin+1 > NMemItg-1)){
			std::cerr << "ERROR in LCGWSpinBBHHHarm1::Computehpc : bin = " << bin << " (x = " << x << ") and bin+1 " << bin+1 << " are not included in [0," << NMemItg-1 << "]." << Endl;  
			throw std::invalid_argument("ERROR in LCGWSpinBBHHHarm1::Computehpc : The required bin does not exist !");
		}
		
		LinIO = tr - bin;
		LinIY = 1. - LinIO;
		om		= LinIO * coordI[bin+1][0] + LinIY * coordI[bin][0];
		phi		= LinIO * coordI[bin+1][12] + LinIY * coordI[bin][12];
		iota	= LinIO * coordI[bin+1][1] + LinIY * coordI[bin][1];
		alpha	= LinIO * coordI[bin+1][2] + LinIY * coordI[bin][2];
		Psi = phi - 2.*M*om*log(om/om0);
		S1x	= LinIO * coordI[bin+1][6] + LinIY * coordI[bin][6];
		S1y = LinIO * coordI[bin+1][7] + LinIY * coordI[bin][7];
		S1z = LinIO * coordI[bin+1][8] + LinIY * coordI[bin][8];
		S2x = LinIO * coordI[bin+1][9] + LinIY * coordI[bin][9];
		S2y = LinIO * coordI[bin+1][10] + LinIY * coordI[bin][10];
		S2z = LinIO * coordI[bin+1][11] + LinIY * coordI[bin][11];
		//tt = LinIO * (tIc[iItg]-dt/2.) + LinIY * (tIc[iItg]+dt/2.);
		
		
		
		
		ci = cos(iota);
		si = sin(iota);
		ciby2 = cos(0.5*iota);
		siby2 = sin(0.5*iota);
		v = pow(M*om, 1./3.);
		
		/*
		Computeh0PN(h0PNp, h0PNc);
		h05PNp = 0.; h05PNc = 0.;
		Computeh05PN(h05PNp, h05PNc);
		h1PNp = 0.; h1PNc = 0.;         
		Computeh1PN(h1PNp, h1PNc);
		*/
		
		h05PNp = 0.; h05PNc = 0.;
		h1PNp = 0.; h1PNc = 0.;         
		ComputehAllPN(h0PNp, h0PNc, h05PNp, h05PNc, h1PNp, h1PNc);
		
		ApplyTaper = true;
		if(ApplyTaper)
			wk = Amp*halfhann(t, tend, tTaper);
		else
			wk = Amp;
		hp = wk*v*v*(h0PNp + v*h05PNp + v*v*h1PNp);
		hc = wk*v*v*(h0PNc + v*h05PNc + v*v*h1PNc);
		
		/*	
		 Lnx = sin(iota)*cos(alpha);
		 Lny = sin(iota)*sin(alpha);
		 Lnz = cos(iota);
		 
		 LnBx = Lnx*ex[0] + Lny*ey[0] + Lnz*ez[0];
		 LnBy = Lnx*ex[1] + Lny*ey[1] + Lnz*ez[1];
		 LnBz = Lnx*ex[2] + Lny*ey[2] + Lnz*ez[2];
		 PhLnB = atan2(LnBy,LnBx);
		 ThLnB = atan2(sqrt(LnBx*LnBx+LnBy*LnBy),LnBz);
		 PolB = atan2( cThS*cos(phiS-PhLnB)*sin(ThLnB)-cos(ThLnB)*sThS , sin(ThLnB)*sin(phiS-PhLnB) );
		 c2PB = cos(2.*PolB);
		 s2PB = sin(2.*PolB);
		 */
		
		hBpLast = -c2PB*hp - s2PB*hc;
		hBcLast =  s2PB*hp - c2PB*hc;
		//hBp = hp;
		//hBp = hc;
		
		/*
		if((t>6.7e7)){
			Cout.precision(12);
			Cout << t << " "  << om << " "  << phi << " "  << iota << " "  << alpha << " "  << S1x << " "  << S1y << " "  << S1z << " "  << S2x << " " << S2y << " " << S2z;
			Cout << " " << hp << " " << hc << " " << Amp << " " << wk << " " << h0PNp << " " << h0PNc << " " << h1PNp << " " << h1PNc << " " << v;
			//Cout << " " << c2PB << " " << s2PB << " " << hBpLast << " " << hBcLast 
			Cout << Endl;
		}
		*/
		
        //Cout << t << " " << hBpLast << " " << hBcLast << Endl;
        
	}else{
		hBpLast = 0.0;
		hBcLast = 0.0;
	}
}




// *******************
// *  Other methods  *
// *******************

void LCGWSpinBBHHHarm1::DispInfo(char * BTab)
{
	MT->o->precision(12);
	
	if(MT->Disp()){
		DispInfoBase(BTab);
	}
    
    Cout << BTab << "\t- Polarization using initial data    = " << PolB << " rad   (tan = " << tan(PolB) << ")" << Endl;
    Cout << BTab << "\t- thetad    = " << theta << " rad "<< Endl;
	Cout << BTab << "\t- Beta    = " << Beta << " rad "<< Endl;
	Cout << BTab << "\t- Lambda  = " << Lambda << " rad "<< Endl;
	Cout << BTab << "\t- m1      = " << mass1 << " MSun "<< Endl;
	Cout << BTab << "\t- m2    = " << mass2 << " MSun "<< Endl;
	Cout << BTab << "\t- tc    = " << tc << " s "<< Endl;
	Cout << BTab << "\t- DL    = " << DL << " kpc "<< Endl;
	Cout << BTab << "\t- chi1  = " << chi1 << Endl;
	Cout << BTab << "\t- chi2  = " << chi2 << Endl;
	Cout << BTab << "\t- PolarAngleOfSpin1      = " << PolarAngleOfSpin1 << " rad "<< Endl;
	Cout << BTab << "\t- PolarAngleOfSpin2      = " << PolarAngleOfSpin2 << " rad "<< Endl;
	Cout << BTab << "\t- AzimuthalAngleOfSpin1  = " << AzimuthalAngleOfSpin1 << " rad "<< Endl;
	Cout << BTab << "\t- AzimuthalAngleOfSpin2  = " << AzimuthalAngleOfSpin2 << " rad "<< Endl;
	Cout << BTab << "\t- phic                   = " << phic << " rad "<< Endl;
	Cout << BTab << "\t- phi0                   = " << Phi0 << " rad "<< Endl;
	Cout << BTab << "\t- InitialPolarAngleL      = " << InitialPolarAngleL << " rad "<< Endl;
	Cout << BTab << "\t- InitialAzimuthalAngleL  = " << InitialAzimuthalAngleL << " rad "<< Endl;
	Cout << BTab << "\t- Mchirp  = " << Mchirp << " MSun "<< Endl;
	Cout << BTab << "\t- mu      = " << mum << " MSun "<< Endl;
	Cout << BTab << "\t- eta     = " << eta << Endl;
	Cout << BTab << "\t- Mtot  = " << Mtot << " MSun " << Endl;
	Cout << BTab << "\t- PolNR = " << PolNR << " rad (tan = " << tan(PolNR) << " )" << Endl;
	Cout << BTab << "\t- ThdNR = " << ThdNR << " rad" << Endl;
	Cout << BTab << "\t- PhdNR = " << PhdNR << " rad" << Endl;
	Cout << BTab << "\t- AmpS1 = " << AmpS1 << " " << Endl;
	Cout << BTab << "\t- chi1 = " << chi1 << " " << Endl;
	Cout << BTab << "\tExtra infos :" << Endl;
	//Cout << BTab << "\t+ Start integration at t = " << getSpecialParam(0) << " s  = " << getSpecialParam(1) << " M" << Endl;
	Cout << BTab << "\t- Ln in SSB = " << LnB << Endl;
    Cout << BTab << "\t- S1 in SSB = " << S1B << Endl;
    Cout << BTab << "\t- S2 in SSB = " << S2B << Endl;
	Cout << BTab << "\t- J in NRF  = " << JN << Endl; 
	Cout << BTab << "\t- Ln in NRF = " << LnN << Endl;
    Cout << BTab << "\t- S1 in NRF = " << S1N << Endl;
    Cout << BTab << "\t- S2 in NRF = " << S2N << Endl;
    Cout << BTab << "\t- TypeExtraParamCompute = " << TypeExtraParamCompute << Endl;
    
	
}


void LCGWSpinBBHHHarm1::DispAllParam(std::ostream * out)
{
	for(int iP=0; iP<NParams; iP++)
		(*out) << " " << getParam(iP);
}


void LCGWSpinBBHHHarm1::DispAllParamName(std::ostream * out)
{
	(*out) << " Bet Lam m1 m2 tc DL chi1 chi2 thS1 thS2 phS1 phS2 Phic Phi0 thL phL Mc mu eta";
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

void LCGWSpinBBHHHarm1::getTaperInfo(double & tTaper_s, double & tendTaper_s, double & FreqMaxTaper_s)
{
	tTaper_s = tTaper;
	tendTaper_s = tend;
	FreqMaxTaper_s = FreqMax;
}


// ***************************************
// * Linking and initialization methods  *
// ***************************************



// *********************
// *  Running methods  *
// *********************




// *******************
// *  Other methods  *
// *******************




void LCGWSpinBBHHHarm1::ComputeInspiral()
{
	
	
	//! ***** Declaration and initilaization
	double coord0[13];
	double coordn[13];
	double rhs[13];
	
	coord0[0] = om;
	coord0[1] = iota;
	coord0[2] = alpha;
	coord0[3] = Lnx;
	coord0[4] = Lny;
	coord0[5] = Lnz;
	coord0[6] = S1x;
	coord0[7] = S1y;
	coord0[8] = S1z;
	coord0[9] = S2x;
	coord0[10] = S2y;
	coord0[11] = S2z;
	coord0[12] = Phi0;
	
	
	om0 = om;
	double iota0 = iota;
	double alpha0 = alpha;
	double Lnx0 = Lnx;
	double Lny0 = Lny;
	double Lnz0 = Lnz;
	double S1x0 = S1x;
	double S1y0 = S1y;
	double S1z0 = S1z;
	double S2x0 = S2x; 
	double S2y0 = S2y;
	double S2z0 = S2z;
    
	double newstep = dt;
	
	double t=t0Intg;
	int ind0 = floor((t0Intg-t0Start)/dt);
	int ind = ind0;
	double fnyby4 = 0.125/dt;
	double om_term = 2.*M_PI*fnyby4;
	
	
	//! *** Store the intial point of integration which is the starting point if there is no backward integration
	tIc = t;
	for (int i=0; i<13; i++)
		coordI[0][i] = coord0[i];
	
	
	
	//CheckData = true;
    /*
	LS1 = Lnx*S1x + Lny*S1y + Lnz*S1z;
	LS2 = Lnx*S2x + Lny*S2y + Lnz*S2z;
	S1S2 = S1x*S2x + S1y*S2y + S1z*S2z;
    */
	
	Cout << "Initial condition : " << Endl;
	Cout << t0Intg << " " << t0Intg/M << " " << coord0[0] << " " << M*coord0[0] << " " << LS1 << " " << LS2 << " " << S1S2 ;
	for(int i=1; i<13; i++)
		Cout << " " << coord0[i];
	Cout << " 0." << Endl;
	
	
	//! **** Backward integration : obtain values for the first time to be computed
	
	if ( t0Intg > t0Start && ind0 > 0){ 
		back = true;
		
		while (ind > 0){
			om_prev = om;
			Derivs(t, coord0, rhs, 13);
			newstep =  Integrator(coord0, rhs, coordn, t, dt, 13); 
			//for(int i=0; i<13; i++)
			//	Cout << coord0[i] << " ";
			//Cout << Endl;
			om = coordn[0];
			t = t-dt;
			//Cout << t << " " << om << Endl;
			if (om > om_prev){
				std::cerr << " ERROR: frequency condition  reached in back ward integration: " \
				<< om  << "   " << om_prev << "  " << t << Endl;
				break;
			}
			iota = coordn[1];
			alpha = coordn[2];
			Lnx = coordn[3];
			Lny = coordn[4];
			Lnz = coordn[5];
			S1x = coordn[6];
			S1y = coordn[7];
			S1z = coordn[8];
			S2x = coordn[9];
			S2y = coordn[10];
			S2z = coordn[11];
			phi = coordn[12];
			ind --;
			for (int i=0; i<13; i++){
				coord0[i] = coordn[i];
			}
			/*
			 tm[ind] = t;
			 Phase[ind] = phi;
			 freq[ind] = om;
			 io[ind] = iota;
			 al[ind] = alpha;
			 s1x[ind] = S1x;
			 s1y[ind] = S1y;
			 s1z[ind] = S1z;
			 s2x[ind] = S2x;
			 s2y[ind] = S2y;
			 s2z[ind] = S2z;
			 */
			
			/*
             // =====================  Checking display  =====================
             if((t<1000.)||((t>=t0Intg-1000.*dt)&&(t<t0Intg+1000.*dt))){
             //CheckData = true;
             LS1 = Lnx*S1x + Lny*S1y + Lnz*S1z;
             LS2 = Lnx*S2x + Lny*S2y + Lnz*S2z;
             S1S2 = S1x*S2x + S1y*S2y + S1z*S2z;
             Cout << t << " " << t/M << " " << coordn[0] << " " << M*coordn[0] << " " << LS1 << " " << LS2 << " " << S1S2 ;
             for(int i=1; i<13; i++)
             Cout << " " << coordn[i];
             Cout << " 0. " << Endl;
             }else{
             CheckData = false;
             }
             // ==========================================
             */
		}
		
		//! *** Store the point reach at the end of backward integration which is the starting point of time 
		tIc = t;
		for (int i=0; i<13; i++)
			coordI[0][i] = coord0[i];
		
		
		//! *** Set the value at the inital condition for preparing the forward integration */
		t= t0Intg;
		ind = ind0;
		coord0[0] = om0;
		coord0[1] = iota0;
		coord0[2] = alpha0;
		coord0[3] = sin(iota0)*cos(alpha0);
		coord0[4] = sin(iota0)*sin(alpha0);
		coord0[5] = cos(iota0);
		coord0[6] = S1x0;
		coord0[7] = S1y0;
		coord0[8] = S1z0;
		coord0[9] = S2x0;
		coord0[10] = S2y0;
		coord0[11] = S2z0;
		coord0[12] = Phi0;
		
		om = om0;
		iota = iota0;
		alpha = alpha0;
		Lnx = sin(iota0)*cos(alpha0);
		Lny = sin(iota0)*sin(alpha0);
		Lnz = cos(iota0);
		S1x = S1x0;
		S1y = S1y0;
		S1z = S1z0;
		S2x = S2x0;
		S2y = S2y0;
		S2z = S2z0;
		phi = Phi0;
	}   
	FreqMin = om/(4.*M_PI);
	
	//! **** Forward integration : for reaching the termination condition : store end time and max frequency
	om_prev = om;
	double Ener = ComputeEnergy();
	En_prev = Ener;
	back = false;
	ind = ind0;
	bool ReachEnd(true);
	while ( t<= maxDur){
		/*
		 tm[ind] = t;
		 Phase[ind] = phi;
		 freq[ind] = om;
		 io[ind] = iota;
		 al[ind] = alpha;
		 s1x[ind] = S1x;
		 s1y[ind] = S1y;
		 s1z[ind] = S1z;
		 s2x[ind] = S2x;
		 s2y[ind] = S2y;
		 s2z[ind] = S2z;
		 */
		Derivs(t, coord0, rhs, 13);
		newstep =  Integrator(coord0, rhs, coordn, t, dt, 13);
		ind++;
		t += dt;
		om = coordn[0];
		iota = coordn[1];
		alpha = coordn[2];
		Lnx = coordn[3];
		Lny = coordn[4];
		Lnz = coordn[5];
		S1x = coordn[6];
		S1y = coordn[7];
		S1z = coordn[8];
		S2x = coordn[9];
		S2y = coordn[10];
		S2z = coordn[11];
		phi = coordn[12];
		Ener = ComputeEnergy();
		if (Ener > En_prev){
			if(MT->Disp()){
				Cout << "Instability condition is met at MECO: dE/domega = 0: orb. ang. freq. = " << om << "  " << om_prev << Endl;
				Cout << std::setprecision(15) << "*** we reached merger condition at t = " << t << Endl;
			}
			tend = t-dt;
			FreqMax = 2.*om_prev/M_PI;
			ReachEnd = false;
			break;
		}
		if (om < om_prev){
			if(MT->Disp()){
				Cout << "Instability condition is met at domega/dt = 0: orb. ang. freq. = " << om << "   " << om_prev << Endl;
				Cout << std::setprecision(15) << "*** we reached merger condition at t = " << t << Endl;
			}
			tend = t-dt;
			FreqMax = 2.*om_prev/M_PI;
			ReachEnd = false;
			break;
		}
		if (om >= om_term){
			if(MT->Disp()){
				Cout << "Reached termination condition f<f_ny/4: orb. ang. freq. = " << om << "   " << om_term << Endl;
				Cout << std::setprecision(15) << "*** we reached merger condition at t = " << t << Endl;
			}
			tend = t-dt;
			FreqMax = 2.*om_prev/M_PI;
			ReachEnd = false;
			break;
		}
		om_prev = om; 
		En_prev = Ener; 
		for (int i=0; i<13; i++){
			coord0[i] = coordn[i];
		}
	}// while loop
	
	if(ReachEnd){
		tend = t-dt;
		FreqMax = 2.*om_prev/M_PI;
		ApplyTaper = false;
		if(MT->Disp()){
			Cout << "*** we reached end of observation (without merger)" << Endl;
		}
	}
	
    //tend = tm[ind_end -1];
	
	
	for (int iM=1; iM<NMemItg; iM++)
		for (int i=0; i<13; i++)
			coordI[iM][i] = coordI[0][i];
	
    
	/*
	 std::ofstream fout28("TesThSpinInsp.dat");
	 for (int i=0; i<ind_end; i++){
	 fout28 << std::setprecision(10) << tm[i] <<"   " <<  freq[i] << "   " << Phase[i];
	 fout28 << " " << io[i] << " " << al[i] ;
	 fout28 << " " << s1x[i] << " " << s1y[i] << " " << s1z[i] ;
	 fout28 << " " << s2x[i] << " " << s2y[i] << " " << s2z[i] ;
	 fout28 << Endl;
	 }
	 fout28.close();
	 */
	
	CheckData = false;
}




double LCGWSpinBBHHHarm1::ComputeApprOrbFreq(double t, double Tc)
{
	
	//! *** WARNING : THIS FUNCTION HAVE TO BE EXACTLY THE SAME AS THE ONE DEFINE IN LCGWSpinBBHNR1 TO KEEP COHERENT WAVEFORM 
	
	/*! Here I use 2PN approximate expression  for the frequency as function of t 
	 * neglecting the precession
	 */
	
    /*LS1 = Lnx*S1x + Lny*S1y + Lnz*S1z;
	 LS2 = Lnx*S2x + Lny*S2y + Lnz*S2z;
	 S1S2 = S1x*S2x + S1y*S2y + S1z*S2z;  computed in the constructor*/
	
	beta = (1.0/12.0)*( x1*LS1*(113.0*m1M2 + 75.0*eta) +  x2*LS2*(113.0*m2M2 + 75.0*eta) );
    sigma = (1./48.0)*eta*x1*x2*( -247.0*S1S2 + 721.0*LS1*LS2 );
    
	//Cout << t << " " << Tc << " " << x1 << " " << x2 << " " << LS1 << " " << LS2 << " " << S1S2 << " , " << m1M2 << " " << m2M2 << " " << beta << " " << sigma << Endl;
	
    
    tau = 0.2*eta*(Tc-t)/M;
    tau38 = pow(tau, -0.375);
    tau14 = pow(tau, -0.25);
    
	//   Cout << "stas check: " << beta << "  " << sigma << "   " << tau << Endl;
    
    x = 0.125*tau38*(1.0 + (743./2688. + 11.*eta/32.)*tau14 - 0.3*(M_PI - 0.25*beta)*tau38 + (1855099./14450688. + 56975./258048.*eta + 371./2048.*eta*eta - 3./64.*sigma)*tau14*tau14);
    
	//Cout << eta << " , " << tau << " , " << x << " / " << M << Endl;
	
    return(x/M);	
}


double LCGWSpinBBHHHarm1::ComputeApprOrbPhase(double t, double Tc){
	
	double xPhase;
	/*! Here I use 2PN approximate expression  for the orbital phase as function of t
	 * neglecting the precession
	 */
	
    /*LS1 = Lnx*S1x + Lny*S1y + Lnz*S1z;
	 LS2 = Lnx*S2x + Lny*S2y + Lnz*S2z;
	 S1S2 = S1x*S2x + S1y*S2y + S1z*S2z;  computed in the constructor*/
	
    beta = (1.0/12.0)*( x1*LS1*(113.0*m1M2 + 75.0*eta) +  x2*LS2*(113.0*m2M2 + 75.0*eta) );
    sigma = (1./48.0)*eta*x1*x2*( -247.0*S1S2 + 721.0*LS1*LS2 );
    
    
    tau = 0.2*eta*(Tc-t)/M;
    tau38 = pow(tau, -0.375);
    tau14 = pow(tau, -0.25);
	tau58 = pow(tau, 0.625);
    
	//   Cout << "stas check: " << beta << "  " << sigma << "   " << tau << Endl;
    
	xPhase = phic -(tau58/eta)*(1. + (3715./8064.+55./96.*eta)*tau14 - (3./16.)*(4.*M_PI - beta)*tau38 + (9275495./14450688.+284875./258048.*eta+1855./2048.*eta*eta-15.*sigma/64.)*tau14*tau14);
    
    return(xPhase);
	
}


double LCGWSpinBBHHHarm1::Integrator(double* ystart, double* dydx, double* yend, double tstart, double dt, int n){
	
	int i;
	double newstep;
	double aa[7];
	double bb[7][6];
	double cc[7];
	double dc[7];
	for (i=0; i<7; i++){
		aa[i] = 0.;
		cc[i] = 0.;
		dc[i] = 0.;
		for (int ii=0; ii<6; ii++){
			bb[i][ii] = 0.0;
		}
	}
	
	aa[2]=0.2; aa[3]=0.3; aa[4]=0.6; aa[5]=1.0; aa[6]=0.875;
	bb[2][1]=0.2; bb[3][1]= 3.0/40.0; bb[3][2]= 9.0/40.0;
	bb[4][1]=0.3; bb[4][2] = -0.9; bb[4][3] = 1.2;
	bb[5][1]=-11.0/54.0; bb[5][2]=2.5; bb[5][3]=-70.0/27.0; bb[5][4]=35.0/27.0;
	bb[6][1]=1631.0/55296.0; bb[6][2]=175.0/512.0; bb[6][3]=575.0/13824.0;
	bb[6][4]=44275.0/110592.0; bb[6][5]=253.0/4096.0;
	cc[1]=37.0/378.0; cc[3]=250.0/621.0; cc[4]=125.0/594.0; cc[6]=512.0/1771.0;
	dc[5]=-277.0/14336.0; dc[1]=cc[1]-2825.0/27648.0; dc[3]=cc[3]-18575.0/48384.0;
	dc[4]=cc[4]-13525.0/55296.0; dc[6]=cc[6]-0.25;
	
	
	double ak2[n];
	double ak3[n];
	double ak4[n];
	double ak5[n];
	double ak6[n];
	double ytemp[n];
	double yerr[n];
	
	for(i=0; i<n; i++){                                     // first step
		ytemp[i] = ystart[i] + bb[2][1]*dt*dydx[i];
	}
	Derivs(tstart+aa[2]*dt, ytemp, ak2, n);
	for(i=0; i<n; i++){                                   // second step
		ytemp[i] = ystart[i] + dt*( bb[3][1]*dydx[i] + bb[3][2]*ak2[i] );
	}
	Derivs(tstart+aa[3]*dt, ytemp, ak3, n);
	for(i=0; i<n; i++){                                   // third step
		ytemp[i] = ystart[i] + dt*( bb[4][1]*dydx[i] + bb[4][2]*ak2[i] + bb[4][3]*ak3[i] );
	}
	Derivs(tstart+aa[4]*dt, ytemp, ak4, n);
	for(i=0; i<n; i++){                                   // fourth step
		ytemp[i] = ystart[i] + dt*( bb[5][1]*dydx[i] + bb[5][2]*ak2[i] + bb[5][3]*ak3[i] +
								   bb[5][4]*ak4[i] );
	}
	Derivs(tstart+aa[5]*dt, ytemp, ak5, n);
	for(i=0; i<n; i++){                                   // fifth step
		ytemp[i] = ystart[i] + dt*( bb[6][1]*dydx[i] + bb[6][2]*ak2[i] + bb[6][3]*ak3[i] +
								   bb[6][4]*ak4[i] + bb[6][5]*ak5[i] );
	}
	Derivs(tstart+aa[6]*dt, ytemp, ak6, n);
	for(i=0; i<n; i++){                                   // sixth step
		yend[i] = ystart[i] + dt*( cc[1]*dydx[i] + cc[3]*ak3[i] + cc[4]*ak4[i] + cc[6]*ak6[i] );
	}
	for(i=0; i<n; i++)                                    // estimate error
		yerr[i] = dt*( dc[1]*dydx[i] + dc[3]*ak3[i] + dc[4]*ak4[i] + dc[5]*ak5[i] + dc[6]*ak6[i] );
	
	
	//! Compute improved step size as a suggestion.             
	double yscal;                  
	for(i=0; i<n; i++){
		yscal = fabs(ystart[i]) + fabs(dt*dydx[i]) + 1.0e-30;
		yerr[i] = yerr[i]/yscal;
	}
	double tollerance = 1.0e-6;
	double maxerr = fabs(yerr[0]);
	
	for(i=1; i<n; i++){
		if(fabs(yerr[i])>maxerr)
            maxerr = fabs(yerr[i]);
	}
	newstep = dt;
	double ratio = tollerance/maxerr;
	if (ratio < 1.0)
		newstep = 0.9*dt*pow(ratio,0.25);
	else
		newstep = 0.9*dt*pow(ratio,0.2);
	
	
	return(newstep);
}




//   Old Derivs by Antoine

void LCGWSpinBBHHHarm1::Derivs(double x, double* y, double* dydx, int n){
	
	double Mom, Mom13, Mom23, Mom2, etaByMom13, Mom2ByMBy2;	
	double SpinOrb, SpinSpin, PN0, PN2, PN3, PN4, PN5, PN6, PN7, LS1, LS2, S1S2;
	double VecS1x, VecS1y, VecS1z, VecS2x, VecS2y, VecS2z, VecLx, VecLy, VecLz;
	
	if(CheckData){
		Cout << "y[i]:";
		for(int i=0; i<13; i++)
			Cout << " " << y[i] ;
		Cout << Endl;
	}
	
	
	om = y[0];
	iota = y[1];
	alpha = y[2];
	Lnx = y[3];
	Lny = y[4];
	Lnz = y[5];
	S1x = y[6];
	S1y = y[7];
	S1z = y[8];
	S2x = y[9];
	S2y = y[10];
	S2z = y[11];
	phi = y[12];
	
	Mom = om*M;
	Mom13 = pow(Mom, 1./3.);
	Mom23 = Mom13*Mom13;
	Mom2  = Mom*Mom;
	etaByMom13 = eta/Mom13;
	Mom2ByMBy2 = 0.5*Mom*Mom/M;
	
	
	LS1 = Lnx*S1x + Lny*S1y + Lnz*S1z;
	LS2 = Lnx*S2x + Lny*S2y + Lnz*S2z;
	S1S2 = S1x*S2x + S1y*S2y + S1z*S2z;
	double Sl = x1*m1*m1*LS1 + x2*m2*m2*LS2;
	double Sigl = M*(x2*m2*LS2 - x1*m1*LS1);
	
	if(nonspin){
		SpinOrb = 0.0;
		SpinSpin = 0.0;
	}else{
		SpinOrb = -(1.0/12.0)*( x1*LS1*(113.0*m1*m1/(M*M) + 75.0*eta) + \
							   x2*LS2*(113.0*m2*m2/(M*M) + 75.0*eta) )*Mom + \
		( (-31811./1008. + 5039.*eta/84.)*Sl + (-473./84. + 1231./56.*eta)*dm/M*Sigl )*Mom*Mom13*Mom13/(M*M);
		SpinSpin = -(1./48.0)*eta*x1*x2*( 247.0*S1S2 - 721.0*LS1*LS2 )*Mom*Mom13;
	}
	
	PN0 = (96.0/5.0)*(eta/(M*M))*Mom23*Mom2*Mom;
	
	PN2 = -( (743.0 + 924.0*eta)/336.0 )*Mom13*Mom13;
	
	PN3 = 4.0*M_PI*Mom;
	
	PN4 = ( 34103.0/18144.0 + (13661.0/2016.0)*eta + (59.0/18.0)*eta*eta )*Mom*Mom13;
	
	PN5 = -(1.0/672.0)*( 4159.0 + 15876.*eta )*M_PI*Mom*Mom13*Mom13;
	
	PN6 = ( (16447322263.0/139708800.0 - (1712.0/105.0)*LC::PN_GAMMA + (16.0/3.0)*M_PI*M_PI) + \
		   (-273811877.0/1088640.0 + (451.0/48.0)*M_PI*M_PI - (88.0/3.0)*(1039.0/4620.0))*eta + \
		   (541.0/896.0)*eta*eta - (5605.0/2592.0)*eta*eta*eta - \
		   (856.0/105.0)*log(16.0*Mom13*Mom13) )*Mom*Mom;
	
	PN7 = ( -4415.0/4032.0 + 358675./6048.*eta + 91495./1512.*eta*eta) * M_PI*Mom13*Mom*Mom;
	//(661775.0/12096.0)*eta + (149789.0/3024.0)*eta*eta )*M_PI*Mom13*Mom*Mom;
	
	dydx[0] = PN0*( 1.0 + PN2 + PN3 + SpinOrb + PN4 + SpinSpin + PN5 + PN6 + PN7 ); //frequency (omega)
	//dydx[0] = PN0*( 1.0 + PN2 + PN3 + SpinOrb + PN4 + SpinSpin + PN5 ); //frequency (omega)
	
	if (!nonspin){  
		//! Modulation of S1
        VecS1x = Mom2ByMBy2 * ( (etaByMom13)*(4.0 + 3.0*m2/m1)*Lnx + x2*m2M2*(S2x - 3.0*LS2*Lnx) );
		
        VecS1y = Mom2ByMBy2 * ( (etaByMom13)*(4.0 + 3.0*m2/m1)*Lny + x2*m2M2*(S2y - 3.0*LS2*Lny) );
		
        VecS1z = Mom2ByMBy2 * ( (etaByMom13)*(4.0 + 3.0*m2/m1)*Lnz + x2*m2M2*(S2z - 3.0*LS2*Lnz) );
		
		//! Modulation of S2
		
        VecS2x = Mom*Mom/(2.0*M) * ( (etaByMom13)*(4.0 + 3.0*m1/m2)*Lnx + x1*m1M2*(S1x - 3.0*LS1*Lnx) );
		
        VecS2y = Mom*Mom/(2.0*M) * ( (etaByMom13)*(4.0 + 3.0*m1/m2)*Lny + x1*m1M2*(S1y - 3.0*LS1*Lny) );
		
        VecS2z = Mom*Mom/(2.0*M) * ( (etaByMom13)*(4.0 + 3.0*m1/m2)*Lnz + x1*m1M2*(S1z - 3.0*LS1*Lnz) );
		
		//! Modulation of Ln
		
        VecLx = Mom*Mom/(2.0*M)* ( (4.0 + 3.0*m2/m1)*x1*m1M2*S1x + (4.0 + 3.0*m1/m2)*x2*m2M2*S2x  -\
								  3.0*Mom13*eta*x1*x2*( LS2*S1x + LS1*S2x) );
		
        VecLy = Mom*Mom/(2.0*M)* ( (4.0 + 3.0*m2/m1)*x1*m1M2*S1y + (4.0 + 3.0*m1/m2)*x2*m2M2*S2y  -\
								  3.0*Mom13*eta*x1*x2*(LS2*S1y + LS1*S2y) );
		
        VecLz = Mom*Mom/(2.0*M)* ( (4.0 + 3.0*m2/m1)*x1*m1M2*S1z + (4.0 + 3.0*m1/m2)*x2*m2M2*S2z  -\
								  3.0*Mom13*eta*x1*x2*(LS2*S1z + LS1*S2z) );
		
        dydx[1] = VecLy*cos(alpha) - VecLx*sin(alpha); // iota
		
        dydx[2] = VecLz - (cos(iota)/sin(iota))*( VecLx*cos(alpha) + VecLy*sin(alpha) ); // alpha
		
        //! the order see above (the same as y[i])
        dydx[3] = VecLy*Lnz - VecLz*Lny;
        dydx[4] = VecLz*Lnx - VecLx*Lnz;
        dydx[5] = VecLx*Lny - VecLy*Lnx;
		
        dydx[6] = S1z*VecS1y - S1y*VecS1z;
        dydx[7] = S1x*VecS1z - S1z*VecS1x;
        dydx[8] = S1y*VecS1x - S1x*VecS1y;
		
        dydx[9] = S2z*VecS2y - S2y*VecS2z;
        dydx[10] = S2x*VecS2z - S2z*VecS2x;
        dydx[11] = S2y*VecS2x - S2x*VecS2y;
		
        dydx[12] = om - dydx[2]*cos(iota); //phase
		
	} else{
		
        dydx[1] = 0.0;
        dydx[2] = 0.0;
        dydx[3] = 0.0;
        dydx[4] = 0.0;
        dydx[5] = 0.0;
        dydx[6] = 0.0;
        dydx[7] = 0.0;
        dydx[8] = 0.0;
        dydx[9] = 0.0;
        dydx[10] = 0.0;
        dydx[11] = 0.0;
        dydx[12] = om;
	}
	
	if(back)
		for(int i=0; i<13; i++)
			dydx[i] *= (-1.0);
	
	if(CheckData){
		Cout << "dydx[i]:";
		for(int i=0; i<13; i++)
			Cout << " " << dydx[i] ;
		Cout << Endl;
	}
	
	
}


void LCGWSpinBBHHHarm1::Computeh0PN(double& h0PNp, double& h0PNc){
	
	double Psi2 = 2.*Psi;
	double fact = (-1.5 - 0.5*c2th);
	h0PNp = fact*pow(ciby2, 4.)*cos(2.*alpha + Psi2) - 2.*pow(ciby2, 3.)*s2th*siby2*cos(alpha + Psi2) \
	+ 2.*ciby2*s2th*pow(siby2, 3.)*cos(alpha - Psi2) + fact*pow(siby2, 4.)*cos(2.*alpha-Psi2) -\
	1.5*sth*sth*si*si*cos(Psi2);
	
	h0PNc = 4.*ciby2*sth*pow(siby2,3.)*sin(alpha - Psi2) - 2.*cth*pow(siby2, 4.)*sin(2.*alpha - Psi2) - \
	4.*pow(ciby2, 3.)*sth*siby2*sin(alpha + Psi2) - 2.*cth*pow(ciby2, 4.)*sin(2.*alpha + Psi2); 
	/*
	 //Check: (Kidder)
	 double sa = sin(alpha);
	 double ca = cos(alpha);
	 double C = 0.5*cth*cth*( sa*sa - ci*ci*ca*ca ) + 0.5*( ci*ci*sa*sa - ca*ca ) - 0.5*sth*sth*si*si - \
	 0.25*s2th*sin(2.*iota)*ca;
	 double S = 0.5*(1.0 + cth*cth)*ci*sin(2.*alpha) + 0.5*s2th*si*sa;
	 h0PNp = 2.*(C*cos(Psi2) + S*sin(Psi2));
	 C = -0.5*cth*sin(2.*alpha)*(1.+ci*ci) - 0.5*sth*sin(2.*iota)*sa;
	 S = -cth*ci*cos(2.*alpha) - sth*si*ca;
	 
	 h0PNc = 2.*(C*cos(Psi2) + S*sin(Psi2));
	 */
	
}

void LCGWSpinBBHHHarm1::Computeh05PN(double& h05PNp, double& h05PNc)
{
	
	double Psi3 = 3.*Psi;
	double alph2 = 2.*alpha;
	double alph3 = 3.*alpha;
	double del = dm/M;
	double c2i = cos(2.*iota);
	double s2i = sin(2.*iota); 
	double s3i = sin(3.*iota);
	double c3i = cos(3.*iota);
	
	h05PNp = pow(ciby2, 6.)*cos(alph3 + Psi3)*(-45.*sth - 9.*s3th)/32. + ciby2*ciby2*cos(alpha + Psi)*(-175.*sth +\
																									   4.*ci*(87.*sth - 5.*s3th) + c2i*(-5.*sth + 15.*s3th) + 13.*s3th )/256. + siby2*siby2*cos(alpha - Psi)*(175.*sth +\
																																																			  4.*ci*(87.*sth - 5.*s3th) + c2i*(5.*sth - 15.*s3th) - 13.*s3th )/256. + pow(ciby2, 4.)*cos(alph3 + Psi)*\
	siby2*siby2*(-5.*sth - s3th)/32. + pow(ciby2, 4.)*cos(alpha + Psi3)*siby2*siby2*(-45.*sth + 135.*s3th)/32. +\
	ciby2*ciby2*pow(siby2,4.)*cos(alpha - Psi3)*(45.*sth - 135.*s3th)/32. + ciby2*ciby2*pow(siby2, 4.)*cos(alph3 - Psi)*\
	(5.*sth + s3th)/32. + sth*pow(siby2, 6.)*cos(alph3 - Psi3)*(27. + 9.*c2th)/16. + cth*sth*sth*pow(si,3.)*cos(Psi3)*\
	45./16. + cos(alph2 + Psi)*cth*( si*(-85./256. - c2th*(1. + 4.*ci + 3.*c2i)/128.) - 11.*s2i/64. - s3i/256. ) + \
	cos(alph2 + Psi3)*cth*( si*(45./256. + c2th*(81. + 108.*ci + 27.*c2i)/128.) + 9.*s2i/64. + 9.*s3i/256. ) + \
	cos(alph2 - Psi)*cth*(si*(-85. + c2th)/256.  + s2i*(11. + c2th)/64. + s3i*(-1. - 3.*c2th)/256. ) + \
	cos(alph2 - Psi3)*cth*(si*(45. + 135.*c2th)/256. + s2i*(-9. - 27.*c2th)/64. + s3i*(9. + 27.*c2th)/256. ) + \
	cos(Psi)*cth*sth*sth*(si + 5.*s3i)/64.;
	
	h05PNp *= del;
	
	h05PNc = -5.625*ciby2*ciby2*s2th*pow(siby2, 4.)*sin(alpha - Psi3) + 4.5*c2th*ciby2*pow(siby2, 5.)*sin(alph2 - Psi3) +\
	1.125*s2th*pow(siby2, 6.)*sin(alph3 - Psi3) + s2th*sin(alpha - Psi)*(-2. + 43.*ci - 46.*c2i + 5.*c3i)/256. +\
	(-1. - 0.25*c2th*(1. - ci))*ciby2*pow(siby2, 3.)*sin(alph2 - Psi) + 0.125*ciby2*ciby2*s2th*pow(siby2, 4.)*sin(alph3 - Psi)\
	+ 0.5*sth*sth*s2i*sin(Psi) + sin(alpha + Psi)*s2th*(2. + 43.*ci + 46.*c2i + 5.*c3i)/256. +\
	(-1. - 0.25*c2th*(1. + ci))*pow(ciby2, 3.)*siby2*sin(alph2 + Psi) - 0.125*pow(ciby2, 4.)*s2th*siby2*siby2*sin(alph3 + Psi)\
	+ 5.625*pow(ciby2, 4.)*s2th*siby2*siby2*sin(alpha + Psi3) + 4.5*c2th*pow(ciby2, 5.)*siby2*sin(alph2 + Psi3) -\
	1.125*pow(ciby2, 6.)*s2th*sin(alph3 + Psi3);
	
	h05PNc *= del;
	
	// Check (Kidder):
	
	/*   double sa = sin(alpha);
	 double ca = cos(alpha);
	 double  Ka = -sth*sa;
	 double  Kb = cth*si - sth*ci*ca;
	 double C = 0.5*cth*cth*( sa*sa - ci*ci*ca*ca ) + 0.5*( ci*ci*sa*sa - ca*ca ) - 0.5*sth*sth*si*si - \
	 0.25*s2th*sin(2.*iota)*ca;
	 double S = 0.5*(1.0 + cth*cth)*ci*sin(2.*alpha) + 0.5*s2th*si*sa;
	 double K = 0.5*cth*cth*(sa*sa + ci*ci*ca*ca) - 0.5*( ci*ci*sa*sa + ca*ca ) + 0.5*sth*sth*si*si + \
	 0.25*s2th*s2i*ca;
	 
	 double P05p = 0.25*(dm/M)*( 9.0*(Ka*S + Kb*C)*cos(Psi3) + 9.0*(Kb*S - Ka*C)*sin(Psi3) +\
	 (3.0*Ka*S - 3.0*Kb*C - 2.0*Kb*K)*cos(Psi) - (3.0*Kb*S + 3.0*Ka*C - 2.0*Ka*K)*sin(Psi) );
	 h05PNp = -P05p;
	 
	 K = -0.5*cth*sin(alph2)*si*si + 0.5*sth*sin(2.*iota)*sa;
	 C = -0.5*cth*sin(2.*alpha)*(1.+ci*ci) - 0.5*sth*sin(2.*iota)*sa;
	 S = -cth*ci*cos(2.*alpha) - sth*si*ca;
	 
	 double P05c = 0.25*(dm/M)*( 9.0*(Ka*S + Kb*C)*cos(Psi3) + 9.0*(Kb*S - Ka*C)*sin(Psi3) +\
	 (3.0*Ka*S - 3.0*Kb*C - 2.0*Kb*K)*cos(Psi) - (3.0*Kb*S + 3.0*Ka*C - 2.0*Ka*K)*sin(Psi) );
	 h05PNc = -P05c;
	 */
	
	
}

void LCGWSpinBBHHHarm1::Computeh1PN(double& h1PNp, double& h1PNc){
	
	double alph2 = 2.*alpha;
	double Psi2 = 2.*Psi;
	double alph3 = 3.*alpha;
	double alph4 = 4.*alpha;
	double Psi4 = 4.*Psi;
	double c2i = cos(2.*iota);
	double s2i = sin(2.*iota);
	double ci3by2 = cos(1.5*iota);
	double ci5by2 = cos(2.5*iota);
	double si3by2 = sin(1.5*iota);
	double si5by2 = sin(2.5*iota); 
	double fct1 = (1./3. - eta);
	
	
	
	h1PNp = cos(alph2 + Psi2)*pow(ciby2, 4.)*( (59. + 40.*c2th - 3.*c4th + 2.*(5. - 44.*c2th + 7.*c4th)*ci/3. + \
												(-5. - 4.*c2th - 7.*c4th)*c2i/3.) + eta*(-25. - 208.*c2th/3. + 9.*c4th + (-10. + 88.*c2th - 14.*c4th)*ci +\
																						 (5. + 4.*c2th + 7.*c4th)*c2i)  )/16. + cos(alph4 + Psi4)*pow(ciby2, 8.)*sth*sth*(-fct1)*(6.+2.*c2th) +\
	cos(alph3 + Psi4)*pow(cth, 3.)*pow(ciby2, 7.)*sth*siby2*32.*fct1 + cos(alph3 + Psi2)*pow(ciby2, 5.)*s2th*\
	siby2*(-fct1)*0.5*(5. - c2th + 4.*cth*cth*ci) + cos(alph2 + Psi4)*pow(ciby2, 6.)*siby2*siby2*(-fct1)*\
	(10. + 8.*c2th + 14.*c4th) + cos(alph4 + Psi2)*pow(ciby2, 6.)*sth*sth*siby2*siby2*(-fct1)*0.5*(3.+c2th) +\
	cos(alpha + Psi4)*pow(ciby2, 5.)*s2th*pow(siby2, 3.)*(8.-56.*c2th)*(fct1) + cos(alpha - Psi2)*pow(siby2, 3.)*\
	( -ciby2*(72. + 31.*ci + c2i)*s2th/12. + s4th*(19.*ciby2 + 14.*ci3by2 + 7.*ci5by2)*(fct1)/16. + eta*(ciby2*(64./3. +\
																												31.*ci + c2i)*s2th/4.) ) + cos(alph2 - Psi2)*pow(siby2, 4.)*( (59.+ 40.*c2th - 3.*c4th)/16.  + \
																																											 eta*(-25. - 52./3.*c2th + 9.*c4th)/16. +  (fct1)*((-5. + 44.*c2th - 7.*c4th)*ci/8. +(-5.-4.*c2th-7.*c4th)*c2i/16.) )\
	+ cos(alpha - Psi4)*pow(ciby2, 3.)*s2th*pow(siby2, 5.)*(-fct1)*(8.-56.*c2th) + cos(alph3 - Psi2)*pow(siby2, 5.)*s2th*\
	ciby2*(fct1)*(5. - c2th - 4.*cth*cth*ci)*0.5 + cos(alph2 - Psi4)*ciby2*ciby2*pow(siby2, 6.)*(-fct1)*\
	(10. + 8*c2th + 14.*c4th) + cos(alph4 - Psi2)*ciby2*ciby2*sth*sth*pow(siby2, 6.)*(-fct1)*(3.+c2th)*0.5 + \
	cos(alph3-Psi4)*32.*ciby2*sth*pow(cth, 3.)*pow(siby2, 7.)*(-fct1) + cos(alph4 - Psi4)*sth*sth*pow(siby2, 8.)*\
	(-fct1)*(6.+2.*c2th) + cos(Psi2)*sth*sth*si*si*( (349. - eta*135.)/96.  + (-fct1)*(25.*c2th + \
																					   (25.+ 35.*c2th)*c2i)/32.) + cos(Psi4)*sth*sth*pow(si, 4)*(-fct1)*0.25*(25.+35.*c2th) + cos(alpha+Psi2)*pow(ciby2, 3.)*\
	(s2th*siby2*( 6. - 16.*eta/3. + (-fct1)*(31.*ci-c2i)*0.25 ) + s4th*(19.*siby2 - 14*si3by2 + 7.*si5by2)*(-fct1)/16. );
	
	double chia_x = 0.5*(x1*S1x - x2*S2x);
	double chia_y = 0.5*(x1*S1y - x2*S2y);
	double chia_z = 0.5*(x1*S1z - x2*S2z);
	double chis_x = 0.5*(x1*S1x + x2*S2x);
	double chis_y = 0.5*(x1*S1y + x2*S2y);
	double chis_z = 0.5*(x1*S1z + x2*S2z);
	double fct2 = chia_x*cth - chia_z*sth;
	double fct3 = chis_x*cth - chis_z*sth;
	
	h1PNp = h1PNp + cos(alpha + Psi)*ciby2*ciby2*fct2  + cos(alpha - Psi)*(chia_x*0.5*cth*(1.-ci) - \
																		   chia_z*sth*siby2*siby2) - chia_y*( (sin(alpha - Psi)*siby2*siby2 + sin(alpha + Psi)*ciby2*ciby2)*cth \
																											 + sin(Psi)*sth*si ) + (dm/M)*(cos(alpha + Psi)*ciby2*ciby2*fct3 + cos(alpha - Psi)*\
																																		   (cth*0.5*chis_x*(1.-ci) - chis_z*sth*siby2*siby2) - chis_y*(cth*(siby2*siby2*sin(alpha-Psi) + \
																																																			ciby2*ciby2*sin(alpha+Psi)) + sth*si*sin(Psi)) );       
	
	h1PNc = sin(alpha - Psi4)*pow(ciby2, 3.)*pow(siby2, 5.)*(12.*sth + 28.*s3th)*fct1 + sin(alph2 - Psi4)*ciby2*ciby2*\
	pow(siby2, 6.)*(4.*cth + 28.*c3th)*(-fct1) + sin(alph3 - Psi4)*ciby2*pow(siby2, 7.)*(4.*sth - 12.*s3th)*fct1 \
	+ sin(alph4 - Psi4)*8.*cth*sth*sth*pow(siby2, 8.)*(-fct1)  + sin(alpha - Psi2)*ciby2*pow(siby2, 3.)*\
	( (-79./8. + eta*103./24.)*sth + (ci*(9.*sth - 19.*s3th)*0.25 + c2i*(3.*sth + 7.*s3th)*0.125 - 9.*s3th*0.125)*fct1 ) +\
	sin(alph2 - Psi2)*pow(siby2, 4.)*( cth*(47.*0.125  - 119.*eta/24.) + (s3th*3.*0.125 + ci*0.5*(7.*cth + c3th) - \
																		  c2i*0.125*(cth + 7.*c3th))*fct1 ) + sin(alph3 - Psi2)*pow(siby2,5.)*ciby2*sth*( 4. - (1.+3.*c2th)*ci )*fct1\
	+ sin(alph4 - Psi2)*2.*cth*ciby2*ciby2*sth*sth*pow(siby2, 6.)*(-fct1) + 7.5*sin(Psi2)*cth*ci*sth*sth*si*si*(-fct1)\
	+ sin(alpha+Psi2)*siby2*pow(ciby2, 3.)*( sth*(237. - 103.*eta)/24. + (ci*(9.*sth -19.*s3th)*0.25 - \
																		  c2i*(3.*sth + 7.*s3th)*0.125)*fct1 ) + sin(alph2+Psi2)*pow(ciby2, 4.)*(cth*(141.- eta*119.)/24. + (3.*c3th*0.125 -\
																																											 ci*(7.*cth + c3th)*0.5 - c2i*(cth + 7.*c3th)*0.125)*fct1 ) + sin(alph3+Psi2)*pow(ciby2, 5.)*sth*siby2*(4. + \
																																																																					(1.+3.*c2th)*ci)*(-fct1) + sin(alph4+Psi2)*2.*cth*pow(ciby2,6.)*sth*sth*siby2*siby2*(-fct1) + \
	sin(alpha+Psi4)*pow(ciby2, 5.)*pow(siby2,3.)*(12.*sth + 28.*s3th)*(-fct1) + \
	sin(alph2 + Psi4)*pow(ciby2, 6.)*siby2*siby2*(4.*cth + 28.*c3th)*(-fct1) +\
	sin(alph3+Psi4)*pow(ciby2, 7.)*sth*siby2*(8.+24.*c2th)*fct1 + 8.*sin(alph4+Psi4)*cth*pow(ciby2, 8.)*sth*sth*(-fct1);
	
	
	
	h1PNc = h1PNc + chia_y*0.5*(1.+ci)*cos(alpha+Psi) + chia_y*siby2*siby2*cos(alpha-Psi) + sin(alpha-Psi)*0.5*cth*fct2*(1.-ci) \
	+ sin(Psi)*sth*si*fct2 + sin(alpha+Psi)*0.5*cth*fct2*(1.+ci) + (dm/M)*( 0.5*cos(alpha+Psi)*chis_y*(1.+ci) + \
																		   cos(alpha-Psi)*chis_y*siby2*siby2 + sin(alpha-Psi)*0.5*cth*fct3*(1.-ci) + sin(Psi)*sth*si*fct3 +\
																		   sin(alpha+Psi)*cth*0.5*fct3*(1.+ci) );
	/*         
	 // Check:(Kidder)
	 double sa = sin(alpha);
	 double ca = cos(alpha);
	 double  Ka = -sth*sa;
	 double  Kb = cth*si - sth*ci*ca;
	 double  Kc = cth*ci*ca + si*sth;
	 double  Dx = x2*(m2/M)*S2x - x1*(m1/M)*S1x; 
	 double  Dy = x2*(m2/M)*S2y - x1*(m1/M)*S1y; 
	 double  Dz = x2*(m2/M)*S2z - x1*(m1/M)*S1z;
	 double d = Dz*sth - Dx*cth;
	 
	 double C = 0.5*cth*cth*( sa*sa - ci*ci*ca*ca ) + 0.5*( ci*ci*sa*sa - ca*ca ) - 0.5*sth*sth*si*si - \
	 0.25*s2th*sin(2.*iota)*ca;
	 double S = 0.5*(1.0 + cth*cth)*ci*sin(2.*alpha) + 0.5*s2th*si*sa;
	 double K = 0.5*cth*cth*(sa*sa + ci*ci*ca*ca) - 0.5*( ci*ci*sa*sa + ca*ca ) + 0.5*sth*sth*si*si + \
	 0.25*s2th*s2i*ca;
	 
	 double DC = -( Dy*sa*cth + d*ca);
	 double DS = -( Kc*Dy - d*ci*sa );
	 
	 double Q = 2.*(C*cos(Psi2) + S*sin(Psi2));
	 double P1p =  (8./3.)*(1.-3.*eta)*( ((Ka*Ka - Kb*Kb)*C - 2.*Ka*Kb*S)*cos(Psi4) + \
	 ((Ka*Ka - Kb*Kb)*S + 2.*Ka*Kb*C)*sin(Psi4)  ) + \
	 (1./6.)*( (4.*(1.- 3.*eta)*(Kb*Kb + Ka*Ka) + (13. - eta))*Q - 4.*(1.- 3.*eta)*\
	 ((Ka*Ka - Kb*Kb)*cos(Psi2) + 2.*Ka*Kb*sin(Psi2))*K );
	 
	 h1PNp = -P1p;             
     
	 K = -0.5*cth*sin(alph2)*si*si + 0.5*sth*sin(2.*iota)*sa;
	 C = -0.5*cth*sin(2.*alpha)*(1.+ci*ci) - 0.5*sth*sin(2.*iota)*sa;
	 S = -cth*ci*cos(2.*alpha) - sth*si*ca;    
	 
	 DC = Dy*ca - d*cth*sa;
	 DS = -Dy*ci*sa - Kc*d;
	 Q = 2.*(C*cos(Psi2) + S*sin(Psi2));
	 
	 double P1c =  (8./3.)*(1.-3.*eta)*( ((Ka*Ka - Kb*Kb)*C - 2.*Ka*Kb*S)*cos(Psi4) + \
	 ((Ka*Ka - Kb*Kb)*S + 2.*Ka*Kb*C)*sin(Psi4)  ) + \
	 (1./6.)*( (4.*(1.- 3.*eta)*(Kb*Kb + Ka*Ka) + (13. - eta))*Q - 4.*(1.- 3.*eta)*\
	 ((Ka*Ka - Kb*Kb)*cos(Psi2) + 2.*Ka*Kb*sin(Psi2))*K );
     
	 h1PNc = -P1c; 
	 */           
	
}


void LCGWSpinBBHHHarm1::ComputehAllPN(double& h0PNp, double& h0PNc, double& h05PNp, double& h05PNc, double& h1PNp, double& h1PNc)
{
	
	double Psi2 = 2.*Psi;
	double fact = (-1.5 - 0.5*c2th);
	double Psi3 = 3.*Psi;
	double alph2 = 2.*alpha;
	double alph3 = 3.*alpha;
	double del = dm/M;
	double c2i = cos(2.*iota);
	double s2i = sin(2.*iota); 
	double s3i = sin(3.*iota);
	double c3i = cos(3.*iota);
	double alph4 = 4.*alpha;
	double Psi4 = 4.*Psi;
	double ci3by2 = cos(1.5*iota);
	double ci5by2 = cos(2.5*iota);
	double si3by2 = sin(1.5*iota);
	double si5by2 = sin(2.5*iota); 
	double fct1 = (1./3. - eta);
	
	double chia_x = 0.5*(x1*S1x - x2*S2x);
	double chia_y = 0.5*(x1*S1y - x2*S2y);
	double chia_z = 0.5*(x1*S1z - x2*S2z);
	double chis_x = 0.5*(x1*S1x + x2*S2x);
	double chis_y = 0.5*(x1*S1y + x2*S2y);
	double chis_z = 0.5*(x1*S1z + x2*S2z);
	double fct2 = chia_x*cth - chia_z*sth;
	double fct3 = chis_x*cth - chis_z*sth;
	
	double ciby2p2 = ciby2*ciby2;
	double ciby2p3 = ciby2p2*ciby2;
	double ciby2p4 = ciby2p3*ciby2;
	double ciby2p5 = ciby2p4*ciby2;
	double ciby2p6 = ciby2p5*ciby2;
	double ciby2p7 = ciby2p6*ciby2;
	double ciby2p8 = ciby2p7*ciby2;
	
	double siby2p2 = siby2*siby2;
	double siby2p3 = siby2p2*siby2;
	double siby2p4 = siby2p3*siby2;
	double siby2p5 = siby2p4*siby2;
	double siby2p6 = siby2p5*siby2;
	double siby2p7 = siby2p6*siby2;
	double siby2p8 = siby2p7*siby2;
	
	
	
	/*
	 double cPs = cos(Psi);
	 double cPs2 = cos(Psi2);
	 double cPs3 = cos(Psi3);
	 double cPs4 = cos(Psi4);
	 double sPs = sin(Psi);
	 double sPs2 = sin(Psi2);
	 */ 
	
	double cPs = cos(Psi);
	double sPs = sin(Psi);
	double cPs2 = 2.*cPs*cPs - 1.;
	double cPs3 = cPs*(4.*cPs*cPs - 3.);
	double cPs4 = 2.*cPs2*cPs2 - 1.;
	double sPs2 = 2.*sPs*cPs;
	double sPs3 = sPs*(3. - 4.*sPs*sPs);
	double sPs4 = 2.*sPs2*cPs2;
	
	double cal = cos(alpha);
	double sal = sin(alpha);
	double cal2 = 2.*cal*cal - 1.;
	double cal3 = cal*(4.*cal*cal - 3.);
	double cal4 = 2.*cal2*cal2 - 1.;
	double sal2 = 2.*sal*cal;
	double sal3 = sal*(3. - 4.*sal*sal);
	double sal4 = 2.*sal2*cal2;
	
	double calpP   = cal *cPs  - sal *sPs ;
	double calmP   = cal *cPs  + sal *sPs ;
	double calpP2  = cal *cPs2 - sal *sPs2;
	double calmP2  = cal *cPs2 + sal *sPs2;
	double calpP3  = cal *cPs3 - sal *sPs3;
	double calmP3  = cal *cPs3 + sal *sPs3;
	double calpP4  = cal *cPs4 - sal *sPs4;
	double calmP4  = cal *cPs4 + sal *sPs4;
	double cal2pP  = cal2*cPs  - sal2*sPs ;
	double cal2mP  = cal2*cPs  + sal2*sPs ;
	double cal2pP2 = cal2*cPs2 - sal2*sPs2;
	double cal2mP2 = cal2*cPs2 + sal2*sPs2;
	double cal2pP3 = cal2*cPs3 - sal2*sPs3;
	double cal2mP3 = cal2*cPs3 + sal2*sPs3;
	double cal2pP4 = cal2*cPs4 - sal2*sPs4;
	double cal2mP4 = cal2*cPs4 + sal2*sPs4;
	double cal3pP  = cal3*cPs  - sal3*sPs ;
	double cal3mP  = cal3*cPs  + sal3*sPs ;
	double cal3pP2 = cal3*cPs2 - sal3*sPs2;
	double cal3mP2 = cal3*cPs2 + sal3*sPs2;
	double cal3pP3 = cal3*cPs3 - sal3*sPs3;
	double cal3mP3 = cal3*cPs3 + sal3*sPs3;
	double cal3pP4 = cal3*cPs4 - sal3*sPs4;
	double cal3mP4 = cal3*cPs4 + sal3*sPs4;
	double cal4pP2 = cal4*cPs2 - sal4*sPs2;
	double cal4mP2 = cal4*cPs2 + sal4*sPs2;
	double cal4pP4 = cal4*cPs4 - sal4*sPs4;
	double cal4mP4 = cal4*cPs4 + sal4*sPs4;
	
	/*
	 double calpP = cos(alpha + Psi);
	 double calmP = cos(alpha - Psi);
	 double calpP2 = cos(alpha + Psi2);
	 double calmP2 = cos(alpha - Psi2);
	 double calpP3 = cos(alpha + Psi3);
	 double calmP3 = cos(alpha - Psi3);
	 double calpP4 = cos(alpha + Psi4);
	 double calmP4 = cos(alpha - Psi4);
	 double cal2pP = cos(alph2 + Psi);
	 double cal2mP = cos(alph2 - Psi);
	 double cal2pP2 = cos(alph2 + Psi2);
	 double cal2mP2 = cos(alph2 - Psi2);
	 double cal2pP3 = cos(alph2 + Psi3);
	 double cal2mP3 = cos(alph2 - Psi3);
	 double cal2pP4 = cos(alph2 + Psi4);
	 double cal2mP4 = cos(alph2 - Psi4);
	 double cal3mP = cos(alph3 - Psi);
	 double cal3pP = cos(alph3 + Psi);
	 double cal3pP2 = cos(alph3 + Psi2);
	 double cal3mP2 = cos(alph3 - Psi2);
	 double cal3pP3 = cos(alph3 + Psi3);
	 double cal3mP3 = cos(alph3 - Psi3);
	 double cal3pP4 = cos(alph3 + Psi4);
	 double cal3mP4 = cos(alph3 - Psi4);
	 double cal4pP2 = cos(alph4 + Psi2);
	 double cal4mP2 = cos(alph4 - Psi2);
	 double cal4pP4 = cos(alph4 + Psi4);
	 double cal4mP4 = cos(alph4 - Psi4);
	 */
	
	
	double salpP   = sal *cPs  + cal *sPs ;
	double salmP   = sal *cPs  - cal *sPs ;
	double salpP2  = sal *cPs2 + cal *sPs2;
	double salmP2  = sal *cPs2 - cal *sPs2;
	double salpP3  = sal *cPs3 + cal *sPs3;
	double salmP3  = sal *cPs3 - cal *sPs3;
	double salpP4  = sal *cPs4 + cal *sPs4;
	double salmP4  = sal *cPs4 - cal *sPs4;
	double sal2pP  = sal2*cPs  + cal2*sPs ;
	double sal2mP  = sal2*cPs  - cal2*sPs ;
	double sal2pP2 = sal2*cPs2 + cal2*sPs2;
	double sal2mP2 = sal2*cPs2 - cal2*sPs2;
	double sal2pP3 = sal2*cPs3 + cal2*sPs3;
	double sal2mP3 = sal2*cPs3 - cal2*sPs3;
	double sal2pP4 = sal2*cPs4 + cal2*sPs4;
	double sal2mP4 = sal2*cPs4 - cal2*sPs4;
	double sal3pP  = sal3*cPs  + cal3*sPs ;
	double sal3mP  = sal3*cPs  - cal3*sPs ;
	double sal3pP2 = sal3*cPs2 + cal3*sPs2;
	double sal3mP2 = sal3*cPs2 - cal3*sPs2;
	double sal3pP3 = sal3*cPs3 + cal3*sPs3;
	double sal3mP3 = sal3*cPs3 - cal3*sPs3;
	double sal3pP4 = sal3*cPs4 + cal3*sPs4;
	double sal3mP4 = sal3*cPs4 - cal3*sPs4;
	double sal4pP2 = sal4*cPs2 + cal4*sPs2;
	double sal4mP2 = sal4*cPs2 - cal4*sPs2;
	double sal4pP4 = sal4*cPs4 + cal4*sPs4;
	double sal4mP4 = sal4*cPs4 - cal4*sPs4;
	
	/*
	 double salpP = sin(alpha + Psi);
	 double salmP = sin(alpha - Psi);
	 double salpP2 = sin(alpha + Psi2);
	 double salmP2 = sin(alpha - Psi2);
	 double salpP3 = sin(alpha + Psi3);
	 double salmP3 = sin(alpha - Psi3);
	 double salpP4 = sin(alpha + Psi4);
	 double salmP4 = sin(alpha - Psi4);
	 double sal2pP = sin(alph2 + Psi);
	 double sal2mP = sin(alph2 - Psi);
	 double sal2pP2 = sin(alph2 + Psi2);
	 double sal2mP2 = sin(alph2 - Psi2);
	 double sal2pP3 = sin(alph2 + Psi3);
	 double sal2mP3 = sin(alph2 - Psi3);
	 double sal2pP4 = sin(alph2 + Psi4);
	 double sal2mP4 = sin(alph2 - Psi4);
	 double sal3pP = sin(alph3 + Psi);
	 double sal3mP = sin(alph3 - Psi);
	 double sal3pP2 = sin(alph3 + Psi2);
	 double sal3mP2 = sin(alph3 - Psi2);
	 double sal3pP3 = sin(alph3 + Psi3);
	 double sal3mP3 = sin(alph3 - Psi3);
	 double sal3pP4 = sin(alph3 + Psi4);
	 double sal3mP4 = sin(alph3 - Psi4);
	 double sal4pP2 = sin(alph4 + Psi2);
	 double sal4mP2 = sin(alph4 - Psi2);
	 double sal4pP4 = sin(alph4 + Psi4);
	 double sal4mP4 = sin(alph4 - Psi4);
	 */
	
	
	
	
	//! ********** Compute 0 PN
	
	h0PNp = fact*ciby2p4*cal2pP2 - 2.*ciby2p3*s2th*siby2*calpP2 \
	+ 2.*ciby2*s2th*siby2p3*calmP2 + fact*siby2p4*cal2mP2 -\
	1.5*sth2*si*si*cPs2;
	
	h0PNc = 4.*ciby2*sth*siby2p3*salmP2 - 2.*cth*siby2p4*sal2mP2 - \
	4.*ciby2p3*sth*siby2*salpP2 - 2.*cth*ciby2p4*sal2pP2; 
	
	
	
	
	//! ********** Compute 0.5 PN
	
	h05PNp = ciby2p6*cal3pP3*(-45.*sth - 9.*s3th)/32. + ciby2p2*calpP*(-175.*sth +\
																				 4.*ci*(87.*sth - 5.*s3th) + c2i*(-5.*sth + 15.*s3th) + 13.*s3th )/256. + siby2p2*calmP*(175.*sth +\
																																										 4.*ci*(87.*sth - 5.*s3th) + c2i*(5.*sth - 15.*s3th) - 13.*s3th )/256. + ciby2p4*cal3pP*\
	siby2p2*(-5.*sth - s3th)/32. + ciby2p4*calpP3*siby2p2*(-45.*sth + 135.*s3th)/32. +\
	ciby2p2*siby2p4*calmP3*(45.*sth - 135.*s3th)/32. + ciby2p2*siby2p4*cal3mP*\
	(5.*sth + s3th)/32. + sth*siby2p6*cal3mP3*(27. + 9.*c2th)/16. + cth*sth2*pow(si,3.)*cPs3*\
	45./16. + cal2pP*cth*( si*(-85./256. - c2th*(1. + 4.*ci + 3.*c2i)/128.) - 11.*s2i/64. - s3i/256. ) + \
	cal2pP3*cth*( si*(45./256. + c2th*(81. + 108.*ci + 27.*c2i)/128.) + 9.*s2i/64. + 9.*s3i/256. ) + \
	cal2mP*cth*(si*(-85. + c2th)/256.  + s2i*(11. + c2th)/64. + s3i*(-1. - 3.*c2th)/256. ) + \
	cal2mP3*cth*(si*(45. + 135.*c2th)/256. + s2i*(-9. - 27.*c2th)/64. + s3i*(9. + 27.*c2th)/256. ) + \
	cPs*cthsth2*(si + 5.*s3i)/64.;
	
	h05PNp *= del;
	
	h05PNc = -5.625*ciby2p2*s2th*siby2p4*salmP3 + 4.5*c2th*ciby2*siby2p5*sal2mP3 +\
	1.125*s2th*siby2p6*sal3mP3 + s2th*salmP*(-2. + 43.*ci - 46.*c2i + 5.*c3i)/256. +\
	(-1. - 0.25*c2th*(1. - ci))*ciby2*siby2p3*sal2mP + 0.125*ciby2p2*s2th*siby2p4*sal3mP\
	+ 0.5*sth2*s2i*sPs + salpP*s2th*(2. + 43.*ci + 46.*c2i + 5.*c3i)/256. +\
	(-1. - 0.25*c2th*(1. + ci))*ciby2p3*siby2*sal2pP - 0.125*ciby2p4*s2th*siby2p2*sal3pP\
	+ 5.625*ciby2p4*s2th*siby2p2*salpP3 + 4.5*c2th*ciby2p5*siby2*sal2pP3 -\
	1.125*ciby2p6*s2th*sal3pP3;
	
	h05PNc *= del;
	
	
	
	//! ********** Compute 1 PN
	
	h1PNp = cal2pP2*ciby2p4*( (59. + 40.*c2th - 3.*c4th + 2.*(5. - 44.*c2th + 7.*c4th)*ci/3. + \
										 (-5. - 4.*c2th - 7.*c4th)*c2i/3.) + eta*(-25. - 208.*c2th/3. + 9.*c4th + (-10. + 88.*c2th - 14.*c4th)*ci +\
																				  (5. + 4.*c2th + 7.*c4th)*c2i)  )/16. + cal4pP4*ciby2p8*sth2*(-fct1)*(6.+2.*c2th) +\
	cal3pP4*cth3*ciby2p7*sth*siby2*32.*fct1 + cal3pP2*ciby2p5*s2th*\
	siby2*(-fct1)*0.5*(5. - c2th + 4.*cth2*ci) + cal2pP4*ciby2p6*siby2p2*(-fct1)*\
	(10. + 8.*c2th + 14.*c4th) + cal4pP2*ciby2p6*sth2*siby2p2*(-fct1)*0.5*(3.+c2th) +\
	calpP4*ciby2p5*s2th*siby2p3*(8.-56.*c2th)*(fct1) + calmP2*siby2p3*\
	( -ciby2*(72. + 31.*ci + c2i)*s2th/12. + s4th*(19.*ciby2 + 14.*ci3by2 + 7.*ci5by2)*(fct1)/16. + eta*(ciby2*(64./3. +\
																												31.*ci + c2i)*s2th/4.) ) + cal2mP2*siby2p4*( (59.+ 40.*c2th - 3.*c4th)/16.  + \
																																									  eta*(-25. - 52./3.*c2th + 9.*c4th)/16. +  (fct1)*((-5. + 44.*c2th - 7.*c4th)*ci/8. +(-5.-4.*c2th-7.*c4th)*c2i/16.) )\
	+ calmP4*ciby2p3*s2th*siby2p5*(-fct1)*(8.-56.*c2th) + cal3mP2*siby2p5*s2th*\
	ciby2*(fct1)*(5. - c2th - 4.*cth2*ci)*0.5 + cal2mP4*ciby2p2*siby2p6*(-fct1)*\
	(10. + 8*c2th + 14.*c4th) + cal4mP2*ciby2p2*sth2*siby2p6*(-fct1)*(3.+c2th)*0.5 + \
	cal3mP4*32.*ciby2*sthcth3*siby2p7*(-fct1) + cal4mP4*sth2*siby2p8*\
	(-fct1)*(6.+2.*c2th) + cPs2*sth2*si*si*( (349. - eta*135.)/96.  + (-fct1)*(25.*c2th + \
																				  (25.+ 35.*c2th)*c2i)/32.) + cPs4*sth2*pow(si, 4)*(-fct1)*0.25*(25.+35.*c2th) + calpP2*ciby2p3*\
	(s2th*siby2*( 6. - 16.*eta/3. + (-fct1)*(31.*ci-c2i)*0.25 ) + s4th*(19.*siby2 - 14*si3by2 + 7.*si5by2)*(-fct1)/16. );
	
	
	
	h1PNp = h1PNp + calpP*ciby2p2*fct2  + calmP*(chia_x*0.5*cth*(1.-ci) - \
												 chia_z*sth*siby2p2) - chia_y*( (salmP*siby2p2 + salpP*ciby2p2)*cth \
																			   + sPs*sth*si ) + (dm/M)*(calpP*ciby2p2*fct3 + calmP*\
																										(cth*0.5*chis_x*(1.-ci) - chis_z*sth*siby2p2) - chis_y*(cth*(siby2p2*salmP + \
																																									 ciby2p2*salpP) + sth*si*sPs) );       
	
	h1PNc = salmP4*ciby2p3*siby2p5*(12.*sth + 28.*s3th)*fct1 + sal2mP4*ciby2p2*\
	siby2p6*(4.*cth + 28.*c3th)*(-fct1) + sal3mP4*ciby2*siby2p7*(4.*sth - 12.*s3th)*fct1 \
	+ sal4mP4*8.*cth*sth2*siby2p8*(-fct1)  + salmP2*ciby2*siby2p3*\
	( (-79./8. + eta*103./24.)*sth + (ci*(9.*sth - 19.*s3th)*0.25 + c2i*(3.*sth + 7.*s3th)*0.125 - 9.*s3th*0.125)*fct1 ) +\
	sal2mP2*siby2p4*( cth*(47.*0.125  - 119.*eta/24.) + (s3th*3.*0.125 + ci*0.5*(7.*cth + c3th) - \
																   c2i*0.125*(cth + 7.*c3th))*fct1 ) + sal3mP2*siby2p5*ciby2*sth*( 4. - (1.+3.*c2th)*ci )*fct1\
	+ sal4mP2*2.*cth*ciby2p2*sth2*siby2p6*(-fct1) + 7.5*sPs2*cth*ci*sth2*si*si*(-fct1)\
	+ salpP2*siby2*ciby2p3*( sth*(237. - 103.*eta)/24. + (ci*(9.*sth -19.*s3th)*0.25 - \
																   c2i*(3.*sth + 7.*s3th)*0.125)*fct1 ) + sal2pP2*ciby2p4*(cth*(141.- eta*119.)/24. + (3.*c3th*0.125 -\
																																							   ci*(7.*cth + c3th)*0.5 - c2i*(cth + 7.*c3th)*0.125)*fct1 ) + sal3pP2*ciby2p5*sth*siby2*(4. + \
																																																															   (1.+3.*c2th)*ci)*(-fct1) + sal4pP2*2.*cth*pow(ciby2,6.)*sth2*siby2p2*(-fct1) + \
	salpP4*ciby2p5*siby2p3*(12.*sth + 28.*s3th)*(-fct1) + \
	sal2pP4*ciby2p6*siby2p2*(4.*cth + 28.*c3th)*(-fct1) +\
	sal3pP4*ciby2p7*sth*siby2*(8.+24.*c2th)*fct1 + 8.*sal4pP4*cth*ciby2p8*sth2*(-fct1);
	
	
	
	h1PNc = h1PNc + chia_y*0.5*(1.+ci)*calpP + chia_y*siby2p2*calmP + salmP*0.5*cth*fct2*(1.-ci) \
	+ sPs*sth*si*fct2 + salpP*0.5*cth*fct2*(1.+ci) + (dm/M)*( 0.5*calpP*chis_y*(1.+ci) + \
																	  calmP*chis_y*siby2p2 + salmP*0.5*cth*fct3*(1.-ci) + sPs*sth*si*fct3 +\
																	  salpP*cth*0.5*fct3*(1.+ci) );	
	
}



void LCGWSpinBBHHHarm1::ComputehAllPN2(double& h0PNp, double& h0PNc, double& h05PNp, double& h05PNc, double& h1PNp, double& h1PNc)
{
	double h1PNSOp = 0.0;
    double h1PNSOc = 0.0;
	
	double del = dm/M;
	double chia_x = 0.5*(x1*S1x - x2*S2x);
	double chia_y = 0.5*(x1*S1y - x2*S2y);
	double chia_z = 0.5*(x1*S1z - x2*S2z);
	double chis_x = 0.5*(x1*S1x + x2*S2x);
	double chis_y = 0.5*(x1*S1y + x2*S2y);
	double chis_z = 0.5*(x1*S1z + x2*S2z);
	double fct2 = chia_x*cth - chia_z*sth;
	double fct3 = chis_x*cth - chis_z*sth;
	
	
	//! **  iota and  mixed
	double  nu = .25*(1.-del*del);
	double  ChpM1P3nu = -1. + 3. * nu;
	double  ChpSHio = sin(iota / 0.2e1);
	double  ChpS1io = sin(iota);  
	double  ChpCHio = cos(iota / 0.2e1);
	double  ChpC1io = cos(iota);         
	double  ChpSHioE2 = ChpSHio * ChpSHio;
	double  ChpCHioE2 = ChpCHio * ChpCHio;
	double  ChpSHioE3 = ChpSHioE2 * ChpSHio;
	double  ChpSHioE4 = ChpSHioE3 * ChpSHio;
	double  ChpSHioE5 = ChpSHioE4 * ChpSHio;
	double  ChpSHioE6 = ChpSHioE5 * ChpSHio;
	double  ChpSHioE7 = ChpSHioE6 * ChpSHio;
	double  ChpSHioE8 = ChpSHioE7 * ChpSHio;    
	double  ChpCHioE3 = ChpCHioE2 * ChpCHio;
	double  ChpCHioE4 = ChpCHioE3 * ChpCHio;
	double  ChpCHioE5 = ChpCHioE4 * ChpCHio;
	double  ChpCHioE6 = ChpCHioE5 * ChpCHio;
	double  ChpCHioE7 = ChpCHioE6 * ChpCHio;
	double  ChpCHioE8 = ChpCHioE7 * ChpCHio; 
	double  ChpS1ioE2 = ChpS1io * ChpS1io;       
	double  ChpS1ioE3 = ChpS1ioE2 * ChpS1io;
	double  ChpS1ioE4 = ChpS1ioE3 * ChpS1io;
	double  ChpC1ioE2 = ChpC1io * ChpC1io;  
	double  ChpC1ioE3 = ChpC1ioE2 * ChpC1io;  
	double  ChpS2io = 2. * ChpS1io * ChpC1io ;
	double  ChpS3io = ChpS1io * (4. * ChpC1ioE2 - 1.);
	double  ChpC2io = 2. * ChpC1ioE2 - 1.;
	double  ChpC3io = 4. * ChpC1ioE3 - 3. * ChpC1io;    
	double  ChpC3by2io = ChpCHio * ChpC1io - ChpSHio * ChpS1io;
	double  ChpS3by2io = ChpSHio * ChpC1io + ChpCHio * ChpS1io;
	double  ChpC5by2io = ChpCHio * ChpC2io - ChpSHio * ChpS2io;
	double  ChpS5by2io = ChpSHio * ChpC2io + ChpCHio * ChpS2io;      
	double  ChpC2ioXS3th = ChpC2io * ChpS3th;
	double  ChpC2thXS1io = ChpC2th * ChpS1io;       
	double  ChpC1ioXS3th = ChpC1io * ChpS3th;
	double  ChpC1ioXC4th = ChpC1io * ChpC4th;
	double  ChpS2ioXC2th = ChpS2io * ChpC2th;
	double  ChpC2ioXS1th = ChpC2io * ChpS1th;
	double  ChpC1ioXS1th = ChpC1io * ChpS1th;
	double  ChpC1ioXC2th = ChpC1io * ChpC2th;
	double  ChpSHioXS2th = ChpSHio * ChpS2th;
	double  ChpSHioXS4th = ChpSHio * ChpS4th;
	double  ChpCHioXS2th = ChpCHio * ChpS2th;
	double  ChpCHioXS4th = ChpCHio * ChpS4th;
	double  ChpC2ioXC2th = ChpC2io * ChpC2th;
	double  ChpC2ioXC4th = ChpC2io * ChpC4th;
	double  ChpC2ioXS2th = ChpC2io * ChpS2th;
	double  ChpC3ioXS2th = ChpC3io * ChpS2th;
	double  ChpC1ioXC1th = ChpC1th * ChpC1io;  
	double  ChpC2ioXC3th = ChpC2io * ChpC3th;
	double  ChpC1ioXC3th = ChpC1io * ChpC3th;
	double  ChpC1thXC2io = ChpC1th * ChpC2io;    
	double  ChpC1thXS1io = ChpC1th * ChpS1io;
	double  ChpC1thXS2io = ChpC1th * ChpS2io;
	double  ChpC1thXS3io = ChpC1th * ChpS3io;
	double  ChpS3thXC1th = ChpS3th * ChpC1th;
	double  ChpS1ioXC2thXC1io = ChpC2thXS1io * ChpC1io;
	double  ChpS1ioXC2thXC2io = ChpC2thXS1io * ChpC2io;
	double  ChpCHioXC2ioXS2th = ChpCHioXS2th * ChpC2io;
	double  ChpCHioXC1ioXS2th = ChpCHioXS2th *  ChpC1io;
	double  ChpSHioXC1ioXS2th = ChpCHioXS2th  * ChpC1io;
	double  ChpSHioXC2ioXS2th = ChpSHio * ChpC2ioXS2th;
	double  ChpC1thXC1ioXS1th = ChpC1ioXC1th * ChpS1th;
	double  ChpC1ioXS1thXC2th = ChpC1ioXS1th * ChpC2th;
	double  Chp3PC2thXCHioE4 = Chp3PC2th * ChpCHioE4;
	double  ChpCHioE3XSHioXS2th = ChpCHioE3 * ChpSHioXS2th;
	double  ChpS1thE2XS1ioE2 = ChpS1thE2 * ChpS1ioE2;
	double  Chp3PC2thXSHioE4 = Chp3PC2th * ChpSHioE4;
	double  ChpSHioE3XCHioXS2th = ChpSHioE3 * ChpCHioXS2th;
	double  ChpS3ioX3PC2th = ChpS3io * Chp3PC2th;
	double  ChpS3thXS1ioE2 = ChpS3th * ChpS1ioE2;
	double  ChpS1ioE2X5S1thPS3thE2 = ChpS1ioE2 * Chp5S1thPS3th * Chp5S1thPS3th;
	double  ChpM1P3nuXS3thXS1ioE2 = ChpM1P3nu * ChpS3thXS1ioE2;
	double  ChpS1thE2XM1P3nuX3PC2th = ChpS1thE2 * ChpM1P3nu * Chp3PC2th;
	double  ChpSHioXS4thXM1P3nu = ChpSHioXS4th * ChpM1P3nu;
	double  ChpSHioE3XCHioXC1ioXS2thXM1P3nu = ChpSHioE3 * ChpCHioXC1ioXS2th * ChpM1P3nu;
	double  ChpS1ioE2X5S1thPS3thXM1P3nu = ChpS1ioE2 * Chp5S1thPS3th * ChpM1P3nu;
	double  ChpSHioE3XM1P3nu = ChpSHioE3 * ChpM1P3nu;
	double  ChpS2thXM1P3nu = ChpS2th * ChpM1P3nu;
	double  ChpC1ioXC1thE2 = ChpC1io * ChpC1thE2;
	double  Chp5S1thPS3thXC1ioXC1th = Chp5S1thPS3th * ChpC1ioXC1th;
	double  ChpCHioE4XSHioE2XC1th = ChpCHioE4 * ChpSHioE2 * ChpC1th;
	double  ChpC1thX5S1thPS3th = ChpC1th * Chp5S1thPS3th;
	double  ChpC1thXC1ioXS3th = ChpC1th * ChpC1ioXS3th;
	double  ChpC1thX5S1thPS3thXC1ioE2 = ChpC1thX5S1thPS3th * ChpC1ioE2;
	double  Chp5S1thPS3thXC1ioXC1thXC1ioE2 = Chp5S1thPS3thXC1ioXC1th * ChpC1ioE2;
	double  ChpC1thXS1thE2XM1P3nu = ChpC1th * ChpS1thE2 * ChpM1P3nu;
	double  ChpSHioXM1Pnu = ChpSHio * ChpM1P3nu;
	double  ChpCHioXM1P3nu = ChpCHio * ChpM1P3nu;
	double  ChpCHioE6XSHioE2 = ChpCHioE6 * ChpSHioE2;
	double  ChpCHioE2XSHioE6 = ChpCHioE2 * ChpSHioE6;
	double  ChpS3thXS1io = ChpS3th * ChpS1io;
	double  ChpS3thXSHioE2 = ChpS3th * ChpSHioE2;
	double  ChpSHioE2X5S1thPS3th = ChpSHioE2 * Chp5S1thPS3th;
	double  ChpC1ioXC1thXS1io = ChpC1ioXC1th * ChpS1io;
	double  ChpC1thXC2thXS1io = ChpC1th * ChpC2thXS1io;
	
	
	//! ** Others coefficients
	
	double  ChpCalp  = cos(alpha);
	double  ChpSalp  = sin(alpha);
	double  ChpCPsi  = cos(Psi);
	double  ChpSPsi  = sin(Psi);     //  erreur 
	double  ChpCalpE2  = ChpCalp * ChpCalp;
	double  ChpCalpE3 = ChpCalpE2 * ChpCalp ;
	double  ChpCalpE4 = ChpCalpE3 * ChpCalp;
	double  ChpCPsiE2  = ChpCPsi * ChpCPsi;
	double  ChpCPsiE3 = ChpCPsiE2 * ChpCPsi;
	double  ChpCPsiE4 = ChpCPsiE3 * ChpCPsi;      
	double  ChpC2alp = 2. * ChpCalpE2 - 1.;
	double  ChpC3alp = 4.  * ChpCalpE3 - 3. * ChpCalp;   // erreur
	double  ChpC4alp = 8. * ChpCalpE4 - 8. * ChpCalpE2 + 1.;
	double  ChpC2Psi = 2. * ChpCPsiE2 - 1.;
	double  ChpC3Psi = ChpCPsi * (4. *  ChpCPsiE2 - 3.);  // erreur
	double  ChpC4Psi = 8. * ChpCPsiE4 - 8. * ChpCPsiE2 + 1.;
	double  ChpS2alp = 2. * ChpSalp * ChpCalp;
	double  ChpS3alp =  ChpSalp * (4. * ChpCalpE2 -1.); // erreur
	double  ChpS4alp = ChpSalp * (8. * ChpCalpE3 - 4. * ChpCalp);
	double  ChpS2Psi = 2. * ChpSPsi * ChpCPsi ;
	double  ChpS3Psi = ChpSPsi * (4. * ChpCPsiE2 -1.);   // erreur
	double  ChpS4Psi = ChpSPsi * (8. * ChpCPsiE3 - 4. * ChpCPsi);
	double  cP4alpP4Psi  =  ChpC4alp * ChpC4Psi - ChpS4alp * ChpS4Psi       ;
	double  cP3alpP4Psi  =  ChpC3alp * ChpC4Psi - ChpS3alp * ChpS4Psi        ;
	double  cP2alpP4Psi  =  ChpC2alp * ChpC4Psi - ChpS2alp * ChpS4Psi      ;
	double  cP1alpP4Psi  =  ChpCalp * ChpC4Psi - ChpSalp * ChpS4Psi    ;
	double  cP0alpP4Psi  = ChpC4Psi     ;
	double  cM1alpP4Psi  =  ChpCalp * ChpC4Psi + ChpSalp * ChpS4Psi        ;
	double  cM2alpP4Psi  =  ChpC2alp * ChpC4Psi + ChpS2alp * ChpS4Psi       ;
	double  cM3alpP4Psi  =  ChpC3alp * ChpC4Psi + ChpS3alp * ChpS4Psi     ;
	double  cM4alpP4Psi  =  ChpC4alp * ChpC4Psi + ChpS4alp * ChpS4Psi   ;
	double  sP4alpP4Psi  =    ChpS4alp * ChpC4Psi + ChpC4alp * ChpS4Psi     ;
	double  sP3alpP4Psi  =   ChpS3alp * ChpC4Psi + ChpC3alp * ChpS4Psi      ;
	double  sP2alpP4Psi  =    ChpS2alp * ChpC4Psi + ChpC2alp * ChpS4Psi     ;
	double  sP1alpP4Psi  =   ChpSalp * ChpC4Psi + ChpCalp * ChpS4Psi      ;
	double  sP0alpP4Psi  =    ChpS4Psi     ;
	double  sM1alpP4Psi  =  - ChpSalp * ChpC4Psi + ChpCalp * ChpS4Psi      ;
	double  sM2alpP4Psi  =  - ChpS2alp * ChpC4Psi + ChpC2alp * ChpS4Psi        ;
	double  sM3alpP4Psi  =   - ChpS3alp * ChpC4Psi + ChpC3alp * ChpS4Psi       ;
	double  sM4alpP4Psi  =   - ChpS4alp * ChpC4Psi + ChpC4alp * ChpS4Psi       ;
	double  cP4alpP3Psi  = ChpC4alp * ChpC3Psi - ChpS4alp * ChpS3Psi      ;
	double  cP3alpP3Psi  = ChpC3alp * ChpC3Psi - ChpS3alp * ChpS3Psi        ;
	double  cP2alpP3Psi  =   ChpC2alp * ChpC3Psi - ChpS2alp * ChpS3Psi      ;
	double  cP1alpP3Psi  =     ChpCalp * ChpC3Psi - ChpSalp * ChpS3Psi    ;
	double  cP0alpP3Psi  =       ChpC3Psi   ;
	double  cM1alpP3Psi  =      ChpCalp * ChpC3Psi + ChpSalp * ChpS3Psi   ;
	double  cM2alpP3Psi  =   ChpC2alp * ChpC3Psi + ChpS2alp * ChpS3Psi      ;
	double  cM3alpP3Psi  =     ChpC3alp * ChpC3Psi + ChpS3alp * ChpS3Psi    ;
	double  cM4alpP3Psi  =    ChpC4alp * ChpC3Psi + ChpS4alp * ChpS3Psi     ;
	double  sP4alpP3Psi  =   ChpS4alp * ChpC3Psi + ChpC4alp * ChpS3Psi        ;
	double  sP3alpP3Psi  =   ChpS3alp * ChpC3Psi + ChpC3alp * ChpS3Psi       ;
	double  sP2alpP3Psi  =    ChpS2alp * ChpC3Psi + ChpC2alp * ChpS3Psi      ;
	double  sP1alpP3Psi  =    ChpSalp * ChpC3Psi + ChpCalp * ChpS3Psi      ;
	double  sP0alpP3Psi  =     ChpS3Psi     ;
	double  sM1alpP3Psi  =    - ChpSalp * ChpC3Psi + ChpCalp * ChpS3Psi     ;
	double  sM2alpP3Psi  =   - ChpS2alp * ChpC3Psi + ChpC2alp * ChpS3Psi      ;
	double  sM3alpP3Psi  =   - ChpS3alp * ChpC3Psi + ChpC3alp * ChpS3Psi      ;
	double  sM4alpP3Psi  =    - ChpS4alp * ChpC3Psi + ChpC4alp * ChpS3Psi     ;
	double  cP4alpP2Psi  =  ChpC4alp * ChpC2Psi - ChpS4alp * ChpS2Psi       ;
	double  cP3alpP2Psi  =  ChpC3alp * ChpC2Psi - ChpS3alp * ChpS2Psi       ;
	double  cP2alpP2Psi  =  ChpC2alp * ChpC2Psi - ChpS2alp * ChpS2Psi       ;
	double  cP1alpP2Psi  =  ChpCalp * ChpC2Psi - ChpSalp * ChpS2Psi       ;
	double  cP0alpP2Psi  =  ChpC2Psi     ;
	double  cM1alpP2Psi  =    ChpCalp * ChpC2Psi + ChpSalp * ChpS2Psi     ;
	double  cM2alpP2Psi  = ChpC2alp * ChpC2Psi + ChpS2alp * ChpS2Psi         ;
	double  cM3alpP2Psi  =   ChpC3alp * ChpC2Psi + ChpS3alp * ChpS2Psi       ;
	double  cM4alpP2Psi  =   ChpC4alp * ChpC2Psi + ChpS4alp * ChpS2Psi       ;
	double  sP4alpP2Psi  =    ChpS4alp * ChpC2Psi + ChpC4alp * ChpS2Psi      ;
	double  sP3alpP2Psi  =   ChpS3alp * ChpC2Psi + ChpC3alp * ChpS2Psi        ;
	double  sP2alpP2Psi  =   ChpS2alp * ChpC2Psi + ChpC2alp * ChpS2Psi        ;
	double  sP1alpP2Psi  =     ChpSalp * ChpC2Psi + ChpCalp * ChpS2Psi      ;
	double  sP0alpP2Psi  =    ChpS2Psi      ;
	double  sM1alpP2Psi  =    - ChpSalp * ChpC2Psi + ChpCalp * ChpS2Psi      ;
	double  sM2alpP2Psi  =  - ChpS2alp * ChpC2Psi + ChpC2alp * ChpS2Psi       ;
	double  sM3alpP2Psi  =   - ChpS3alp * ChpC2Psi + ChpC3alp * ChpS2Psi      ;
	double  sM4alpP2Psi  =  - ChpS4alp * ChpC2Psi + ChpC4alp * ChpS2Psi       ;
	double  cP4alpP1Psi  = ChpC4alp * ChpCPsi - ChpS4alp * ChpSPsi         ;
	double  cP3alpP1Psi  = ChpC3alp * ChpCPsi - ChpS3alp * ChpSPsi        ;
	double  cP2alpP1Psi  =   ChpC2alp * ChpCPsi - ChpS2alp * ChpSPsi      ;
	double  cP1alpP1Psi  =    ChpCalp * ChpCPsi - ChpSalp * ChpSPsi     ;
	double  cP0alpP1Psi  =   ChpCPsi       ;
	double  cM1alpP1Psi  =  ChpCalp * ChpCPsi + ChpSalp * ChpSPsi       ;
	double  cM2alpP1Psi  =  ChpC2alp * ChpCPsi + ChpS2alp * ChpSPsi       ;
	double  cM3alpP1Psi  =    ChpC3alp * ChpCPsi + ChpS3alp * ChpSPsi     ;
	double  cM4alpP1Psi  =    ChpC4alp * ChpCPsi + ChpS4alp * ChpSPsi     ;
	double  sP4alpP1Psi  =   ChpS4alp * ChpCPsi + ChpC4alp * ChpSPsi       ;
	double  sP3alpP1Psi  =   ChpS3alp * ChpCPsi + ChpC3alp * ChpSPsi        ;
	double  sP2alpP1Psi  =    ChpS2alp * ChpCPsi + ChpC2alp * ChpSPsi       ;
	double  sP1alpP1Psi  =    ChpSalp * ChpCPsi + ChpCalp * ChpSPsi       ;
	double  sP0alpP1Psi  =    ChpSPsi      ;
	double  sM1alpP1Psi  =  - ChpSalp * ChpCPsi + ChpCalp * ChpSPsi       ;
	double  sM2alpP1Psi  =  - ChpS2alp * ChpCPsi + ChpC2alp * ChpSPsi       ;
	double  sM3alpP1Psi  =   - ChpS3alp * ChpCPsi + ChpC3alp * ChpSPsi      ;
	double  sM4alpP1Psi  =   - ChpS4alp * ChpCPsi + ChpC4alp * ChpSPsi      ;
	double  cP4alpP0Psi  =    ChpC4alp     ;
	double  cP3alpP0Psi  =    ChpC3alp     ;
	double  cP2alpP0Psi  =     ChpC2alp    ;
	double  cP1alpP0Psi  =      ChpCalp   ;
	double  sP4alpP0Psi  =   ChpS4alp      ;
	double  sP3alpP0Psi  =   ChpS3alp      ;
	double  sP2alpP0Psi  =   ChpS2alp      ;
	double  sP1alpP0Psi  =    ChpSalp     ;

	
	
	h0PNp = -0.5000000000e0 * Chp3PC2thXCHioE4 * cP2alpP2Psi - 0.5000000000e0 * Chp3PC2thXSHioE4 *\
    cM2alpP2Psi - 0.2000000000e1 * ChpCHioE3XSHioXS2th * cP1alpP2Psi + 0.2000000000e1 * ChpSHioE3XCHioXS2th *\
    cM1alpP2Psi - 0.1500000000e1 * ChpS1thE2XS1ioE2 * cP0alpP2Psi;
    
	
    
    h05PNp = (-0.1406250000e1 * cP1alpP3Psi * ChpSHioE2 * ChpS1thM3S3th * ChpCHioE4 + 0.1406250000e1 * cM1alpP3Psi * ChpSHioE4 * ChpS1thM3S3th * ChpCHioE2 + 0.1359375000e1 * ChpC1ioXS1th * ChpSHioE2 * cM1alpP1Psi + 0.1562500000e-1 * cM2alpP1Psi * ChpC1th * ChpS2ioXC2th - 0.1171875000e-1 * cM2alpP1Psi * ChpC1th * ChpS3ioX3PC2th + 0.8437500000e0 * cP2alpP3Psi * ChpC1th * ChpS1ioXC2thXC1io + 0.2109375000e0 * cP2alpP3Psi * ChpC1th * ChpS1ioXC2thXC2io - 0.3125000000e-1 * cP2alpP1Psi * ChpC1th * ChpS1ioXC2thXC1io - 0.2343750000e-1 * cP2alpP1Psi * ChpC1th * ChpS1ioXC2thXC2io + 0.7812500001e-1 * ChpS1thE2 * cP0alpP1Psi * ChpC1thXS3io + 0.1562500000e-1 * ChpS1thE2 * cP0alpP1Psi * ChpC1thXS1io + 0.1406250000e0 * cP2alpP3Psi * ChpC1thXS2io + 0.3515625000e-1 * cP2alpP3Psi * ChpC1thXS3io + 0.6328125000e0 * cP2alpP3Psi * ChpC1thXC2thXS1io + 0.1757812500e0 * cP2alpP3Psi * ChpC1thXS1io - 0.1718750000e0 * cP2alpP1Psi * ChpC1thXS2io - 0.3906250000e-2 * cP2alpP1Psi * ChpC1thXS3io - 0.7812500001e-2 * cP2alpP1Psi * ChpC1thXC2thXS1io - 0.3320312501e0 * cP2alpP1Psi * ChpC1thXS1io + 0.3125000000e-1 * cM2alpP1Psi * ChpC1thXS3io + 0.1718750000e0 * cM2alpP1Psi * ChpC1thXS2io + 0.3906250000e-2 * cM2alpP1Psi * ChpC1thXC2thXS1io - 0.3320312501e0 * cM2alpP1Psi * ChpC1thXS1io + 0.2812500000e1 * ChpC1th * ChpS1thE2 * ChpS1ioE3 * cP0alpP3Psi - 0.1406250000e0 * ChpC1thXS2io * cM2alpP3Psi * Chp1P3C2th + 0.3515625000e-1 * ChpC1thXS3io * cM2alpP3Psi * Chp1P3C2th + 0.1757812500e0 * ChpC1thXS1io * cM2alpP3Psi * Chp1P3C2th + 0.1875000000e0 * ChpS3th * cP1alpP1Psi * ChpCHioE2 - 0.7812500001e-1 * ChpC1ioXS3th * cP1alpP1Psi * ChpCHioE2 - 0.1367187500e0 * Chp5S1thPS3th * cP1alpP1Psi * ChpCHioE2 + 0.5859375002e-1 * ChpC2ioXS3th * cP1alpP1Psi * ChpCHioE2 + 0.1359375000e1 * ChpC1ioXS1th * cP1alpP1Psi * ChpCHioE2 - 0.1953125000e-1 * ChpC2ioXS1th * cP1alpP1Psi * ChpCHioE2 - 0.5859375002e-1 * ChpC2ioXS3th * ChpSHioE2 * cM1alpP1Psi - 0.7812500001e-1 * ChpC1ioXS3th * ChpSHioE2 * cM1alpP1Psi + 0.1953125000e-1 * ChpC2ioXS1th * ChpSHioE2 * cM1alpP1Psi - 0.2812500000e0 * ChpCHioE6 * Chp5S1thPS3th * cP3alpP3Psi - 0.3125000000e-1 * ChpCHioE4 * ChpSHioE2X5S1thPS3th * cP3alpP1Psi - 0.1125000000e0 * Chp3PC2th * ChpSHioE6 * cM3alpP3Psi * ChpS3th + 0.1125000000e0 * Chp3PC2th * ChpSHioE6 * cM3alpP3Psi * Chp5S1thPS3th + 0.3125000000e-1 * ChpCHioE2 * Chp5S1thPS3th * ChpSHioE4 * cM3alpP1Psi + 0.1367187500e0 * ChpSHioE2X5S1thPS3th * cM1alpP1Psi - 0.1875000000e0 * ChpS3thXSHioE2 * cM1alpP1Psi) * del;
    
	
	h1PNp = -0.1333333333e0 * ChpSHioE5 * ChpM1P3nu * cM3alpP2Psi * ChpC1thE3 * ChpC3by2io * ChpS3th -\
    0.1444444444e1 * cM2alpP2Psi * ChpM1P3nu * Chp3PC2thXSHioE4 + 0.5866666668e2 * ChpCHioE3 * ChpSHioE5 *\
    ChpS2thXM1P3nu * cM1alpP4Psi + 0.2133333333e1 * ChpC1thE3 * ChpCHioE7 * ChpSHioXM1Pnu * cP3alpP4Psi *\
    ChpS3th - 0.2133333333e1 * ChpC1thE3 * ChpCHioE7 * ChpSHioXM1Pnu * cP3alpP4Psi * Chp5S1thPS3th -\
    0.2133333333e1 * ChpC1thE3 * ChpSHioE7 * ChpCHioXM1P3nu * cM3alpP4Psi * ChpS3th + 0.2133333333e1 *\
    ChpC1thE3 * ChpSHioE7 * ChpCHioXM1P3nu * cM3alpP4Psi * Chp5S1thPS3th - 0.5866666668e2 * ChpS2th *\
    ChpCHioE5 * ChpSHioE3XM1P3nu * cP1alpP4Psi + 0.1866666667e2 * ChpS2th * ChpCHioE5 * ChpSHioE3XM1P3nu *\
    cP1alpP4Psi * Chp3PC2th + 0.1458333333e0 * cP1alpP2Psi * ChpCHioE3 * ChpS4th * ChpS5by2io * ChpM1P3nu -\
    0.1866666667e2 * ChpCHioE3 * ChpSHioE5 * ChpS2thXM1P3nu * cM1alpP4Psi * Chp3PC2th + 0.1777777778e1 *\
    cM1alpP2Psi * ChpM1P3nu * ChpSHioE3XCHioXS2th - 0.2916666666e0 * cM1alpP2Psi * ChpC3by2io * ChpS4th *\
    ChpSHioE3XM1P3nu - 0.1458333333e0 * cM1alpP2Psi * ChpC5by2io * ChpS4th * ChpSHioE3XM1P3nu -\
    0.6666666668e1 * ChpS1ioE4 * cP0alpP4Psi * ChpS1thE2 * ChpM1P3nu + 0.3812499999e1 * cP2alpP2Psi *\
    ChpCHioE4 * ChpM1P3nu - 0.1444444444e1 * cP2alpP2Psi * ChpM1P3nu * Chp3PC2thXCHioE4 + 0.3812499999e1 *\
    cM2alpP2Psi * ChpSHioE4 * ChpM1P3nu - 0.1777777778e1 * cP1alpP2Psi * ChpM1P3nu * ChpCHioE3XSHioXS2th +\
    0.3958333333e0 * cP1alpP2Psi * ChpCHioE3 * ChpSHioXS4thXM1P3nu + 0.8333333334e-1 * cM1alpP2Psi *\
    ChpCHioXC2ioXS2th * ChpSHioE3XM1P3nu + 0.2916666666e1 * ChpS1ioE4 * cP0alpP4Psi *\
    ChpS1thE2XM1P3nuX3PC2th - 0.2533333333e0 * cP0alpP2Psi * Chp5S1thPS3th * ChpS3thXS1ioE2 + 0.1266666667e0 *\
    cP0alpP2Psi * ChpS3th * ChpS3thXS1ioE2 + 0.1562500000e-1 * cP0alpP2Psi * ChpC2ioXS1th *\
    ChpS1ioE2X5S1thPS3thXM1P3nu - 0.3645833333e-1 * cP0alpP2Psi * ChpC2ioXS3th * ChpM1P3nuXS3thXS1ioE2 -\
    0.2395833333e-1 * cP0alpP2Psi * ChpM1P3nu * ChpS1ioE2X5S1thPS3thE2 - 0.1562500000e-1 * cP0alpP2Psi *\
    ChpC2ioXS1th * ChpM1P3nuXS3thXS1ioE2 + 0.3645833333e-1 * cP0alpP2Psi * ChpC2ioXS3th *\
    ChpS1ioE2X5S1thPS3thXM1P3nu - 0.4999999999e-1 * cP0alpP2Psi * ChpS3th * ChpM1P3nuXS3thXS1ioE2 +\
    0.6666666668e0 * ChpCHioE8 * ChpS1thE2XM1P3nuX3PC2th * cP4alpP4Psi + 0.1666666667e0 * ChpS1thE2XM1P3nuX3PC2th *\
    ChpCHioE6XSHioE2 * cP4alpP2Psi + 0.1666666667e0 * ChpS1thE2XM1P3nuX3PC2th * ChpCHioE2XSHioE6 *\
    cM4alpP2Psi + 0.6666666668e0 * ChpSHioE8 * ChpS1thE2XM1P3nuX3PC2th * cM4alpP4Psi - 0.1666666667e0 *\
    ChpCHioE5 * cP3alpP2Psi * ChpSHioXS4thXM1P3nu + 0.7395833334e-1 * cP0alpP2Psi * Chp5S1thPS3th *\
    ChpM1P3nuXS3thXS1ioE2 - 0.2916666666e0 * cP1alpP2Psi * ChpCHioE3 * ChpS4th * ChpS3by2io * ChpM1P3nu +\
    0.2583333333e1 * cM1alpP2Psi * ChpSHioE3XCHioXC1ioXS2thXM1P3nu + 0.1266666667e0 * cP0alpP2Psi *\
    ChpS1ioE2X5S1thPS3thE2 - 0.3958333333e0 * cM1alpP2Psi * ChpCHioXS4th * ChpSHioE3XM1P3nu -\
    0.1333333333e0 * ChpCHioE5 * cP3alpP2Psi * ChpM1P3nu * ChpS3by2io * ChpC1thE3 * ChpS3th + 0.1333333333e0 *\
    ChpCHioE5 * cP3alpP2Psi * ChpM1P3nu * ChpS3by2io * ChpC1thE3 * Chp5S1thPS3th + 0.1333333333e0 * ChpSHioE5 *\
    ChpM1P3nu * cM3alpP2Psi * ChpC1thE3 * ChpC3by2io * Chp5S1thPS3th + 0.1055555556e1 * Chp3PC2thXCHioE4 *\
    cP2alpP2Psi + 0.1055555556e1 * Chp3PC2thXSHioE4 * cM2alpP2Psi + 0.4222222222e1 * ChpCHioE3XSHioXS2th *\
    cP1alpP2Psi - 0.4222222222e1 * ChpSHioE3XCHioXS2th * cM1alpP2Psi + 0.6666666668e0 * ChpCHioE5 * cP3alpP2Psi *\
    ChpM1P3nu * ChpSHioXS2th + 0.1666666667e0 * ChpSHioE5 * ChpM1P3nu * cM3alpP2Psi * ChpCHioXS4th - 0.6666666668e0 *\
    ChpSHioE5 * ChpM1P3nu * cM3alpP2Psi * ChpCHioXS2th + 0.6666666668e0 * ChpM1P3nu * Chp5p4C2thP7C4th *\
    ChpCHioE6XSHioE2 * cP2alpP4Psi - 0.2916666666e0 * cP2alpP2Psi * ChpCHioE4 * ChpC1ioXC4th * ChpM1P3nu +\
    0.1458333333e0 * cP2alpP2Psi * ChpCHioE4 * ChpC2ioXC4th * ChpM1P3nu - 0.2083333333e0 * cP2alpP2Psi *\
    ChpCHioE4 * ChpC1io * ChpM1P3nu + 0.1833333333e1 * cP2alpP2Psi * ChpCHioE4 * ChpC1ioXC2th * ChpM1P3nu +\
    0.1875000000e0 * cP2alpP2Psi * ChpCHioE4 * ChpC4th * ChpM1P3nu + 0.8333333334e-1 * cP2alpP2Psi * ChpCHioE4 *\
    ChpC2ioXC2th * ChpM1P3nu + 0.1041666667e0 * cP2alpP2Psi * ChpCHioE4 * ChpC2io * ChpM1P3nu + 0.2083333333e0 *\
    cM2alpP2Psi * ChpSHioE4 * ChpC1io * ChpM1P3nu + 0.1041666667e0 * cM2alpP2Psi * ChpSHioE4 * ChpC2io *\
    ChpM1P3nu + 0.2916666666e0 * cM2alpP2Psi * ChpSHioE4 * ChpC1ioXC4th * ChpM1P3nu - 0.1833333333e1 *\
    cM2alpP2Psi * ChpSHioE4 * ChpC1ioXC2th * ChpM1P3nu + 0.8333333334e-1 * cM2alpP2Psi * ChpSHioE4 *\
    ChpC2ioXC2th * ChpM1P3nu + 0.1875000000e0 * cM2alpP2Psi * ChpSHioE4 * ChpC4th * ChpM1P3nu + 0.1458333333e0 *\
    cM2alpP2Psi * ChpSHioE4 * ChpC2ioXC4th * ChpM1P3nu + 0.6666666668e0 * ChpM1P3nu * Chp5p4C2thP7C4th *\
    ChpCHioE2XSHioE6 * cM2alpP4Psi + 0.2583333333e1 * cP1alpP2Psi * ChpCHioE3 * ChpSHioXC1ioXS2th *\
    ChpM1P3nu - 0.8333333334e-1 * cP1alpP2Psi * ChpCHioE3 * ChpSHioXC2ioXS2th * ChpM1P3nu;
    
	
    
	h1PNSOp = ((-0.5000000000e0 * cM1alpP1Psi * ChpC1ioXC1th + ChpCHioE2 * cP1alpP1Psi * ChpC1th + 0.5000000000e0 *\
				cM1alpP1Psi * ChpC1th) * chis_x + (-0.1e1 * ChpCHioE2 * ChpC1th * sP1alpP1Psi + 0.2000000000e0 *\
												   ChpS3thXS1io * sP0alpP1Psi + ChpSHioE2 * ChpC1th * sM1alpP1Psi - 0.2000000000e0 * Chp5S1thPS3th * ChpS1io *\
												   sP0alpP1Psi) * chis_y + (-0.2000000000e0 * Chp5S1thPS3th * cP1alpP1Psi * ChpCHioE2 +\
																			0.2000000000e0 * ChpS3thXSHioE2 * cM1alpP1Psi - 0.2000000000e0 * ChpSHioE2X5S1thPS3th * cM1alpP1Psi +\
																			0.2000000000e0 * ChpS3th * cP1alpP1Psi * ChpCHioE2) * chis_z) * del + (-0.5000000000e0 * cM1alpP1Psi *\
																																				   ChpC1ioXC1th + ChpCHioE2 * cP1alpP1Psi * ChpC1th + 0.5000000000e0 * cM1alpP1Psi * ChpC1th) * chia_x +\
    (-0.1e1 * ChpCHioE2 * ChpC1th * sP1alpP1Psi + 0.2000000000e0 * ChpS3thXS1io * sP0alpP1Psi + ChpSHioE2 *\
     ChpC1th * sM1alpP1Psi - 0.2000000000e0 * Chp5S1thPS3th * ChpS1io * sP0alpP1Psi) * chia_y + (-0.2000000000e0 *\
																								 Chp5S1thPS3th * cP1alpP1Psi * ChpCHioE2 + 0.2000000000e0 * ChpS3thXSHioE2 * cM1alpP1Psi -\
																								 0.2000000000e0 * ChpSHioE2X5S1thPS3th * cM1alpP1Psi + 0.2000000000e0 * ChpS3th *\
																								 cP1alpP1Psi * ChpCHioE2) * chia_z;
    
    
	
    
	//   h1PNp +=  h1PNSOp;
    
    
    
	h0PNc = -0.2000000000e1 * ChpC1th * ChpCHioE4 * sP2alpP2Psi + 0.2000000000e1 * ChpC1th * ChpSHioE4 *\
    sM2alpP2Psi + 0.8000000000e0 * ChpCHioE3 * ChpSHio * sP1alpP2Psi * ChpS3th - 0.8000000000e0 * ChpCHioE3 *\
    ChpSHio * sP1alpP2Psi * Chp5S1thPS3th + 0.8000000000e0 * ChpCHio * ChpSHioE3 * sM1alpP2Psi *\
    ChpS3th - 0.8000000000e0 * ChpCHio * ChpSHioE3 * sM1alpP2Psi * Chp5S1thPS3th;
    
	
    
	
    
    
    h05PNc = (0.4500000000e1 * Chp3PC2th * sP2alpP3Psi * ChpSHio * ChpCHioE5 - 0.2500000000e0 * Chp3PC2th * sP2alpP1Psi * ChpCHioE3 * ChpSHio - 0.2500000000e0 * ChpC1ioXC2th * sP2alpP1Psi * ChpCHioE3 * ChpSHio + 0.5625000000e1 * sP1alpP3Psi * ChpSHioE2 * ChpS2th * ChpCHioE4 + 0.5625000000e1 * sM1alpP3Psi * ChpSHioE4 * ChpS2th * ChpCHioE2 - 0.1250000000e0 * ChpCHioE4 * ChpS2th * ChpSHioE2 * sP3alpP1Psi - 0.1250000000e0 * ChpCHioE2 * ChpS2th * ChpSHioE4 * sM3alpP1Psi + 0.2500000000e0 * ChpCHio * ChpSHioE3 * sM2alpP1Psi * Chp3PC2th - 0.2500000000e0 * ChpCHio * ChpSHioE3 * sM2alpP1Psi * ChpC1ioXC2th - 0.4500000000e1 * ChpCHio * ChpSHioE5 * sM2alpP3Psi * Chp3PC2th + 0.5000000000e0 * sP0alpP1Psi * ChpS2io * ChpS1thE2 + 0.1562500000e-1 * ChpC1thXS1th * sP1alpP1Psi + 0.3359375000e0 * ChpC1thXC1ioXS1th * sP1alpP1Psi + 0.1796875000e0 * ChpC2ioXS2th * sP1alpP1Psi + 0.1953125000e-1 * ChpC3ioXS2th * sP1alpP1Psi + 0.1562500000e-1 * ChpC1thXS1th * sM1alpP1Psi - 0.3359375000e0 * ChpC1thXC1ioXS1th * sM1alpP1Psi + 0.1796875000e0 * ChpC2ioXS2th * sM1alpP1Psi - 0.1953125000e-1 * ChpC3ioXS2th * sM1alpP1Psi - 0.1350000000e2 * sP2alpP3Psi * ChpSHio * ChpCHioE5 - 0.2500000000e0 * sP2alpP1Psi * ChpCHioE3 * ChpSHio - 0.1125000000e1 * ChpCHioE6 * ChpS2th * sP3alpP3Psi - 0.1125000000e1 * ChpS2th * ChpSHioE6 * sM3alpP3Psi + 0.2500000000e0 * ChpCHio * ChpSHioE3 * sM2alpP1Psi + 0.1350000000e2 * ChpCHio * ChpSHioE5 * sM2alpP3Psi) * del;
	
	
    
    
    h1PNc =   0.2666666667e0 * ChpCHioE5 * ChpSHioXM1Pnu * sP3alpP2Psi * Chp5S1thPS3th - 0.1666666666e0 * ChpCHioE5 * ChpSHioXM1Pnu * sP3alpP2Psi * ChpC1ioXS1th - 0.6666666667e0 * ChpC1thXS1thE2XM1P3nu * ChpCHioE2XSHioE6 * sM4alpP2Psi - 0.2916666666e0 * ChpSHioE4 * sM2alpP2Psi * ChpC2ioXC3th * ChpM1P3nu + 0.1666666666e0 * ChpSHioE4 * sM2alpP2Psi * ChpC1ioXC3th * ChpM1P3nu + 0.1652777778e1 * ChpSHioE4 * sM2alpP2Psi * ChpC1th * ChpM1P3nu - 0.4166666667e-1 * ChpSHioE4 * sM2alpP2Psi * ChpC1thXC2io * ChpM1P3nu + 0.1166666667e1 * ChpCHioE4 * sP2alpP2Psi * ChpC1ioXC1th * ChpM1P3nu - 0.1250000000e0 * ChpCHioE4 * sP2alpP2Psi * ChpM1P3nu * ChpC3th + 0.2916666666e0 * ChpCHioE4 * sP2alpP2Psi * ChpC2ioXC3th * ChpM1P3nu + 0.1666666666e0 * ChpCHioE4 * sP2alpP2Psi * ChpC1ioXC3th * ChpM1P3nu - 0.1652777778e1 * ChpCHioE4 * sP2alpP2Psi * ChpC1th * ChpM1P3nu + 0.9333333332e1 * ChpM1P3nu * ChpCHioE6XSHioE2 * sP2alpP4Psi * ChpC3th + 0.1333333333e1 * ChpM1P3nu * ChpCHioE6XSHioE2 * sP2alpP4Psi * ChpC1th + 0.1333333333e1 * ChpSHioE7 * ChpS1thM3S3th * ChpCHioXM1P3nu * sM3alpP4Psi + 0.2666666667e0 * ChpSHioE5 * ChpCHioXM1P3nu * sM3alpP2Psi * Chp5S1thPS3th + 0.1666666666e0 * ChpSHioE5 * ChpCHioXM1P3nu * sM3alpP2Psi * ChpC1ioXS1th - 0.2666666667e0 * ChpSHioE5 * ChpCHioXM1P3nu * sM3alpP2Psi * ChpS3th - 0.5000000001e0 * ChpSHioE5 * ChpCHioXM1P3nu * sM3alpP2Psi * ChpC1ioXS3th - 0.2666666667e1 * ChpSHioE8 * ChpC1thXS1thE2XM1P3nu * sM4alpP4Psi + 0.1166666667e1 * ChpSHioE4 * sM2alpP2Psi * ChpC1ioXC1th * ChpM1P3nu + 0.1250000000e0 * ChpSHioE4 * sM2alpP2Psi * ChpM1P3nu * ChpC3th + 0.4166666667e-1 * ChpCHioE4 * sP2alpP2Psi * ChpC1thXC2io * ChpM1P3nu - 0.9333333332e1 * ChpM1P3nu * ChpCHioE2XSHioE6 * sM2alpP4Psi * ChpC3th - 0.1333333333e1 * ChpM1P3nu * ChpCHioE2XSHioE6 * sM2alpP4Psi * ChpC1th + 0.8533333333e1 * ChpCHioE5 * ChpSHioE3XM1P3nu * sP1alpP4Psi * ChpS3th + 0.8000000000e0 * ChpCHioE5 * ChpSHioE3XM1P3nu * sP1alpP4Psi * Chp5S1thPS3th + 0.2916666666e0 * ChpCHioE3 * sP1alpP2Psi * ChpC2ioXS3th * ChpSHioXM1Pnu + 0.1250000000e0 * ChpCHioE3 * sP1alpP2Psi * ChpC2ioXS1th * ChpSHioXM1Pnu + 0.8333333333e0 * ChpCHioE3 * sP1alpP2Psi * ChpC1ioXS1th * ChpSHioXM1Pnu + 0.3166666666e1 * ChpCHioE3 * sP1alpP2Psi * ChpC1ioXS1thXC2th * ChpSHioXM1Pnu - 0.8888888889e-1 * ChpCHioE3 * sP1alpP2Psi * ChpS3th * ChpSHioXM1Pnu - 0.2861111112e0 * ChpCHioE3 * sP1alpP2Psi * Chp5S1thPS3th * ChpSHioXM1Pnu + 0.1250000000e0 * ChpCHio * sM1alpP2Psi * ChpC2ioXS1th * ChpSHioE3XM1P3nu - 0.8888888889e-1 * ChpCHio * sM1alpP2Psi * ChpS3th * ChpSHioE3XM1P3nu - 0.3166666666e1 * ChpCHio * sM1alpP2Psi * ChpC1ioXS1thXC2th * ChpSHioE3XM1P3nu + 0.2916666666e0 * ChpCHio * sM1alpP2Psi * ChpC2ioXS3th * ChpSHioE3XM1P3nu - 0.2861111112e0 * ChpCHio * sM1alpP2Psi * Chp5S1thPS3th * ChpSHioE3XM1P3nu - 0.8333333333e0 * ChpCHio * sM1alpP2Psi * ChpC1ioXS1th * ChpSHioE3XM1P3nu + 0.5000000001e0 * ChpCHioE5 * ChpSHioXM1Pnu * sP3alpP2Psi * ChpC1ioXS3th - 0.2666666667e0 * ChpCHioE5 * ChpSHioXM1Pnu * sP3alpP2Psi * ChpS3th + 0.2500000000e1 * ChpC1ioXC1th * ChpM1P3nu * ChpS1thE2XS1ioE2 * sP0alpP2Psi - 0.5333333334e0 * ChpCHioE7 * Chp1P3C2th * ChpSHioXM1Pnu * sP3alpP4Psi * Chp5S1thPS3th + 0.6666666667e0 * ChpC1thXS1thE2XM1P3nu * ChpCHioE6XSHioE2 * sP4alpP2Psi - 0.1688888888e1 * ChpCHioE3 * ChpSHio * sP1alpP2Psi * ChpS3th + 0.1688888888e1 * ChpCHioE3 * ChpSHio * sP1alpP2Psi * Chp5S1thPS3th - 0.1688888888e1 * ChpCHio * ChpSHioE3 * sM1alpP2Psi * ChpS3th + 0.1688888888e1 * ChpCHio * ChpSHioE3 * sM1alpP2Psi * Chp5S1thPS3th + 0.8000000000e0 * ChpCHioE3 * ChpSHioE5 * ChpM1P3nu * sM1alpP4Psi * Chp5S1thPS3th + 0.8533333333e1 * ChpCHioE3 * ChpSHioE5 * ChpM1P3nu * sM1alpP4Psi * ChpS3th + 0.2666666667e1 * ChpCHioE8 * ChpC1thXS1thE2XM1P3nu * sP4alpP4Psi + 0.5333333334e0 * ChpCHioE7 * Chp1P3C2th * ChpSHioXM1Pnu * sP3alpP4Psi * ChpS3th - 0.4222222222e1 * ChpC1th * ChpSHioE4 * sM2alpP2Psi + 0.4222222222e1 * ChpC1th * ChpCHioE4 * sP2alpP2Psi;
	
	
    
    
    h1PNSOc = ((-0.5000000000e0 * sM1alpP1Psi * ChpC1thE2 + 0.5000000000e0 * sP1alpP1Psi * ChpC1thE2 + 0.5000000000e0 * sM1alpP1Psi * ChpC1ioXC1thE2 + 0.5000000000e0 * sP1alpP1Psi * ChpC1ioXC1thE2 + 0.2000000000e0 * sP0alpP1Psi * Chp5S1thPS3th * ChpC1thXS1io - 0.2000000000e0 * sP0alpP1Psi * ChpS3th * ChpC1thXS1io) * chis_x + (0.5000000000e0 * cP1alpP1Psi * ChpC1io + 0.5000000000e0 * cP1alpP1Psi + 0.1000000000e1 * ChpSHioE2 * cM1alpP1Psi) * chis_y + (-0.1000000000e0 * sP1alpP1Psi * ChpC1thX5S1thPS3th - 0.1000000000e0 * sP1alpP1Psi * Chp5S1thPS3thXC1ioXC1th + 0.1000000000e0 * sM1alpP1Psi * ChpC1thX5S1thPS3th - 0.1000000000e0 * sM1alpP1Psi * Chp5S1thPS3thXC1ioXC1th + 0.1000000000e0 * sM1alpP1Psi * ChpC1thXC1ioXS3th + 0.1000000000e0 * sP1alpP1Psi * ChpC1thXC1ioXS3th - 0.1000000000e0 * sM1alpP1Psi * ChpS3thXC1th + 0.1000000000e0 * sP1alpP1Psi * ChpS3thXC1th - 0.4000000000e-1 * sP0alpP1Psi * ChpS1io * Chp5S1thPS3th * Chp5S1thPS3th + 0.8000000000e-1 * sP0alpP1Psi * Chp5S1thPS3th * ChpS3thXS1io - 0.4000000000e-1 * sP0alpP1Psi * ChpS3th * ChpS3thXS1io) * chis_z) * del + (-0.5000000000e0 * sM1alpP1Psi * ChpC1thE2 + 0.5000000000e0 * sP1alpP1Psi * ChpC1thE2 + 0.5000000000e0 * sM1alpP1Psi * ChpC1ioXC1thE2 + 0.5000000000e0 * sP1alpP1Psi * ChpC1ioXC1thE2 + 0.2000000000e0 * sP0alpP1Psi * Chp5S1thPS3th * ChpC1thXS1io - 0.2000000000e0 * sP0alpP1Psi * ChpS3th * ChpC1thXS1io) * chia_x + (0.5000000000e0 * cP1alpP1Psi * ChpC1io + 0.5000000000e0 * cP1alpP1Psi + 0.1000000000e1 * ChpSHioE2 * cM1alpP1Psi) * chia_y + (-0.1000000000e0 * sP1alpP1Psi * ChpC1thX5S1thPS3th - 0.1000000000e0 * sP1alpP1Psi * Chp5S1thPS3thXC1ioXC1th + 0.1000000000e0 * sM1alpP1Psi * ChpC1thX5S1thPS3th - 0.1000000000e0 * sM1alpP1Psi * Chp5S1thPS3thXC1ioXC1th + 0.1000000000e0 * sM1alpP1Psi * ChpC1thXC1ioXS3th + 0.1000000000e0 * sP1alpP1Psi * ChpC1thXC1ioXS3th - 0.1000000000e0 * sM1alpP1Psi * ChpS3thXC1th + 0.1000000000e0 * sP1alpP1Psi * ChpS3thXC1th - 0.4000000000e-1 * sP0alpP1Psi * ChpS1io * Chp5S1thPS3th * Chp5S1thPS3th + 0.8000000000e-1 * sP0alpP1Psi * Chp5S1thPS3th * ChpS3thXS1io - 0.4000000000e-1 * sP0alpP1Psi * ChpS3th * ChpS3thXS1io) * chia_z;
	
    
	//  h1PNc +=  h1PNSOc;
    
    
	
	
	/*   
	 // test
	 h0PNp=0.0;
	 h0PNc=0.0;
	 h05PNp=0.0;
	 h05PNc=0.0; 
	 h1PNp = 0.0; 
	 h1PNc = 0.0;
	 */
    
}


/*
double LCGWSpinBBHHHarm1::halfhann(double t, double t0, double t1){
	
 	double w;
	
 	if(t > t0 && t > t1){
 		if(t0 > t1) return 0;
 		else return 1;
 	} else if(t < t0 && t < t1){
 		if(t0 < t1) return 0;
 		else return 1;
 	} else {
 		w = cos(M_PI*(t-t1)/(2.0*(t0-t1)));
 		return w*w;
 	}
	
}
*/


double LCGWSpinBBHHHarm1::ComputeEnergy(){
	
	/*! We compute the orbial energy to 3PN order in non-spinnin part 
	 * 1PN beyond leadin in s-o interaction
	 */
	
	double Mom = om*M;
	double Mom23 = pow(Mom, 2./3.);
	double Mom43 = Mom23*Mom23;
	double Mom53 = Mom23 *Mom;
	double S1S2 = S1x*S2x + S1y*S2y + S1z*S2z;
	double LS1 = Lnx*S1x + Lny*S1y + Lnz*S1z;
	double LS2 = Lnx*S2x + Lny*S2y + Lnz*S2z;
	
	double SOr = 2.*( x1*( 4./3.*m1M2 + eta )*LS1 + x2*( 4./3.*m2M2 + eta )*LS2 );
	
	double SOr1PN = x1*( m1M2*(8.-31./9.*eta) + (3.-10./3.*eta)*eta )+ \
	x2*( m2M2*(8.-31./9.*eta) + (3.-10./3.*eta)*eta );
	
	double SSp = (1./eta)*(x1*x2*m1M2*m2M2)*(S1S2 - 3.0*LS1*LS2);
	
	double Energy = -0.5*mu*Mom23*( 1.0 - Mom23*(9.0+eta)/12.0 + SOr*Mom \
								   - Mom43*(81.0-57.0*eta+eta*eta)/24.0 + Mom43*SSp +  Mom53*SOr1PN - \
								   Mom*Mom*(675./64. -(34445./576. - M_PI*M_PI*205./96.)*eta +\
											155./96.*eta*eta - 35./5184.*eta*eta*eta ) );
	
	return(Energy);
}






// end of LISACODE-GWSpinBBHHHarm1.cpp