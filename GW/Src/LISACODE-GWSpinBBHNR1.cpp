/*
 *  LISACODE-GWSpinBBHNR1.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 16/08/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-GWSpinBBHNR1.h"


// *********************
// ***  Constructor  ***
// *********************

LCGWSpinBBHNR1::LCGWSpinBBHNR1()
: LCGW()
{  

	LnN.setTools(MT);
	S1N.setTools(MT);
	S2N.setTools(MT);
	JN.setTools(MT);
	LnB.setTools(MT);
	S1B.setTools(MT);
	S2B.setTools(MT);

    
	initNULL(false);
	LoadNRdata();
}


LCGWSpinBBHNR1::LCGWSpinBBHNR1(LCTools * MT_n)
: LCGW(MT_n)
{

	LnN.setTools(MT);
	S1N.setTools(MT);
	S2N.setTools(MT);
	JN.setTools(MT);
	LnB.setTools(MT);
	S1B.setTools(MT);
	S2B.setTools(MT);
  
	
	initNULL(false);
	LoadNRdata();
}


LCGWSpinBBHNR1::~LCGWSpinBBHNR1()
{
	initNULL(true);
}




// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCGWSpinBBHNR1::initNULL(bool CleanMem)
{
	initNULLBase(CleanMem);
	
	if(CleanMem){
		
		if(GWSBHHH != NULL)
			delete GWSBHHH;
		
		if(NRdata != NULL){
			for(int iH=0; iH<NHarm; iH++){
				if(NRdata[iH] != NULL){
					for(int iC=0; iC<2; iC++){
						if(NRdata[iH][iC] != NULL){
							MT->Free(NRdata[iH][iC], NRdataN*sizeof(double));
						}
					}
					MT->Free(NRdata[iH], 2*sizeof(double*));
				}
			}
			MT->Free(NRdata, NHarm*sizeof(double**));
		}
		
		if(lHarm != NULL)
			MT->Free(lHarm, NHarm*sizeof(int));
		
		if(mHarm != NULL)
			MT->Free(mHarm, NHarm*sizeof(int));
		
		if(SpherHarmVal != NULL)
			MT->Free(SpherHarmVal, NHarm*sizeof(dcomplex));
		
		if(NRFileName != NULL)
			MT->Free(NRFileName, (strlen(NRFileName)+1) * sizeof(char));
	}else{
		NRFileName = NULL;
		MT->stripcopy("NRdata.txt", NRFileName);
	}
	
	NParams = 22;
	
	
	
	Mtot = 1.;
	tc = 0.;
	DL = 0.;
	Thd = 0.;
	Phd = 0.;
    NRPolarization = 0.;
    RandPolarization = 0.;
	
	mrat = 1.;
	
	chi1 = 0.;
	chi2 = 0.;
	thS1 = 0.;
	thS2 = 0.;
	phS1 = 0.;
	phS2 = 0.;
	Phi0 = 0.;
	thL = 0.;
	phL = 0.;
	
	m1 = 0.5;
	m2 = 0.5;
	
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
	
	M = 1.;
	dist = 1.;
	thBS1 = 0.;
	thBS2 = 0.;
	phBS1 = 0.;
	phBS2 = 0.;
	thBL = 0.;
	phBL = 0.;
	
	
	
	NHarm = 0;
	lHarm = NULL;
	mHarm = NULL;
	NRdata = NULL;
	NRdataN = 0;
	NRdatat0 = 0.;
	NRdatatend = 0.;
	NRdatadt = 0.;
	tInitObs = 0.;
	toMInitHyb = 0.;
	omMInit = 0.;
	tShiftHybObs = 0.;
	SpherHarmVal = NULL;
	FreqMaxM = 0.5;
	
	GWSBHHH = NULL;
	tTaper = 0.;
	tendTaper = 0.;
	ApplyTaper = false;
	PSIN = 0.;
	tInfot0 = 0.;
	tInfotMaxDiff = 0.;
	
    dhidtau = NULL;
	
	tStartWave = 0.;
	tEndWave = 0.;
	
	t0Start = 0.;
	maxDur = 0.;
	dt = 1.;
	
	
	Amp = 0.;
	thetaS = 0.;
	phiS = 0.;
	cThS = 0.;
	sThS = 1.;
	c2psi = 0.;
	s2psi = 1.;
	
	
}


// ***************************
// *  Configuration methods  *
// ***************************


void LCGWSpinBBHNR1::config(ezxml_t xmlbloc)
{
	//! *** Read sky position, polarization and name
	configBase(xmlbloc);
	
	//! ** Reset declination to compute quantities related to sky position 
	setParam(0, Beta);
	
	
	ezxml_t param;
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		
		//! *** Read ecliptic lattitude 
		if(MT->wcmp(ezxml_attr(param,"Name"),"EclipticLatitude"))
			setParam(0, MT->gXMLAngle(param));
		
		//! *** Read ecliptic colattitude 
		if(MT->wcmp(ezxml_attr(param,"Name"),"EclipticColatitude"))
			setParam(0, M_PI/2. - MT->gXMLAngle(param));
		
		//! *** Read ecliptic longitude
		if(MT->wcmp(ezxml_attr(param,"Name"),"EclipticLongitude"))
			setParam(1, MT->gXMLAngle(param));
		
		//! *** Read polarization
		if(MT->wcmp(ezxml_attr(param,"Name"),"Polarization"))
			setParam(2, MT->gXMLAngle(param));
		
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"TotalMass"))
			setParam(3, MT->gXMLAstroMass(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"CoalescenceTime"))
			setParam(4, MT->gXMLTime(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Distance"))
			setParam(5, MT->gXMLAstroDistance(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"PolarAngleTotalMomentuminSBB"))
			setParam(6, MT->gXMLAngle(param) );
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"AzimuthalAngleTotalMomentuminSBB"))
			setParam(7, MT->gXMLAngle(param) );
		
		//if(MT->wcmp(ezxml_attr(param,"Name"),"FileName"))
		//	MT->stripcopy((*param).txt,NRFileName);
		
	}
}


void LCGWSpinBBHNR1::config(int iParam, double ParamVal)
{
	
}


// ********************
// *  Access methods  *
// ********************


double LCGWSpinBBHNR1::getParam(int iP)
{
	ComputeExtraParam();
	
	switch (iP) {
			//! *** Free parameters
		case 0:
			return (Beta);
			break;
		case 1:
			return (Lambda);
			break;
		case 2:
			return (Polarization);  // ToDO we have to cancel it, but it doesn't affect the code even if we decide to let it here
			break;
		case 3:
			return (Mtot);
			break;
		case 4:
			return (tc);
			break;
		case 5:
			return (DL);
			break;
		case 6:
			return (thetaJ);    // return (Thd);   // changed by Sofiane
			break;
		case 7:
			return (phiJ);  // return (Phd);    // changed by Sofiane
			break;
			
			//! *** Fixed parameters
		case 8:
			return (mrat);
			break;
		case 9:
			return (chi1);
			break;
		case 10:
			return (chi2);
			break;
		case 11:
			return (thS1);
			break;
		case 12:
			return (thS2);
			break;
		case 13:
			return (phS1);
			break;
		case 14:
			return (phS2);
			break;
		case 15:
			return (Phi0);
			break;
		case 16:
			return (thBS1);
			break;
		case 17:
			return (thBS2);
			break;
		case 18:
			return (phBS1);
			break;
		case 19:
			return (phBS2);
			break;
		case 20:
			return (thBL);
			break;
		case 21:
			return (phBL);
			break;
			
		default:
			Cout << "ERROR in LCGWSpinBBHNR1::getParam : The parameter " << iP << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWSpinBBHNR1::getParam : The parameter is unknown !");
	}
	return(0.);
}


void LCGWSpinBBHNR1::setParam(int iP, double Param_n)
{
	
	switch (iP) {
		case 0:
			Beta = Param_n;
			break;
		case 1:
			Lambda = Param_n;
			break;
		case 2:
			RandPolarization = Param_n;
			//RandPolarization = 4.0;
			break;
		case 3:
			Mtot = Param_n;
			M = Mtot * LC::TSUN;
			break;
		case 4:
			tc = Param_n;
			break;
		case 5:
			DL = Param_n;
			dist = DL*LC::kpc_s;
			break;
		case 6:
			 thetaJ = Param_n; //    Thd = Param_n;  // changed by Sofiane
			break;
		case 7:
			phiJ = Param_n;     //   Phd = Param_n;  // changed by Sofiane
			break;
        case 8:
			mrat = Param_n;     //   mrat = Param_n;  // Added by Sofiane
			break;
        case 9:
			chi1 = Param_n;     //   Phd = Param_n;  // Added by Sofiane
			break;
        case 10:
			chi2 = Param_n;     //   Phd = Param_n;  // Added by Sofiane
			break;
			
		case 16:
			thBS1 = Param_n;
			break;
		case 17:
			thBS2 = Param_n;
			break;
		case 18:
			phBS1 = Param_n;
			break;
		case 19:
			phBS2 = Param_n;
			break;
		case 20:
			thBL = Param_n;
			break;
		case 21:
			phBL = Param_n;
			break;
			
		default:
			Cout << "ERROR in LCGWSpinBBHNR1::setParam : The parameter " << iP << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWSpinBBHNR1::setParam  : The parameter is unknown !");
	}
	
	
	
	if((iP==0)||(iP==1)){
		thetaS = M_PI/2 - Beta;
		phiS = Lambda;
		cThS = cos(thetaS);
		sThS = sin(thetaS);
	}
	
	if((iP==16)||(iP==17))
		S1B.xyz(sin(thBS1)*cos(phBS1), sin(thBS1)*sin(phBS1), cos(thBS1));
	
	if((iP==18)||(iP==19))
		S2B.xyz(sin(thBS2)*cos(phBS2), sin(thBS2)*sin(phBS2), cos(thBS2));
	
	if((iP==20)||(iP==21))
		LnB.xyz(sin(thBL)*cos(phBL), sin(thBL)*sin(phBL), cos(thBL));
	
	
	if((iP==0)||(iP==1)||(iP==2)||(iP==3)||(iP==6)||(iP==7)||((iP>=16)&&(iP<=21)))
		NeedExtraParamCompute = true;
	
	if((iP==2)||(iP==6)||(iP==7))
		TypeExtraParamCompute = 1;  
	
	if((iP>=16)&&(iP<=21))
		TypeExtraParamCompute = 2;  
	
	Cout << "TypeExtraParamCompute for param " << iP << " = " << TypeExtraParamCompute << Endl;
}


void LCGWSpinBBHNR1::getRange(int iP, double &Pmin, double &Pmax)
{
	double MtotMin(1.e5);
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
			Pmin = MtotMin;
			Pmax = 7.e7;   // changed by Sofiane
			//B//Pmin = 1.e6;
			//B//Pmax = 1.e8;
			break;
		case 4:
			Pmin = 0.0;
			//B//Pmin = 7.8e6;
			Pmax = (NRdatatend - NRdatat0 ) * (MtotMin*LC::TSUN);
			Cout << "Pmax = " << Pmax << Endl;
			break;
		case 5:
			Pmin = 1.0e3;
			//B//Pmax = 1.e7;
			Pmax = 2.5e8;
			break;
		case 6:
			Pmin = 0.0;
			Pmax = M_PI;
			break;
		case 7:
			Pmin = 0.0;
			Pmax = 2.0*M_PI;
			break;
		default:
			Cout << "ERROR in LCGWSpinBBHNR1::getRange : The parameter " << iP << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWSpinBBHNR1::getRange : The parameter is unknown !");
			break;
	}
}


double LCGWSpinBBHNR1::getDelta(int iP)
{
	double tmpMin, tmpMax, DRes;
	getRange(iP, tmpMin, tmpMax);
	switch (iP){
		case 0:
			DRes = 2.e-4;   // changed by Sofiane
          //  DRes = 2.e-6;
			break;
		case 1:
			DRes = 2.e-4;   // changed by Sofiane
			break;
		case 2:
			DRes = 1.e-5;
			break;
		case 3:
			DRes = 3.e-1;
			break;
		case 4:
			DRes = 1.e-2;
			break;
		case 6:
			DRes = 2.e-4;   // changed by Sofiane
			break;
		case 7:
			DRes = 2.e-4;   // changed by Sofiane
			break;
			
		default:
			DRes = (tmpMax-tmpMin)*5.0e-8;
	};
	//Cout << iP << " --> Delta = " << DRes<< Endl;
	return (DRes);
}



double LCGWSpinBBHNR1::getSpecialParam(int iPS)
{
	switch (iPS){
			
		default:
			Cout << "ERROR in LCGWSpinBBHNR1::getSpecialParam : The parameter " << iPS << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWSpinBBHNR1::getSpecialParam : The parameter is unknown !");
			break;
	}
	return(0.);
}


void LCGWSpinBBHNR1::setSpecialParam(int iPS, double SpecParam_n)
{
	switch (iPS){
		case 0 :
			if (SpecParam_n > 0)
				ApplyTaper = true;
			else
				ApplyTaper = false;			
			break;
            
        case 1:
            OutputType = 0; //! 0--> Standard : \f$ h_{+,\times} \f$ 
            if((SpecParam_n>=0.5)&&(SpecParam_n<=1.5))
                OutputType = 1; //! 1 --> For tc derivative : \f$ \partial h_{+,\times} / \partial \tau \f$ 
            if((SpecParam_n>=1.5)&&(SpecParam_n<=2.5))
                OutputType = 2; //! 2 --> For total mass derivative : \f$ \partial h_{+,\times} / \partial \tau \f$             
            break;
            
		default:
			Cout << "ERROR in LCGWSpinBBHNR1::setSpecialParam : The parameter " << iPS << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWSpinBBHNR1::setSpecialParam : The parameter is unknown !");
			break;
	}
}



void LCGWSpinBBHNR1::RandParam(int iP)
{
	double Pmin(0.0);
	double Pmax(1.0);
	
	//! ** Only for free parameters 
	if(iP<8){
		getRange(iP, Pmin, Pmax);
		setParam(iP,MT->RandUniform(Pmin, Pmax));
		
		// ** Specifical case : non-uniform
		if(iP == 0)
			setParam(0, M_PI/2.0 - acos(MT->RandUniform(-1.0, 1.0)));
		if(iP == 3)
			setParam(3, pow(10.,(MT->RandUniform(log10(Pmin), log10(Pmax)))) );
	}
	
}

void LCGWSpinBBHNR1::setTimeInfo(double t0_n, double dt_n, double TObs_n, double tAskMin_n, double tDOrbMax_n,  double tMaxDiff)
{
    t0Real = t0_n; 
	tAskMin = tAskMin_n;
	tDOrbMax = tDOrbMax_n;
	
	//! *** Define the integration step and the number of steps to store
	dt = dt_n;
	
	/*! ** Set time for starting integration to be sure that we have enough data and be able to manage any kind of GW direction relative to detector position :
	 *	Actually the reference time for GW is defined at the solar system barycenter and the time of the simulation at the detector */
	t0Start = tAskMin_n;
	maxDur = TObs_n;
	
	tInfot0 = t0_n;
	tInfotMaxDiff = tMaxDiff;
	
}


// ***************************************
// * Linking and initialization methods  *
// ***************************************

int LCGWSpinBBHNR1::init()
{	
	initBase();
#ifdef _DEBUG_GW_    
	char fNCheck[512];
    int iRCheck(MT->ifloor(MT->RandUniform(0, 10000.)));
    sprintf(fNCheck,"CheckGWSpinBBHNR1_%d.txt",iRCheck);
    std::cerr << "DEBUG:GW:LCGWSpinBBHNR1 : File = " << fNCheck << Endl;
    DEBUGfCheck = new std::ofstream(fNCheck);
    DEBUGDispAll=true;
#endif
	double t0tcM(9.);
	
	ComputeExtraParam();
    
	if(GWSBHHH == NULL)
		GWSBHHH = new LCGWSpinBBHHHarm1(MT);
	GWSBHHH->setParam(0, Beta);
	GWSBHHH->setParam(1, Lambda);
	GWSBHHH->setParam(4, tc);
	GWSBHHH->setParam(5, DL);
	GWSBHHH->setParam(6, chi1);
	GWSBHHH->setParam(7, chi2);
	GWSBHHH->setParam(8, thBS1);
	GWSBHHH->setParam(9, thBS2);
	GWSBHHH->setParam(10, phBS1);
	GWSBHHH->setParam(11, phBS2);
	GWSBHHH->setParam(12, Phi0);
	GWSBHHH->setParam(13, thBL);
	GWSBHHH->setParam(14, phBL);
	GWSBHHH->setParam(15, Mtot*pow(m1*m2/(M*M), 0.6));
	GWSBHHH->setParam(17, m1*m2/(M*M));
    GWSBHHH->setParam(21, thetaJ);    // added by Sofiane
	GWSBHHH->setParam(22, phiJ);      // added by Sofiane
	GWSBHHH->setSpecialParam(3, omMInit);
	GWSBHHH->setTimeInfo(tInfot0, dt, maxDur, tAskMin, tDOrbMax, tInfotMaxDiff );
	GWSBHHH->init();
	GWSBHHH->DispInfo(" >>> SpinBBH in NR <<< ");
	
	//! *** Copy the informations about the taper from inspiral waveform
	GWSBHHH->getTaperInfo(tTaper,tendTaper,FreqMaxTaper);
   
    
	
	
	Amp =  M / dist;
    
	
	//tShifttc = 5. * t0tcM*t0tcM*t0tcM*t0tcM  / 256. * M / eta;
	
	
	/*! *** Determine the relation between the hybrid time and observation time using the following condition : 
	 *	At the initial data point, the hybrid orbital frequency should match the one based on PN 
	 */
	
	//double LS1, LS2, S1S2, m1M2, m2M2, eta;
	double DomLow, DomHigh, DomMean, tLow(0.), tHigh(tc), tMean, omMin, omMax;
	
	//m1M2 = pow( (1.+1./mrat) ,-2.);
	//m2M2 = pow( (1.+mrat) ,-2.);
	
	//LS1 = LnB.scalar(S1B);
	//LS2 = LnB.scalar(S2B);
	//S1S2 = S1B.scalar(S2B);
	
	
	
	tLow = 0.;
	tHigh = tc - 5. / 256. * M / eta;
	
	omMin = GWSBHHH->ComputeApprOrbFreq(tLow, tc);
	//Cout << "omMin = " << omMin << " vs " << ComputeApprOrbFreq(tLow, tc, chi1, chi2, LS1, LS2, S1S2, m1M2, m2M2, eta) << Endl;
	
	omMax = GWSBHHH->ComputeApprOrbFreq(tHigh, tc);
	//Cout << "omMax = " << omMax << " vs " << ComputeApprOrbFreq(tHigh, tc, chi1, chi2, LS1, LS2, S1S2, m1M2, m2M2, eta) << Endl;
	
	DomLow  = omMin - omMInit/M;
	DomHigh = omMax - omMInit/M;
	
	MT->o->precision(12);
	
	if( DomLow*DomHigh <= 0.0 ){
		Cout << "Looking for the time corresponding to initial omega " << omMInit/M << " in the range [ " << omMin << " , " << omMax << " ] corresponding to the time range [ " << tLow << " , " << tHigh << " ] ... " << Endl;
		do{
			tMean = (tHigh+tLow)/2.0;
			DomMean = GWSBHHH->ComputeApprOrbFreq(tMean, tc) - omMInit/M;
			if(DomLow*DomMean>0){
				tLow = tMean;
				DomLow = DomMean;
			}else{
				tHigh = tMean;
				DomHigh = DomMean;
			}
			//Cout << "DomLow ( " << tLow << " ) = " << DomLow << "   -   DomHigh ( " << tHigh << " ) = " << DomHigh << Endl; 
		}while(fabs(DomMean/(omMInit/M))>1.e-6);
		Cout << "\t ==> Found initial omega " << omMInit/M << " at t = " << tMean << " s ( precision = " << DomMean << " )" << Endl; 
		tInitObs = tMean;
	}else {
		Cout << "WARNING : The required value of omega " << omMInit/M << " is not in the waveform possible value, i.e. in the range [ " << omMin << " , " << omMax << " ] corresponding to the time range [ " << tLow << " , " << tHigh << " ]" << Endl;
		tInitObs = tc - 5. * 9.*9.*9.*9.  / 256. * M / eta ; // Default value corresponding at 9 M
	}
	
	
	
	
	
	tShiftHybObs = tInitObs - toMInitHyb*M;
	
	Cout << "Time relation between hybrid and observation (matching of initial orbital frequency " << omMInit <<  " with PN orb. freq. ) :" << Endl;
	Cout << "\t + Initial hybrid time      = " << toMInitHyb*M << " s = " << toMInitHyb << " M" << Endl;
	Cout << "\t + Initial observation time = " << tInitObs << " s = " << tInitObs/M << " M" << Endl;
	Cout << "\t + Time between initial and merger (observation time) = " << (tc-tInitObs) << " s = " << (tc-tInitObs)/M << " M" << Endl;
	Cout << "\t + Shift between observation and hybrid = " << tShiftHybObs << " s = " << tShiftHybObs/M << " M" << Endl;
    Cout << " NRdatadt*M  =   " << NRdatadt*M << " sec" << Endl;  // added by Sofiane
 	
	
	
	//! ** Find the range of time when the waveform is define
	tStartWave = 0.;
	tEndWave = 0.;
	if(NRdata[0]!=NULL){
		//	tStartWave = tc  + (NRdatat0+10.*NRdatadt)*M;
        tStartWave = tShiftHybObs + (NRdatat0+10.*NRdatadt)*M;
        
        
		//Cout << tc << Endl;
		//Cout << NRdata[0]->getxend() << "  ==>  " << NRdata[0]->getxend()*M <<  Endl;
		//Cout << NRdata[0]->getdx() << "  ==>  " << NRdata[0]->getdx()*M << Endl;
		
		//    tEndWave = tc  + (NRdatatend-10.*NRdatadt)*M;
		tEndWave = tShiftHybObs + (NRdatatend-10.*NRdatadt)*M;
	}
	Cout << "Waveform defined in time range = [ " << tStartWave << " : " << tEndWave << " ] sec" << Endl;
	
	//! ** Compute and store spherical harmonics
  
	for(int iH=0; iH<NHarm; iH++){
		if(NRdata[iH] != NULL)
			SpherHarmVal[iH] = ComputeSpherHarm(lHarm[iH], mHarm[iH], Thd, Phd);
		else
			SpherHarmVal[iH] = 0.;
        
	}
	
	FreqMin = 0.125*pow(0.2*(mrat/((1+mrat)*(1+mrat)))*(tc)/M, -0.375)/(4.*M*M_PI);
	
	if(ApplyTaper){
		FreqMax = FreqMaxTaper;
		Cout << "\t Taper (in LCGWSpinBBHNR1::init) : applied with parameters : tTaper = " << tTaper << " , tend = " << tendTaper <<  " , FreqMax = " << FreqMaxTaper << Endl;
	}else{ 
		FreqMax = FreqMaxM/M;
	}
	
	//FreqMax = 1./(M_PI*M*pow(10.,3./2.));
	
	//Cout << "Mtot = " << M << " s : t = " << tc - pow(10.,4.) * (5.*M)/((mrat/((1+mrat)*(1+mrat)))*256) <<  " ==> FreqMax = " << FreqMax << Endl;
	
	return 0;
}


// *********************
// *  Running methods  *
// *********************

void LCGWSpinBBHNR1::Computehpc(double t)
{
	double hpl, hcr;
	double AmpT;
	dcomplex Yp;
	dcomplex Yn;
	dcomplex img(0.0, 1.0);
	dcomplex h, dhdtau;
	
	hpl = 0.0;
	hcr = 0.0;
	
#ifdef _DEBUG_GW_  
    if((DEBUGDispAll)&&(t>1000.)){
        DispAllParam(MT->o);
        DEBUGDispAll=false;
    }
#endif
    

	
	if ((t>tStartWave)&&(t<=tEndWave)){
		
		
		for(int iH=2; iH<lmax+1; iH++){
			for(int jH=1; jH<iH+1; jH++){
				Yp = SpherHarm(iH, jH);
				Yn = SpherHarm(iH, -jH);
                
                if((OutputType==1)||(OutputType==2))
                    dhdtau = (dhpharmdtau(iH, jH, t) - img * dhcharmdtau(iH, jH, t))*Yp + (dhpharmdtau(iH, -jH, t) - img * dhcharmdtau(iH, -jH, t))*Yn;
                
                if((OutputType==0)||(OutputType==2))
                    h = (hpharm(iH, jH, t) - img * hcharm(iH, jH, t))*Yp + (hpharm(iH, -jH, t) - img * hcharm(iH, -jH, t))*Yn;
                
                //if((t>1000.)&&(t<1100.))
                //    Cout << t << " " << t/M << " "  << OutputType << " " << iH << " " << jH << " " << real(h) << " " << imag(h) << " " << real(dhdtau) << " " << imag(dhdtau) << Endl; 
                
                
                
                if(OutputType==1){
                    h = dhdtau/M;
                }
                
                if(OutputType==2){
                    h = h/M - dhdtau*(t-tShiftHybObs)/(M*M);
                    h *= LC::TSUN;
                }
                    
                //if((OutputType==2)&&(t>100.)&&(jH==2))
                //    Cout << " " << real(h) << " " << imag(h) << Endl; 
                
				hpl += real(h);
				hcr += -imag(h);

				//Cout << iH << " " << jH << " : " << h << " " << hpharm(iH, jH, t) << " " << hcharm(iH, jH, t) << " " << hpl << " " << hcr << " " << hpharm(iH, -jH, t) << " " << hcharm(iH, -jH, t)  << Endl;
			} 
			
		}
		
		/*	
		 if(lmax < 3){
		 Yp = SpherHarm(2, 2);
		 Yn = SpherHarm(2, -2);
		 h = (hpharm(2, 2, t) - img * hcharm(2, 2, t))*Yp + (hpharm(2, -2, t) - img * hcharm(2, -2, t))*Yn;
		 hpl += real(h);
		 hcr += -imag(h);
		 
		 Yp = SpherHarm(2, 1);
		 Yn = SpherHarm(2, -1);
		 h = (hpharm(2, 1, t) - img * hcharm(2, 1, t))*Yp + (hpharm(2, -1, t) - img * hcharm(2, -1, t))*Yn;
		 hpl += real(h);
		 hcr += -imag(h);
		 }
		 
		 
		 if(lmax < 4){
		 Yp = SpherHarm(3, 3);
		 Yn = SpherHarm(3, -3);
		 h = (hpharm(3, 3, t) - img * hcharm(3, 3, t))*Yp + (hpharm(3, -3, t) - img * hcharm(3, -3, t))*Yn;
		 hpl += real(h);
		 hcr += -imag(h);
		 }
		 
		 
		 if(lmax < 5){
		 Yp = SpherHarm(4, 4);
		 Yn = SpherHarm(4, -4);
		 h = (hpharm(4, 4, t) - img * hcharm(4, 4, t))*Yp + (hpharm(4, -4, t) - img * hcharm(4, -4, t))*Yn;
		 hpl += real(h);
		 hcr += -imag(h);
		 }
		 
		 
		 if(lmax < 6){
		 Yp = SpherHarm(5, 5);
		 Yn = SpherHarm(5, -5);
		 h = (hpharm(5, 5, t) - img * hcharm(5, 5, t))*Yp + (hpharm(5, -5, t) - img * hcharm(5, -5, t))*Yn;
		 hpl += real(h);
		 hcr += -imag(h);
		 }
		 
		 
		 if(lmax < 7){
		 Yp = SpherHarm(6, 6);
		 Yn = SpherHarm(6, -6);
		 h = (hpharm(6, 6, t) - img * hcharm(6, 6, t))*Yp + (hpharm(6, -6, t) - img * hcharm(6, -6, t))*Yn;
		 hpl += real(h);
		 hcr += -imag(h);
		 }
		 
		 
		 if(lmax < 8){
		 Yp = SpherHarm(7, 7);
		 Yn = SpherHarm(7, -7);
		 h = (hpharm(7, 7, t) - img * hcharm(7, 7, t))*Yp + (hpharm(7, -7, t) - img * hcharm(7, -7, t))*Yn;
		 hpl += real(h);
		 hcr += -imag(h);
		 }
		 
		 */
		
		//!  h+ and hx in SSB are compute in LCGW::ComputehBpc(double t) of LISACODE-GW.cpp
		
		//hBpLast = - Amp * (hpl * c2psi + hcr * s2psi);   // Are they defined before ?
		//hBcLast = - Amp * (-hpl * s2psi + hcr * c2psi);  // Are they defined before ?
		
		if(ApplyTaper)
			AmpT = Amp*GWSBHHH->halfhann(t, tendTaper, tTaper);
		else
			AmpT = Amp;
		
	//	hBpLast = - AmpT * hpl;
	//	hBcLast = - AmpT * hcr;

        
        hBpLast = - Amp * (hpl * c2psi + hcr * s2psi);   // Are they defined before ?
		hBcLast = - Amp * (-hpl * s2psi + hcr * c2psi);  // Are they defined before ?
		
#ifdef _DEBUG_GW_     
        (*DEBUGfCheck) << " " << hpl << " " << hcr << " " << hBpLast << " " << hBcLast  ;
#endif
        
		
	}else{
		hBpLast = 0.0;
		hBcLast = 0.0;
	}
    
	
}


// *******************
// *  Other methods  *
// *******************

void LCGWSpinBBHNR1::DispInfo(char * BTab)
{
	ComputeExtraParam();
	
	MT->o->precision(12);
	
	if(MT->Disp()){
		DispInfoBase(BTab);
	}
    
    
    
	
	Cout << BTab << "\tFree parameters :" << Endl;
	Cout << BTab << "\t- Beta    = " << Beta << " rad "<< Endl;
	Cout << BTab << "\t- Lambda  = " << Lambda << " rad "<< Endl;
    Cout << BTab << "\t- phiJ    = " << phiJ << " rad "<< Endl;
     Cout << BTab << "\t- thetaJ    = " << thetaJ << " rad "<< Endl;
	Cout << BTab << "\t- Polarization = " << NRPolarization << " rad  (tan = " << tan(NRPolarization) << ")" << Endl;
	Cout << BTab << "\t- Mtot  = " << Mtot << " MSun  ==> m1 = " << Mtot/(1.+1./mrat) << " MSun , m2 = " << Mtot/(mrat+1.) << " MSun" << Endl;  
	Cout << BTab << "\t- tc    = " << tc << " s "<< Endl;
	Cout << BTab << "\t- DL    = " << DL << " kpc "<< Endl;
	Cout << BTab << "\t- Pol. ang. of detector in src = " << Thd << " rad "<< Endl;
	Cout << BTab << "\t- Az.  ang. of detector in src = " << Phd << " rad "<< Endl;
	Cout << BTab << "\tFixed parameters read/computed from file describing  numerical relativity data : " << NRFileName << Endl;
	Cout << BTab << "\t- q    = " << mrat << Endl;
	Cout << BTab << "\t- chi1 = " << chi1 << Endl;
	Cout << BTab << "\t- chi2 = " << chi2 << Endl;
	Cout << BTab << "\t- thS1 in src frame = " << thS1 << " rad " << Endl;
	Cout << BTab << "\t- thS2 in src frame = " << thS2 << " rad " << Endl;
	Cout << BTab << "\t- phS1 in src frame = " << phS1 << " rad " << Endl;
	Cout << BTab << "\t- phS2 in src frame = " << phS2 << " rad " << Endl;
	Cout << BTab << "\t- Inital phase = " << Phi0 << " rad " << Endl;
	Cout << BTab << "\t- thS1 in SSB = " << thBS1 << " rad " << Endl;
	Cout << BTab << "\t- thS2 in SSB = " << thBS2 << " rad " << Endl;
	Cout << BTab << "\t- phS1 in SSB = " << phBS1 << " rad " << Endl;
	Cout << BTab << "\t- phS2 in SSB = " << phBS2 << " rad " << Endl;
	Cout << BTab << "\t- thBL in SSB  = " << thBL << " rad " << Endl;
	Cout << BTab << "\t- phBL in SSB  = " << phBL << " rad " << Endl;
	Cout << BTab << "\tExtra infos :" << Endl;
	Cout << BTab << "\t- Ln in SSB = " << LnB << Endl;
    Cout << BTab << "\t- S1 in SSB = " << S1B << Endl;
    Cout << BTab << "\t- S2 in SSB = " << S2B << Endl;
	Cout << BTab << "\t- J in NRF  = " << JN << Endl; 
	Cout << BTab << "\t- Ln in NRF = " << LnN << Endl;
    Cout << BTab << "\t- S1 in NRF = " << S1N << Endl;
    Cout << BTab << "\t- S2 in NRF = " << S2N << Endl;
    
	Cout << BTab << "\t- harmonics 2,2    = " << SpherHarm(2,2) << "  "<< Endl;
	Cout << BTab << "\t- harmonics 2,-2    = " << SpherHarm(2,-2) << "  "<< Endl;
    
    Cout << BTab << "\t- PolB recalcule    = " << PSIN << " rad "<< Endl;
	
}


void LCGWSpinBBHNR1::DispAllParam(std::ostream * out)
{
	for(int iP=0; iP<NParams; iP++)
		(*out) << " " << getParam(iP);
}


void LCGWSpinBBHNR1::DispAllParamName(std::ostream * out)
{
	(*out) << " Bet Lam Pol Mtot tc DL Phi0 Inc0";
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

void LCGWSpinBBHNR1::LoadNRdata()
{
	char tmpLine[16384];
	std::vector<char*> Words(0);
	int l, m;
	char fNInHarm[16384];
	LCDataFileRead * tmpDataRead;
	int tmpNRec, tmpNRDataN;
	double tmpNRdatat0, tmpNRdatatend, tmpNRdatadt;
    int NHarmReq;
	

	
	if(ReadNRParam(NRFileName, mrat, eta, chi1, chi2, Phi0, FreqMax, omMInit, toMInitHyb, S1N, S2N, LnN)){
		Cout << "ERROR on LCGWSpinBBHNR1::LoadNRdata : Problem in openning NR main data file " << NRFileName << " ." << Endl;
		throw std::invalid_argument("ERROR on LCGWSpinBBHNR1::LoadNRdata : Problem in openning NR main data file.");
	}
		
	thS1 = S1N.th();
	phS1 = S1N.ph();
	thS2 = S2N.th();
	phS2 = S2N.ph();
	thL = LnN.th();
	phL = LnN.ph();
	
	
	std::ifstream fInMain(NRFileName);
	
	//! **** Pass through fixed parameters 
	int iLine(0);
	do{
		iLine++;
		fInMain.getline(tmpLine,16384,'\n');
		MT->wextract(tmpLine, Words);
	}while((!MT->wcmp(Words[0],"Harmonics"))&&(!fInMain.eof()));
	
	
    NHarmReq=0;
    if((MT->wcmp(Words[0],"Harmonics"))&&(Words.size()>=2)){
        NHarmReq = atoi(Words[1]);
    }else{
        Cout << "ERROR in LCGWSpinBBHNR1::LoadNRdata : We need to know the number of harmonics : 'Harmonics 6' (NEW format !!!)" << Endl;
        throw std::invalid_argument("ERROR in LCGWSpinBBHNR1::LoadNRdata : We need to know the number of harmonics : 'Harmonics 6' (NEW format !!!)" );
    }
    Cout << "Number of harmonics required = " << NHarmReq << Endl;        
    
	//! *** Find the number of harmonics
	lmax = 0;
	NHarm = 0;
	do{
		iLine++;
		fInMain.getline(tmpLine,16384,'\n');
		MT->wextract(tmpLine, Words);		
		l = atoi(Words[0]);
		if(l>lmax)
			lmax = l;
		NHarm++;
	}while((Words.size()>=3)&&(!fInMain.eof())&&(NHarm<NHarmReq));
    Cout << "Number of harmonics read = " << NHarm << Endl; 
 
    if(NHarmReq!=NHarm){
        Cout << "ERROR in LCGWSpinBBHNR1::LoadNRdata :  We do not have the required number of harmonics, " << NHarmReq << " . We found " << NHarm << "harmonics." << Endl;
        throw std::invalid_argument("ERROR in LCGWSpinBBHNR1::LoadNRdata : We do not have the required number of harmonics !" );
    }
	
	fInMain.close();
	fInMain.clear();
	
	
	//! *** Allow memory
	NRdata = (double***) MT->AllocMemory(NHarm*sizeof(double**));
	lHarm = (int*) MT->AllocMemory(NHarm*sizeof(int));
	mHarm = (int*) MT->AllocMemory(NHarm*sizeof(int));
	SpherHarmVal = (dcomplex*) MT->AllocMemory(NHarm*sizeof(dcomplex));
	
	
	//! *** Read harmonic files and load them 
	fInMain.open(NRFileName);
	do{
		fInMain.getline(tmpLine,16384,'\n');
		MT->wextract(tmpLine, Words);	
	}while((!MT->wcmp(Words[0],"Harmonics"))&&(!fInMain.eof()));
    
	for(int iHarm=0; iHarm<NHarm; iHarm++){

		fInMain.getline(tmpLine,16384,'\n');
		MT->wextract(tmpLine, Words);
		
		//Cout << tmpLine << " ==> " << Words[0] << " " << Words[1] << " " << Words[2] << Endl;
		
		lHarm[iHarm] = atoi(Words[0]);
		mHarm[iHarm] = atoi(Words[1]);
		
		tmpDataRead = new LCDataFileRead(MT, Words[2]);
		tmpDataRead->init(NRdata[iHarm], tmpNRec, tmpNRDataN);
		tmpNRdatat0 = tmpDataRead->getx0();
		tmpNRdatatend = tmpDataRead->getxend();
		tmpNRdatadt = tmpDataRead->getdx();
		delete tmpDataRead;
		
		/*
		for(int i=0; i<tmpNRDataN; i++){
			Cout << i ;
			for(int j=0; j<tmpNRec; j++){
				Cout << " " << NRdata[iHarm][j][i];
			}
			Cout << Endl;
		}
		*/
		
		
		
		
		if(tmpNRec != 2){
			Cout << "LCGWSpinBBHNR1::LoadNRdata : The data file should contain only 2 records." << Endl;
			throw std::invalid_argument("LCGWSpinBBHNR1::LoadNRdata : The data file should contain only 2 records.");
		}
		if(iHarm==0){
			NRdataN = tmpNRDataN;
			NRdatat0 = tmpNRdatat0;
			NRdatatend = tmpNRdatatend;
			NRdatadt = tmpNRdatadt;
		}else{
			if(NRdataN != tmpNRDataN){
				Cout << "LCGWSpinBBHNR1::LoadNRdata : The harmonic " <<iHarm<<" (l="<<lHarm[iHarm]<<",m="<<mHarm[iHarm]<<") don't have the same number of data " << tmpNRDataN << " than the previous one " <<iHarm<<" (l="<<lHarm[iHarm]<<",m="<<mHarm[iHarm]<<")  " << NRdataN << "  !" << Endl;
				throw std::invalid_argument("LCGWSpinBBHNR1::LoadNRdata : The harmonic don't have the same number of data than the previous one !");
			}
			if(fabs((NRdatat0-tmpNRdatat0)/(NRdatat0+tmpNRdatat0))>1.e-12){
				MT->o->precision(15);
				Cout << "LCGWSpinBBHNR1::LoadNRdata : The harmonic " <<iHarm<<" (l="<<lHarm[iHarm]<<",m="<<mHarm[iHarm]<<") don't have the same initial time " << tmpNRdatat0 << " than the previous one " <<iHarm<<" (l="<<lHarm[iHarm]<<",m="<<mHarm[iHarm]<<")  " << NRdatat0 << "  !" << Endl;
				throw std::invalid_argument("LCGWSpinBBHNR1::LoadNRdata : The harmonic don't have the same number of data than the previous one !");
			}
			if(fabs((NRdatatend-tmpNRdatatend)/(NRdatatend+tmpNRdatatend))>1.e-12){
				MT->o->precision(15);
				Cout << "LCGWSpinBBHNR1::LoadNRdata : The harmonic " <<iHarm<<" (l="<<lHarm[iHarm]<<",m="<<mHarm[iHarm]<<") don't have the same final time " << tmpNRdatatend << " than the previous one " <<iHarm<<" (l="<<lHarm[iHarm]<<",m="<<mHarm[iHarm]<<")  " << NRdatatend << "  !" << Endl;
				throw std::invalid_argument("LCGWSpinBBHNR1::LoadNRdata : The harmonic don't have the same number of data than the previous one !");
			}
			if(fabs((NRdatadt-tmpNRdatadt)/(NRdatadt+tmpNRdatadt))>1.e-12){
				MT->o->precision(15);
				Cout << "LCGWSpinBBHNR1::LoadNRdata : The harmonic " <<iHarm<<" (l="<<lHarm[iHarm]<<",m="<<mHarm[iHarm]<<") don't have the same time step " << tmpNRdatadt << " than the previous one " <<iHarm<<" (l="<<lHarm[iHarm]<<",m="<<mHarm[iHarm]<<")  " << NRdatadt << "  !" << Endl;
				throw std::invalid_argument("LCGWSpinBBHNR1::LoadNRdata : The harmonic don't have the same number of data than the previous one !");
			}
		}
		MT->MemDisplay();
		
		//for(int i=0; i<10; i++)
		//	Cout << NRdata[iHarm][0][i] << " " << NRdata[iHarm][1][i] << Endl;
		
	}
	
	
	MemUse += NHarm * (sizeof(double**) + 2*sizeof(int) + sizeof(dcomplex) + 2*( NRdataN*sizeof(double) + sizeof(double*)) );
	
	
	for(int iW=0; iW<Words.size(); iW++)
		MT->Free(Words[iW],256*sizeof(char));
}



double LCGWSpinBBHNR1::hpharm(int l, int m, double t)
{
	int iH(IndHarm(l,m));
	if(iH != -1)
        //	return(gHarmDataLIN(iH, 1, (t-tc)/M));
        return(InterLagrangeNR(iH, 0, (t-tShiftHybObs)/M, LCGWSpinBBHNR1_ORDERLAGINTER));  
	else
		return 0.;
}


double LCGWSpinBBHNR1::hcharm(int l, int m, double t)
{
	int iH(IndHarm(l,m));
	if(iH != -1)
		//	return(gHarmDataLIN(iH, 1, (t-tc)/M));
        return(InterLagrangeNR(iH, 1, (t-tShiftHybObs)/M, LCGWSpinBBHNR1_ORDERLAGINTER));  
	else
		return 0.;
}


/*

double LCGWSpinBBHNR1::dhpharmdtau(int l, int m, double t)
{
    return (dhharmdtau(0, l, m, t));
}


double LCGWSpinBBHNR1::dhcharmdtau(int l, int m, double t)
{
    return (dhharmdtau(1, l, m, t));
}

*/




/*
 Added by Sofiane
 */  

double LCGWSpinBBHNR1::dhpharmdtau(int l, int m, double t)
{
    double epsilon(.1);
    double alpha(epsilon/M);
    double hpDer;
    
    hpDer = 0. ;
    hpDer += (3.) * hpharm(l, m, t - 4. * epsilon );
    hpDer += (-32.) * hpharm(l, m, t - 3. * epsilon );
    hpDer += (168.) * hpharm(l, m, t - 2. * epsilon );
    hpDer += (-672.) * hpharm(l, m, t - epsilon );
    hpDer += (672.) * hpharm(l, m, t + epsilon );
    hpDer += (-168.) * hpharm(l, m, t + 2. * epsilon );
    hpDer += (32.) * hpharm(l, m, t + 3. * epsilon );
    hpDer += (-3.) * hpharm(l, m, t + 4. * epsilon );
    hpDer /= 840. ;
   // hpDer /= epsilon ;
    hpDer /= alpha ;
    
    return (hpDer);
}


double LCGWSpinBBHNR1::dhcharmdtau(int l, int m, double t)
{
    double epsilon(.1);
    double alpha(epsilon/M);
    double hcDer;
    
    hcDer = 0. ;
    hcDer += (3.) * hcharm(l, m, t - 4. * epsilon );
    hcDer += (-32.) * hcharm(l, m, t - 3. * epsilon );
    hcDer += (168.) * hcharm(l, m, t - 2. * epsilon );
    hcDer += (-672.) * hcharm(l, m, t - epsilon );
    hcDer += (672.) * hcharm(l, m, t + epsilon );
    hcDer += (-168.) * hcharm(l, m, t + 2. * epsilon );
    hcDer += (32.) * hcharm(l, m, t + 3. * epsilon );
    hcDer += (-3.) * hcharm(l, m, t + 4. * epsilon );
    hcDer /= 840. ;
   // hcDer /= epsilon ;
    hcDer /= alpha ;
    
    return (hcDer);
}


// end Added by Sofiane 







double LCGWSpinBBHNR1::dhharmdtau(int ihphc, int l, int m, double t)
{
    int ib;
	int iH(IndHarm(l,m));
	if(dhidtau==NULL)
        dhidtau = (double*)MT->AllocMemory((LCGWSpinBBHNR1_ORDERLAGINTER+4)*sizeof(double));
    
	if(iH != -1){
        double toM((t-tShiftHybObs)/M);
        
        //! ** Computing derivative with finite coefficients
        double xr((toM-NRdatat0)/NRdatadt);        
        int ibin(MT->ifloor(xr));
        int No2(MT->iceil(LCGWSpinBBHNR1_ORDERLAGINTER/2.+1));
        
        for(int ii=0; ii<(LCGWSpinBBHNR1_ORDERLAGINTER+4) ; ii++){
            ib = ii + ibin - No2;
            dhidtau[ii] = 0.;
            dhidtau[ii] += ( 1./12.)*NRdata[iH][ihphc][ib-2] ;
            dhidtau[ii] += (-2./3. )*NRdata[iH][ihphc][ib-1] ;
            dhidtau[ii] += ( 2./3. )*NRdata[iH][ihphc][ib+1] ;
            dhidtau[ii] += (-1./12.)*NRdata[iH][ihphc][ib+2] ;
            dhidtau[ii] /= NRdatadt;
        }
        
        return(MT->InterLagrange(dhidtau, LCGWSpinBBHNR1_ORDERLAGINTER+3, No2+(xr-ibin), LCGWSpinBBHNR1_ORDERLAGINTER));  
	}else{
		return 0.;
    }
}


int LCGWSpinBBHNR1::IndHarm(int l, int m)
{
	int iH(0);
	// ** Find the index of harmonics
	while ((iH<NHarm)&&((lHarm[iH]!=l)||(mHarm[iH]!=m)))
		iH++;
	if(iH==NHarm)
		return -1;
	else
		return(iH);
}

double LCGWSpinBBHNR1::gHarmDataLIN(int iH, int ihphc, double toM)
{
	double xr;
	int bin;
	//! ***** At the moment just linear interpolation
	xr = (toM-NRdatat0)/NRdatadt ;
	//xr = (toM-NRdatat0 - 0.32805e5 / 0.256e3 * M / eta )/NRdatadt ;
	bin = MT->ifloor(xr);
	if((bin<0)||(bin+1 > NRdataN-1)){
		std::cerr << "ERROR in LCGWSpinBBHNR1::gHarmData : bin = " << bin << " (x = " << toM << ") is not included in [0," << NRdataN-1 << "]." << Endl;  
		throw std::invalid_argument("LCGWSpinBBHNR1::gHarmData(): The required bin does not exist !");
	}
	return( (bin + 1.0 - xr) * (NRdata[iH][ihphc][bin]) + (xr-bin) * (NRdata[iH][ihphc][bin+1]) );
	
}



double LCGWSpinBBHNR1::InterLagrangeNR(int iH, int ihphc, double toM, int order)
{
    double xr;
	xr = (toM-NRdatat0)/NRdatadt ;
    
    return ( MT->InterLagrange(NRdata[iH][ihphc], NRdataN, xr, order) );
}




dcomplex LCGWSpinBBHNR1::ComputeSpherHarm( int l, int m, double theta, double phi)
{   
	
	int k1;
	int k2;
	int s = 2;    // a corriger
	dcomplex img(0.0, 1.0);
	double dnumer;
	double ddnom;
	double dlms;
	dcomplex Y;
	
	k1 =  (m-s < 0) ? 0 : m-s ; 
	k2 = (l+m < l-s )? l+m : l-s ;
	dlms = 0.0;
	
	for (int k =k1; k<k2+1 ; k++){
		
		dnumer = pow(-1, k) * sqrt(factorial(l+m) * factorial(l-m) * factorial(l+s) * factorial(l-s)) ;
		ddnom = factorial(l+m-k) * factorial(l-s-k) * factorial(k) * factorial(k+s-m) ;
		
		dlms += dnumer / ddnom * pow(cos(.5*theta), 2 * l + m - s - 2 * k) * pow(sin(.5*theta), 2 * k + s - m) ;
	}
	
	Y = pow(-1.,s) * sqrt(.25*(2.*l + 1.)/M_PI) * dlms * exp(img * ((double)(m)) * phi);
	
	return Y;
}


dcomplex LCGWSpinBBHNR1::SpherHarm(int l, int m)
{
	int iH(IndHarm(l,m));
	if(iH != -1)
		return(SpherHarmVal[iH]);
	else
		return 0.;
	
}


int LCGWSpinBBHNR1::factorial (int num)
{
	if (num==0) return 1;
	if (num==1) return 1;
	return factorial(num-1)*num; 
}



void LCGWSpinBBHNR1::ComputeExtraParam()
{
    
    //*************************************
    // Added by Sofiane
    //*************************************
    
    double stJ = sin(thetaJ);
    double ctJ = cos(thetaJ);
    
    //double down = cThS*stJ*cos(phiS-phiJ) - ctJ*sThS;
    //double up = stJ*sin(phiS - phiJ);
    //psi = atan2(-up, down);
    //  double up = cThS*stJ*cos(phiS-phiJ) - ctJ*sThS;  // commented by Sofiane
    //   double down = stJ*sin(phiS - phiJ);    // commented by Sofiane
    
  
    
    
    Thd = acos(-ctJ*cThS - stJ*sThS*cos(phiS - phiJ));
    
     ConvertLS1S2dirNR2SSBSofVer(Phd, thetaS, phiS, thetaJ, phiJ, LnN, S1N, S2N, AmpL, AmpS1, AmpS2, LnB, S1B, S2B);  // see how to use it in more elegant way
    
    
    double up =  cos(Phd - phiS);    // introduced by Sofiane
    double down = cos(thetaS)*sin(Phd - phiS);  // introduced by Sofiane
    
    //PolB = atan2(-up, down);   // there is perhapse a problem of sign with up, it depends in the convention used I added the sign but I am not sure to be checked
    double PolB = atan2(up, down);
    if(PolB<0.)
        PolB += M_PI;  // changed by Sofiane *2
    NRPolarization = PolB; 
    Polarization = NRPolarization;   // added by Sofiane
    c2psi = cos(2.*PolB);
	s2psi = sin(2.*PolB);
    
    
    //*************************************
    // End added by Sofiane
    //*************************************
    
    
    
    
	if(NeedExtraParamCompute){
		ComputeTotalAngMom();
		if(TypeExtraParamCompute==1){
			//! *** Compute Ln, S1 and S2 in SSB
		//	ConvertLS1S2dirNR2SSB(Thd, Phd, Polarization, JN, LnN, S1N, S2N, AmpL, AmpS1, AmpS2, LnB, S1B, S2B, PSIN);
            // Use Sofiane function
            ConvertLS1S2dirNR2SSBSofVer(Phd, thetaS, phiS, thetaJ, phiJ, LnN, S1N, S2N, AmpL, AmpS1, AmpS2, LnB, S1B, S2B);  // introduced by Sofiane
                                        
            double up, down;
            
            
            up =  cos(Phd - phiS);    // introduced by Sofiane
            down = cos(thetaS)*sin(Phd - phiS);  // introduced by Sofiane
            
            
            NRPolarization = atan2(up, down);
            if(NRPolarization<0.)
                NRPolarization += M_PI;
            Polarization = NRPolarization; 
            PSIN = NRPolarization;
            
            
            
            Thd = acos(-cos(thetaJ)*cos(thetaS) - sin(thetaJ)*sin(thetaS)*cos(phiS - phiJ));                          
                                    
            
            
			//! *** Compute angles in SSB
			thBL = LnB.th();
			phBL = LnB.ph();
			thBS1 = S1B.th();
			phBS1 = S1B.ph();
			thBS2 = S2B.th();
			phBS2 = S2B.ph();
		}
		if(TypeExtraParamCompute==2) { 
	//		ComputeThdPhdPol(LnB, S1B, S2B, LnN, S1N, S2N, AmpL, AmpS1, AmpS2, Thd, Phd, Polarization);    // canceled by Sofiane
        ComputeThdPhdPolSofVer(LnB, S1B, S2B, LnN, S1N, S2N, AmpL, AmpS1, AmpS2, Thd, Phd, NRPolarization, thetaJ, phiJ);   // added by Sofiane
             Polarization = NRPolarization;   // added by Sofiane
            PSIN = NRPolarization;
            c2psi = cos(2.*NRPolarization);
            s2psi = sin(2.*NRPolarization);
        }
	}
	NeedExtraParamCompute = false;
}


void LCGWSpinBBHNR1::ComputeTotalAngMom()
{
	//! *** Compute m1 and m2
	m1 =  M/(1.+1./mrat);
	m2 = M/(mrat+1.);
	
	//! *** Amplitude of L, S1 and S2
	AmpL  = m1*m2*pow(omMInit, -1./3);
	AmpS1 = chi1*m1*m1;
	AmpS2 = chi2*m2*m2;
	Cout << "AmpL = " << AmpL << Endl;
	Cout << "omMInit = " << omMInit << Endl;
	
	//! *** Components of J in NR frame
	JN = LnN*AmpL + S1N*AmpS1 + S2N*AmpS2; 
	Cout << "Sofiane Total Angular Momentum = " << JN  << Endl;  // added by Sofiane
	
}



void LCGWSpinBBHNR1::ConvertSpinDirToSSB(int iSpin)
{
	double thS, phS, th, ph;
	double iota, stheta, nez;
	double thBS, phBS;
	double phi0(0.);
	LCVector ex(MT), ey(MT), spin(MT), nW(MT), ez(MT);
	
	
	if(iSpin==1){
		thS = thS1;
		phS = phS1;
	}
	if(iSpin==2){
		thS = thS2;
		phS = phS2;
	}
	
	th = M_PI/2. - Beta;
	ph = Lambda;
	
	
	//! **** Compute \f$ \theta_L \f$ and \f$ \phi_L \f$ :
	//! ** Declaration of coordinates of ^k, ^theta, ^phi and ^Ln in SSB frame and cos / sin of some angles
	double x_k, y_k, z_k;
	double x_th, y_th, z_th;
	double x_ph, y_ph, z_ph;
	// double LnBx, LnBy, LnBz;
	double sPhS, cPhS, sPsi, cPsi, sThd, cThd; 
	
	sPhS = sin(phiS);
	cPhS = cos(phiS);
	sPsi = sin(2.*M_PI-NRPolarization);
	cPsi = cos(2.*M_PI-NRPolarization);
	sThd = sin(Thd);
	cThd = cos(Thd);
	
	x_k = - sThS * cPhS;
	y_k = - sThS * sPhS;
	z_k = - cThS;
	
	x_th = cThS * cPhS;
	y_th = cThS * sPhS;
	z_th = - sThS;
    
	x_ph = -y_k * z_th + z_k * y_th;
	y_ph = -z_k * x_th + x_k * z_th;
	z_ph = -x_k * y_th + y_k * x_th;
    
	
    
	
	
	//! ** Following convention of Sofiane for the angle psi (^p = cos psi ^theta + sin psi ^phi) :
	//	 LnBx = (- sThd * sPsi * x_th + sThd * cPsi * x_ph + cThd * x_k  );
	//	 LnBy = (- sThd * sPsi * y_th + sThd * cPsi * y_ph + cThd * y_k  );
	//	 LnBz = (- sThd * sPsi * z_th + sThd * cPsi * z_ph + cThd * z_k  );
	
	
	/*  
	 // new
	 
	 LnBx = (- sThd * sPsi * x_th - sThd * cPsi * x_ph + cThd * x_k  );
	 LnBy = (- sThd * sPsi * y_th - sThd * cPsi * y_ph + cThd * y_k  );
	 LnBz = (- sThd * sPsi * z_th - sThd * cPsi * z_ph + cThd * z_k  );
	 */
    
    /*
	 //! ** Following convention of Antoine for the angle psi (^p = cos psi ^theta - sin psi ^phi) :
	 LnBx = -( sThd * cPsi * x_th - sThd * sPsi * x_ph - cThd * x_k  );
	 LnBy = -( sThd * cPsi * y_th - sThd * sPsi * y_ph - cThd * y_k  );
	 LnBz = -( sThd * cPsi * z_th - sThd * sPsi * z_ph - cThd * z_k  );   
	 */
  	
	// with all angles
	
	
    LnB.x( -LnN.y()*cPsi*x_th+LnN.y()*sPsi*x_ph+LnN.x()*cThd*sPsi*x_th+LnN.x()*cThd*cPsi*x_ph-LnN.z()*sThd*sPsi*x_th-LnN.z()*sThd*cPsi*x_ph+x_k*LnN.x()*sThd+x_k*LnN.z()*cThd );
    LnB.y( -LnN.y()*cPsi*y_th+LnN.y()*sPsi*y_ph+LnN.x()*cThd*sPsi*y_th+LnN.x()*cThd*cPsi*y_ph-LnN.z()*sThd*sPsi*y_th-LnN.z()*sThd*cPsi*y_ph+y_k*LnN.x()*sThd+y_k*LnN.z()*cThd );
    LnB.z( -LnN.y()*cPsi*z_th+LnN.y()*sPsi*z_ph+LnN.x()*cThd*sPsi*z_th+LnN.x()*cThd*cPsi*z_ph-LnN.z()*sThd*sPsi*z_th-LnN.z()*sThd*cPsi*z_ph+z_k*LnN.x()*sThd+z_k*LnN.z()*cThd );
    
    
    double n_LnB = LnB.norm();
    
	//	thBL =   atan2(sqrt(LnBx*LnBx + LnBy*LnBy),LnBz);  // acos(LnBz/n_LnB); corrected by Sofiane
    thBL =  acos(LnB.z()/n_LnB);
    
	//    if(thBL<0.)
	//        thBL +-M_PI;
    
	//	phBL =M_PI*(1-LnBy/fabs(LnBy)) + LnBy/fabs(LnBy) * acos(LnBx / n_LnB / sqrt(1. - LnBz*LnBz / (n_LnB*n_LnB)) );  //  atan2(LnBy,LnBx); corrected by Sofiane
	
    phBL = LnB.ph();
    
    
    if(phBL<0.)
  		phBL += 2.*M_PI;
	
	//  double up = cThS*sin(thBL)*cos(phiS-phBL) - cos(thBL)*sThS;   /: commented by Sofiane
	//  double down = sin(thBL)*sin(phiS - phBL);    // commented by Sofiane
    
    ConvertLS1S2dirNR2SSBSofVer(Phd, thetaS, phiS, thetaJ, phiJ, LnN, S1N, S2N, AmpL, AmpS1, AmpS2, LnB, S1B, S2B);  // see how to use it in more elegant way
    
    
    double up =  cos(Phd - phiS);    // introduced by Sofiane
    double down = cos(thetaS)*sin(Phd - phiS);  // introduced by Sofiane
    
    PSIN =atan2(up,down);
    
    if(PSIN<0.)  
        PSIN +=2*M_PI;
	
	
	iota = acos(sin(th)*sin(thBL)*cos(ph - phBL) + cos(th)*cos(thBL)); 
	
    
	stheta = fabs(sin(iota));
	if (iota == 0.0)
		Cout << "WARNING in LCGWSpinBBHNR1::ConvertSpinDirToSSB : degenerate case need to consider separately: L colinear with n !" << Endl;
	nW.p[0] = sin(th)*cos(ph) ;
	nW.p[1] = sin(th)*sin(ph) ;
	nW.p[2] = cos(th) ;
	ez.p[0] = sin(thBL)*cos(phBL) ;  
	ez.p[1] = sin(thBL)*sin(phBL) ;
	ez.p[2] = cos(thBL) ;
	ey.p[0] = (nW.p[1]*ez.p[2] - nW.p[2]*ez.p[1])/stheta ;
	ey.p[1] = (nW.p[2]*ez.p[0] - nW.p[0]*ez.p[2])/stheta ;
	ey.p[2] = (nW.p[0]*ez.p[1] - nW.p[1]*ez.p[0])/stheta ;
	
	nez = nW.p[0]*ez.p[0] + nW.p[1]*ez.p[1] + nW.p[2]*ez.p[2] ;
	ex.p[0] = (ez.p[0]*nez - nW.p[0])/stheta ;
	ex.p[1] = (ez.p[1]*nez - nW.p[1])/stheta ;
	ex.p[2] = (ez.p[2]*nez - nW.p[2])/stheta ;
	
	//! ** if it is eccentric orbit, we need to rotate ex, ey by -phi0
	
	for(int i=0; i<3; i++)
		spin.p[i] = sin(thS)*cos(phS)*( cos(phi0)*ex.p[i] - sin(phi0)*ey.p[i] ) + sin(thS)*sin(phS)*( sin(phi0)*ex.p[i] + cos(phi0)*ey.p[i] ) + cos(thS)*ez.p[i] ;
	
	thBS = acos(spin.p[2]);
	phBS = atan2(spin.p[1], spin.p[0]);
	if (phBS < 0.)
		phBS = phBS + 2.*M_PI;
	
	if(iSpin==1){
		thBS1 = thBS;
		phBS1 = phBS;
	}
	if(iSpin==2){
		thBS2 = thBS;
		phBS2 = phBS;
	}	
	
	S1B.x( sin(thBS1)*cos(phBS1) );
	S1B.y( sin(thBS1)*sin(phBS1) );
	S1B.z( cos(thBS1) );
	
	S2B.x( sin(thBS2)*cos(phBS2) );
	S2B.y( sin(thBS2)*sin(phBS2) );
	S2B.z( cos(thBS2) );
	
	
}






  /*   
 double LCGWSpinBBHNR1::ComputeApprOrbFreq(double t, 
 double Tc,
 double x1,
 double x2,
 double LS1, 
 double LS2,
 double S1S2,
 double m1M2,
 double m2M2,
 double eta)
 {
 double beta, sigma, tau, tau38, tau14, x;
 
 //! *** WARNING : THIS FUNCTION HAVE TO BE EXACTLY THE SAME AS THE ONE DEFINE IN LCGWSpinBBHNR1 TO KEEP COHERENT WAVEFORM 
 
 /*! Here I use 2PN approximate expression  for the frequency as function of t 
 * neglecting the precession
 /
 
 /*LS1 = Lnx*S1x + Lny*S1y + Lnz*S1z;
 LS2 = Lnx*S2x + Lny*S2y + Lnz*S2z;
 S1S2 = S1x*S2x + S1y*S2y + S1z*S2z;  computed in the constructor/
 
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
  */   



// end of LISACODE-GWSpinBBHNR1.cpp