/*
 *  LISACODE-Detector.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 12/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


#include "LISACODE-Detector.h"


// *********************
// ***  Constructor  ***
// *********************

LCDetector::LCDetector()
{
	MT = new LCTools;
	MT->LocTools = true;
	initNULLBase(false);
}


LCDetector::LCDetector(LCTools * MT_n)
{
	MT = MT_n;
	initNULLBase(false);
}



LCDetector::~LCDetector()
{
	initNULLBase(true);
}


void LCDetector::initNULLBase(bool CleanMem)
{
	if(CleanMem){
		
		if(ArmGW != NULL)
			delete ArmGW;
		
		if(Noises!=NULL){
			for(int i=0; i<NNoises; i++)
				if(Noises[i] != NULL)
					delete Noises[i];
			MT->Free(Noises, NNoises*sizeof(LCNoise*));
		}
		
		if(Measures!=NULL){
			for(int i=0; i<NMeasures; i++)
				if(Measures[i] != NULL)
					delete Measures[i];
			MT->Free(Measures, NMeasures*sizeof(LCOBPM*));
		}
		
		if(Clocks!=NULL){
			for(int i=0; i<NClocks; i++)
				if(Clocks[i] != NULL)
					delete Clocks[i];
			MT->Free(Clocks, NClocks*sizeof(LCUSO*));
		}
		
		if(Lasers!=NULL){
			for(int i=0; i<NLasers; i++)
				if(Lasers[i]!=NULL)
					MT->Free(Lasers[i], sizeof(LCLaser));
			MT->Free(Lasers, NLasers*sizeof(LCLaser*));
		}
		
		if(Telesc!=NULL){
			for(int i=0; i<NTelesc; i++)
				if(Telesc[i] != NULL)
					MT->Free(Telesc[i], sizeof(LCTelescope));
			MT->Free(Telesc, NTelesc*sizeof(LCTelescope*));
		}
		
		
		
	}
	ArmGW = NULL;
	Noises = NULL;
	NNoises = 0;
	Measures = NULL;
	NMeasures = 0;
	Clocks = NULL;
	NClocks = 0;
	Lasers = NULL;
	NLasers = 0;
	Telesc = NULL;
	NTelesc = 0;
	NoiseInterp = UND;
	NoiseInterpVal = 0;
	
	
	Orb = NULL;
	
	dtPhy = 1.0;
	dtMes = 1.0;
	tStoreDelay = 30.;
	NDataPha = 30;
	
	UpdateShotNoise = true;
}



// ***************************
// *  Configuration methods  *
// ***************************


void LCDetector::configDetector(ezxml_t xmlbloc)
{
	ezxml_t param, section;
	//! **** Configuration of detector parameters
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"StepPhysic"))
			dtPhy = MT->gXMLTime(param);
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"UpdateShotNoise")){
			if(MT->wcmp((*param).txt,"On"))
				UpdateShotNoise = true;
			else
				UpdateShotNoise = false;
		}
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"InterpolationNoises")){
			if(MT->wcmp(ezxml_attr(param,"Type"),"Lagrange"))
				NoiseInterp = LAG;
			if(MT->wcmp(ezxml_attr(param,"Type"),"Linear"))
				NoiseInterp = LIN;
			if(MT->wcmp(ezxml_attr(param,"Type"),"Trunc"))
				NoiseInterp = TRU;
			NoiseInterpVal = atoi((*param).txt);
			
		}
	}
	
	for (section = ezxml_child(xmlbloc, "XSIL"); section; section = section->next) {
		//! **** Configuration of optical bench + phasemeter
		if(MT->wcmp(ezxml_attr(section,"Type"),"OpticalBenchPhasemeter"))
			configMeasure(section);
		
		//! **** Configuration of USO
		if(MT->wcmp(ezxml_attr(section,"Type"),"USO"))
			configClock(section);
		
		//! **** Configuration of Laser
		if(MT->wcmp(ezxml_attr(section,"Type"),"Laser"))
			configLaser(section);
		
		//! **** Configuration of Telescope
		if(MT->wcmp(ezxml_attr(section,"Type"),"Telescope"))
			configTelescope(section);
		
	}
	
}


void LCDetector::configNoise(ezxml_t section)
{
	ezxml_t noisedata;
	for (noisedata = ezxml_child(section, "XSIL"); noisedata; noisedata = noisedata->next) {
		ezxml_t param;
		char * SourceType(NULL);
        char * SpectralType(NULL);
        
        
        for(param = ezxml_child(noisedata,"Param"); param; param = param->next){
            if(MT->wcmp(ezxml_attr(param,"Name"),"SourceType"))
                MT->stripcopy((*param).txt, SourceType);
            if(MT->wcmp(ezxml_attr(param,"Name"),"SpectralType"))
                MT->stripcopy((*param).txt, SpectralType);
        }
        
        if((MT->wcmp(SourceType, "PseudoRandomNoise"))&&(MT->wcmp(SpectralType, "White"))){
            NNoises++;
            Noises = (LCNoise**) MT->ReAllocMemory(Noises, (NNoises-1)*sizeof(LCNoise*), NNoises*sizeof(LCNoise*));
            Noises[NNoises-1] = new LCNoiseWhite(MT);
            Noises[NNoises-1]->config(noisedata);
        }
        
        if((MT->wcmp(SourceType, "PseudoRandomNoise"))&&((MT->wcmp(SpectralType, "Filter_f"))
                                                         ||(MT->wcmp(SpectralType, "Filter_1of"))
                                                         ||(MT->wcmp(SpectralType, "WhitePhase"))
                                                         ||(MT->wcmp(SpectralType, "BlueFrequency"))
                                                         ||(MT->wcmp(SpectralType, "WhiteFrequency"))
                                                         ||(MT->wcmp(SpectralType, "RedFrequency"))
                                                         ||(MT->wcmp(SpectralType, "PinkFrequency"))
                                                         ||(MT->wcmp(SpectralType, "PinkAcceleration"))
                                                         ||(MT->wcmp(SpectralType, "Filter_1of_1of2"))
                                                         ||(MT->wcmp(SpectralType, "Filter_1of_1of32"))
                                                         ||(MT->wcmp(SpectralType, "PreStabLaserNoiseFreq")))){
            NNoises++;
            Noises = (LCNoise**) MT->ReAllocMemory(Noises, (NNoises-1)*sizeof(LCNoise*), NNoises*sizeof(LCNoise*));
            Noises[NNoises-1] = new LCNoiseFilter(MT);
            Noises[NNoises-1]->config(noisedata);
        }
        
        if(MT->wcmp(SourceType, "DataFile")){
            NNoises++;
            Noises = (LCNoise**) MT->ReAllocMemory(Noises, (NNoises-1)*sizeof(LCNoise*), NNoises*sizeof(LCNoise*));
            Noises[NNoises-1] = new LCNoiseFile(MT);
            Noises[NNoises-1]->config(noisedata);
        }
        
        //! ** HERE, Add new type of noise
		
        if(SpectralType != NULL)
            MT->Free(SpectralType, (strlen(SpectralType)+1) * sizeof(char));
        if(SourceType != NULL)
            MT->Free(SourceType, (strlen(SourceType)+1) * sizeof(char));
    }
	
}



void LCDetector::configClock(ezxml_t xmlbloc)
{
	ezxml_t param;
	char * ClockType(NULL);
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		if(MT->wcmp(ezxml_attr(param,"Name"),"ClockType")){
			MT->stripcopy((*param).txt, ClockType);
			if(MT->wcmp(ClockType, "NoiseDeriv")){
				NClocks++;
				Clocks = (LCUSO**) MT->ReAllocMemory(Clocks, (NClocks-1)*sizeof(LCUSO*), NClocks*sizeof(LCUSO*));
				Clocks[NClocks-1] = new LCUSONoiseDeriv(MT);
				Clocks[NClocks-1]->config(xmlbloc);
			}
		}
	}
	if(ClockType != NULL)
		MT->Free(ClockType, (strlen(ClockType)+1) * sizeof(char));
	
}


void LCDetector::configMeasure(ezxml_t xmlbloc)
{
	ezxml_t param;
	char * OpticalBenchType(NULL);
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		if(MT->wcmp(ezxml_attr(param,"Name"),"OpticalBenchType")){
			MT->stripcopy((*param).txt, OpticalBenchType);
			if(MT->wcmp(OpticalBenchType, "Std2002")){
				NMeasures++;
				Measures = (LCOBPM**) MT->ReAllocMemory(Measures, (NMeasures-1)*sizeof(LCOBPM*), NMeasures*sizeof(LCOBPM*));
				Measures[NMeasures-1] = new LCOBPMStd2002(MT);
				Measures[NMeasures-1]->config(xmlbloc);
			}
			if(MT->wcmp(OpticalBenchType, "Std2010")){
				NMeasures++;
				Measures = (LCOBPM**) MT->ReAllocMemory(Measures, (NMeasures-1)*sizeof(LCOBPM*), NMeasures*sizeof(LCOBPM*));
				Measures[NMeasures-1] = new LCOBPMStd2010(MT);
				Measures[NMeasures-1]->config(xmlbloc);
			}
			if(MT->wcmp(OpticalBenchType, "Octa1")){
				NMeasures++;
				Measures = (LCOBPM**) MT->ReAllocMemory(Measures, (NMeasures-1)*sizeof(LCOBPM*), NMeasures*sizeof(LCOBPM*));
				Measures[NMeasures-1] = new LCOBPMOcta1(MT);
				Measures[NMeasures-1]->config(xmlbloc);
			}
			//! ** HERE, Add new type of optical bench
		}
	}
	if(OpticalBenchType != NULL)
		MT->Free(OpticalBenchType, (strlen(OpticalBenchType)+1) * sizeof(char));
}


void LCDetector::configLaser(ezxml_t xmlbloc)
{
	ezxml_t param;
	
	//! **** Create the laser
	NLasers++;
	Lasers = (LCLaser**) MT->ReAllocMemory(Lasers, (NLasers-1)*sizeof(LCLaser*), NLasers*sizeof(LCLaser*));
	Lasers[NLasers-1] = (LCLaser*) MT->AllocMemory(sizeof(LCLaser));
	Lasers[NLasers-1]->iSC = -1;
	Lasers[NLasers-1]->iSC = -1;
	Lasers[NLasers-1]->power = 1.;
	Lasers[NLasers-1]->freq = 1064.e-9;
	
	//! **** Set the position
	const char * TmpName;
	TmpName = ezxml_attr(xmlbloc, "Name");
	int ipos(strlen(TmpName)-1);
	
	//! ** Start from the end : detect if the direction of optical bench
	if(TmpName[ipos] == 's'){
		Lasers[NLasers-1]->IndirectDir = 1;
		ipos--;
	}else{
		Lasers[NLasers-1]->IndirectDir = 0;
	}
	
	//! *** Then : Spacecraft index
	char TmpiSC[2];
	TmpiSC[0] = TmpName[ipos];
	TmpiSC[1] = '\0';
	Lasers[NLasers-1]->iSC = atoi(TmpiSC);
	
	//if((Lasers[NLasers-1]->iSC<1)||(Lasers[NLasers-1]->iSC>3))
	//	throw std::invalid_argument("ERROR in LCDetector::configLaser : Index of spacecraft should be in 1,2 or 3 !" );
	
	
	
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		if(MT->wcmp(ezxml_attr(param,"Name"),"Power"))
			Lasers[NLasers-1]->power = atof((*param).txt);
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"Wavelength"))
			Lasers[NLasers-1]->freq = LC::c_SI/MT->gXMLLength(param);
		
	}
	
}



void LCDetector::configTelescope(ezxml_t xmlbloc)
{
	ezxml_t param;
	char TmpiSC[2];
	
	//! **** Create the laser
	NTelesc++;
	Telesc = (LCTelescope**) MT->ReAllocMemory(Telesc, (NTelesc-1)*sizeof(LCTelescope*), NTelesc*sizeof(LCTelescope*));
	Telesc[NTelesc-1] = (LCTelescope*) MT->AllocMemory(sizeof(LCTelescope));
	Telesc[NTelesc-1]->iSC = -1;
	Telesc[NTelesc-1]->iSC = -1;
	Telesc[NTelesc-1]->diameter = .4;
	
	//! **** Set the position
	const char * TmpName;
	TmpName = ezxml_attr(xmlbloc, "Name");
	int ipos(strlen(TmpName)-1);
	
	//! ** Start from the end : detect if the direction of optical bench
	if(TmpName[ipos] == 's'){
		Telesc[NTelesc-1]->IndirectDir = 1;
		ipos--;
	}else{
		if(isdigit(TmpName[ipos])&&(!(isdigit(TmpName[ipos-1])))){	
			Telesc[NTelesc-1]->IndirectDir = 0;
		}else{
			//! ** There are 2 digits at the end of the name , so it's not a LISA like detector (not 3 SCs)
			TmpiSC[0] = TmpName[ipos];
			TmpiSC[1] = '\0';
			Telesc[NTelesc-1]->iSCLink = atoi(TmpiSC);
			Telesc[NTelesc-1]->IndirectDir = -1*atoi(TmpiSC); 
			ipos--;
		}
	}
	
	//! *** Then : Spacecraft index
	TmpiSC[0] = TmpName[ipos];
	TmpiSC[1] = '\0';
	Telesc[NTelesc-1]->iSC = atoi(TmpiSC);
	
	
	if(Telesc[NTelesc-1]->IndirectDir >= 0){
		//! ** Computation of index of emitter
		if(Telesc[NTelesc-1]->IndirectDir)
			Telesc[NTelesc-1]->iSCLink = (Telesc[NTelesc-1]->iSC+1)%3 + 1;
		else
			Telesc[NTelesc-1]->iSCLink = (Telesc[NTelesc-1]->iSC)%3 + 1;
	}
	
	//Cout << Telesc[NTelesc-1]->iSC << " , " << Telesc[NTelesc-1]->IndirectDir << " : link " << Telesc[NTelesc-1]->iSCLink << Endl;
	//if((Telesc[NTelesc-1]->iSC<1)||(Telesc[NTelesc-1]->iSC>3))
	//	throw std::invalid_argument("ERROR in LCDetector::configTelescope : Index of spacecraft should be in 1,2 or 3 !" );
	
	
	
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		if(MT->wcmp(ezxml_attr(param,"Name"),"Diameter")){
			Telesc[NTelesc-1]->diameter = MT->gXMLLength(param);
		}
	}
}


// ***************************************
// * Linking and initialization methods  *
// ***************************************

void LCDetector::LinkGWs( LCGW ** LCGWs, int NGWs )
{
	if(Orb == NULL)
		throw std::invalid_argument("ERROR in LCDetector::LinkGWs : The orbits should be link to the detector before the gravitational waves ! ");
	if(ArmGW != NULL)
		delete ArmGW;
	ArmGW = new LCArmResp(MT);
	ArmGW->LinkOrb(Orb);
	ArmGW->LinkGWs(LCGWs, NGWs);
}


void LCDetector::init()
{
	//! ***** Pre computation of shot noise factor in telescope
	PreComputeFactShotNoise();
	
	//! ***** Arm response
	ArmGW->init();
	
	//! ***** Noises :
	for(int i=0; i<NNoises; i++){
		//! *** Set interpolation if needed
		if(NoiseInterp != UND)
			Noises[i]->setInterpolation(NoiseInterp, NoiseInterpVal);
		//! *** Set time informations
		Noises[i]->setTimeInfo(dtMes, dtPhy, 0., tStoreDelay-(MIN(t0,0.)));
		//! *** Initialization
		Noises[i]->init();
			
	}
	
	
	//! ***** Clocks :
	for(int i=0; i<NClocks; i++){
		//! *** Set time informations
		Clocks[i]->setTimeInfo(dtMes, dtPhy, 0., tStoreDelay);
		//! *** Initialization
		Clocks[i]->init();
	}
	
	
	//! ***** Measurement systems : optical benches + photodiodes + phsemeters
	for (int i=0; i<NMeasures; i++) {
		
		//! *** Set time informations
		Measures[i]->setTimeInfo(t0, dtMes, dtPhy, 0., tStoreDelay, NDataPha);
		
		/*! *** Link the orbits */
		Measures[i]->LinkOrbits(Orb);
		
		/*! *** Link the arm response */
		Measures[i]->LinkArmResp(ArmGW);
		
		//! *** Link noises
		Measures[i]->LinkNoise(Noises, NNoises);
		
		//! *** Link clocks
		Measures[i]->LinkClock(Clocks, NClocks);
		
		//! *** Link factor on shot noise
		for(int iTel=0; iTel<NTelesc; iTel++)
			if(Telesc[iTel] != NULL)
				if( (Telesc[iTel]->iSC == Measures[i]->getiSC()) && (Telesc[iTel]->iSCLink == Measures[i]->getiSCLink()) )
					Measures[i]->LinkFactShotNoise(&(Telesc[iTel]->FactShotNoise));
		
		//! *** Initialization
		Measures[i]->init();
		
	}
	
}


bool LCDetector::LinkMes2Output(double * pOut, char * ObsName)
{
	bool NoFoundOut(true);
	for(int i=0; i<NMeasures; i++){
		if(Measures[i]->MatchName(ObsName)){
			NoFoundOut = false;
			Measures[i]->LinkLastMesOutput(pOut);
		}
	}
	return(NoFoundOut);
}


// ********************
// *  Access methods  *
// ********************

void LCDetector::setTimeInfo( double t0_n, double dtMes_n, double tStoreDelay_n, int NDataPha_n)
{
	t0 = t0_n;
	dtMes = dtMes_n;
	tStoreDelay = tStoreDelay_n;
	NDataPha = NDataPha_n;
}


LCSerie2 * LCDetector::LinkPhaMes(int iSC_cmp, int IndirectDir_cmp, const char * Name_cmp, int nchar)
{
	LCSerie2 * res(NULL);
	for(int iM=0; iM<NMeasures; iM++){
		if(Measures[iM]->CmpPosName(iSC_cmp, IndirectDir_cmp, Name_cmp, nchar))
			res = Measures[iM]->getPhaMes();
		
	}
	return(res);
}


// *********************
// *  Running methods  *
// *********************

void LCDetector::RunStep(double t)
{
	//! ***** Run one step of the noises (add noise)
	for(int i=0; i<NNoises; i++){
		Noises[i]->RunStep();
	}
	
	
	//! ***** Run one step of the clock (compute new time)
	for(int i=0; i<NClocks; i++){
		Clocks[i]->RunStep(t);
	}
	
	
	//! ***** Computation of Shot noise factors
	if(UpdateShotNoise){
		for(int iR=0; iR<NTelesc; iR++){
			ComputeFactShotNoise(iR, t);
		}
	}
	
	
	//! *****  Run one step of measurement systems : optical benches + photodiodes + phasemeters
	for (int i=0; i<NMeasures; i++) {
		Measures[i]->MeasurePha(t);
		
	}
	
}


void LCDetector::ComputeFactShotNoise(int iR, double tGlobal)
{	
	if(Telesc[iR] != NULL)
		Telesc[iR]->FactShotNoise = Telesc[iR]->PreFactShotNoise * Orb->Arm(Telesc[iR]->iSCLink, Telesc[iR]->iSC, tGlobal);
}




// ********************
// *  Others methods  *
// ********************

void LCDetector::PreComputeFactShotNoise()
{
	/*! ***** For each telescope compute the factor on shot noise which is :
	 *	\f$ Fact = PreFact \times t_{travel} \f$ 
	 *	with \f$ PreFact = 8740942.01607775 {\lambda_laser^{3/2} \over \sqrt{P_{laser}} D^2_telescope } =  8740942.01607775 { (c/f_laser)^{3/2} \over \sqrt{P_{laser}} D^2_telescope } \f$
	 */
	
	//! ** Computation of pre-computed factor which does not change during the simulation
	for(int iR=0; iR<NTelesc; iR++){
		if(Telesc[iR] != NULL){
			
			//Cout << iR << " : " << Telesc[iR]->iSC << " <- " << Telesc[iR]->iSCLink ;  
			
			//! ** Look for the laser emitting the beam 
			int iLaser(-1);
			for(int iL=0; iL<NLasers; iL++)
				if( (Telesc[iR]->iSCLink == Lasers[iL]->iSC) && (Telesc[iR]->IndirectDir != Lasers[iL]->IndirectDir) )
					iLaser = iL;
			
			/*! ** Computation of pre-computed factor which does not change during the simulation 
			 *	\f[ \delta_{SN} \propto { L \over P D^2 } \lambda^{3/2} \f]
			 *	\f[ \delta_{SN} = \delta_{SN,0} \left( L \over 5 \times 10^9 \textrm{m}\right) \left( 1 \textrm{m} \over P \right) \left( 0.4 \textrm{m} \over D \right)^2  \left( \lambda \over 1064 \textrm{nm} \right)^{3/2} \f]
			 */
			//Telesc[iR]->PreFactShotNoise = ( (0.4*0.4) / ( pow(1064e-9,3./2.) * 5.e9 / LC::c_SI) ) * ( pow(LC::c_SI/Lasers[iLaser]->freq,3./2.) / ( Lasers[iLaser]->power * Telesc[iR]->diameter * Telesc[iR]->diameter ) );
			
			/* ** Computation of pre-computed factor which does not change during the simulation 
			 *	\f[ \delta_{SN} \propto { L \over P D^2 } \lambda^{3/2} \f]
			 *	\f[ \delta_{SN} = \delta_{SN,0} \left( L \over 5 \times 10^9 \textrm{m}\right) {\left( 1 \textrm{m} \over P \right)}^{1/2} \left( 0.4 \textrm{m} \over D \right)^2  \left( \lambda \over 1064 \textrm{nm} \right)^{3/2} \f]
			 */
			
			//Cout << " ==> iL = " << iLaser << Endl; 
			
			Telesc[iR]->PreFactShotNoise = ( (0.4*0.4) / ( pow(1064e-9,3./2.) * 5.e9 / LC::c_SI) ) * ( pow(LC::c_SI/Lasers[iLaser]->freq,3./2.) / ( sqrt(Lasers[iLaser]->power) * Telesc[iR]->diameter * Telesc[iR]->diameter ) );
			
			
			//! ** First computation of Factor 
			Telesc[iR]->FactShotNoise = Telesc[iR]->PreFactShotNoise * Orb->getNominalArm();
			//Cout << Telesc[iR]->PreFactShotNoise << " " <<  Orb->getNominalArm() <<  " " << Telesc[iR]->FactShotNoise << Endl;
		}
	}	
}


void LCDetector::DispInfo(char * BTab)
{
	if(MT->Disp()){
		char BTab2[1024];
		strcpy(BTab2,BTab);
		strcat(BTab2, "\t");
		
		Cout << BTab << "Detector : " << Endl;
		
		Cout << BTab2 << "Time step physical = " << dtPhy << Endl;
		
		//! ** Display the noises
		if(Noises != NULL)
			for(int i=0; i<NNoises; i++)
				if(Noises[i] != NULL)
					Noises[i]->DispInfo(BTab2);
		
		//! ** Display the clocks
		if(Clocks != NULL)
			for(int i=0; i<NClocks; i++)
				if(Clocks[i] != NULL)
					Clocks[i]->DispInfo(BTab2);
		
		//! ** Display the measures
		if(Measures != NULL)
			for(int i=0; i<NMeasures; i++)
				if(Measures[i] != NULL)
					Measures[i]->DispInfo(BTab2,true);
		
		//! ** Display the beams
		if(Lasers != NULL){
			for(int i=0; i<NLasers; i++){
				if(Lasers[i] != NULL){
					Cout << BTab2 << "Laser : spacecraft = " << Lasers[i]->iSC << "  , direction = " << Lasers[i]->IndirectDir << Endl;
					if(MT->DispDet()){
						Cout << BTab2 << "\t- Power = "  << Lasers[i]->power << " W" << Endl;
						Cout << BTab2 << "\t- Frequency = "  << Lasers[i]->freq << " Hz" << Endl;
					}
				}
			}
		}
		
		
		
		//! ** Display the telescope
		if(Telesc != NULL){
			for(int i=0; i<NTelesc; i++){
				if(Telesc[i] != NULL){
					Cout << BTab2 << "Telescope : spacecraft = " << Telesc[i]->iSC << "  , direction = " << Telesc[i]->IndirectDir << Endl;
					if(MT->DispDet())
						Cout << BTab2 << "\t- Diameter = "  << Telesc[i]->diameter << " m" << Endl;
				}
			}
		}
	}
}

// end of LISACODE-Detector.cpp
																						  




