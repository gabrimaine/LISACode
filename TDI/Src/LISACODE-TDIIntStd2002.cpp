/*
 *  LISACODE-LCTDIIntStd2002.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 20/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


#include "LISACODE-TDIIntStd2002.h"



// *********************
// ***  Constructor  ***
// *********************


LCTDIIntStd2002::LCTDIIntStd2002(LCTools * MT_n)
: LCTDIInt(MT_n)
{
	initNULL(false);
}


LCTDIIntStd2002::~LCTDIIntStd2002()
{
	initNULL(true);
}



// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCTDIIntStd2002::initNULL(bool CleanMem)
{
	initNULLBase(CleanMem);
}


void LCTDIIntStd2002::config(ezxml_t xmlbloc)
{
	ezxml_t param;
	char * TDIIntType(NULL);
	
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		
		//! **** Checking that it's the right configuration of optical bench + phasemeter
		if(MT->wcmp(ezxml_attr(param,"Name"),"TDIIntermediateType")){
			MT->stripcopy((*param).txt, TDIIntType);
			if(!MT->wcmp(TDIIntType, "Std2002"))
				throw std::invalid_argument("ERROR in LCTDIIntStd2002::config : The xml bloc does not correspond to a type of the optical bench + phasemeter 'Std2002' .");
		}
	
		//! **** Checking that it's the right configuration of optical bench + phasemeter
		if(MT->wcmp(ezxml_attr(param,"Name"),"Interpolation")){
			if(MT->wcmp(ezxml_attr(param,"Type"),"Lagrange"))
				InterpType = LAG;
			if(MT->wcmp(ezxml_attr(param,"Type"),"Linear"))
				InterpType = LIN;
			if(MT->wcmp(ezxml_attr(param,"Type"),"Trunc"))
				InterpType = TRU;
			InterpUtilValue = atoi((*param).txt);
		}
	}
	
	//! *** Set the name and the position
	setNameEqPos(ezxml_attr(xmlbloc,"Name"));
	
	if(TDIIntType != NULL)
		MT->Free(TDIIntType, (strlen(TDIIntType)+1) * sizeof(char) );
}


void LCTDIIntStd2002::LinkDetector(LCDetector * LISA)
{
	//! *** Delete previous allocation
	if(Sig!=NULL)
		MT->Free(Sig, NSig*sizeof(LCSerie2*));
	Sig = NULL;
	
	if(!IndirectDir){
		/*! *** If it's in direct direction
		 *	Sig[0]	: \f$ s_i \f$  : scientific meausrement  [local]
		 *	Sig[1]	: \f$ \tau_{i+1} \f$  : extern back measurement  [tau..]
		 *	Sig[2]	: \f$ \tau'_{i+1} \f$  : extern join back measurement [tau..]
		 */
		NSig = 3;
		Sig = (LCSerie2**) MT->AllocMemory(NSig*sizeof(LCSerie2*));
		for(int i=0; i<NSig; i++)
			Sig[i] = NULL;
		int iSCp1(1+(iSC)%3);
		Sig[0] = LISA->LinkPhaMes(iSC, 0, "sci", 3);
		Sig[1] = LISA->LinkPhaMes(iSCp1, 0, "tau", 3);
		Sig[2] = LISA->LinkPhaMes(iSCp1, 1, "tau", 3);
		
	}else{
		/*! *** If it's in indirect direction
		 *	Sig[0]	: \f$ s'_i \f$  : scientific meausrement  [local]
		 *	Sig[1]	: \f$ \tau_i \f$  : extern back measurement  [tau..]
		 *	Sig[2]	: \f$ \tau'_i \f$  : extern join back measurement [tau..]
		 */
		NSig = 3;
		Sig = (LCSerie2**) MT->AllocMemory(NSig*sizeof(LCSerie2*));
		for(int i=0; i<NSig; i++)
			Sig[i] = NULL;
		Sig[0] = LISA->LinkPhaMes(iSC, 1, "sci", 3);
		Sig[1] = LISA->LinkPhaMes(iSC, 0, "tau", 3);
		Sig[2] = LISA->LinkPhaMes(iSC, 1, "tau", 3);
	}
	
}


void LCTDIIntStd2002::LinkDelays(LCSerie2** AllDelays, int NAllDelays)
{
	//! *** Delete previous allocation
	if(Delays!=NULL)
		MT->Free(Delays, NDelays*sizeof(LCSerie2*));
	Delays = NULL;
	
	if(!IndirectDir){
		/*! *** If it's in direct direction
		 *	Delays[0] : 
		 */
		NDelays = 1;
		Delays = (LCSerie2**) MT->AllocMemory(NDelays*sizeof(LCSerie2*));
		Delays[0] = NULL;
		if(NAllDelays>(iSC+1)%3)
			Delays[0] = AllDelays[(iSC+1)%3];
		
	}
	//! *** If it's in indirect direction : no delays
		 
}


void LCTDIIntStd2002::init()
{
	initBase();
}


int LCTDIIntStd2002::getNMaxDelay()
{
	return(1);
}


void LCTDIIntStd2002::RunStep(double t)
{
	double TmpV(0.);
	//! **** If it's indirect direction 
	if(!IndirectDir){
		//! *** Compute \f$ \eta_i = s_i - { D_{i+2} \tau_{i+1} - D_{i+2} \tau'_{i+1}  \over 2 } \f$
		double tDel(tShift - Delays[0]->gData(tShift, InterpTypeDelay, InterpUtilValueDelay));
		TmpV = gN(0,tShift) - 0.5*( gN(1,tDel) - gN(2,tDel) );
	}else{
		//! *** Compute \f$ \eta'_i = s'_i + { \tau_{i} - \tau'_{i}  \over 2 } \f$
		TmpV = gN(0,tShift) + 0.5*( gN(1,tShift) - gN(2,tShift) ) ;
	}
	
	Dat->addData(TmpV);
}


void LCTDIIntStd2002::DispInfo(char * BTab)
{
	if(MT->Disp()){
		DispInfoBase(BTab);
	}
	
}



// ***********************
// ***  Local mehtods  ***
// ***********************



// *******************
// *  Other methods  *
// *******************



// end of LISACODE-InterTDI.cpp
