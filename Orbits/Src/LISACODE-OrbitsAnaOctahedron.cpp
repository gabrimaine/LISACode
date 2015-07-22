/*
 *  LISACODE-OrbitsOctaAnalytic.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 10/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-OrbitsAnaOctahedron.h"

// *********************
// ***  Constructor  ***
// *********************

/**
 *
 */
LCOrbitsAnaOctahedron::LCOrbitsAnaOctahedron()
: LCOrbits()
{
	MT = new LCTools;
	MT->LocTools = true;
	initNULL(false);
}


LCOrbitsAnaOctahedron::LCOrbitsAnaOctahedron(LCTools * MT_n)
: LCOrbits(MT_n)
{
	MT = MT_n;
	initNULL(false);
}


LCOrbitsAnaOctahedron::~LCOrbitsAnaOctahedron()
{
	initNULL(true);
}




// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCOrbitsAnaOctahedron::initNULL(bool CleanMem)
{	
	NSC = 6;
	Rgc = LC::au_m - 1.5e9;
	omega = 2.*M_PI/LC::Yr_SI;
	Lo2 = L0m / 2.;
	Lsq2o2 =  Lo2 * sqrt(2.) ;
	for(int i=0; i<6; i++){
		pSCOff[i] = 0.;
		for(int j=0; j<3; j++)
			posRL1[i].p[j] = 0.;
	}
	
}


void LCOrbitsAnaOctahedron::config(ezxml_t xmlbloc)
{
	configBase(xmlbloc);
	ezxml_t param;
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		if(strcmp(ezxml_attr(param,"Name"),"Armlength")==0){
			L0m = MT->gXMLLength(param);
		}
		if(strcmp(ezxml_attr(param,"Name"),"OrderDelay")==0){
			OrderArm = atoi(ezxml_txt(param));
		}
		if(strcmp(ezxml_attr(param,"Name"),"OrbitApproximation")==0){
			char * tmpStr(NULL);
			MT->stripcopy((*param).txt, tmpStr);
			
			if(MT->wcmp(tmpStr,"Static"))
				move = 0;
			if(MT->wcmp(tmpStr,"Moving"))
				move = 1;
			if(tmpStr != NULL)
				MT->Free(tmpStr, (strlen(tmpStr)+1) * sizeof(char));
		}
		if(strcmp(ezxml_attr(param,"Name"),"OffsetSC1")==0){
			pSCOff[0] = atof(ezxml_txt(param));
		}
		if(strcmp(ezxml_attr(param,"Name"),"OffsetSC2")==0){
			pSCOff[1] = atof(ezxml_txt(param));
		}
		if(strcmp(ezxml_attr(param,"Name"),"OffsetSC3")==0){
			pSCOff[2] = atof(ezxml_txt(param));
		}
		if(strcmp(ezxml_attr(param,"Name"),"OffsetSC4")==0){
			pSCOff[3] = atof(ezxml_txt(param));
		}
		if(strcmp(ezxml_attr(param,"Name"),"OffsetSC5")==0){
			pSCOff[4] = atof(ezxml_txt(param));
		}
		if(strcmp(ezxml_attr(param,"Name"),"OffsetSC6")==0){
			pSCOff[5] = atof(ezxml_txt(param));
		}
	}
}


void LCOrbitsAnaOctahedron::config(int iParam, double ParamVal)
{
	
}


void LCOrbitsAnaOctahedron::init()
{
	initBase();
	
	Lo2 = L0m / 2.;
	Lsq2o2 =  Lo2 * sqrt(2.) ;
	/*
	RxySC[0] = Rgc;
	RxySC[1] = sqrt((Rgc+L0m/2.)*(Rgc+L0m/2.) + L0m*L0m/4.);
	RxySC[2] = sqrt((Rgc+L0m/2.)*(Rgc+L0m/2.) + L0m*L0m/4.);
	RxySC[3] = Rgc;
	RxySC[4] = sqrt((Rgc-L0m/2.)*(Rgc-L0m/2.) + L0m*L0m/4.);
	RxySC[5] = sqrt((Rgc-L0m/2.)*(Rgc-L0m/2.) + L0m*L0m/4.);
	
	phixySC[0] = 0.; 
	phixySC[1] = atan2(-L0m/2, RxySC[1]);
	phixySC[2] = atan2( L0m/2, RxySC[2]);
	phixySC[3] = 0.; 
	phixySC[4] = atan2( L0m/2, RxySC[3]);
	phixySC[5] = atan2(-L0m/2, RxySC[4]);
	 
	 
	MT->o->precision(15);
	for(int i=0; i<6; i++)
		Cout << i << " : RxySC = " << RxySC[i] << " ,  phixySC = " << phixySC[i] << Endl;
	*/
	//! *** Nominal position of spacecrafts
	posRL1[0].x( 0. );
	posRL1[0].y( 0. );
	posRL1[0].z( Lsq2o2 );
	
	posRL1[1].x(  Lo2 );
	posRL1[1].y( -Lo2 );
	posRL1[1].z(  0. );
	
	posRL1[2].x( Lo2 );
	posRL1[2].y( Lo2 );
	posRL1[2].z(  0. );
	
	posRL1[3].x(  0. );
	posRL1[3].y(  0. );
	posRL1[3].z( -Lsq2o2 );
	
	posRL1[4].x( -Lo2 );
	posRL1[4].y(  Lo2 );
	posRL1[4].z(  0. );
	
	posRL1[5].x( -Lo2 );
	posRL1[5].y( -Lo2 );
	posRL1[5].z(  0. );
	
	//! *** Add offset
	for(int iSC=1; iSC<=6; iSC++){
		if((iSC==1)||(iSC==4))
			posRL1[iSC-1] *= (1.+pSCOff[iSC-1]);
		else
			posRL1[iSC-1] *= (1.+pSCOff[iSC-1]/sqrt(2));
	}
	
	
	/*! *** Store the value corresponding at t = 0 */
	
	tStorePos = 0.;
	for(int i=1; i<=NSC; i++)
		SCposStore[i-1] = position(i,0.0);
	
	tStoreArm = 0.;
	int iArm(0);
	for(int iem=1; iem<=NSC; iem++)
		for(int ire=1; ire<=NSC; ire++)
			if(iem!=ire)
				ArmStore[iArm++] = ArmCompute(iem, ire, 0.0);
}


double LCOrbitsAnaOctahedron::ArmSpecific(int em, int rec, double trec)
{
	//*! Not used
	return(0.0);
}


LCVector LCOrbitsAnaOctahedron::position(int iSC, double t)
{
	LCVector r(MT);

	
	//r.p[0] = RxySC[iSC-1] * cos(omega*t + phixySC[iSC-1]);
	//r.p[1] = RxySC[iSC-1] * sin(omega*t + phixySC[iSC-1]);
	
	//! *** Position of L1
	r.p[0] = Rgc * cos(omega*t);
	r.p[1] = Rgc * sin(omega*t);
	r.p[2] = 0.;
	
	//! *** Add spacecraft position relative to L1
	r += posRL1[iSC-1];
		
	return r;
}


LCVector LCOrbitsAnaOctahedron::velocity(int iSC, double t)
{
	LCVector v(MT);
	
	
	/*! *** Computation of velocity :
	 */
	
	
	v.p[0] = 0.;
	v.p[1] = 0.;
	v.p[2] = 0.;
	
	
	return v;
}

void LCOrbitsAnaOctahedron::DispInfo(char * BTab)
{
	if(MT->Disp()){
		Cout << BTab << "Analytic orbit :" << Endl;
		DispInfoBase(BTab);
		Cout << "    - Armlength             = " << L0m << " m" << Endl;
	}
}	


// ***********************
// ***  Local mehtods  ***
// ***********************




// ********************
// *  Access methods  *
// ********************




// ********************
// *  Others methods  *
// ********************



// end of LISACODE-OrbitsOctaLCOrbitsAnaOctahedronAnalytic.cpp