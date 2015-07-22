/*
 *  LISACODE-LCTDIGen.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 20/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


#include "LISACODE-TDIGen.h"



// *********************
// ***  Constructor  ***
// *********************


LCTDIGen::LCTDIGen(LCTools * MT_n)
{
	MT = MT_n;
	initNULL(false);
}


LCTDIGen::~LCTDIGen()
{
	initNULL(true);
}



void LCTDIGen::initNULL(bool CleanMem)
{
	if(CleanMem){
		
		if(p != NULL){
			for(int i=0; i<Np; i++)
				if(p[i].Del != NULL)
					MT->Free(p[i].Del, p[i].NDel*sizeof(LCSerie2*));
			MT->Free(p, Np*sizeof(TDIPack));
		}
		
		if((Result != NULL)&&(ResultLocalAlloc))
			MT->Free(Result, sizeof(double));
	}
	
	p = NULL;
	Np = 0;
	Result = NULL;
	ResultLocalAlloc = true;
	
	tShiftDelay = 0.;
	tShiftSig = 0.;
	DtShiftSigDelay = 0.;
	
	InterpType = LAG;
	InterpUtilValue = 20;
	InterpTypeDelay = LAG;
	InterpUtilValueDelay = 6;
	
	powDigit = 10;
	
	NbMaxDelays = 0;
	strcpy(Name,"");
	
}



// ***************************
// *  Configuration methods  *
// ***************************


void LCTDIGen::config(ezxml_t xmlbloc)
{
	ezxml_t param;
	char * TDIGenType(NULL);
	
	//! *** Set the name
	int NMin( MIN(7, strlen(ezxml_attr(xmlbloc,"Name"))) );
	//Cout << ezxml_attr(xmlbloc,"Name") << " : " << NMin << Endl;
	strncpy(Name, ezxml_attr(xmlbloc,"Name"), NMin );
	Name[NMin] = '\0';
	
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		
		//! **** Checking that it's the right configuration of optical bench + phasemeter
		if(MT->wcmp(ezxml_attr(param,"Name"),"GeneratorType")){
			MT->stripcopy((*param).txt, TDIGenType);
			if(MT->wcmp(TDIGenType, "Preregistred"))
				PreRegistred(Name);
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
		if(MT->wcmp(ezxml_attr(param,"Name"),"InterpolationDelay")){
			if(MT->wcmp(ezxml_attr(param,"Type"),"Lagrange"))
				InterpTypeDelay = LAG;
			if(MT->wcmp(ezxml_attr(param,"Type"),"Linear"))
				InterpTypeDelay = LIN;
			if(MT->wcmp(ezxml_attr(param,"Type"),"Trunc"))
				InterpTypeDelay = TRU;
			InterpUtilValueDelay = atoi((*param).txt);
			
		}
	}
	
	if(TDIGenType != NULL)
		MT->Free(TDIGenType, (strlen(TDIGenType)+1) * sizeof(char) );
	
}


void LCTDIGen::addPack(int CodePack, double Fact)
{
	Np++;
	p = (TDIPack*) MT->ReAllocMemory(p, (Np-1)*sizeof(TDIPack), Np*sizeof(TDIPack) );
	p[Np-1].Code = abs(CodePack);
	p[Np-1].Fact = Fact;
	p[Np-1].Sig = NULL;
	p[Np-1].Del = NULL;
	p[Np-1].NDel = 0;
	if(CodePack<0)
		p[Np-1].Fact = -p[Np-1].Fact;
	
}


void LCTDIGen::PreRegistred(const char * generatorname)
{
	bool UnKnowTDI(true);
	powDigit = 10;
	
	//! ***** TDI generator without factor applied on packs
	if((strcmp(generatorname,"Alpha")==0)||(strcmp(generatorname,"alpha")==0)){
		UnKnowTDI = false;
		addPack(-1,1.);
		addPack(-32,1.);
		addPack(-133,1.);
		addPack(4,1.);
		addPack(455,1.);
		addPack(56,1.);
		if(2>NbMaxDelays)
			NbMaxDelays = 2;
	}
	if((strcmp(generatorname,"Beta")==0)||(strcmp(generatorname,"beta")==0)){
		UnKnowTDI = false;
		addPack(-121,1.);
		addPack(-2,1.);
		addPack(-13,1.);
		addPack(64,1.);
		addPack(5,1.);
		addPack(566,1.);
		if(2>NbMaxDelays)
			NbMaxDelays = 2;
	}
	if((strcmp(generatorname,"Gamma")==0)||(strcmp(generatorname,"gamma")==0)){
		UnKnowTDI = false;
		addPack(-21,1.);
		addPack(-232,1.);
		addPack(-3,1.);
		addPack(464,1.);
		addPack(45,1.);
		addPack(6,1.);
		if(2>NbMaxDelays)
			NbMaxDelays = 2;
	}
	if((strcmp(generatorname,"Zeta")==0)||(strcmp(generatorname,"zeta")==0)){
		UnKnowTDI = false;
		addPack(-11,1.);
		addPack(-22,1.);
		addPack(-33,1.);
		addPack(44,1.);
		addPack(55,1.);
		addPack(66,1.);
		if(1>NbMaxDelays)
			NbMaxDelays = 1;
	}
	if((strcmp(generatorname,"X1s1")==0)||(strcmp(generatorname,"X")==0)||(strcmp(generatorname,"Xf")==0)){
		UnKnowTDI = false;
		addPack(1,1.);
		addPack(35,1.);
		addPack(364,1.);
		addPack(3653,1.);
		addPack(-4,1.);
		addPack(-53,1.);
		addPack(-521,1.);
		addPack(-5235,1.);
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if((strcmp(generatorname,"X1s2")==0)||(strcmp(generatorname,"Y")==0)||(strcmp(generatorname,"Yf")==0)){
		UnKnowTDI = false;
		addPack(2,1.);
		addPack(16,1.);
		addPack(145,1.);
		addPack(1461,1.);
		addPack(-5,1.);
		addPack(-61,1.);
		addPack(-632,1.);
		addPack(-6316,1.);
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if((strcmp(generatorname,"X1s3")==0)||(strcmp(generatorname,"Z")==0)||(strcmp(generatorname,"Zf")==0)){
		UnKnowTDI = false;
		addPack(3,1.);
		addPack(24,1.);
		addPack(256,1.);
		addPack(2542,1.);
		addPack(-6,1.);
		addPack(-42,1.);
		addPack(-413,1.);
		addPack(-4124,1.);
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if((strcmp(generatorname,"X2s1")==0)||(strcmp(generatorname,"X2")==0)||(strcmp(generatorname,"Xf2")==0)){
		UnKnowTDI = false;
		addPack(1,1.);
		addPack(35,1.);
		addPack(364,1.);
		addPack(3653,1.);
		addPack(36524,1.);
		addPack(365253,1.);
		addPack(3652521,1.);
		addPack(36525235,1.);
		addPack(-4,1.);
		addPack(-53,1.);
		addPack(-521,1.);
		addPack(-5235,1.);
		addPack(-52361,1.);
		addPack(-523635,1.);
		addPack(-5236364,1.);
		addPack(-52363653,1.);
		if(7>NbMaxDelays)
			NbMaxDelays = 7;
	}
	if((strcmp(generatorname,"X2s2")==0)||(strcmp(generatorname,"Y2")==0)||(strcmp(generatorname,"Yf2")==0)){
		UnKnowTDI = false;
		addPack(2,1.);
		addPack(16,1.);
		addPack(145,1.);
		addPack(1461,1.);
		addPack(14635,1.);
		addPack(146361,1.);
		addPack(1463632,1.);
		addPack(14636316,1.);
		addPack(-5,1.);
		addPack(-61,1.);
		addPack(-632,1.);
		addPack(-6316,1.);
		addPack(-63142,1.);
		addPack(-631416,1.);
		addPack(-6314145,1.);
		addPack(-63141461,1.);
		if(7>NbMaxDelays)
			NbMaxDelays = 7;
	}
	if((strcmp(generatorname,"X2s3")==0)||(strcmp(generatorname,"Z2")==0)||(strcmp(generatorname,"Zf2")==0)){
		UnKnowTDI = false;
		addPack(3,1.);
		addPack(24,1.);
		addPack(256,1.);
		addPack(2542,1.);
		addPack(25416,1.);
		addPack(254142,1.);
		addPack(2541413,1.);
		addPack(25414124,1.);
		addPack(-6,1.);
		addPack(-42,1.);
		addPack(-413,1.);
		addPack(-4124,1.);
		addPack(-41253,1.);
		addPack(-412524,1.);
		addPack(-4125256,1.);
		addPack(-41252542,1.);
		if(7>NbMaxDelays)
			NbMaxDelays = 7;
	}
	// Beacon
	if(strcmp(generatorname,"P1")==0){
		UnKnowTDI = false;
		addPack(25,1.);
		addPack(-63,1.);
		addPack(-22,1.);
		addPack(66,1.);
		addPack(642,1.);
		addPack(-216,1.);
		addPack(1463,1.);
		addPack(-1425,1.);
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if(strcmp(generatorname,"Q1")==0){
		UnKnowTDI = false;
		addPack(36,1.);
		addPack(-41,1.);
		addPack(-33,1.);
		addPack(44,1.);
		addPack(453,1.);
		addPack(-324,1.);
		addPack(2541,1.);
		addPack(-2536,1.);
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if(strcmp(generatorname,"R1")==0){
		UnKnowTDI = false;
		addPack(14,1.);
		addPack(-52,1.);
		addPack(-11,1.);
		addPack(55,1.);
		addPack(561,1.);
		addPack(-135,1.);
		addPack(3652,1.);
		addPack(-3614,1.);
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	// Monitor
	if(strcmp(generatorname,"E1")==0){
		UnKnowTDI = false;
		addPack(542,1.);
		addPack(56,1.);
		addPack(-316,1.);
		addPack(-32,1.);
		addPack(-144,1.);
		addPack(141,1.);
		addPack(4,1.);
		addPack(-1,1.);
		if(2>NbMaxDelays)
			NbMaxDelays = 2;
	}
	if(strcmp(generatorname,"F1")==0){
		UnKnowTDI = false;
		addPack(653,1.);
		addPack(64,1.);
		addPack(-124,1.);
		addPack(-13,1.);
		addPack(-255,1.);
		addPack(252,1.);
		addPack(5,1.);
		addPack(-2,1.);
		if(2>NbMaxDelays)
			NbMaxDelays = 2;
	}
	if(strcmp(generatorname,"G1")==0){
		UnKnowTDI = false;
		addPack(461,1.);
		addPack(45,1.);
		addPack(-235,1.);
		addPack(-21,1.);
		addPack(-366,1.);
		addPack(363,1.);
		addPack(6,1.);
		addPack(-3,1.);
		if(2>NbMaxDelays)
			NbMaxDelays = 2;
	}
	// Relay
	if(strcmp(generatorname,"U1")==0){
		UnKnowTDI = false;
		addPack(145,1.);
		addPack(1464,1.);
		addPack(-5,1.);
		addPack(-64,1.);
		addPack(16,1.);
		addPack(2,1.);
		addPack(-6542,1.);
		addPack(-656,1.);
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if(strcmp(generatorname,"V1")==0){
		UnKnowTDI = false;
		addPack(256,1.);
		addPack(2545,1.);
		addPack(-6,1.);
		addPack(-45,1.);
		addPack(24,1.);
		addPack(3,1.);
		addPack(-4653,1.);
		addPack(-464,1.);
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if(strcmp(generatorname,"W1")==0){
		UnKnowTDI = false;
		addPack(364,1.);
		addPack(3656,1.);
		addPack(-4,1.);
		addPack(-56,1.);
		addPack(35,1.);
		addPack(1,1.);
		addPack(-5461,1.);
		addPack(-545,1.);
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if(strcmp(generatorname,"Zeta1")==0){
		UnKnowTDI = false;
		addPack(3214,1.);
		addPack(-414,1.);
		addPack(-3252,1.);
		addPack(452,1.);
		addPack(3255,1.);
		addPack(-455,1.);
		addPack(-5633,1.);
		addPack(133,1.);
		addPack(5636,1.);
		addPack(-136,1.);
		addPack(-5641,1.);
		addPack(141,1.);
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if(strcmp(generatorname,"Zeta2")==0){
		UnKnowTDI = false;
		addPack(1325,1.);
		addPack(-525,1.);
		addPack(-1363,1.);
		addPack(563,1.);
		addPack(1366,1.);
		addPack(-566,1.);
		addPack(-6411,1.);
		addPack(211,1.);
		addPack(6414,1.);
		addPack(-214,1.);
		addPack(-6452,1.);
		addPack(252,1.);
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if(strcmp(generatorname,"Zeta3")==0){
		UnKnowTDI = false;
		addPack(2136,1.);
		addPack(-636,1.);
		addPack(-2141,1.);
		addPack(641,1.);
		addPack(2144,1.);
		addPack(-644,1.);
		addPack(-4522,1.);
		addPack(322,1.);
		addPack(4525,1.);
		addPack(-325,1.);
		addPack(-4563,1.);
		addPack(363,1.);
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	
	//! ***** TDI generator with factor applied on packs
	if((strcmp(generatorname,"A")==0)||(strcmp(generatorname,"Aa")==0)||(strcmp(generatorname,"AX")==0)){
		UnKnowTDI = false;
		// Packs of X * -1/sqrt(2)
		addPack(1, -1.0/sqrt(2.));
		addPack(35, -1.0/sqrt(2.));
		addPack(364, -1.0/sqrt(2.));
		addPack(3653, -1.0/sqrt(2.));
		addPack(-4, -1.0/sqrt(2.));
		addPack(-53, -1.0/sqrt(2.));
		addPack(-521, -1.0/sqrt(2.));
		addPack(-5235, -1.0/sqrt(2.));
		//! ** Pack of Z * -1/3
		addPack(3, -1.0/sqrt(2.));
		addPack(24, -1.0/sqrt(2.));
		addPack(256, -1.0/sqrt(2.));
		addPack(2542, -1.0/sqrt(2.));
		addPack(-6, -1.0/sqrt(2.));
		addPack(-42, -1.0/sqrt(2.));
		addPack(-413, -1.0/sqrt(2.));
		addPack(-4124, -1.0/sqrt(2.));
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if((strcmp(generatorname,"E")==0)||(strcmp(generatorname,"Ea")==0)||(strcmp(generatorname,"EX")==0)){
		UnKnowTDI = false;
		// Packs of X * 1./sqrt(6.)
		addPack(1, 1./sqrt(6.) );
		addPack(35, 1./sqrt(6.) );
		addPack(364, 1./sqrt(6.) );
		addPack(3653, 1./sqrt(6.) );
		addPack(-4, 1./sqrt(6.) );
		addPack(-53, 1./sqrt(6.) );
		addPack(-521, 1./sqrt(6.) );
		addPack(-5235, 1./sqrt(6.) );
		// Packs of Y * -2.0/sqrt(6.0)
		addPack(2, -2.0/sqrt(6.0) );
		addPack(16, -2.0/sqrt(6.0) );
		addPack(145, -2.0/sqrt(6.0) );
		addPack(1461, -2.0/sqrt(6.0) );
		addPack(-5, -2.0/sqrt(6.0) );
		addPack(-61, -2.0/sqrt(6.0) );
		addPack(-632, -2.0/sqrt(6.0) );
		addPack(-6316, -2.0/sqrt(6.0) );
		//! ** Pack of Z * 1./sqrt(6.)
		addPack(3, 1./sqrt(6.) );
		addPack(24, 1./sqrt(6.) );
		addPack(256, 1./sqrt(6.) );
		addPack(2542, 1./sqrt(6.) );
		addPack(-6, 1./sqrt(6.) );
		addPack(-42, 1./sqrt(6.) );
		addPack(-413, 1./sqrt(6.) );
		addPack(-4124, 1./sqrt(6.) );
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if((strcmp(generatorname,"T")==0)||(strcmp(generatorname,"Ta")==0)||(strcmp(generatorname,"TX")==0)){
		UnKnowTDI = false;
		// Packs of X * 1./sqrt(3.)
		addPack(1, 1./sqrt(3.) );
		addPack(35, 1./sqrt(3.) );
		addPack(364, 1./sqrt(3.) );
		addPack(3653, 1./sqrt(3.) );
		addPack(-4, 1./sqrt(3.) );
		addPack(-53, 1./sqrt(3.) );
		addPack(-521, 1./sqrt(3.) );
		addPack(-5235, 1./sqrt(3.) );
		// Packs of Y * 1./sqrt(3.)
		addPack(2, 1./sqrt(3.) );
		addPack(16, 1./sqrt(3.) );
		addPack(145, 1./sqrt(3.) );
		addPack(1461, 1./sqrt(3.) );
		addPack(-5, 1./sqrt(3.) );
		addPack(-61, 1./sqrt(3.) );
		addPack(-632, 1./sqrt(3.) );
		addPack(-6316, 1./sqrt(3.) );
		//! ** Pack of Z * 1./sqrt(3.)
		addPack(3, 1./sqrt(3.) );
		addPack(24, 1./sqrt(3.) );
		addPack(256, 1./sqrt(3.) );
		addPack(2542, 1./sqrt(3.) );
		addPack(-6, 1./sqrt(3.) );
		addPack(-42, 1./sqrt(3.) );
		addPack(-413, 1./sqrt(3.) );
		addPack(-4124, 1./sqrt(3.) );
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if((strcmp(generatorname,"Am")==0)||strcmp(generatorname,"AMLDC")==0){
		UnKnowTDI = false;
		// Packs of X * 2/3
		addPack(1, 2./3. );
		addPack(35, 2./3. );
		addPack(364, 2./3. );
		addPack(3653, 2./3. );
		addPack(-4, 2./3. );
		addPack(-53, 2./3. );
		addPack(-521, 2./3. );
		addPack(-5235, 2./3. );
		// Packs of Y * -1/3
		addPack(2, -1./3. );
		addPack(16, -1./3. );
		addPack(145, -1./3. );
		addPack(1461, -1./3. );
		addPack(-5, -1./3. );
		addPack(-61, -1./3. );
		addPack(-632, -1./3. );
		addPack(-6316, -1./3. );
		//! ** Pack of Z * -1/3
		addPack(3, -1./3. );
		addPack(24, -1./3. );
		addPack(256, -1./3. );
		addPack(2542, -1./3. );
		addPack(-6, -1./3. );
		addPack(-42, -1./3. );
		addPack(-413, -1./3. );
		addPack(-4124, -1./3. );
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if((strcmp(generatorname,"Em")==0)||(strcmp(generatorname,"EMLDC")==0)){
		UnKnowTDI = false;
		// Packs of Y * -1/sqrt(3)
		addPack(2, -1./sqrt(3.) );
		addPack(16, -1./sqrt(3.) );
		addPack(145, -1./sqrt(3.) );
		addPack(1461, -1./sqrt(3.) );
		addPack(-5, -1./sqrt(3.) );
		addPack(-61, -1./sqrt(3.) );
		addPack(-632, -1./sqrt(3.) );
		addPack(-6316, -1./sqrt(3.) );
		//! ** Pack of Z * 1/sqrt(3)
		addPack(3, 1./sqrt(3.) );
		addPack(24, 1./sqrt(3.) );
		addPack(256, 1./sqrt(3.) );
		addPack(2542, 1./sqrt(3.) );
		addPack(-6, 1./sqrt(3.) );
		addPack(-42, 1./sqrt(3.) );
		addPack(-413, 1./sqrt(3.) );
		addPack(-4124, 1./sqrt(3.) );
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	if((strcmp(generatorname,"Tm")==0)||(strcmp(generatorname,"TMLDC")==0)){
		UnKnowTDI = false;
		// Packs of X * sqrt(2)/3
		addPack(1, sqrt(2.0)/3.0 );
		addPack(35, sqrt(2.0)/3.0 );
		addPack(364, sqrt(2.0)/3.0 );
		addPack(3653, sqrt(2.0)/3.0 );
		addPack(-4, sqrt(2.0)/3.0 );
		addPack(-53, sqrt(2.0)/3.0 );
		addPack(-521, sqrt(2.0)/3.0 );
		addPack(-5235, sqrt(2.0)/3.0 );
		// Packs of Y * sqrt(2)/3
		addPack(2, sqrt(2.0)/3.0 );
		addPack(16, sqrt(2.0)/3.0 );
		addPack(145, sqrt(2.0)/3.0 );
		addPack(1461, sqrt(2.0)/3.0 );
		addPack(-5, sqrt(2.0)/3.0 );
		addPack(-61, sqrt(2.0)/3.0 );
		addPack(-632, sqrt(2.0)/3.0 );
		addPack(-6316, sqrt(2.0)/3.0 );
		//! ** Pack of Z * sqrt(2)/3
		addPack(3, sqrt(2.0)/3.0 );
		addPack(24, sqrt(2.0)/3.0 );
		addPack(256, sqrt(2.0)/3.0 );
		addPack(2542, sqrt(2.0)/3.0 );
		addPack(-6, sqrt(2.0)/3.0 );
		addPack(-42, sqrt(2.0)/3.0 );
		addPack(-413, sqrt(2.0)/3.0 );
		addPack(-4124, sqrt(2.0)/3.0 );
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	
	//! ******************* TDI Octahedron
	
	if((strcmp(generatorname,"oY1")==0)){
		UnKnowTDI = false;
		powDigit = 100;
		
		addPack( 12, 1.);
		addPack( 15, 1.);
		addPack(-13, 1.);
		addPack(-16, 1.);
		
		addPack(-42, 1.);
		addPack(-45, 1.);
		addPack( 43, 1.);
		addPack( 46, 1.);
		
		addPack(-2421, 1.);
		addPack( 2124, 1.);
		////addPack(23, 1.);
		////addPack(26, 1.);
		
		addPack(-5451, 1.);
		addPack( 5154, 1.);
		////addPack(53, 1.);
		////addPack(56, 1.);
		
		addPack( 3431, 1.);
		addPack(-3134, 1.);
		////addPack(32, 1.);
		////addPack(35, 1.);
		
		addPack( 6461, 1.);
		addPack(-6164, 1.);
		////addPack(62, 1.);
		////addPack(65, 1.);
	}
	
	if((strcmp(generatorname,"oY2")==0)){
		UnKnowTDI = false;
		powDigit = 100;
		
		addPack(-1512, 1.);
		addPack( 1215, 1.);
		////addPack( 13, 1.);
		////addPack( 16, 1.);
		
		addPack(-4542, 1.);
		addPack( 4245, 1.);
		////addPack( 43, 1.);
		////addPack( 46, 1.);
		
		addPack( 21, 1.);
		addPack( 24, 1.);
		addPack(-23, 1.);
		addPack(-26, 1.);
		
		addPack(-51, 1.);
		addPack(-54, 1.);
		addPack( 53, 1.);
		addPack( 56, 1.);
		
		////addPack( 31, 1.);
		////addPack( 34, 1.);
		addPack( 3532, 1.);
		addPack(-3235, 1.);
		
		////addPack( 61, 1.);
		////addPack( 64, 1.);
		addPack( 6562, 1.);
		addPack(-6265, 1.);
		
		
		if(2>NbMaxDelays)
			NbMaxDelays = 2;
	}
	
	if((strcmp(generatorname,"oY3")==0)){
		UnKnowTDI = false;
		powDigit = 100;
		
		////addPack( 12, 1.);
		////addPack( 15, 1.);
		addPack( 1613, 1.);
		addPack(-1316, 1.);
		
		////addPack( 42, 1.);
		////addPack( 45, 1.);
		addPack( 4643, 1.);
		addPack(-4346, 1.);
		
		////addPack( 21, 1.);
		////addPack( 24, 1.);
		addPack(-2623, 1.);
		addPack( 2326, 1.);
		
		////addPack( 51, 1.);
		////addPack( 54, 1.);
		addPack(-5653, 1.);
		addPack( 5356, 1.);
		
		addPack(-31, 1.);
		addPack(-34, 1.);
		addPack( 32, 1.);
		addPack( 35, 1.);
		
		addPack( 61, 1.);
		addPack( 64, 1.);
		addPack(-62, 1.);
		addPack(-65, 1.);
		
		
		if(2>NbMaxDelays)
			NbMaxDelays = 2;
	}
	
	if((strcmp(generatorname,"oS10")==0)){
		UnKnowTDI = false;
		powDigit = 100;
		
		////addPack( 12, 1.);
		addPack( 1215, 1.);
		addPack( 121215, 1.);
		////addPack( 13, 1.);
		addPack(-1216, 1.);
		addPack(-121216, 1.);
		
		addPack( 42, 1.);
		addPack(-1242, 1.);
		addPack( 45, 1.);
		addPack( 121245, 1.);
		addPack(-43, 1.);
		addPack( 1243, 1.);
		addPack(-46, 1.);
		addPack(-121246, 1.);
		
		addPack( 1221, 1.);
		addPack(-121221, 1.);
		////addPack( 24, 1.);
		addPack(-1223, 1.);
		addPack( 121226, 1.);
		
		addPack(-51, 1.);
		addPack(-121251, 1.);
		addPack(-54, 1.);
		addPack(-1254, 1.);
		addPack( 53, 1.);
		addPack( 56, 1.);
		addPack( 1256, 1.);
		addPack( 121256, 1.);
		
		addPack(-1231, 1.);
		addPack( 121231, 1.);
		////addPack( 34, 1.);
		addPack( 1232, 1.);
		addPack(-121235, 1.);
		
		addPack( 61, 1.);
		addPack( 121261, 1.);
		addPack( 64, 1.);
		addPack( 1264, 1.);
		addPack(-62, 1.);
		addPack(-65, 1.);
		addPack(-1265, 1.);
		addPack(-121265, 1.);
		
		
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	
	
	if((strcmp(generatorname,"oS20")==0)){
		UnKnowTDI = false;
		powDigit = 100;
		
		////addPack( 12, 1.);
		addPack( 1215, 1.);
		addPack(-1213, 1.);
		////addPack( 16, 1.);
		
		addPack(-42, 1.);
		addPack(-45, 1.);
		addPack( 1245, 1.);
		addPack( 43, 1.);
		addPack(-1243, 1.);
		addPack( 46, 1.);
		
		addPack( 21, 1.);
		addPack(-1221, 1.);
		addPack( 24, 1.);
		addPack(-23, 1.);
		addPack( 1223, 1.);
		addPack(-26, 1.);
		
		addPack(-1251, 1.);
		////addPack( 54, 1.);
		addPack( 1253, 1.);
		////addPack( 56, 1.);
		
		addPack( 1231, 1.);
		////addPack( 34, 1.);
		////addPack( 32, 1.);
		addPack(-1235, 1.);
		
		addPack(-61, 1.);
		addPack( 1261, 1.);
		addPack(-64, 1.);
		addPack( 62, 1.);
		addPack( 65, 1.);
		addPack(-1265, 1.);
		
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	
	
	if((strcmp(generatorname,"oS30")==0)){
		UnKnowTDI = false;
		powDigit = 100;
		
		////addPack( 12, 1.);
		addPack( 1215, 1.);
		addPack( 121215, 1.);
		addPack(-121213, 1.);
		addPack(-1216, 1.);
		
		addPack( 42, 1.);
		addPack(-1242, 1.);
		addPack( 45, 1.);
		addPack( 121245, 1.);
		addPack(-43, 1.);
		addPack( 1243, 1.);
		addPack(-121243, 1.);
		addPack(-46, 1.);
		
		addPack( 1221, 1.);
		addPack(-121221, 1.);
		////addPack( 24, 1.);
		addPack(-1223, 1.);
		addPack( 121223, 1.);
		////addPack( 26, 1.);
		
		addPack(-51, 1.);
		addPack(-121251, 1.);
		addPack(-54, 1.);
		addPack(-1254, 1.);
		addPack( 53, 1.);
		addPack( 121253, 1.);
		addPack( 56, 1.);
		addPack( 1256, 1.);
		
		addPack( 121231, 1.);
		addPack( 1234, 1.);
		////addPack( 32, 1.);
		addPack(-1235, 1.);
		addPack(-121235, 1.);
		
		addPack( 61, 1.);
		addPack(-1261, 1.);
		addPack( 121261, 1.);
		addPack( 64, 1.);
		addPack(-62, 1.);
		addPack( 1262, 1.);
		addPack(-65, 1.);
		addPack(-121265, 1.);
		
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	
	
	if((strcmp(generatorname,"oS4q1mD")==0)){
		UnKnowTDI = false;
		powDigit = 100;
		
		addPack( 12, 1.);
		addPack( -1212, 1.);
		addPack( 15, 1.);
		addPack( -121215, 1.);
		addPack( -13, 1.);
		addPack( 1213, 1.);
		addPack( -16, 1.);
		addPack( 121216, 1.);
		
		
		////addPack( 42, 1.);
		addPack( 1245, 1.);
		addPack( -121245, 1.);
		////addPack( 43, 1.);
		addPack( -1246, 1.);
		addPack( 121246, 1.);
		
		addPack( -1221, 1.);
		addPack( 121221, 1.);
		////addPack( 24, 1.);
		////addPack( 23, 1.);
		addPack( 1226, 1.);
		addPack( -121226, 1.);
		
		addPack( -51, 1.);
		addPack( 121251, 1.);
		addPack( -54, 1.);
		addPack( 1254, 1.);
		addPack( 53, 1.);
		addPack( -1253, 1.);
		addPack( 56, 1.);
		addPack( -121256, 1.);
		
		addPack( 1231, 1.);
		addPack( -121231, 1.);
		////addPack( 34, 1.);
		////addPack( 32, 1.);
		addPack( -1235, 1.);
		addPack( 121235, 1.);
		
		addPack( 61, 1.);
		addPack( -121261, 1.);
		addPack( 64, 1.);
		addPack( -1264, 1.);
		addPack( -62, 1.);
		addPack( 1262, 1.);
		addPack( -65, 1.);
		addPack( 121265, 1.);
		
		 
		if(4>NbMaxDelays)
			NbMaxDelays = 4;
	}
	
	
	/*
	if((strcmp(generatorname,"oNew")==0)){
		UnKnowTDI = false;
		powDigit = 100;
		
		addPack( 12, 1.);
		addPack( 15, 1.);
		addPack( 13, 1.);
		addPack( 16, 1.);
		
		addPack( 42, 1.);
		addPack( 45, 1.);
		addPack( 43, 1.);
		addPack( 46, 1.);
		
		addPack( 21, 1.);
		addPack( 24, 1.);
		addPack( 23, 1.);
		addPack( 26, 1.);
		
		addPack( 51, 1.);
		addPack( 54, 1.);
		addPack( 53, 1.);
		addPack( 56, 1.);
		
		addPack( 31, 1.);
		addPack( 34, 1.);
		addPack( 32, 1.);
		addPack( 35, 1.);
		
		addPack( 61, 1.);
		addPack( 64, 1.);
		addPack( 62, 1.);
		addPack( 65, 1.);
		
		if(3>NbMaxDelays)
			NbMaxDelays = 3;
	}
	*/
	
	if(UnKnowTDI){
		std::cerr << "ERROR in LCTDIGen::PreRegistred : The TDI generator " << generatorname << " is unknown !" << Endl;
		throw std::invalid_argument("ERROR in LCTDIGen::PreRegistred : The TDI generator  is unknown !");
	}
}
	
	
	
	// ***************************************
	// * Linking and initialization methods  *
	// ***************************************
	
	int LCTDIGen::iSigCodePack(int iP)
	{
		double pow10NDigit(1.*powDigit);
		if(iP>=Np)
			throw std::invalid_argument("ERROR in LCTDIGen::getSigPack : The pack doesn't exist");
		//Cout << p[iP].Code << "  :  " << (abs(p[iP].Code)/10.0) << "  -  " << floor(abs(p[iP].Code)/10.0) << "  ==>  " << MT->ifloor(10.0*((abs(p[iP].Code)/10.0)-floor(abs(p[iP].Code)/10.0))+1.e-6) << Endl;
		
		return(MT->ifloor( pow10NDigit*((abs(p[iP].Code)/pow10NDigit) - floor(abs(p[iP].Code)/pow10NDigit)) + 1.e-6 ) );
	}
	
	
	void LCTDIGen::LinkSigDetector(LCDetector * LISA)
	{
		//! **** Loop over the packs
		for(int iP=0; iP<Np; iP++){
			
			//! *** Extract the index of the signal or receiver
			int iS = iSigCodePack(iP);
			
			//if(MT->DispDet())
			//	Cout << " LinkSigDetector :: " << p[iP].Code << " ==> Sig " << iS << Endl;
			
			if(powDigit==10){
				if(iS<=3){
					//! ** Link in direct direction 
					p[iP].Sig = LISA->LinkPhaMes(iS, 0, "sci", 3);
				}else{
					//! ** Link in indirect direction 
					p[iP].Sig = LISA->LinkPhaMes(iS-3, 1, "sci", 3);
				}
			}else{
				p[iP].Sig = LISA->LinkPhaMes( (iS/10), -((iS-10*((int)(iS/10)))), "sci", 3);
			}
		}
		
	}
	
	
	void LCTDIGen::LinkSigTDIInt(LCTDIInt ** TDII, int NTDII)
	{
		//! **** Loop over the packs
		for(int iP=0; iP<Np; iP++){
			p[iP].Sig = NULL;
			
			//! *** Extract the index of the signal
			int iS = iSigCodePack(iP);
			
			//if(MT->DispDet())
			//	Cout << " LinkSigTDIInt :: " << p[iP].Code <<  " ==> Eta " << iS << Endl;
			
			for(int i=0; i<NTDII; i++){
				if(TDII[i] != NULL){
					if(powDigit==10){
						if(iS<=3){
							//! ** Link in direct direction
							if((!TDII[i]->getIndirectDir())&&(TDII[i]->getiSC()==iS))
								p[iP].Sig = TDII[i]->getDat();
						}else{
							//! ** Link in indirect direction
							if((TDII[i]->getIndirectDir())&&(TDII[i]->getiSC()==iS-3))
								p[iP].Sig = TDII[i]->getDat();
						}
					}else{
						if( (TDII[i]->getiSC()==(iS/10)) && (TDII[i]->getIndirectDir()==-((iS-10*((int)(iS/10))))) )
							p[iP].Sig = TDII[i]->getDat();
					}
					
				}
			}
			
			if(p[iP].Sig == NULL){
				Cout << "ERROR in LCTDIGen::LinkSigTDIInt : We don't find any signal corresponding to index " << iS << " !" << Endl;
				throw std::invalid_argument("ERROR in LCTDIGen::LinkSigTDIInt : We don't find any signal corresponding to index : it should be in [1,6] !" );
			}
		}
	}
	
	
	void LCTDIGen::LinkDelays(LCSerie2** AllDelays, int NAllDelays, char ** NameDelays)
	{
		int TmpInfo, TmpCodeDelay, TmpIndexDelay;
		char TmpC[10];
		
		//! **** Loop over the packs
		for(int iP=0; iP<Np; iP++){
			p[iP].NDel = 0;
			p[iP].Del = NULL;
			
			//\\if(MT->DispDet())
			//\\	Cout << "LinkDelays :: " << p[iP].Code << " : "; 
			TmpInfo = MT->ifloor(abs(p[iP].Code)/powDigit);
			//\\Cout << iP << " : ";
			while(TmpInfo != 0){
				//\\if(MT->DispDet())
				//\\	Cout << " (" << TmpInfo << ") ";
				TmpCodeDelay = MT->ifloor( powDigit * ((TmpInfo/((double)powDigit))-floor(TmpInfo/((double)powDigit))+1e-6)) ;
				TmpC[3]='\0'; 
				if(powDigit==10){
					if(TmpCodeDelay<=3){
						sprintf(TmpC,"D%d",(TmpCodeDelay-1)%3+1);
						TmpC[2]='\0'; 
					}else{ 
						sprintf(TmpC,"D%dp",(TmpCodeDelay-1)%3+1);
					}
				}else{
					sprintf(TmpC,"D%d%d", TmpCodeDelay/10, (TmpCodeDelay-10*((int)(TmpCodeDelay/10))) ); 
				}
				TmpIndexDelay = -1;
				for(int iD=0; iD<NAllDelays; iD++)
					if((TmpC[1]==NameDelays[iD][1])&&(TmpC[2]==NameDelays[iD][2]))
						TmpIndexDelay = iD;
				//\\if(MT->DispDet())
				//\\	Cout << " " << TmpC << "/" << NameDelays[TmpIndexDelay] << " ";
				if(TmpIndexDelay==-1){
					Cout << "ERROR in LCTDIGen::LinkDelay : Cannot find the delay corresponding to " << TmpC[1] << TmpC[2] << " !" << Endl;
					throw std::invalid_argument("ERROR in LCTDIGen::LinkDelay : Cannot find the delay !" );
				}
				
				//\\if(MT->DispDet())
				//\\	Cout  << TmpIndexDelay << " , ";
				p[iP].NDel++;
				p[iP].Del = (LCSerie2**) MT->ReAllocMemory(p[iP].Del, (p[iP].NDel-1)*sizeof(LCSerie2*), p[iP].NDel*sizeof(LCSerie2*));
				p[iP].Del[p[iP].NDel-1] = NULL;
				if( (TmpIndexDelay>=0) && (TmpIndexDelay<NAllDelays) ){
					p[iP].Del[p[iP].NDel-1] = AllDelays[TmpIndexDelay];
				}else{
					Cout << "ERROR in LCTDIGen::LinkDelay : TmpIndexDelay = " << TmpIndexDelay << " is a wrong index of delay : it should be in [0," << NAllDelays-1 << "] !" << Endl;
					throw std::invalid_argument("ERROR in LCTDIGen::LinkDelay : Wrong index of delay : it should be in [0,NArm] !" );
				}
				TmpInfo = MT->ifloor(TmpInfo/((double)powDigit));
			};
			//\\if(MT->DispDet())
			//\\	Cout << Endl;
		}
		
	}
	
	
	void LCTDIGen::init()
	{
		//! *** Allocation of memory for the last measurement 
		if(Result == NULL)
			Result = (double*)MT->AllocMemory(sizeof(double));
		
	}
	
	bool LCTDIGen::MatchName(char * ObsName)
	{
		return(MT->wcmp(ObsName,Name));
	}
	
	
	void LCTDIGen::LinkResult2Output(double * pOutData)
	{
		
		if((Result != NULL)&&(ResultLocalAlloc))
			MT->Free(Result, sizeof(double));
		ResultLocalAlloc = false;
		Result = pOutData;
		
	}
	
	
	
	// ********************
	// *  Access methods  *
	// ********************
	
	
	int LCTDIGen::getNDatInterp()
	{
		return(MT->iceil(InterpUtilValue));
	}
	
	
	int LCTDIGen::getNMaxDelay()
	{
		return(NbMaxDelays);
	}
	
	
	void LCTDIGen::setTimeInfo(double tShiftSig_n, double tShiftDelay_n)
	{
		tShiftSig = tShiftSig_n;
		tShiftDelay = tShiftDelay_n;
		DtShiftSigDelay = tShiftSig - tShiftDelay ;
	}
	
	
	// *********************
	// *  Running methods  *
	// *********************
	
	double LCTDIGen::RunStep(double t)
	{
		//! ***** Apply TDI
		(*Result) = 0.0;
		double TotalDelay(0.0), InterResult(0.0);
		
		//! *** For each pack ...
		for(int iP=0; iP<Np; iP++){
			TotalDelay = tShiftDelay;
			InterResult = 0.0;
			
			//! *** ... compute the total combined delay ...
			for(int iD=p[iP].NDel-1; iD>=0; iD--){
				TotalDelay -= p[iP].Del[iD]->gData(TotalDelay, InterpTypeDelay , InterpUtilValueDelay);
			}
			
			/*! *** Correct the total delay to be relative to the signal : \n
			 * The total delay is done relative to the delay for easy compute of combine delay but it's not the same reference for the signal so we have 
			 * to correct it by the difference DtShiftSigDelay = tShiftDelay - tShiftSig  
			 */
			TotalDelay += DtShiftSigDelay;
			
			//! *** ... make the interpolation ...
			InterResult = p[iP].Fact * p[iP].Sig->gData(TotalDelay, InterpType, InterpUtilValue);
			
			//! *** ... and add to the result (easy no ?)
			(*Result) += InterResult;
			
		}
		return((*Result));
	}
	
	
	// *******************
	// *  Other methods  *
	// *******************
	
	void LCTDIGen::DispInfo(char * BTab)
	{
		if(MT->Disp()){
			Cout << "TDI generator : " << Name << Endl;
			if(MT->DispDet()){
				for(int i=0; i<Np; i++){
					Cout << BTab << "\t- Pack " << i << " : " << p[i].Fact << " x " << p[i].Code;
					Cout << Endl;
				}
			}
		}
	}
	
	
	
	
	
	
	// end of LISACODE-InterTDI.cpp
