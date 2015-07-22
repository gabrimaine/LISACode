/*
 *  LISACODE-OrbitsData.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 10/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-OrbitsData.h"

// *********************
// ***  Constructor  ***
// *********************

/**
 *
 */
LCOrbitsData::LCOrbitsData()
: LCOrbits()
{
	initNULL(false);
}


LCOrbitsData::LCOrbitsData(LCTools * MT_n)
: LCOrbits(MT_n)
{
	initNULL(false);
}


LCOrbitsData::~LCOrbitsData()
{
	initNULL(true);
}




// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCOrbitsData::initNULL(bool CleanMem)
{	
	initNULLBase(CleanMem);
	
	if(CleanMem){
		if(fDat != NULL)
			delete fDat;
	}
	
	NSC = 3;
	
	fDat = NULL;
	InterpType = LAG;
	InterpUtilValue = 3;
	
	
	for(int iSC=0; iSC<3; iSC++){
		for(int i=0; i<3; i++){
			iColPos[iSC][i] = -1;
			iColVel[iSC][i] = -1;
		}
	}
	
	VelocityInFile = false;
	ConvTime = 1.;
	ConvPos = 1.;
	ConvVel = 1.;
	//! ** Default shift of one week
	TimeShift = 7.*LC::Dy_SI; 
	
		
}


void LCOrbitsData::config(ezxml_t orbxmlbloc)
{
	
	
	bool NoOrbFound(true);
	bool FirstTimeRef(true);
	bool Stop(false);
	char NameRec[8192];
	
	
	ezxml_t param, section;
	for(param = ezxml_child(orbxmlbloc,"Param"); param; param = param->next){
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"TimeConversion"))
			ConvTime = MT->gXMLTime(param);
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"PositionConversion"))
			ConvPos = MT->gXMLLength(param);
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"VelocityConversion"))
			ConvVel = atof(ezxml_txt(param));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"NominalArmlength"))
			L0m = MT->gXMLLength(param);
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"TimeShift"))
			TimeShift = MT->gXMLTime(param);
		
	}
	
	if(!MT->wcmp(ezxml_attr(orbxmlbloc,"Type"),"OrbitsInfo")){
		
		if(fDat != NULL)
			delete fDat;
		strcpy(NameRec,"None");
		
		for (section = ezxml_child(orbxmlbloc, "XSIL"); section; section = section->next) {
			if((NoOrbFound)&&(MT->wcmp(ezxml_attr(section,"Type"),"TimeSeries"))){
				strcpy(NameRec,ezxml_attr(section,"Name"));
				
				NoOrbFound = false;
				fDat = new LCDataFileRead(MT);
				fDat->config(section);
			}
		}
		
		
		
		if(NoOrbFound)
			throw std::invalid_argument("ERROR in LCOrbitsData::config : Don't find informations about orbits' file !");
		
		//! *** Extract the columns index from the name
		if(!MT->wcmp(NameRec,"None")){
			
			int iCol(-1);
			
			int ipIA(0);
			if((NameRec[0]='t')&&(NameRec[1]=',')){
				FirstTimeRef = true;
				ipIA = 2;
			}
			
			
			while (!Stop) {
				
				char tmpWord[100];
				int ipL(0);
				
				iCol++;
				
				//! * Copy the word
				while ((NameRec[ipIA]!=',')&&(NameRec[ipIA]!='\0')&&(NameRec[ipIA]!='\n')&&(NameRec[ipIA]!=' ')&&(ipL<99)){
					tmpWord[ipL++] = NameRec[ipIA++];  
				}
				//! * Finish the word
				tmpWord[ipL] = '\0';
				
				//Cout << "fOutObsNames[" << NfOut-1 << "][" << NOutName-1 << "] = " << OutName[NOutName-1] << Endl;
				
				//! * Stop reding condition
				if(NameRec[ipIA]!=',')
					Stop = true;
				else
					ipIA++;
				
				if(MT->wcmp(tmpWord,"x1"))
					iColPos[0][0] = iCol;
				if(MT->wcmp(tmpWord,"y1"))
					iColPos[0][1] = iCol;
				if(MT->wcmp(tmpWord,"z1"))
					iColPos[0][2] = iCol;
				
				if(MT->wcmp(tmpWord,"x2"))
					iColPos[1][0] = iCol;
				if(MT->wcmp(tmpWord,"y2"))
					iColPos[1][1] = iCol;
				if(MT->wcmp(tmpWord,"z2"))
					iColPos[1][2] = iCol;
				
				if(MT->wcmp(tmpWord,"x3"))
					iColPos[2][0] = iCol;
				if(MT->wcmp(tmpWord,"y3"))
					iColPos[2][1] = iCol;
				if(MT->wcmp(tmpWord,"z3"))
					iColPos[2][2] = iCol;
				
				if(tmpWord[0] == 'v'){
					
					VelocityInFile = true;
					
					if(MT->wcmp(tmpWord,"vx1"))
						iColVel[0][0] = iCol;
					if(MT->wcmp(tmpWord,"vy1"))
						iColVel[0][1] = iCol;
					if(MT->wcmp(tmpWord,"vz1"))
						iColVel[0][2] = iCol;
					
					if(MT->wcmp(tmpWord,"vx2"))
						iColVel[1][0] = iCol;
					if(MT->wcmp(tmpWord,"vy2"))
						iColVel[1][1] = iCol;
					if(MT->wcmp(tmpWord,"vz2"))
						iColVel[1][2] = iCol;
					
					if(MT->wcmp(tmpWord,"vx3"))
						iColVel[2][0] = iCol;
					if(MT->wcmp(tmpWord,"vy3"))
						iColVel[2][1] = iCol;
					if(MT->wcmp(tmpWord,"vz3"))
						iColVel[2][2] = iCol;
				}
				
			}
		}
		
		// ** Check that all the component are associated to a columns
		for(int iSC=0; iSC<3; iSC++){
			for(int i=0; i<3; i++){
				if(iColPos[iSC][i] == -1){
					Cout << "ERROR in LCOrbitsData::config : The component " << i+1 << " of the position of the spacecraft " << iSC << " is not defined ! " << Endl; 
					throw std::invalid_argument("ERROR in LCOrbitsData::config : One component of the position of a spacecraft is not defined !");
				}
				if((VelocityInFile)&&(iColVel[iSC][i] == -1)){
					Cout << "ERROR in LCOrbitsData::config : The component " << i+1 << " of the velocity of the spacecraft " << iSC << " is not defined ! " << Endl; 
					throw std::invalid_argument("ERROR in LCOrbitsData::config : One component of the velocity of a spacecraft is not defined !");
				}
			}
		}
	}
}
	   

void LCOrbitsData::config(int iParam, double ParamVal)
{
	
}


void LCOrbitsData::init()
{
	fDat->init();
}


double LCOrbitsData::ArmSpecific(int em, int rec, double trec)
{
	return(0.0);
}


LCVector LCOrbitsData::position(int iSC, double t)
{
	LCVector p(MT);
	//p( 0, fDat->gData(3*(iSC-1)+0, t, InterpType, InterpUtilValue) );
	//p( 1, fDat->gData(3*(iSC-1)+1, t, InterpType, InterpUtilValue) );
	//p( 2, fDat->gData(3*(iSC-1)+2, t, InterpType, InterpUtilValue) );
	
	for(int k=0; k<3; k++)
		p.p[k] = ConvPos * fDat->gData(iColPos[iSC-1][k], (t+TimeShift)/ConvTime+fDat->getx0(), InterpType, InterpUtilValue);
	
	return p;
}


LCVector LCOrbitsData::velocity(int iSC, double t)
{
	LCVector p0(MT), p0dt(MT);
	LCVector v(MT);
	
	
	if(VelocityInFile){
		for(int k=0; k<3; k++)
			v.p[k] = ConvVel * fDat->gData(iColVel[iSC-1][k], (t+TimeShift)/ConvTime+fDat->getx0(), InterpType, InterpUtilValue);
		
	}else{
		/*! *** Compute velocity as derivative of the time serie of position. Sign -1 because we are considering a 
		 * and not a delay as a usual time serie (to be checked).
		 */
		for(int k=0; k<3; k++)
			v.p[k] = - fDat->DerivBackOrder2Spe(k, (t+TimeShift)/ConvTime+fDat->getx0(), fDat->getdx(), InterpType, InterpUtilValue);
		//double dt(1.0);
		//p0   = position(iSC, t);
		//p0dt = position(iSC, t+dt);
		//v = ( p0dt - p0 ) / dt;
	}
	
	return v;
	
}


void LCOrbitsData::DispInfo(char * BTab)
{
	if(MT->Disp()){
		Cout << BTab << "Orbit read in file :" << Endl;
		DispInfoBase(BTab);
		if(fDat != NULL)
			fDat->ControlDisplay();
		Cout << BTab << "Conversion of time = " << ConvTime << " s" << Endl;
		Cout << BTab << "Conversion of position = " << ConvPos << " m" << Endl;
		Cout << BTab << "Conversion of velocity = " << ConvVel << " m/s" << Endl;
		Cout << BTab << "Time shift = " << TimeShift << " s" << Endl;
		for(int iSC=0; iSC<3; iSC++){
			Cout << BTab << "\t + Spacecraft " << iSC+1 << " : position =  ";
			for(int i=0; i<3; i++)
				Cout << iColPos[iSC][i] << " , ";
			if(VelocityInFile){
				Cout << "   velocity = ";
				for(int i=0; i<3; i++)
					Cout << iColVel[iSC][i] << " , ";
			}
			Cout << Endl;
		}
	}
	
}	


// ***********************
// ***  Local mehtods  ***
// ***********************




// ********************
// *  Access methods  *
// ********************

void LCOrbitsData::setFile(char * NewFileName)
{
	if(fDat != NULL)
		delete fDat;
	fDat = new LCDataFileRead(MT);
	fDat->setFile(NewFileName);
}




// ********************
// *  Others methods  *
// ********************




// end of LISACODE-OrbitsData.cpp