/*
 *  LISACODE-OBPMOcta1.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 13/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


#include "LISACODE-OBPMOcta1.h"


// *********************
// ***  Constructor  ***
// *********************

LCOBPMOcta1::LCOBPMOcta1()
: LCOBPM()
{
	initNULL(false);	
}


LCOBPMOcta1::LCOBPMOcta1(LCTools * MT_n)
: LCOBPM(MT_n)
{	
	initNULL(false);
}



LCOBPMOcta1::~LCOBPMOcta1()
{
	initNULL(true);
}




// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCOBPMOcta1::initNULL(bool CleanMem)
{
	initNULLBase(CleanMem);
	
	TypeOBPM = 0;
}


void LCOBPMOcta1::config(ezxml_t xmlbloc)
{
	ezxml_t param;
	char * OpticalBenchType(NULL);
	
	
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		
		//! **** Checking that it's the right configuration of optical bench + phasemeter
		if(MT->wcmp(ezxml_attr(param,"Name"),"OpticalBenchType")){
			MT->stripcopy((*param).txt, OpticalBenchType);
			if(!MT->wcmp(OpticalBenchType, "Octa1"))
				throw std::invalid_argument("ERROR in LCOBPMOcta1::config : The xml bloc does not correspond to a type of the optical bench + phasemeter 'Octa1' .");
		}
		
	}
	if(OpticalBenchType != NULL)
		MT->Free(OpticalBenchType, (strlen(OpticalBenchType)+1) * sizeof(char) );
	
	//! *** Set the name and the position
	setNamePos(ezxml_attr(xmlbloc,"Name"));
	
	
	if((Name[0]=='s')&&(Name[1]=='c')&&(Name[2]=='i')){
		TypeOBPM = 0;
		IndirectDirLink = 0;
		iSCLink = -1*IndirectDir;
	}
	
	
	
	
	//! *** Configuration of filter
	configFilter(xmlbloc);
	
}


void LCOBPMOcta1::config(int iParam, double ParamVal)
{
	
}


void LCOBPMOcta1::LinkNoise(LCNoise ** AllNoises, int NAllNoises)
{
	//! ***** For 'sci' optical bench path + phasemeter : 
	if(TypeOBPM == 0){
		
		//! *** Initialization of noises if there not yet allocated
		if(pn == NULL){
			Npn = 9;
			pn = (LCNoise**) MT->AllocMemory(Npn*sizeof(LCNoise*));
			for (int i=0; i<Npn; i++)
				pn[i] = NULL;
		}
		
		/*! **** Loop on all noises and place them using their name. Emitter I and receiver J
		 *  List of noises for 'sci' phasemeter 
		 *	pn[0] : \f$ s_{IJ}^{SN}		\f$	: shot noise [sn..]
		 *	pn[1] : \f$ p_J				\f$	: local laser noise [c..]
		 *	pn[2] : \f$ ax_J^{DRS}		\f$	: local acceleration noise along x [ax..]
		 *	pn[3] : \f$ ay_J^{DRS}		\f$	: local acceleration noise along y [ay..]
		 *	pn[4] : \f$ az_J^{DRS}		\f$	: local acceleration noise along z [az..]
		 *	pn[5] : \f$ p_I				\f$	: distant laser noise [c..]
		 *	pn[6] : \f$ ax_I^{DRS}		\f$	: distant acceleration noise along x [ax..]
		 *	pn[7] : \f$ ay_I^{DRS}		\f$	: distant acceleration noise along y [ay..]
		 *	pn[8] : \f$ az_I^{DRS}		\f$	: distant acceleration noise along z [az..]
		 */
		for(int i=0; i<NAllNoises; i++){
			//! *** Local noise 
			if(AllNoises[i]->getiSC()==iSC){
			   if(AllNoises[i]->CmpName("c",1))
					pn[1] = AllNoises[i];
				if(AllNoises[i]->CmpName("ax",2))
					pn[2] = AllNoises[i];
				if(AllNoises[i]->CmpName("ay",2))
					pn[3] = AllNoises[i];
				if(AllNoises[i]->CmpName("az",2))
					pn[4] = AllNoises[i];
				//! ** Local noise of the optical bench
				if(AllNoises[i]->getIndirectDir()==IndirectDir){
					if(AllNoises[i]->CmpName("sn",2))
						pn[0] = AllNoises[i];
				}
				
			}
			//! ** Linked spacecraft
			if(AllNoises[i]->getiSC()==iSCLink){
				if(AllNoises[i]->CmpName("c",1))
					pn[5] = AllNoises[i];
				if(AllNoises[i]->CmpName("ax",2))
					pn[6] = AllNoises[i];
				if(AllNoises[i]->CmpName("ay",2))
					pn[7] = AllNoises[i];
				if(AllNoises[i]->CmpName("az",2))
					pn[8] = AllNoises[i];
			}
		}
	}
}


void LCOBPMOcta1::LinkFactShotNoise(double * FactShotNoise_n)
{
	if(TypeOBPM == 0)
		FactShotNoise = FactShotNoise_n;
	if(TypeOBPM == 1)
		FactShotNoise = NULL;
}


void LCOBPMOcta1::init()
{
	initBase();
	
	
	
	//char tmpN[8000];
	//sprintf(tmpN, "Check_%s.txt", Name);
	//fCheck.open(tmpN);
}


double LCOBPMOcta1::MeasurePho(double TimeLocal, double TShiftEff)
{
	double TmpBeamExt, TmpBeamLoc;
	double TShiftDist;
	double  ResMes;
	LCVector n(MT);
	
	n = ( orb->position(iSC,TimeLocal) - orb->position(iSCLink, TimeLocal) ).unit(); 
	
	TShiftDist = TShiftEff - orb->Arm(iSCLink, iSC, TimeLocal);
	
	//if(TimeLocal>10.)
	//	Cout << Name << " : delay = " << TShiftDist << Endl;
	
	if(TypeOBPM == 0){
		//! *** Computing the GWs signal here is more precise but slower when the physical time step is small
		GWSignal = sGW(TimeLocal);
		//fCheck << TimeLocal << " " << GWSignal << Endl;
		//Cout << Name << " : t = " << TimeLocal << " GWSignal = " << GWSignal << Endl;
		
		/*! ** Compute the incoming beam path throught the distant optical bench and local optical bench :
		 *	\f$ b_{ext}(t_{USO}) = D_{IJ}(t_{USO}) (p_{I} + n_{IJ,x} a_{I,x} + n_{IJ,y} a_{I,y} + n_{IJ,z} a_{I,z})  - n_{IJ,x} a_{J,x} - n_{IJ,y} a_{J,y} - n_{IJ,z} a_{J,z} \f$ 
		 */
		TmpBeamExt = GWSignal ;
		TmpBeamExt += gN(5, TShiftDist );
		TmpBeamExt += n.x() * gN(6, TShiftDist );
		TmpBeamExt += n.y() * gN(7, TShiftDist );
		TmpBeamExt += n.z() * gN(8, TShiftDist );
		TmpBeamExt -= n.x() * gN(2, TShiftEff );
		TmpBeamExt -= n.y() * gN(3, TShiftEff );
		TmpBeamExt -= n.z() * gN(4, TShiftEff );
		
		/*! ** Compute the local beam path throught the optical bench :
		 * \f$ b_{loc}(t_{USO}) = b_{J}(t_{USO}) =  s_J^{SN}(t_{USO}) + p_J(t_{USO})  \f$
		 */
		//Cout << gFSN(TimeLocal) << Endl;
		TmpBeamLoc = gFSN(TimeLocal)*gN(0, TShiftEff) + gN(1, TShiftEff) 	;
		
	}
		
	//! ** Interference : \f$ s_i = b_{ext}(t_{USO}) - b_{loc}(t_{USO})  \f$
	ResMes = TmpBeamExt - TmpBeamLoc ;
	//fCheck << TimeLocal << " "  << ResMes << " " << TmpBeamExt << " " << TmpBeamLoc << Endl;
	
	return(ResMes);
}



void LCOBPMOcta1::DispInfo(char * BTab, bool LinkedModule)
{
	if(MT->Disp()){
		if(MT->DispDet()) 
			Cout << Endl << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << Endl;
		Cout << BTab << "Measure (optical bench path + photodiode + phasemeter) : " << Name << Endl;
		if(MT->DispDet()) 
			DispInfoBase(BTab,LinkedModule);
		
		if(MT->DispDet()) 
			Cout << Endl << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << Endl;
	}
}


// ***********************
// ***  Local mehtods  ***
// ***********************



// end of LISACODE-OBPMOcta1.cpp

