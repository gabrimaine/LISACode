/*
 *  LISACODE-OBPMStd2002.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 13/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


#include "LISACODE-OBPMStd2002.h"


// *********************
// ***  Constructor  ***
// *********************

LCOBPMStd2002::LCOBPMStd2002()
: LCOBPM()
{
	initNULL(false);	
}


LCOBPMStd2002::LCOBPMStd2002(LCTools * MT_n)
: LCOBPM(MT_n)
{	
	initNULL(false);
}



LCOBPMStd2002::~LCOBPMStd2002()
{
	initNULL(true);
}




// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCOBPMStd2002::initNULL(bool CleanMem)
{
	initNULLBase(CleanMem);
	
	TypeOBPM = 0;
}


void LCOBPMStd2002::config(ezxml_t xmlbloc)
{
	ezxml_t param;
	char * OpticalBenchType(NULL);
	
	
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		
		//! **** Checking that it's the right configuration of optical bench + phasemeter
		if(MT->wcmp(ezxml_attr(param,"Name"),"OpticalBenchType")){
			MT->stripcopy((*param).txt, OpticalBenchType);
			if(!MT->wcmp(OpticalBenchType, "Std2002"))
				throw std::invalid_argument("ERROR in LCOBPMStd2002::config : The xml bloc does not correspond to a type of the optical bench + phasemeter 'Std2002' .");
		}
		
	}
	if(OpticalBenchType != NULL)
		MT->Free(OpticalBenchType, (strlen(OpticalBenchType)+1) * sizeof(char) );
	
	//! *** Set the name and the position
	setNamePos(ezxml_attr(xmlbloc,"Name"));
	
	if((Name[0]=='s')&&(Name[1]=='c')&&(Name[2]=='i')){
		TypeOBPM = 0;		
		if(IndirectDir){
			iSCLink = (iSC+1)%3 + 1;
			IndirectDirLink = 0;
		}else{
			iSCLink = (iSC)%3 + 1;
			IndirectDirLink = 1;
		}
	}
	
	if((Name[0]=='t')&&(Name[1]=='a')&&(Name[2]=='u')){
		TypeOBPM = 1;
		iSCLink = iSC;
		if(IndirectDir)
			IndirectDirLink = 0;
		else
			IndirectDirLink = 1;
	}
	
	
	
	
	//! *** Configuration of filter
	configFilter(xmlbloc);
	
}


void LCOBPMStd2002::config(int iParam, double ParamVal)
{
	
}


void LCOBPMStd2002::LinkNoise(LCNoise ** AllNoises, int NAllNoises)
{
	//! *** For 'sci' optical bench path + phasemeter : 
	if(TypeOBPM == 0){
		
		//! *** Initialization of noises if there not yet allocated
		if(pn == NULL){
			Npn = 5;
			pn = (LCNoise**) MT->AllocMemory(Npn*sizeof(LCNoise*));
			for (int i=0; i<Npn; i++)
				pn[i] = NULL;
		}
		
		/*! *** Loop on all noises and place them using their name
		 *  List of noises for 'sci' phasemeter 
		 *	pn[0] : \f$ s_i^{SN} \f$		: Shot noise [sn..]
		 *	pn[1] : \f$ s_i^{OOPN} \f$		: Other optical path noises [op..]
		 *	pn[2] : \f$ \delta_i^{DRS} \f$	: Acceleraion noise [pm..]
		 *	pn[3] : \f$ p_i \f$				: local laser noise [c..]
		 *	pn[4] for direct OB		: \f$ p_{i+1} \f$	: distant laser noise [c..]
		 *	pn[4] for indirect OB	: \f$ p_{i+2} \f$	: distant laser noise [c..]
		 */
		for(int i=0; i<NAllNoises; i++){
			//! ** Local optical bench
			if( (AllNoises[i]->getiSC()==iSC) && (AllNoises[i]->getIndirectDir()==IndirectDir) ){
				if(AllNoises[i]->CmpName("sn",2))
					pn[0] = AllNoises[i];
				if(AllNoises[i]->CmpName("op",2))
					pn[1] = AllNoises[i];
				if(AllNoises[i]->CmpName("pm",2))
					pn[2] = AllNoises[i];
				if(AllNoises[i]->CmpName("c",1))
					pn[3] = AllNoises[i];
			}
			//! ** Linked optical bench
			if( (AllNoises[i]->getIndirectDir()==IndirectDirLink) && (AllNoises[i]->getiSC()==iSCLink) ){
				if(AllNoises[i]->CmpName("c",1))
					pn[4] = AllNoises[i];
			}
		}
	}
	
	//! *** For 'tau' optical bench path + phasemeter : 
	if(TypeOBPM == 1){
		
		//! *** Initialization of noises if there not yet allocated
		if(pn == NULL){
			Npn = 3;
			pn = (LCNoise**) MT->AllocMemory(Npn*sizeof(LCNoise*));
			for (int i=0; i<Npn; i++)
				pn[i] = NULL;
		}
		
		/*! *** Loop on all noises and place them using their name
		 *  List of noises for 'sci' phasemeter 
		 *	pn[0] : \f$ p_i \f$					: local laser noise [c..]
		 *	pn[1] : \f$ p_{i+1} \f$				: join laser noise [c..]
		 *	pn[2] : \f$ \delta_{i+1}^{DRS} \f$	: join acceleraion noise [pm..]
		 */
		for(int i=0; i<NAllNoises; i++){
			//! ** Local optical bench
			if( (AllNoises[i]->getiSC()==iSC) && (AllNoises[i]->getIndirectDir()==IndirectDir) ){
				if(AllNoises[i]->CmpName("c",1))
					pn[0] = AllNoises[i];
			}
			//! ** Join optical bench
			if( (AllNoises[i]->getIndirectDir()==IndirectDirLink) && (AllNoises[i]->getiSC()==iSCLink) ){
				if(AllNoises[i]->CmpName("c",1))
					pn[1] = AllNoises[i];
				if(AllNoises[i]->CmpName("pm",2))
					pn[2] = AllNoises[i];
			}
		}
	}
}


void LCOBPMStd2002::LinkFactShotNoise(double * FactShotNoise_n)
{
	if(TypeOBPM == 0)
		FactShotNoise = FactShotNoise_n;
	if(TypeOBPM == 1)
		FactShotNoise = NULL;
}


void LCOBPMStd2002::init()
{
	initBase();
	
	//char tmpN[8000];
	//sprintf(tmpN, "Check_%s.txt", Name);
	//fCheck.open(tmpN);
}


double LCOBPMStd2002::MeasurePho(double TimeLocal, double TShiftEff)
{
	double TmpBeamExt(0.), TmpBeamLoc(0.);
	double  ResMes(0.);
	
	if(TypeOBPM == 0){
		//! *** Computing the GWs signal here is more precise but slower when the physical time step is small
		GWSignal = sGW(TimeLocal);
		//fCheck << TimeLocal << " " << GWSignal << Endl;
		//Cout << Name << " : t = " << TimeLocal << " GWSignal = " << GWSignal << Endl;
		
		/*! ** Compute the incoming beam path trought the optical bench :
		 *	\f$ b_{ext}(t_{USO}) = D_j(t_{USO}) p_{i'} = p_{i'} (t_{USO} -L_j (t_{USO}) ) \f$ 
		 */
		TmpBeamExt = GWSignal + gN(4, TShiftEff - orb->Arm(iSCLink, iSC, TimeLocal) );
		/*! ** Compute the local beam path trought the optical bench :
		 * \f$ b_{loc}(t_{USO}) =  s_i^{SN}(t_{USO}) + s_i^{OOPN}(t_{USO}) - 2 \delta_i^{DRS}(t_{USO}) + p_i(t_{USO})  \f$
		 */
		TmpBeamLoc = gFSN(TimeLocal)*gN(0, TShiftEff) + gN(1, TShiftEff) - 2.* gN(2, TShiftEff) + gN(3, TShiftEff)	;
		
	}
	
	if(TypeOBPM == 1){
		
		/*! ** Compute the incoming beam path trought the optical bench :
		 *	\f$ b_{loc}(t_{USO}) =  p_{i+1} (t_{USO}) + 2 \delta_{i+1}^{DRS} (t_{USO}) \f$ 
		 */
		TmpBeamExt = gN(1, TShiftEff) + 2. * gN(2, TShiftEff) ;
		/*! ** Compute the local beam path trought the optical bench :
		 * \f$ b_{loc}(t_{USO}) =  p_i(t_{USO})  \f$
		 */
		TmpBeamLoc = gN(0, TShiftEff) 	;
	}
	
	//! ** Interference : \f$ s_i = b_{ext}(t_{USO}) - b_{loc}(t_{USO})  \f$
	ResMes = TmpBeamExt - TmpBeamLoc ;
	//fCheck << TimeLocal << " "  << ResMes << " " << TmpBeamExt << " " << TmpBeamLoc << Endl;
	
	return(ResMes);
}



void LCOBPMStd2002::DispInfo(char * BTab, bool LinkedModule)
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



// end of LISACODE-OBPMStd2002.cpp

