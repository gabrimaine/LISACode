/*
 *  LISACODE-LISACode.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 12/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


#include "LISACODE-LISACode.h"


// *********************
// ***  Constructor  ***
// *********************

LISACode::LISACode()
{
	MT = new LCTools;
	MT->LocTools = true;
	initNULL(false);
}


LISACode::LISACode(LCTools * MT_n)
{
	MT = MT_n;
	initNULL(false);
}



LISACode::~LISACode()
{
	initNULL(true);
}


void LISACode::initNULL(bool CleanMem)
{
	
	if(CleanMem){
		
		if(GWLocal){
			if(GWs != NULL){
				for(int iGW=0; iGW<NGWs; iGW++)
					if(GWs[iGW] != NULL)
						delete GWs[iGW];
				MT->Free(GWs, NGWs*sizeof(LCGW*));
			}
		}
		
		if(LISA!=NULL)
			delete LISA;
		
		if(Orb!=NULL)
			delete Orb;
		
		if(DelaysTDI != NULL){
			for(int i=0; i<NDelaysTDI; i++)
				if(DelaysTDI[i] != NULL)
					delete DelaysTDI[i];
			MT->Free(DelaysTDI, NDelaysTDI*sizeof(LCSerie2*));
		}
		
		if(LastDelays != NULL)
			MT->Free(LastDelays, NDelaysTDI*sizeof(double*));
		
		if(NameDelays != NULL){
			for(int i=0; i<NDelaysTDI; i++)
				if(NameDelays[i] != NULL)
					MT->Free(NameDelays[i], 4*sizeof(char));
			MT->Free(NameDelays, NDelaysTDI*sizeof(LCSerie2*));
		}
		
		if(TDII != NULL){
			for(int i=0; i<NTDII; i++)
				if(TDII[i] != NULL)
					delete TDII[i];
			MT->Free(TDII, NTDII*sizeof(LCTDIInt*));
		}
		
		if(TDIG != NULL){
			for(int i=0; i<NTDIG; i++)
				if(TDIG[i] != NULL)
					delete TDIG[i];
			MT->Free(TDIG, NTDIG*sizeof(LCTDIGen*));
		}
		
		
		if(fOut!=NULL){
			for(int i=0; i<NfOut; i++)
				if(fOut[i] != NULL)
					delete fOut[i];
			MT->Free(fOut, NfOut*sizeof(LCDataFileWrite*));
		}
		
		if(fOutObsNames != NULL){
			for(int iF=0; iF<NfOut; iF++){
				if(fOutObsNames[iF] != NULL){
					for(int i=0; i<NfOutObsNames[iF]; i++)
						if(fOutObsNames[iF][i]!=NULL)
							MT->Free(fOutObsNames[iF][i], 16*sizeof(char));
					MT->Free(fOutObsNames[iF], NfOutObsNames[iF]*sizeof(char*));
				}
			}
			MT->Free(fOutObsNames, NfOut*sizeof(char**));
		}
		
		if(NfOutObsNames != NULL)
			MT->Free(NfOutObsNames, NfOut*sizeof(int));
		
		//! Free list of time vector. Note the free of individual time vector is done outside
		if(sOut!=NULL)
			MT->Free(sOut, NsOut*sizeof(double*));
		
		if(sOutDat!=NULL){
			for(int i=0; i<NsOut; i++)
				if(sOutDat[i] != NULL)
					MT->Free(sOutDat[i], sizeof(double));
			MT->Free(sOutDat, NsOut*sizeof(double*));
		}
		
		if(NDatsOut!=NULL)
			MT->Free(NDatsOut, NsOut*sizeof(int));
		if(sOutInd!=NULL)
			MT->Free(sOutInd, NsOut*sizeof(int));
		
	}

	
	GWs = NULL;
	NGWs = 0;
	LISA = NULL;
	Orb = NULL;
	DelaysTDI = NULL;
	LastDelays = NULL;
	NameDelays = NULL;
	NDelaysTDI = 0;
	TDII = NULL;
	NTDII = 0;
	TDIG = NULL;
	NTDIG = 0;
	
	fOut = NULL;
	NfOut = 0;
	fOutType.resize(0);
	fOutObsNames = NULL;
	NfOutObsNames = NULL;
    GlobalHeader = NULL;
	
	
	sOut = NULL;
	NDatsOut = NULL;
	NsOut = 0;
	sOutType.resize(0);
	sOutDat = NULL;
	sOutInd = NULL;
	
	t0 = 0.;
	dtMes = 1.0;
	tDur = 10000.;
	tStoreDelay = 30.;
	tShiftTDII = 0.;
	tShiftTDIG = 0.;
	NDataTDIG = 10000;
	NDataReqTDIG = 0;
	NDataReqTDII = 0;
	
	GWOptimizeTime = false;
	
	DispProg = true;
}



// ***************************
// *  Configuration methods  *
// ***************************

void LISACode::config(char * fNInXML)
{
	
	if(LISA == NULL)
		LISA = new LCDetector(MT);
	
	//! ****** Check that the file exist
	Cout << ">>>>> Read configuration " << fNInXML << " ..." << Endl;
	std::ifstream fIn;
	fIn.open(fNInXML);
	if(fIn == NULL){
		Cout << "ERROR: LISACode::config : Can not open the file " << fNInXML << " !" << Endl;
		throw std::invalid_argument("ERROR: LISACode::config : Can not open the file.");
	}
	fIn.close();
	fIn.clear();
	
	
	//! ****** Read XML
	ezxml_t tree, section, param, subsection, orbitdata;
	tree = ezxml_parse_file(fNInXML);
	for (section = ezxml_child(tree, "XSIL"); section; section = section->next) {
		
		if(MT->wcmp(ezxml_attr(section, "Type"),"Simulate")){
			
			//! **** Configuration of simulation parameters
			for(param = ezxml_child(section,"Param"); param; param = param->next){
				
				if(MT->wcmp(ezxml_attr(param,"Name"),"TimeOffset"))
					t0 = MT->gXMLTime(param);
				
				if(MT->wcmp(ezxml_attr(param,"Name"),"Cadence"))
					dtMes = MT->gXMLTime(param);
				
				if(MT->wcmp(ezxml_attr(param,"Name"),"Duration"))
					tDur = MT->gXMLTime(param);
				
			}
			
			for (subsection = ezxml_child(section, "XSIL"); subsection; subsection = subsection->next) {
				if(MT->wcmp(ezxml_attr(subsection, "Type"),"Output")){
					configOuput(subsection);
				}
			}
		}
		
		//! **** Configuration of gravitational waves
		if(MT->wcmp(ezxml_attr(section, "Type"),"SourceData")){
			ezxml_t gwdata; 
			for (gwdata = ezxml_child(section, "XSIL"); gwdata; gwdata = gwdata->next) {
				ezxml_t param;
				
				//! * Configure GW read in file
				if(MT->wcmp(ezxml_attr(gwdata,"Type"),"SampledPlaneWave")){
					AddGW("SampledPlaneWave");
					GWs[NGWs-1]->config(gwdata);
				}else{
					if(MT->wcmp(ezxml_attr(gwdata,"Type"),"PlaneWave")){
						
						for(param = ezxml_child(gwdata,"Param"); param; param = param->next){
							//! * Configure GW modeled in the code
							
							char * SourceType(NULL);
							if(MT->wcmp(ezxml_attr(param,"Name"),"SourceType")){
								MT->stripcopy((*param).txt, SourceType);
								if(MT->wcmp(SourceType, "Stochastic")){
									AddGW("Stochastic");
									GWs[NGWs-1]->config(gwdata);
									GWOptimizeTime = true;
								}
								if((MT->wcmp(SourceType, "GalacticBinaryV"))||(MT->wcmp(SourceType, "GalacticBinary"))){
									AddGW("GalacticBinary");
									GWs[NGWs-1]->config(gwdata);
								}
								if((MT->wcmp(SourceType, "SpinBBHHighHarm"))||(MT->wcmp(SourceType, "SpinBBHHH"))){
									AddGW("SpinBBHHighHarm");
									GWs[NGWs-1]->config(gwdata);
								}
								if((MT->wcmp(SourceType, "NRwave"))){
									AddGW("NRwave");
									GWs[NGWs-1]->config(gwdata);
								}
							       


								if((MT->wcmp(SourceType, "CosmicString"))||(MT->wcmp(SourceType, "CosmicString"))){
									AddGW("CosmicString");
									GWs[NGWs-1]->config(gwdata);
								}

							}
							if(SourceType != NULL)
								MT->Free(SourceType, (strlen(SourceType)+1) * sizeof(char));
						}
					}
				}
			}
		}
									
		
		
		//! **** Configuration of orbits
		if(MT->wcmp(ezxml_attr(section, "Type"),"LISAData")){
			for(orbitdata = ezxml_child(section, "XSIL"); orbitdata; orbitdata = orbitdata->next) {
				if(MT->wcmp(ezxml_attr(orbitdata,"Type"),"LISACode_Orbits")){
					if (Orb != NULL)
						throw std::invalid_argument("ERROR in LISACode::config : The orbit have already been created !  If you want to had complement informations to the orbits use an xml bloc with Type='OrbitsInfo' . ");
					Orb = new LCOrbitsAnaLISA(MT);
					Orb->config(orbitdata);
				}
                if(MT->wcmp(ezxml_attr(orbitdata,"Type"),"MLDC_Orbits")){
					if (Orb != NULL)
						throw std::invalid_argument("ERROR in LISACode::config : The orbit have already been created !  If you want to had complement informations to the orbits use an xml bloc with Type='OrbitsInfo' . ");
					Orb = new LCOrbitsAnaMLDC(MT);
					Orb->config(orbitdata);
				}
				if(MT->wcmp(ezxml_attr(orbitdata,"Type"),"OrbitsFile")){
					if (Orb != NULL)
						throw std::invalid_argument("ERROR in LISACode::config : The orbit have already been created ! If you want to had complement informations to the orbits use an xml bloc with Type='OrbitsInfo' . ");
					Orb = new LCOrbitsData(MT);
					Orb->config(orbitdata);
				}
				if(MT->wcmp(ezxml_attr(orbitdata,"Type"),"Octahedron_FirstOrbits")){
					if (Orb != NULL)
						throw std::invalid_argument("ERROR in LISACode::config : The orbit have already been created ! If you want to had complement informations to the orbits use an xml bloc with Type='OrbitsInfo' . ");
					Orb = new LCOrbitsAnaOctahedron(MT);
					Orb->config(orbitdata);
				}
				if(MT->wcmp(ezxml_attr(orbitdata,"Type"),"OrbitsInfo")){
					if (Orb == NULL)
						throw std::invalid_argument("ERROR in LISACode::config : You have to define orbits before reading an 'OrbitsInfo' bloc !");
					Orb->config(orbitdata);
				}
			}
		}
									
		
		//! **** Configuration of noises 
		if(MT->wcmp(ezxml_attr(section, "Type"),"NoiseData"))
			LISA->configNoise(section);
		
		
		
		if(MT->wcmp(ezxml_attr(section, "Type"),"LISACode")){
			
			for (subsection = ezxml_child(section, "XSIL"); subsection; subsection = subsection->next) {
				
				//Cout << ezxml_attr(subsection,"Type") << Endl;
				
				//! *** Configuration of the detector
				if(MT->wcmp(ezxml_attr(subsection,"Type"),"Detector"))
					LISA->configDetector(subsection);
				
				//! *** Configuration of intermediate TDI
				if(MT->wcmp(ezxml_attr(subsection,"Type"),"TDIIntermediate"))
					configTDIIntermediate(subsection);
				
				//! *** Configuration of TDI generator
				if(MT->wcmp(ezxml_attr(subsection,"Type"),"TDIGenerator")){
					NTDIG++;
					TDIG = (LCTDIGen**) MT->ReAllocMemory(TDIG, (NTDIG-1)*sizeof(LCTDIGen*), NTDIG*sizeof(LCTDIGen*));
					TDIG[NTDIG-1] = new LCTDIGen(MT);
					TDIG[NTDIG-1]->config(subsection);
				}
			}
			
		}
		
	}
	ezxml_free(tree);
}



void LISACode::configTDIIntermediate(ezxml_t xmlbloc)
{
	ezxml_t param;
	char * TDIIntType(NULL);
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		if(MT->wcmp(ezxml_attr(param,"Name"),"TDIIntermediateType")){
			MT->stripcopy((*param).txt, TDIIntType);
			if(MT->wcmp(TDIIntType, "Std2002")){
				NTDII++;
				TDII = (LCTDIInt**) MT->ReAllocMemory(TDII, (NTDII-1)*sizeof(LCTDIInt*), NTDII*sizeof(LCTDIInt*));
				TDII[NTDII-1] = new LCTDIIntStd2002(MT);
				TDII[NTDII-1]->config(xmlbloc);
			}
		}
	}
	if(TDIIntType != NULL)
		MT->Free(TDIIntType, (strlen(TDIIntType)+1) * sizeof(char));
}



void LISACode::configOuput(ezxml_t xmlbloc)
{
	ezxml_t param;
	
	bool NoTypeFound(true), NoNameObsFound(true);
	char * TmpFileName(NULL);
	char * TmpFileTypeC(NULL);
	char * TmpNamesObs(NULL);
	
	int TmpfOutType;
	int TmpNfOutObsNames;
	char ** TmpfOutObsNames(NULL);
	
	TypeFile TmpFileType(UNDEFINED);
	bool FirstTimeRef(false);
	
	//! **** The type of output is define by the 5 first characters
	if(strncmp(ezxml_attr(xmlbloc,"Name"),"Phasemeter",3)==0){
		TmpfOutType = 0;
		NoTypeFound = false;
	}
	
	if(strncmp(ezxml_attr(xmlbloc,"Name"),"TDI",3)==0){
		TmpfOutType = 1;
		NoTypeFound = false;
	}
	
	if(strncmp(ezxml_attr(xmlbloc,"Name"),"Delays",3)==0){
		TmpfOutType = 2;
		NoTypeFound = false;
	}
	
	
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		
		//! *** Read output filename 
		if(MT->wcmp(ezxml_attr(param,"Name"),"FileName"))
			MT->stripcopy((*param).txt, TmpFileName);
		
		//! *** Read output file type
		if(MT->wcmp(ezxml_attr(param,"Name"),"FileType")){
			MT->stripcopy((*param).txt, TmpFileTypeC);
			if(MT->wcmp(TmpFileTypeC,"ASCII")){
				TmpFileType = ASCII;
            }
			if(MT->wcmp(TmpFileTypeC,"Binary")){
				TmpFileType = BINARY;
            }
			if(MT->wcmp(TmpFileTypeC,"XML")){
				TmpFileType = XML;
            }
		}
		
		//! *** Read names of obsevables in output
		//! I know that it's not reader friEndly ! I will try to clarify that in a not rush time !
		
		int ipA(0);
		std::vector<char * > TmpObsAll;
		if(MT->wcmp(ezxml_attr(param,"Name"),"Observables")){
			MT->stripcopy((*param).txt, TmpNamesObs);
			
			NoNameObsFound = false;
			TmpNfOutObsNames = 0;
			
			//! ** Look if a time reference is required
			if((TmpNamesObs[0]='t')&&(TmpNamesObs[1]=',')){
				FirstTimeRef = true;
				ipA = 2;
			}
			
			MT->wextractcoma(TmpNamesObs, ipA, TmpfOutObsNames, TmpNfOutObsNames,16);
			
		}
		
	}
	
	
	
	//! **** Check that everything is defined
	
	if(NoTypeFound){
		std::cerr<< "ERROR in LISACode::configOuput : Don't find the type of output." << Endl;
		throw std::invalid_argument("ERROR in LISACode::configOuput : Don't find the type of output.");
	}
	if(NoNameObsFound){
		std::cerr << "ERROR in LISACode::configOuput : Don't find the names of output." << Endl;
		throw std::invalid_argument("ERROR in LISACode::configOuput : Don't find the names of output.");
	}
	
	if(TmpFileName == NULL){
		std::cerr << "ERROR in LISACode::configOuput : Don't find the output file name." << Endl;
		throw std::invalid_argument("ERROR in LISACode::configOuput : Don't find the output file name.");
	}
	
	if(TmpFileType == UNDEFINED){
		std::cerr << "ERROR in LISACode::configOuput : Don't find the type of file." << Endl;
		throw std::invalid_argument("ERROR in LISACode::configOuput : Don't find the type of file.");
	}
	
	// *** Create the output
	AddOutputFile(TmpfOutType, TmpFileName, TmpFileType, TmpfOutObsNames, TmpNfOutObsNames);
	 
	
	if(TmpFileName != NULL)
		MT->Free(TmpFileName, (strlen(TmpFileName)+1) * sizeof(char));
	if(TmpFileTypeC != NULL)
		MT->Free(TmpFileTypeC, (strlen(TmpFileTypeC)+1) * sizeof(char));
	if(TmpNamesObs != NULL)
		MT->Free(TmpNamesObs, (strlen(TmpNamesObs)+1) * sizeof(char));
	if(TmpfOutObsNames != NULL){
		for(int i=0; i<TmpNfOutObsNames; i++)
			MT->Free(TmpfOutObsNames[i], 16 * sizeof(char));
		MT->Free(TmpfOutObsNames, TmpNfOutObsNames*sizeof(char*));
	}
	
	
}


void LISACode::AddOutputFile(int TmpfOutType, char * TmpFileName, TypeFile TmpFileType, char ** TmpfOutObsNames, int TmpNfOutObsNames)
{
	NfOut++;
	
	// *** Add memory for the file
	fOutObsNames = (char***) MT->ReAllocMemory(fOutObsNames, (NfOut-1)*sizeof(char**), NfOut*sizeof(char**));
	NfOutObsNames = (int*) MT->ReAllocMemory(NfOutObsNames, (NfOut-1)*sizeof(int), NfOut*sizeof(int));
	fOut = (LCDataFileWrite **) MT->ReAllocMemory(fOut, (NfOut-1)*sizeof(LCDataFileWrite*), NfOut*sizeof(LCDataFileWrite*));
	
	
	fOutType.push_back(TmpfOutType);
	
	// *** Add observable names
	NfOutObsNames[NfOut-1] = TmpNfOutObsNames;
	fOutObsNames[NfOut-1] = (char**) MT->AllocMemory(NfOutObsNames[NfOut-1]*sizeof(char*));
	for(int iW=0; iW<NfOutObsNames[NfOut-1]; iW++){
		fOutObsNames[NfOut-1][iW] = (char*) MT->AllocMemory(16*sizeof(char));
		strncpy(fOutObsNames[NfOut-1][iW],TmpfOutObsNames[iW],16);
		//int ip(0);
		//while((fOutObsNames[NfOut-1][iW][ip]!=' ')&&(ip<15))
		//	ip++;
		//fOutObsNames[NfOut-1][iW][ip] = '\0';
	}
	
	fOut[NfOut-1] = new LCDataFileWrite(MT, TmpFileName, TmpFileType);
	
}


/*
void LISACode::ReadNameObs(const char * InAll, int & ipIA, char ** & OutName, int & NOutName)
{
	int ipL;
	bool Stop(false);
	while (!Stop) {
		
		//! * Add memory in the list of words
		NOutName++;
		OutName = (char**) MT->ReAllocMemory(OutName, (NOutName-1)*sizeof(char*), NOutName*sizeof(char*));
		//Cout << "NfOutObsNames[" << NfOut-1 << "] = " << NOutName << Endl;
		
		//! * Allocate memory for the word 
		OutName[NOutName-1] = (char*) MT->AllocMemory(16*sizeof(char));
		
		//! * Copy the word
		ipL = 0;
		while ((InAll[ipIA]!=',')&&(InAll[ipIA]!='\0')&&(InAll[ipIA]!='\n')&&(InAll[ipIA]!=' ')&&(ipL<16)){
			OutName[NOutName-1][ipL++] = InAll[ipIA++];  
		}
		//! * Finish the word
		OutName[NOutName-1][ipL] = '\0';
		
		//Cout << "fOutObsNames[" << NfOut-1 << "][" << NOutName-1 << "] = " << OutName[NOutName-1] << Endl;
		
		//! * Stop reding condition
		if(InAll[ipIA]!=',')
			Stop = true;
		else
			ipIA++;
		
		
	}
}
*/

// ***************************************
// * Linking and initialization methods  *
// ***************************************



int LISACode::init()
{
	char TmpC[100];
	double t0Real;
	
	//! ***** Initialization of orbits
	if(Orb == NULL)
		throw std::invalid_argument("ERROR in LISACode::init : The orbit are not defined ! Check that the orbit are defined in the configuration !");
	Orb->init();
	
	tStoreDelay = 3.0*Orb->getNominalArm();
	
	
	//! ***** Compute the times schedule
	
	//! *** Compute time shift between intermediate TDI computation time and TDI generators computation time
	tShiftTDIG = 0.;
	if(TDIG != NULL){
		for(int i=0; i<NTDIG; i++){
			if(TDIG[i] != NULL){
				double tmpV(dtMes*(MT->iceil(TDIG[i]->getNDatInterp())/2.+2));
				if(tmpV > tShiftTDIG)
					tShiftTDIG = tmpV;
			}
		}
	}
	
	//! *** Compute time shift between intermediate TDI computation time and TDI generators computation time
	tShiftTDII = 0.;
	if(TDII != NULL){
		for(int i=0; i<NTDII; i++){
			if(TDII[i] != NULL){
				double tmpV(dtMes*(MT->iceil(TDII[i]->getNDatInterp())/2.+2));
				if(tmpV > tShiftTDII)
					tShiftTDII = tmpV;
			}
		}
	}
	
	//! *** Compute number of data in TDI output
	NDataTDIG = MT->iceil(tDur/dtMes);
	
	//! *** Compute maximal number of combined delays in TDI generators and then the number of data to be computed in advanced in intermediate TDI in order to compute TDI generators
	int NMaxDelayTDIG(0);
	if(TDIG != NULL){
		for(int i=0; i<NTDIG; i++)
			if(TDIG[i] != NULL)
				if(TDIG[i]->getNMaxDelay() > NMaxDelayTDIG)
					NMaxDelayTDIG = TDIG[i]->getNMaxDelay();
		
	}
	NDataReqTDIG = MT->iceil((2*tShiftTDIG+tStoreDelay*NMaxDelayTDIG)/dtMes) + 1;
	
	//! *** Compute maximal number of combined delays in intermediate TDI and then the number of data to be computed in advanced in detector in order to compute intermediate TDI
	int NMaxDelayTDII(0);
	if(TDII != NULL){
		for(int i=0; i<NTDII; i++)
			if(TDII[i] != NULL)
				if(TDII[i]->getNMaxDelay() > NMaxDelayTDII)
					NMaxDelayTDII = TDII[i]->getNMaxDelay();
	}

	NDataReqTDII = MT->iceil((2*tShiftTDII+tStoreDelay*NMaxDelayTDII)/dtMes) + 1 ;
	
	
	if(TDII != NULL)
		t0Real = t0 - NDataReqTDIG*dtMes + tShiftTDIG - NDataReqTDII*dtMes + tShiftTDII;
	else
		t0Real = t0 - NDataReqTDIG*dtMes + tShiftTDIG;
	
	//Cout << "t0Real = " << t0Real << Endl;
	
	//! ***** Initilization of gravitational waves
	for(int iGW=0; iGW<NGWs; iGW++){
		double tDOrbMax(4000.);
		double tGWmin(t0-tDOrbMax);
		if(GWOptimizeTime){
			double tOrbMin(0.), tOrbMax(0.);
			Orb->tGWMinMax(GWs[iGW]->getParam(0), GWs[iGW]->getParam(1), t0Real, tDur, dtMes, tOrbMin, tOrbMax);
			tGWmin = tOrbMax;
			tDOrbMax = tOrbMax - tOrbMin;
			if(MT->DispDet())
				Cout << "Optimize GW time : tsmallest = " << tGWmin << " s ,  maximal time difference in orbits = " <<  tDOrbMax << " s." << Endl;
		}
		GWs[iGW]->setTimeInfo(t0Real, dtMes, tDur, tGWmin, tDOrbMax, 2.*Orb->getNominalArm() );
        int GWstatus(GWs[iGW]->init());
		if(GWstatus){
            Cout << "WARNING : the initialisation of GW source "<<iGW<<" returns the error code "<<GWstatus<<", it seems that there is a problem ... !"  << Endl;
            std::cerr << "WARNING : the initialisation of GW source "<<iGW<<" returns the error code "<<GWstatus<<", it seems that there is a problem ... !"  << Endl;
			return 1;
        }
	}
	
	
	//! *** Initialization of the list of delays
	NDelaysTDI = Orb->getNArm();
	DelaysTDI = (LCSerie2**) MT->AllocMemory(NDelaysTDI*sizeof(LCSerie2*));
	for(int i=0; i<NDelaysTDI; i++)
		DelaysTDI[i] = new LCSerie2(MT, 0., dtMes, NDataReqTDII+NDataReqTDIG);
	LastDelays = (double**) MT->AllocMemory(NDelaysTDI*sizeof(double*));
	for(int i=0; i<NDelaysTDI; i++)
		LastDelays[i] = NULL;
	NameDelays = (char**) MT->AllocMemory(NDelaysTDI*sizeof(char*));
	for(int i=0; i<NDelaysTDI; i++){
		NameDelays[i] = (char*)MT->AllocMemory(4*sizeof(char));
		NameDelays[i][0] = 'D';
		NameDelays[i][3] = '\0';
	}
	if(Orb->getNSC()==3){
		//! ** LISA specific name
		for(int i=0; i<NDelaysTDI; i++){
			sprintf(TmpC,"%d",1 + i%3);
			NameDelays[i][1] = TmpC[0] ;
			NameDelays[i][2] = (i<3?'\0':'p');
		}
	}else{
		//! ** General name : 'D' + index receiver + index emitter + '\0'
		int iArm(0);
		for(int iR=1; iR<=Orb->getNSC(); iR++){
			for(int iE=1; iE<=Orb->getNSC(); iE++){
				if(iR!=iE){
					sprintf(TmpC,"%d%d", iR, iE);
					NameDelays[iArm][1] = TmpC[0];
					NameDelays[iArm][2] = TmpC[1];
					iArm++;
				}
			}
		}			
	}
	if(MT->DispDet())
		for(int i=0; i<NDelaysTDI; i++)
			Cout << "Name of delay " << i << " : " << NameDelays[i] << Endl;
	
	
	//! ***** Initializationn of the detector
	if(LISA == NULL)
		throw std::invalid_argument("ERROR in LISACode::init : The detector are not defined ! Check that the detector is defined in the configuration !");
	if(TDII != NULL)
		LISA->setTimeInfo(t0Real, dtMes, tStoreDelay, NDataReqTDII);
	else
		LISA->setTimeInfo(t0Real, dtMes, tStoreDelay, NDataReqTDIG);
	
	
	//! *** Link the orbits to the detector
	LISA->LinkOrbits(Orb);
	//! *** Link GWs
	LISA->LinkGWs(GWs, NGWs);
	//! *** Initialization of detector
	LISA->init();
	
	//! ***** Initialization of intermediate TDI
	for(int iI=0; iI<NTDII; iI++){
		if(TDII[iI] != NULL){
			TDII[iI]->setTimeInfo(dtMes, NDataReqTDIG, tShiftTDII);
			//! *** Link to the signal
			TDII[iI]->LinkDetector(LISA);
			//! *** Link to the delays
			TDII[iI]->LinkDelays(DelaysTDI, NDelaysTDI);
			TDII[iI]->init();
		}
	}
	
	
	for(int iG=0; iG<NTDIG; iG++){
		TDIG[iG]->setTimeInfo(tShiftTDIG, tShiftTDII+tShiftTDIG);
		if(TDIG[iG] != NULL){
			//! *** Link to the signal
			if(TDII == NULL){
				TDIG[iG]->LinkSigDetector(LISA);
			}else{
				TDIG[iG]->LinkSigTDIInt(TDII, NTDII);
			}
			//! *** Link to the delays
			TDIG[iG]->LinkDelays(DelaysTDI, NDelaysTDI, NameDelays);
			TDIG[iG]->init();
		}
	}
	
	
	//! **** Initialization of output 
    //! **** Start the global header if it exists 
    if (GlobalHeader!=NULL){
        time_t rawtime;
        struct tm * timeinfo;
        char asciitime[128];
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        asctime_r(timeinfo,asciitime);
        
        char * MyName;
        MyName = new char[512];
        getlogin_r(MyName, 512*sizeof(char));
        
        (*GlobalHeader) << "<?xml version=\"1.0\"?>" << Endl;
        (*GlobalHeader) << "<!DOCTYPE XSIL SYSTEM \"http://www.vallis.org/lisa-xml.dtd\">" << Endl;
        (*GlobalHeader) << "<?xml-stylesheet type=\"text/xsl\" href=\"lisa-xml.xsl\"?>" << Endl;
        (*GlobalHeader) << "<XSIL>" << Endl;
        (*GlobalHeader) << LC::xmlind<<"<Param Name=\"Author\">"<< MyName<<"</Param>" << Endl;
        (*GlobalHeader) << LC::xmlind<<"<Param Name=\"GenerationDate\" Type=\"C format\">"<<asciitime<<"</Param>" << Endl;
        (*GlobalHeader) << LC::xmlind<<"<Comment>lisaXML 1.2 [A. Petiteau (based on M. Vallisneri), July 2013]</Comment>" << Endl;
        (*GlobalHeader) << LC::xmlind<<"<XSIL Type=\"TDIData\">"<<Endl;
        //(*GlobalHeader) << "\t\t<XSIL Type=\"LISACode\">"<<Endl;
        //(*GlobalHeader) << "\t\t\t<XSIL Type=\"Outputs\">"<< Endl;
    }
    
    
	//! **** Initialization of output     
	for(int iO=0; iO<NfOut; iO++){
		fOut[iO]->sett0(t0);
		fOut[iO]->setdt(dtMes);
		
		//! ** If output of phasemeter data, link it in the detector
		if(fOutType[iO]==0){
			fOut[iO]->setNDatExpect(NDataTDIG+NDataReqTDIG+NDataReqTDII);
			for(int iN=0; iN<NfOutObsNames[iO]; iN++){
				bool NoFoundOutput(false);
				NoFoundOutput = LISA->LinkMes2Output(fOut[iO]->AddRecord(fOutObsNames[iO][iN]), fOutObsNames[iO][iN]);
				if(NoFoundOutput){
					Cout << "LISACode::init : The phasemeter output corresponding to " << fOutObsNames[iO][iN] << " has not been found !" << Endl;
					throw std::invalid_argument("LISACode::init : The phasemeter output has not been found !");
				}
			}
			
		}
		
		//! ** If output of TDI data, link it in the TDI generator
		if(fOutType[iO]==1){
			fOut[iO]->setNDatExpect(NDataTDIG);
			for(int iN=0; iN<NfOutObsNames[iO]; iN++){
				for(int i=0; i<NTDIG; i++)
					if(TDIG[i]->MatchName(fOutObsNames[iO][iN]))
						TDIG[i]->LinkResult2Output(fOut[iO]->AddRecord(fOutObsNames[iO][iN]));
			}
		}
		
		//! ** If output of armlength data, link it in the detector
		if(fOutType[iO]==2){
			fOut[iO]->setNDatExpect(NDataTDIG+NDataReqTDIG+NDataReqTDII);
			for(int iN=0; iN<NfOutObsNames[iO]; iN++){
				for(int iA=0; iA<NDelaysTDI; iA++){
					//Cout << fOutObsNames[iO][iN] << " =?= " << NameDelays[iA];
					if(MT->wcmp(fOutObsNames[iO][iN], NameDelays[iA])){
						//if((fOutObsNames[iO][iN][0]==NameDelays[iA][0]) && (fOutObsNames[iO][iN][1]==NameDelays[iA][1])) && ( ((strlen(fOutObsNames[iO][iN])<2)&&(NameDelays[iA][2]=' ')) || ( (strlen(fOutObsNames[iO][iN])>2)&&(fOutObsNames[iO][iN][2]==NameDelays[iA][2])) ) {
						LastDelays[iA] = fOut[iO]->AddRecord(fOutObsNames[iO][iN]);
						//Cout << " ==> " << fOutObsNames[iO][iN] << " == " << NameDelays[iA];
					}
					//Cout << Endl;
				}
			}
			
		}
		
		fOut[iO]->init(GlobalHeader,2);
	}

	if (GlobalHeader!=NULL){
        //! *** Close global header
        //(*GlobalHeader) << "\t\t</XSIL>" << Endl;
        (*GlobalHeader) << LC::xmlind<<"</XSIL>" << Endl;
        (*GlobalHeader) << "</XSIL>" << Endl;
    }
        
	return 0;
}


int LISACode::NDatExpected(int TypeObs)
{
	double NDatExpect(0);
	
	//! ** If output of phasemeter data, link it in the detector
	if(TypeObs==0)
		NDatExpect = NDataTDIG + NDataReqTDIG + NDataReqTDII;
	
	//! ** If output of TDI data, link it in the TDI generator
	if(TypeObs==1)
		NDatExpect = NDataTDIG;
	
	//! ** If output of armlength data, link it in the detector
	if(TypeObs==2)
		NDatExpect = NDataTDIG + NDataReqTDIG + NDataReqTDII;
	
	return(NDatExpect);
}


void LISACode::AddSerieOut(char * ObsName, int TypeObs, int & NDatExpect, bool Allocate)
{
	bool NoFoundOut(true);
	
	NDatExpect = 0;
	
	NsOut++;
	
	//! ** Add the data associated to the serie
	sOutDat = (double**) MT->ReAllocMemory(sOutDat, (NsOut-1)*sizeof(double*), NsOut*sizeof(double*));
	sOutDat[NsOut-1] = (double*)MT->AllocMemory(sizeof(double));
	
	
	
	//!  *** Associate a data to the pointer
	
	//! ** If output of phasemeter data, link it in the detector
	if(TypeObs==0){
		NoFoundOut = LISA->LinkMes2Output(sOutDat[NsOut-1], ObsName);
		NDatExpect = NDataTDIG + NDataReqTDIG + NDataReqTDII;
	}
	
	
	//! ** If output of TDI data, link it in the TDI generator
	if(TypeObs==1){
		for(int i=0; i<NTDIG; i++){
			if(TDIG[i]->MatchName(ObsName)){
				TDIG[i]->LinkResult2Output(sOutDat[NsOut-1]);
				NoFoundOut = false;
				NDatExpect = NDataTDIG;
			}
		}
	}
	
	//! ** If output of armlength data, link it in the detector
	if(TypeObs==2){
		NDatExpect = NDataTDIG + NDataReqTDIG + NDataReqTDII;
		for(int iA=0; iA<NDelaysTDI; iA++){
			if(MT->wcmp(ObsName,NameDelays[iA]) && ((NameDelays[iA][2]!=' ')&&(strlen(ObsName)>2) )){
				LastDelays[iA] = sOutDat[NsOut-1];
				Cout << ObsName << " == " << NameDelays[iA] << Endl;
			}
		}
		
	}
	
	if(NoFoundOut){
		Cout << "LISACode::AddLinkSerieOut : The output corresponding to " << ObsName << " has not been found !" << Endl;
		throw std::invalid_argument("LISACode::AddLinkSerieOut : The output has not been found !");
	}
	
	//! ** Add a new series and allocate all the memory for them (avoiding allocation during the run)
	
	NDatsOut = (int*) MT->ReAllocMemory(NDatsOut, (NsOut-1)*sizeof(int), NsOut*sizeof(int));
	NDatsOut[NsOut-1] = NDatExpect;
	sOut = (double**) MT->ReAllocMemory(sOut, (NsOut-1)*sizeof(double*), NsOut*sizeof(double*));
	sOut[NsOut-1] = NULL;
	if(Allocate)
		sOut[NsOut-1] = (double*) MT->AllocMemory(NDatsOut[NsOut-1]*sizeof(double));
	
	sOutType.push_back(TypeObs);
	
	sOutInd = (int*) MT->ReAllocMemory(sOutInd, (NsOut-1)*sizeof(int), NsOut*sizeof(int));
	sOutInd[NsOut-1] = 0;
	
}


void LISACode::LinkSerieOut(int iS, double * & psOut)
{
	if(iS > NsOut)
		throw std::invalid_argument("ERROR in LISACode::LinkSerieOut : The index of the serie to link is higher than the number of series !");
	if(sOut[iS] == NULL)
		sOut[iS] = psOut;
	else
		psOut = sOut[iS];
}


// ********************
// *  Access methods  *
// ********************

void LISACode::setTimeInfo(double dtMes_n, double t0_n, double tDur_n)
{
	dtMes = dtMes_n;
	t0 = t0_n;
	tDur = tDur_n;
}

void LISACode::setGlobalXMLHeader(char * GlobXMLHeadName)
{
    GlobalHeader = new std::ofstream(GlobXMLHeadName);
}


// *********************
// *  Running methods  *
// *********************

void LISACode::Run()
{
	double t, tTDII, tTDIG;
	int iDisp(0);
	
	//char fNCheck[2000];
	//sprintf(fNCheck,"CheckDelay%d.txt", (int)time(NULL));
	//std::cerr << "Record delay in " << fNCheck << Endl;
	//ofstream fCheck(fNCheck);
	//fCheck.precision(12);
	
	tTDIG = t0 ;
	tTDII = tTDIG - NDataReqTDIG*dtMes + tShiftTDIG;
	t = tTDII - NDataReqTDII*dtMes + tShiftTDII ;
	
	//! ******** Phase 1 : The detector run and fill the output buffer of phasemeter measurements
	if((DispProg)&&(MT->Disp()))
		Cout << "Phase 1 : The detector run and fill the output buffer of phasemeter measurements ..." << Endl;
	for(int iT=0; iT<NDataReqTDII; iT++){
		t += dtMes;
		
		//Cout << t << Endl;
		
		//! *** Run one step of the detector 
		LISA->RunStep(t);
		
		//! *** Store delays
		StoreDelaysTDI(t);
		
		//! *** Records output
		RecordOutputs(1, t, tTDII, tTDIG);
	}
	
	
	//! ******** Phase 2 : Phase 1 + computation of intermediate TDI (because we have enough data in phasmeter buffers)
	if(DispProg){
		if(MT->Disp())
			Cout << "Phase 2 : Phase 1 + computation of intermediate TDI (because we have enough data in phasmeter buffers) ..." << Endl;
		std::cerr << "[........10........20........30........40........50........60........70........80........90.......100]" << Endl;
		std::cerr << "["; MT->o->flush();
	}
	iDisp = 0;
	for(int iT=0; iT<NDataReqTDIG; iT++){
		t += dtMes;
		tTDII += dtMes;
		
		iDisp++;
		if(iDisp>NDataReqTDIG/100){
			iDisp = 0;
			if(DispProg){
				std::cerr << "=";
				fflush(stderr);
			}
		}
		
		//! *** Run one step of the detector
		LISA->RunStep(t);
		
		//! *** Store delays
		StoreDelaysTDI(t);
		
		//! *** Run one step of intermediate TDI
		//fcheckTDII << t-tShiftTDII;
		if(TDII != NULL){
			for(int i=0; i<NTDII; i++){
				if(TDII[i] != NULL){
					TDII[i]->RunStep(tTDII);
					//fcheckTDII << " " << TDII[i]->getLastRes();
				}
			}
		}
		//fcheckTDII << Endl;
		
		//! *** Records output
		RecordOutputs(2, t, tTDII, tTDIG);; 
	}
	if(DispProg)
		std::cerr << "=]" << Endl;
	
	
	//! ******** Phase 3 : Phase 2 + computation of TDI (because we have enough data intermedaite data buffer)
	if(DispProg){
		if(MT->Disp())
			Cout << "Phase 3 : Full computation : Phase 2 + computation of TDI (because we have enough data intermediate data buffer) ... " << Endl;
		std::cerr << "[........10........20........30........40........50........60........70........80........90.......100]" << Endl;
		std::cerr << "["; MT->o->flush();
	}
	iDisp = 0;
	//Cout << "TIME t = " << t+dtMes << Endl;
	for(int iT=0; iT<NDataTDIG; iT++){
		t += dtMes;
		tTDII += dtMes;
		tTDIG += dtMes;
		
		
		iDisp++;
		if(iDisp>NDataTDIG/100){
			iDisp = 0;
			if(DispProg){
				std::cerr << "=";
				fflush(stderr);
			}
		}
		
		
		//! *** Run one step of the detector
		LISA->RunStep(t);
		
		//! *** Store delays
		StoreDelaysTDI(t);
		
		//! *** Run one step of intermediate TDI
		//fcheckTDII << t-tShiftTDII;
		if(TDII != NULL){
			for(int i=0; i<NTDII; i++){
				if(TDII[i] != NULL){
					TDII[i]->RunStep(tTDII);
					//fcheckTDII << " " << TDII[i]->getLastRes();
				}
			}
		}
		//fcheckTDII << Endl;
		
		//fCheck << t << " " << tTDIG << " " << DelaysTDI[0]->getBinValue(0) << " " << DelaysTDI[1]->getBinValue(0) << " " << DelaysTDI[2]->getBinValue(0) << Endl;
		
		//! *** Run one step of TDI generator
		if(TDIG != NULL)
			for(int i=0; i<NTDIG; i++)
				if(TDIG[i] != NULL)
					TDIG[i]->RunStep(tTDIG);
		
		
		//! *** Records output
		RecordOutputs(3, t, tTDII, tTDIG);
	}
	if(DispProg)
		std::cerr << "=]" << Endl;
	//fcheckTDII.close();
}




void LISACode::RecordOutputs(int RunPhase, double t, double tTDII, double tTDIG)
{
	/*!	
	 *	TDI output are recorded only in phase 3
	 *	Intermediate TDI are recorded only in phase 2 and 3
	 *	All the others (phasemeter measurements, noises, etc) are recorded only in any phase
	 */
	
	//! *** Record the file output
	for(int i=0; i<NfOut; i++){
		//! *** Record phasemeter output for any phase
		if(fOutType[i]==0)
			fOut[i]->RecordData(t);
		//! *** Record TDI output for phase 3 only
		if((fOutType[i]==1)&&(RunPhase==3))
			fOut[i]->RecordData(tTDIG);
		//! *** Record delays for any phase
		if(fOutType[i]==2)
			fOut[i]->RecordData(t);
	}
	
	//! *** Store the series output
	for(int iS=0; iS<NsOut; iS++){
		//! *** Record phasemeter output for any phase
		if(sOutType[iS]==0)
			sOut[iS][sOutInd[iS]++] = (*sOutDat[iS]) ;
		//! *** Record TDI output for phase 3 only
		if((sOutType[iS]==1)&&(RunPhase==3))
			sOut[iS][sOutInd[iS]++] = (*sOutDat[iS]) ;
		//! *** Record delays for any phase
		if(sOutType[iS]==2)
			sOut[iS][sOutInd[iS]++] = (*sOutDat[iS]) ;
		
        //if((t>1000.)&&(t<1100.)){
        //    Cout << t << " sOutType["<< iS <<"] (v " << sOutType[iS] << ") = " << (*sOutDat[iS]) << Endl;
        //}
        
        
		if(sOutInd[iS] > NDatsOut[iS]){
			Cout << "ERROR in LISACode::RecordOutputs : We pass the end of the serie " << iS << " which is " << NDatsOut[iS] << " !" << Endl;
			throw std::invalid_argument("ERROR in LISACode::RecordOutputs : We pass the end of the serie !");
		}
	}
	
}

void LISACode::StoreDelaysTDI(double trec)
{
	double Error(0.);
	//Error = MT->RandUniform(-2.0e-2, 2.0e-2);
	
	if(Orb->getNSC()==3){
		for(int iD=0; iD<NDelaysTDI; iD++){
			int iEm(iD < 3 ? 1+(iD+2)%3 : 1+(iD-2)%3 );
			int iRe(iD < 3 ? 1+(iD+4)%3 : 1+(iD-1)%3 );
			//Cout << Endl << iD << " : " << iD+1 << " : " << iEm << " -> " << iRe << Endl;  
			DelaysTDI[iD]->addData( Orb->Arm(iEm, iRe, trec) + Error );
			if(LastDelays[iD]!=NULL)
				(*LastDelays[iD]) = DelaysTDI[iD]->getBinValue(0);
		}
	}else{
		int iD(0);
		for(int iRe=1; iRe<=Orb->getNSC(); iRe++){
			for(int iEm=1; iEm<=Orb->getNSC(); iEm++){
				if(iEm!=iRe){
					DelaysTDI[iD]->addData( Orb->Arm(iEm, iRe, trec) + Error );
					if(LastDelays[iD]!=NULL)
						(*LastDelays[iD]) = DelaysTDI[iD]->getBinValue(0);
					iD++;
				}
			}
		}
		
	}
}



// *********************
// *  Methods for GWs  *
// *********************

void LISACode::AddGW(char * GWTypeName)
{
	NGWs++;
	GWs = (LCGW**) MT->ReAllocMemory(GWs, (NGWs-1)*sizeof(LCGW*), NGWs*sizeof(LCGW*));
	GWs[NGWs-1] = NULL;
	
	if(MT->wcmp(GWTypeName, "SampledPlaneWave"))
		GWs[NGWs-1] = new LCGWFile(MT);
	
	if(MT->wcmp(GWTypeName, "Stochastic"))
		GWs[NGWs-1] = new LCGWStochastic(MT);
	
	if(MT->wcmp(GWTypeName, "GalacticBinary"))
		GWs[NGWs-1] = new LCGWGalBin(MT);
	
	if(MT->wcmp(GWTypeName, "SpinBBHHighHarm"))
		GWs[NGWs-1] = new LCGWSpinBBHHHarm1(MT);
	
	if(MT->wcmp(GWTypeName, "NRwave"))
		GWs[NGWs-1] = new LCGWSpinBBHNR1(MT);

	if(MT->wcmp(GWTypeName, "CosmicString"))
        GWs[NGWs-1] = new LCGWCosmicString(MT);
	
	if(GWs[NGWs-1] == NULL){
		Cout << "ERROR in LISACode::AddGW : The type of GW is unknow (defined type : SampledPlaneWave , Stochastic , GalacticBinary , SpinBBHHighHarm ) !" << Endl;
		throw std::invalid_argument("ERROR in LISACode::AddGW : The type of GW is unknow (defined type : SampledPlaneWave , Stochastic , GalacticBinary , SpinBBHHighHarm, CosmicString ) !");
	}
}

void LISACode::GWsetParam(int iGW, int iParam, double ParValue)
{
	if(NGWs<=iGW)
		throw std::invalid_argument("ERROR in LISACode::GWsetParam : No GW of required index ");
	GWs[iGW]->setParam(iParam, ParValue);
}

double LISACode::GWgetParam(int iGW, int iParam)
{
	if(NGWs<=iGW)
		throw std::invalid_argument("ERROR in LISACode::GWgetParam : No GW of required index ");
	return( GWs[iGW]->getParam(iParam) );
}

void LISACode::GWDispParam(std::ostream * out)
{
	for (int i=0; i<NGWs; i++)
		GWs[i]->DispAllParam(out);
}


void LISACode::GWDispParamName(std::ostream * out)
{
	for (int i=0; i<NGWs; i++)
		GWs[i]->DispAllParamName(out);
}


void LISACode::GWRandParam(int iParamRand)
{
	for(int i=0; i<NGWs; i++)
		GWs[i]->RandParam(iParamRand);
}

void LISACode::GWRandAllParam()
{
	for(int i=0; i<NGWs; i++)
		GWs[i]->RandAllParams();
}


double LISACode::GWgetDeltaParam(int iP)
{
	if(NGWs==0)
		throw std::invalid_argument("ERROR in LISACode::GWAddDeltaParam : No GW ");
	return(GWs[0]->getDelta(iP));
}

double LISACode::GWAddDeltaParam(int iP, double FactDelta)
{
	if(NGWs==0)
		throw std::invalid_argument("ERROR in LISACode::GWAddDeltaParam : No GW ");
	return(GWs[0]->AddDeltaPar(iP, FactDelta));
}


void LISACode::GWAddExtDeltaParam(int iP, double ExtDelta)
{
	if(NGWs==0)
		throw std::invalid_argument("ERROR in LISACode::GWAddDeltaParam : No GW ");
	return(GWs[0]->setParam(iP, GWs[0]->getParam(iP)+ExtDelta ));
}

double LISACode::GWgetFreqMin(int iGW)
{
	if(NGWs<=iGW)
		throw std::invalid_argument("ERROR in LISACode::GWsetParam : No GW of required index ");
	return(GWs[iGW]->getFreqMin());
}

double LISACode::GWgetFreqMax(int iGW)
{
	if(NGWs<=iGW)
		throw std::invalid_argument("ERROR in LISACode::GWsetParam : No GW of required index ");
	return(GWs[iGW]->getFreqMax());
}


void LISACode::GWsetSpecialParam(int iGW, int iSParam, double SpeParValue)
{
	if(NGWs<=iGW)
		throw std::invalid_argument("ERROR in LISACode::GWsetParam : No GW of required index ");
	GWs[iGW]->setSpecialParam(iSParam, SpeParValue);
}


LCGW * LISACode::GWget()
{
	if(NGWs>1){
		Cout << "ERROR : LISACode::GWget : There are more than one GW but only the fisrt one can be returned ==> the other ones will be lost !" << Endl;
		throw std::invalid_argument("ERROR : LISACode::GWget : There are more than one GW but only the fisrt one can be returned ==> the other ones will be lost !");
	}
	GWLocal = false;
	return(GWs[0]);
}


void LISACode::GWLinkExt(LCGW * ExtGW)
{
	NGWs++;
	GWs = (LCGW**) MT->ReAllocMemory(GWs, (NGWs-1)*sizeof(LCGW*), NGWs*sizeof(LCGW*));
	GWs[NGWs-1] = NULL;
	GWs[NGWs-1] = ExtGW;
	GWLocal = false;
}


// ********************
// *  Others methods  *
// ********************

void LISACode::DispInfo(char * BTab)
{
	if(MT->Disp()){
		MT->o->flush();
		Cout << BTab << "Time offset = " << t0 << " s" << Endl;
		Cout << BTab << "Time step measurement = " << dtMes << " s" << Endl;
		Cout << BTab << "Duration = " << tDur << " s" << Endl;
		Cout << BTab << "Duration of the storage for one delay = " << tStoreDelay << " s" << Endl;
		Cout << BTab << "Time shift between measurements and TDI intermediate = " << tShiftTDII << " s" << Endl;
		Cout << BTab << "Time shift between TDI intermediate and TDI = " << tShiftTDIG << " s" << Endl;
		Cout << BTab << "Number of data to be computed in TDI generator = " << NDataTDIG << Endl;
		Cout << BTab << "Number of data required for TDI generator (to be computed in advanced in TDI intermediate) = " << NDataReqTDIG << Endl;
		Cout << BTab << "Number of data required for intermediate TDI (to be computed in advanced in the detector) = " << NDataReqTDII << Endl;
		
		if(Orb != NULL)
			Orb->DispInfo(BTab);
		if(LISA != NULL)
			LISA->DispInfo(BTab);
		
		for(int i=0; i<NTDII; i++)
			if(TDII[i] != NULL)
				TDII[i]->DispInfo(BTab);
		
		for(int i=0; i<NTDIG; i++)
			if(TDIG[i] != NULL)
				TDIG[i]->DispInfo(BTab);
		
		for(int iF=0; iF<NfOut; iF++){
			Cout << BTab << "Output " << iF << " :" << Endl;
			Cout << BTab << "\t- Type = " << fOutType[iF] << Endl;
			Cout << BTab << "\t- Observables :";  
			if(fOutObsNames[iF]!=NULL){
				for(int i=0; i<NfOutObsNames[iF]; i++)
					if(fOutObsNames[iF][i] != NULL)
						Cout << " " << fOutObsNames[iF][i];
			}
			Cout << Endl;
			fOut[iF]->DispInfo(BTab);
		}
		
		for(int iS=0; iS<NsOut; iS++){
			Cout << BTab << "Output " << iS << " :" << Endl;
			Cout << BTab << "\t- Type = " << sOutType[iS] << Endl;
		}
		
		for(int iGW=0; iGW<NGWs; iGW++){
			if(GWs[iGW] != NULL){
				Cout << BTab << "GW source " << iGW << " :" << Endl;
				GWs[iGW]->DispInfo(BTab);
			}
		}
	}
}





// end of LISACODE-LISACode.cpp

