// $Id:  $
/*
 *  LISACODE-FIM.h
 *  V 2.0
 *
 *  Created on 07/05/2005 by  Antoine Petiteau (AEI)
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


/*! Compilation including GSL :
 
 */


#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include "LISACODE-Constants.h"
#include "LISACODE-Tools.h"
#include "LISACODE-Matrix.h"
#include "LISACODE-LISACode.h"


/** \ingroup main Main 
 * \{
 */


/*! Load confusion noise model in a serie */ 
void LoadConfusionNoise(LCTools * MT, LCSerie2 * & ConfNoise, int KeyModel, double fMin, double fMax, double df);


/** \brief The executable compute the FIM curve and by product the GW and noise response
 * \author A. Petiteau
 * \version 2.0
 * \date 25/02/2011
 * 
 *
 *
 */
int main (int argc, char * const argv[])
{
	try {
		LCTools MT;
		int nOptions(0);
		char cmd[16384];
		
		long SeedRand((long)time(NULL));
		double dtMes(0.3);
		double t0(0.);
		double tDur(63115200.0);
		int NtData(-1), tmpNtData(-1);
		int NfData(-1);
		
		int NRun(1);
		
		// *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution :" << Endl;
			Coutm << "\t\t(./)LC2FIM [Options] ConfigFile1.xml ConfigFile2.xml ..." << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << Endl << "\t\t * ConfigFileI.xml (required) : Configuration files. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -P %diP1,%diP2  : Index of parameters to study ('L' before the index mean that we have to use the log). Example : 0,1,2,3,L4,6,7 [default : None]"  << Endl ;
			Coutm << "\t\t * -cv %diP %fDeltaMin %fDeltaMax %fDeltaLogStep : Test convergence for parameter iP with value of Delta between DeltaMin and DeltaMax with step log step of DeltaLogStep [default : None]"  << Endl ;
			Coutm << "\t\t * -nf %sPSDFileName : File containing the PSD [required or option -nfc]"  << Endl ;
			Coutm << "\t\t * -nfc %sPSDFileName %sCodeColumn : File containing the PSD. Example : -nfc PSDFile.txt f,n,X,Y,Z [required]"  << Endl ;
			Coutm << "\t\t * -ca %dKeyModel : Analytic model of confusion noise  [default : None]"  << Endl ;
			Coutm << "\t\t * -I %dKey  : Method used for inverting FIM : 0-Gauss 1-SVD [default : 1]"  << Endl ;
			Coutm << "\t\t * -o %sOutFile  : Output file containing results [default: FIMResults.txt]"  << Endl ;
			Coutm << "\t\t * -tdi %sObsName : Names of TDI generators [default : X]"  << Endl ;
			Coutm << "\t\t * -t0 %ft0 : Initial time [default : " << t0 <<"]"  << Endl ;
			Coutm << "\t\t * -T %fT   : Duration [default : " << tDur  << " (= "  << tDur/LC::Yr_SI << " yr)]"  << Endl ;
			Coutm << "\t\t * -pr %diP  : Index of parameters to randomize [default : None]"  << Endl ;
			Coutm << "\t\t * -s %dseed : Seed for random gennerator. [default : current time]"  << Endl ;
			Coutm << "\t\t\t  NB: If not specified the random seed is the time machine." << Endl;
			Coutm << "\t\t * -dn \t\t: No screen display. [default: false]"  << Endl ;
			Coutm << "\t\t * -dl %sfile \t: Write standard output in a file. [default: no file]"  << Endl ;
			Coutm << "\t\t * -v \t\t: Verbose mode : display full details. [default: false]"  << Endl ;
			Coutm << "\t\t * -h \t\t: This help."  << Endl ;
			Coutm << "\t\t * -V \t\t: Version."  << Endl ;
			Coutm << Endl << " ----------------" << Endl;
			return 0;
			
		}
		
		// *********** Version *************
		if(((argc>1)&&(strcmp(argv[1],"--version")==0))||((argc>1)&&(strcmp(argv[1],"-v")==0))){
			Coutm << " ----- VERSION -----" << Endl;
			Coutm << " LC2FIM : Compute the Fisher information matrix (FIM) - LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		
		// ***** Main declaration *****
		LISACode * Sim(NULL);
		
		
		//! *** List of names of files describing the GW configuration
		char ** fNGWCfg(NULL);
		int NfNGWCfg(0);
		
		//! *** List of names of PSD files and index of columns
		char ** fNPSD(NULL);
		int NfNPSD(0);
		char *** TitColfPSD(NULL);
		int * NTitColfPSD(NULL);
		int NCharTitCol(16);
		LCDataFileRead ** fPSD(NULL);
		
		//! *** Confusion noise
		int ModelConfNoise(-1);
		LCSerie2 * ConfNoise(NULL);
		
		//! *** Name of TDI observables
		char TDINameAll[16384];
		char ** TDIName(NULL);
		int NTDI(0);
		strcpy(TDINameAll, "A,E");
		
		// *** List of pointer on read file and index of signal in it. Size NTDI
		LCDataFileRead ** fTDIPSD(NULL);
		int * ifTDIPSD(NULL);
		
		
		//! * Derivative of time data of the TDI signal of GW. Size NTDI x NtData
		double ** tSigDer(NULL);
		//! * Time data of the TDI signal of GW for parameter - delta : \f$ s(t ; \lambda - \Delta \lambda) \f$ . Size NTDI x NtData
		double ** tSigMinD(NULL);
		//! * Derivative of frequency data of the TDI signal of GW. Size : NParamToStudy x NTDI x NfData = NParamToStudy x NTDI x (NtData+1)/2 
		dcomplex *** fSigDer(NULL);
		
		//! *** FIM values
		LCMatrix FIM;
		LCMatrix CovM;
		double df;
		int iFmin, iFmax;
		dcomplex tmpRes;
		int InvFIMMethod(1);
		std::vector<double> ErrPar(0);
		
		//! *** Output (results) file
		char bNOut[16384], fNOutFIM[16384], fNOutErr[16384];
		std::ofstream fOutFIM, fOutErr;
		strcpy(bNOut, "FIMRessults");
		
		//! *** Index of parameters to be chosen randomly
		int * iParamRand(NULL);
		int NParamRand(0);
		
		//! *** Index of parameters to be chosen randomly
		int * iParStudy(NULL);
		int NParStudy(0);
		bool * lnParStudy(NULL);
		double * ParStuValue(NULL);
		double * ParStuDelta(NULL);
		
		//! *** For testing convergence
		double ConvDMin(-1.), ConvDMax(-1.), ConvDLogStep(-1.);
		int ConviP(-1);
		
		//! HERE I HARD CODED THE TOTAL NUMBER OF PARAMETERS : GALBIN
		int NParTot(12);
		
		// ***** Read options *****
		for(int iarg=1; iarg<argc; iarg++){
			
			if((argc>1)&&(strcmp(argv[iarg],"-P")==0)){
				//! **** Read list of parameters to study
				char tmpIndPar[1064];
				int iPtmp(0);
				char ** tmpWordInd(NULL);
				
				strcpy(tmpIndPar, argv[iarg+1]);
				MT.wextractcoma(tmpIndPar, iPtmp, tmpWordInd, NParStudy, 10);
				
				iParStudy = (int*) MT.AllocMemory(NParStudy*sizeof(int));
				lnParStudy = (bool*) MT.AllocMemory(NParStudy*sizeof(bool));
				ParStuValue = (double*) MT.AllocMemory(NParStudy*sizeof(double));
				ParStuDelta = (double*) MT.AllocMemory(NParStudy*sizeof(double));
				for(int i=0; i<NParStudy; i++){
					if(tmpWordInd[i][0] == 'L'){
						iParStudy[i] = atoi(tmpWordInd[i]+1);
						lnParStudy[i] = true;
					}else{
						iParStudy[i] = atoi(tmpWordInd[i]);
						lnParStudy[i] = false;
					}
					ParStuValue[i] = 0.;
					ParStuDelta[i] = 0.;
				}
				
				for(int i=0; i<NParStudy; i++)
					MT.Free(tmpWordInd[i], 10*sizeof(char));
				MT.Free(tmpWordInd, NParStudy*sizeof(char*));
				
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-cv")==0)){
				if(argc-iarg<4)
					throw std::invalid_argument("ERROR : Need 4 parameters for testinf the convergence !");
				ConviP = atoi(argv[iarg+1]);
				ConvDMin = atof(argv[iarg+2]);
				ConvDMax = atof(argv[iarg+3]);
				ConvDLogStep = atof(argv[iarg+4]);
				nOptions += 5;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-nf")==0)){
				//! **** Read name of one PSD noise file in option and allocate memory for list of title
				NfNPSD++;
				
				fNPSD = (char**) MT.ReAllocMemory(fNPSD, (NfNPSD-1)*sizeof(char*), NfNPSD*sizeof(char*));
				fNPSD[NfNPSD-1] = (char*)MT.AllocMemory(2048*sizeof(char));
				strcpy(fNPSD[NfNPSD-1], argv[iarg+1]);
				
				TitColfPSD = (char***) MT.ReAllocMemory(TitColfPSD, (NfNPSD-1)*sizeof(char**), NfNPSD*sizeof(char**));
				NTitColfPSD = (int*) MT.ReAllocMemory(NTitColfPSD, (NfNPSD-1)*sizeof(int), NfNPSD*sizeof(int));
				TitColfPSD[NfNPSD-1] = NULL;
				NTitColfPSD[NfNPSD-1] = 0;
				
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-nfc")==0)){
				
				//! **** Read name of one PSD noise file in option and allocate memory for list of title
				NfNPSD++;
				
				fNPSD = (char**) MT.ReAllocMemory(fNPSD, (NfNPSD-1)*sizeof(char*), NfNPSD*sizeof(char*));
				fNPSD[NfNPSD-1] = (char*)MT.AllocMemory(2048*sizeof(char));
				strcpy(fNPSD[NfNPSD-1], argv[iarg+1]);
				
				TitColfPSD = (char***) MT.ReAllocMemory(TitColfPSD, (NfNPSD-1)*sizeof(char**), NfNPSD*sizeof(char**));
				NTitColfPSD = (int*) MT.ReAllocMemory(NTitColfPSD, (NfNPSD-1)*sizeof(int), NfNPSD*sizeof(int));
				TitColfPSD[NfNPSD-1] = NULL;
				NTitColfPSD[NfNPSD-1] = 0;
				int ipIS(0);
				MT.wextractcoma(argv[iarg+2], ipIS, TitColfPSD[NfNPSD-1], NTitColfPSD[NfNPSD-1], NCharTitCol);
				
				nOptions +=3;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-ca")==0)){
				ModelConfNoise = atof(argv[iarg+1]); 
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-tdi")==0)){
				// ***** Read gravitational waves file in option
				strcpy(TDINameAll, argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-o")==0)){
				// ***** Read base name of output file in option
				strcpy(bNOut, argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-I")==0)){
				InvFIMMethod = atoi(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-t0")==0)){
				t0 = atof(argv[iarg+1]); 
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-T")==0)){
				tDur = atof(argv[iarg+1]); 
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-r")==0)){
				NRun = atoi(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-pr")==0)){
				NParamRand++;
				iParamRand = (int*) MT.ReAllocMemory(iParamRand, (NParamRand-1)*sizeof(int), NParamRand*sizeof(int));
				iParamRand[NParamRand-1] = atof(argv[iarg+1]); 
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-s")==0)){
				SeedRand = atoi(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-dn")==0)){
				MT.unsetDisp();
				nOptions++;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-dl")==0)){
				MT.setDispInFile(argv[iarg+1],true);
				nOptions+=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-V")==0)){
				MT.setDispDetails();
				nOptions++;
			}
		} 
		MT.setRandSeed(SeedRand);
		
		
		Coutm << "FIM study for parameters : " ;
		for(int iP=0; iP<NParStudy; iP++){
			if(lnParStudy[iP]){
				Coutm << " log(" << iParStudy[iP] << ")";
			}else{
				Coutm << " " << iParStudy[iP];
			}
		}
		Coutm << Endl;
		
		//! ***** Read configuration files in arguments
		
		for(int ia=1+nOptions; ia<argc; ia++){
			NfNGWCfg++;
			fNGWCfg = (char**)MT.ReAllocMemory(fNGWCfg, (NfNGWCfg-1)*sizeof(char*), NfNGWCfg*sizeof(char*));
			fNGWCfg[NfNGWCfg-1] = (char*)MT.AllocMemory(2048*sizeof(char));
			strcpy(fNGWCfg[NfNGWCfg-1], argv[ia]);
		}
		
		//! *** Display for checking which files we are using
		Coutm << "Configuration files : " << Endl;
		for(int i=0; i<NfNGWCfg; i++){
			Coutm << "\t + " << fNGWCfg[i] ; MT.CheckFile("", fNGWCfg[i]);
		}
		
		
		
		//! ***** Read the bloc name of TDI generators and obtain the individual names and the number of generators
		int ipTDINameAll(0);
		MT.wextractcoma(TDINameAll, ipTDINameAll, TDIName, NTDI, NCharTitCol);
		Coutm << "TDI generators : ";
		for(int i=0; i<NTDI; i++)
			Coutm << " " << TDIName[i];
		Coutm << Endl;
		
		
		//! ***** Allocate memory 
		tSigDer = (double **) MT.AllocMemory(NTDI*sizeof(double*));
		tSigMinD = (double **) MT.AllocMemory(NTDI*sizeof(double*));
		fSigDer = (dcomplex ***) MT.AllocMemory(NParStudy*sizeof(dcomplex**));
		FIM.init(&MT, NParStudy, NParStudy);
		CovM.init(&MT, NParStudy, NParStudy);
		ErrPar.resize(NParStudy,0.);
		fTDIPSD = (LCDataFileRead **) MT.AllocMemory(NTDI*sizeof(LCDataFileRead*));
		ifTDIPSD = (int *) MT.AllocMemory(NTDI*sizeof(int));
		for(int iP=0; iP<NParStudy; iP++){
			fSigDer[iP] = (dcomplex **) MT.AllocMemory(NTDI*sizeof(dcomplex*));
		}
		
		//! ***** Initialization of memory 
		for(int i=0; i<NTDI; i++){
			tSigDer[i] = NULL;
			tSigMinD[i] = NULL;
			fTDIPSD[i] = NULL;
			ifTDIPSD[i] = 0;
			for(int iP=0; iP<NParStudy; iP++)
				fSigDer[iP][i] = NULL; 
		}
		for(int iP=0; iP<NParStudy; iP++){
			for(int iP2=0; iP2<NParStudy; iP2++){
				FIM(iP,iP2,0.);
				CovM(iP,iP2,0.);
			}
		}
		
		
		//! ***** Read the PSD files
		fPSD = (LCDataFileRead **) MT.AllocMemory(NfNPSD*sizeof(LCDataFileRead*));
		for(int iF=0; iF<NfNPSD; iF++){
			fPSD[iF] = new LCDataFileRead(&MT, fNPSD[iF]);
			fPSD[iF]->init();
			
			//! *** If we don't have the name of columns (not specified in option) try to find them in the header
			if(TitColfPSD[iF] == NULL)
				fPSD[iF]->getHeader(TitColfPSD[iF], NTitColfPSD[iF], NCharTitCol);
			
			//! ** Loop over the TDI generator trying to match the name with the column of the file
			for(int iTDI=0; iTDI<NTDI; iTDI++){
				for(int iC=0; iC<NTitColfPSD[iF]; iC++){
					if(MT.wcmp(TDIName[iTDI], TitColfPSD[iF][iC])){
						fTDIPSD[iTDI] = fPSD[iF];
						ifTDIPSD[iTDI] = iC - 1 + fPSD[iF]->ReadFirstCol();
					}
				}
			}
		}
		
		//! *** Display for checking PSD files
		Coutm << "PSD noise files :" << Endl;
		for(int iTDI=0; iTDI<NTDI; iTDI++){
			if(fTDIPSD[iTDI] != NULL){
				Coutm << Endl << "\t - " << TDIName[iTDI] << " : signal " << ifTDIPSD[iTDI] << " of :";
				fTDIPSD[iTDI]->ControlDisplay();
			}else{
				std::cerr << "ERROR : No PSD of noises for TDI " << TDIName << " !" << Endl;
				throw std::invalid_argument("ERROR : A PSD of noise file is missing for one TDI !");
			}
		}
		
		
		//! *** Load Confusion noise if needed
		if(ModelConfNoise>=0)
			LoadConfusionNoise(&MT, ConfNoise, ModelConfNoise, fTDIPSD[0]->getx0(), fTDIPSD[0]->getxend(), fTDIPSD[0]->getdx());		
		
		
		//! **** Open and initialize the FIM output file
		sprintf(fNOutFIM, "%s-FIM.txt", bNOut);
		fOutFIM.open(fNOutFIM);
		
		//! *** Write header of FIM output file
		fOutFIM << "#i";
		for(int iP=0; iP<NParTot; iP++)
			fOutFIM << " Par" << iP;
		if(ConviP>=0)
			fOutFIM << " D_" << ConviP << " P-D_" << ConviP << " P+D_" << ConviP;
		for(int iP1=0; iP1<NParStudy; iP1++)
			for(int iP2=0; iP2<NParStudy; iP2++)
				fOutFIM << " FIM" << iParStudy[iP1] << iParStudy[iP2];
		for(int iP1=0; iP1<NParStudy; iP1++)
			for(int iP2=0; iP2<NParStudy; iP2++)
				fOutFIM << " Cov" << iParStudy[iP1] << iParStudy[iP2];
		fOutFIM << Endl;
		fOutFIM.precision(12);
		
		
		//! **** Open and initialize the Err output file
		sprintf(fNOutErr, "%s-Err.txt", bNOut);
		fOutErr.open(fNOutErr);
		
		//! *** Write header of Err output file
		fOutErr << "#i";
		for(int iP=0; iP<NParTot; iP++)
			fOutErr << " Par" << iP;
		for(int iP=0; iP<NParStudy; iP++)
			fOutErr << " Err" << iParStudy[iP];
		fOutErr << Endl;
		fOutErr.precision(12);
		
		
		if(ConviP>=0)
			NRun = MT.ifloor((log10(ConvDMax)-log10(ConvDMin))/ConvDLogStep);
		
		for(int iRun=0; iRun<NRun; iRun++){
			
			std::cerr << "Run : " << iRun << " / " << NRun << Endl;
			
			//! ********* Compute the derivative of time signal : loop over the parameters to study
			
			for(int iP=0 ; iP<NParStudy; iP++){
				double DeltaPar;
				
				//! ***** Two steps : one for plus Delta and one for minus Delta
				for(int iD=0; iD<2; iD++){
					
					//! *** Allocation of the simulator
					Sim  = new LISACode(&MT);
					Sim->setTimeInfo(dtMes, t0, tDur);
					
					//! *** Configure the simulator using the configuration files
					for(int i=0; i<NfNGWCfg; i++)
						Sim->config(fNGWCfg[i]);
					
					//! *** Update time informations or check that there are the same as the previous run
					if((iRun==0)&&(iP==0)&&(iD==0)){
						t0 = Sim->gett0();
						dtMes = Sim->getdtMes();
					}else{
						if(!MT.deq(t0,Sim->gett0()))
							throw std::invalid_argument("ERROR : The initial time changes !");
						if(!MT.deq(dtMes,Sim->getdtMes()))
							throw std::invalid_argument("ERROR : The time step changes !");
					}
					
					
					//! *** Choose randomly the parameters to randomize
					if(iRun!=0){
						for(int iP=0; iP<NParamRand; iP++){
							//! *** THIS CONDITION IS FOR VERIFICATION BINARIES (if inclination = 60 deg, then randomize it) SO IT SHOULD BE REMOVED FOR ANY OTHER SOURCE !!! 
							if((iParamRand[iP]!=3)||(((iParamRand[iP]==3))&&(MT.deq(Sim->GWgetParam(0,iParamRand[iP]),60.*M_PI/180.))))
								Sim->GWRandParam(iParamRand[iP]);
						}
					}
					
					
					//! *** Write parameters value in the output files
					if((iP==0)&&(iD==0)){
						fOutFIM << iRun;
						Sim->GWDispParam(&fOutFIM);
						fOutErr << iRun;
						Sim->GWDispParam(&fOutErr);
						
					}
					
					//! *** Here add plus or minus Delta to the parameter and get Delta
					ParStuValue[iP] = Sim->GWgetParam(0, iParStudy[iP]);
					if((iParStudy[iP] <= 0.)&&(lnParStudy[iP])){
						lnParStudy[iP] = false;
						Coutm << "WARNING : The LOG consideration of parameter " << lnParStudy[iP] << " is CANCELLED because the value is equal to 0 or negative !!!" << Endl;
						std::cerr << "WARNING : The LOG consideration of parameter " << lnParStudy[iP] << " is CANCELLED because the valueis equal to 0 or negative !!!" << Endl;
					}
					
					if(ConviP==iParStudy[iP]){
						//! *** Testing convergence
						DeltaPar = pow(10., log10(ConvDMin) + iRun*ConvDLogStep);
						if(iD == 0){
							fOutFIM << " " << DeltaPar;
							Sim->GWAddExtDeltaParam(iParStudy[iP], -DeltaPar);
						}else{
							Sim->GWAddExtDeltaParam(iParStudy[iP], +DeltaPar);
						}
						fOutFIM << " " << Sim->GWgetParam(0, iParStudy[iP]);
					}else{
						//! *** Standard FIM
						if(iD == 0)
							DeltaPar = Sim->GWAddDeltaParam(iParStudy[iP], -1.0);
						else
							DeltaPar = Sim->GWAddDeltaParam(iParStudy[iP], +1.0);
					}
					ParStuDelta[iP] = DeltaPar;
					
					
					
					//! *** Initialization of the simulator
					Sim->init();
					
					
					//! *** Add output series with allocation of the series if needed
					for(int i=0; i<NTDI; i++){
						Sim->AddSerieOut(TDIName[i], 1, tmpNtData, ((iRun==0)&&(iP==0)) );
						if(NtData == -1)
							NtData = tmpNtData; 
						if(NtData != tmpNtData){
							Coutm << "ERROR : The numbers of data in TDI time vector are not the same : previously it was " << NtData << "but it's " << tmpNtData << "for the last one (" << TDIName[i] << ") !" << Endl;
							throw std::invalid_argument("ERROR : The numbers of data in TDI time vector are not the same !");
						}
					}
					
					//! *** Link output series
					if(iD==0){
						for(int i=0; i<NTDI; i++)
							Sim->LinkSerieOut(i, tSigMinD[i]);
					}else{
						for(int i=0; i<NTDI; i++)
							Sim->LinkSerieOut(i, tSigDer[i]);
					}
					
					//Sim->DispInfo("");
					
					
					//! *** Compute the time series : Run the simulator
					Sim->Run();
					
					//! *** Delete the simulator
					delete Sim;
					Sim = NULL;
					
				}
				
				//! *** Compute the centered derivative
				for(int iTDI=0; iTDI<NTDI; iTDI++)
					for(int i=0; i<NtData; i++)
						tSigDer[iTDI][i] = ( tSigDer[iTDI][i] - tSigMinD[iTDI][i] ) / DeltaPar; 
				
				Coutm << "Delta of parameter " << iP << " = " << DeltaPar << Endl;
				
				//! *** Compute the Fourier transform of the derivative
				
				//! ****** Compute the Fourrier transform
				//! *** Prepare computation
				if((iRun==0)&&(iP==0)){
					NfData = MT.getNfFTreal(NtData);
					df = 1./(dtMes*NtData);
					iFmin = MAX(MT.iceil(fTDIPSD[0]->getx0()/df), 0)+2;
					iFmax = MIN(MT.iceil(fTDIPSD[0]->getxend()/df), NfData)-2;
				}
				if(iRun==0)
					for(int i=0; i<NTDI; i++)
						fSigDer[iP][i] = (dcomplex*) MT.AllocMemory(NfData*sizeof(dcomplex));
				
				//! *** Compute Fourrier transform
				MT.FTMakeFwdMulti(tSigDer, fSigDer[iP], NtData, NTDI);
				
				//! *** Normalize Fourrier transform
				for(int iF=0; iF<NfData; iF++)
					for(int iS=0; iS<NTDI; iS++)
						fSigDer[iP][iS][iF] *= dtMes;
				
				Coutm << "Frequency frame : df = " << df << " Hz ,  Fmin = " << iFmin*df << " Hz ,  Fmax = " << iFmax*df << " Hz" << Endl; 
				
				//! *** If we have to consider the logartihm instead of the value itself, we need to multiply the derivative by the value
				if(lnParStudy[iP]){
					for(int iF=0; iF<NfData; iF++)
						for(int iS=0; iS<NTDI; iS++)
							fSigDer[iP][iS][iF] *= ParStuValue[iP];
				}
				
			}
			
			
			//! ****** Computation of Fisher Information Matrix 
			for(int iP1=0; iP1<NParStudy; iP1++){
				for(int iP2=0; iP2<NParStudy; iP2++){
					double tmpFIM(0.);
					for(int iS=0; iS<NTDI; iS++){
						tmpRes = 0.;
						for(int iF=iFmin; iF<iFmax; iF++){
							//if(iS==0)
							//	Coutm << iF << " " << iF*df << " " << real(fSig[iS][iF]) << " " << imag(fSig[iS][iF]) << " " << real(fSig[iS][iF] * conj(fSig[iS][iF])) << " " << fTDIPSD[iS]->gData(ifTDIPSD[iS], iF*df, LIN, 0) << " " << real(tmpRes) << " " << imag(tmpRes) << Endl;
							if(ConfNoise != NULL)
							tmpRes += fSigDer[iP1][iS][iF] * conj(fSigDer[iP2][iS][iF]) / (fTDIPSD[iS]->gData(ifTDIPSD[iS], iF*df, LIN, 0) + ConfNoise->gData(iF*df, LIN, 0));
							else
							tmpRes += fSigDer[iP1][iS][iF] * conj(fSigDer[iP2][iS][iF]) / fTDIPSD[iS]->gData(ifTDIPSD[iS], iF*df, LIN, 0);
						}
						tmpFIM += 4.0 * df * tmpRes.real();
					}
					FIM(iP1,iP2,tmpFIM);
				}
			}
			
			
			//! **** Display Fisher matrix
			Coutm << "============================ Fisher Information Matrix ============================" << Endl;
			Coutm << FIM << Endl;
			
			//! ** Mathematica format 
			Coutm << "FIM = ";
			for(int iP1=0; iP1<NParStudy; iP1++){
				if(iP1==0)
					Coutm << "{{";
				else 
					Coutm << ",{";
				for(int iP2=0; iP2<NParStudy; iP2++){
					if(iP2==0)
						Coutm <<  FIM(iP1,iP2);
					else
						Coutm << "," <<  FIM(iP1,iP2);
				}
				Coutm << "}";
			}
			Coutm << "}" << Endl;
			
			//! **** Record FIM
			for(int iP1=0; iP1<NParStudy; iP1++)
				for(int iP2=0; iP2<NParStudy; iP2++)
					fOutFIM << " " <<  FIM(iP1,iP2);
			
			
			
			
			//! ***** Invert Fisher matrix 
			CovM = FIM.Inv(InvFIMMethod);
			
			Coutm << "============================ Covariance Matrix ============================" << Endl;
			Coutm << CovM << Endl;
			
			
			//! ***** Normalization of covariance matrix
			for(int iP=0; iP<NParStudy; iP++){
				for(int jP=0; jP<NParStudy; jP++){
					if(iP!=jP)
						CovM(iP, jP, CovM(iP,jP)/sqrt(CovM(iP,iP)*CovM(jP,jP)) );
				}
			}
			
			Coutm << "============================ Normalized covariance Matrix ============================" << Endl;
			Coutm << CovM << Endl;
			
			//! **** Record Covariance matrix
			for(int iP1=0; iP1<NParStudy; iP1++)
				for(int iP2=0; iP2<NParStudy; iP2++)
					fOutFIM << " " <<  CovM(iP1,iP2);
			fOutFIM << Endl;
			
			
			//! ***** Errors on parameters
			for(int iP=0; iP<NParStudy; iP++){
				if(lnParStudy[iP])
					ErrPar[iP] = sqrt(CovM(iP,iP))*ParStuValue[iP];
				else
					ErrPar[iP] = sqrt(CovM(iP,iP));
			}
			
			//! *** Display errors
			Coutm << "============================ Errors on parameters ============================" << Endl;
			for(int iP=0; iP<NParStudy; iP++)
				Coutm << "Par " << iParStudy[iP] << " = " << ParStuValue[iP] << " @ " << ErrPar[iP] << Endl; 
			
			//! *** Record errors
			for(int iP=0; iP<NParStudy; iP++)
				fOutErr << " " << ErrPar[iP];
			fOutErr << Endl;
		}
		
		fOutFIM.close();
		fOutErr.close();
		
		
		// ***** Free memory *****
		
		if(Sim != NULL)
			delete Sim;
		Sim = NULL;
		
		
		if(fNGWCfg != NULL){
			for(int iF=0; iF<NfNGWCfg; iF++)
				if(fNGWCfg[iF] != NULL)
					MT.Free(fNGWCfg[iF], 2048*sizeof(char));
			MT.Free(fNGWCfg, NfNGWCfg*sizeof(char*));
		}
		
		if(fNPSD != NULL){
			for(int i=0; i<NfNPSD; i++)
				if(fNPSD[i] != NULL)
					MT.Free(fNPSD[i], 2048*sizeof(char));
			MT.Free(fNPSD, NfNPSD*sizeof(char*));
		}
		
		if(TitColfPSD != NULL){
			for(int i=0; i<NfNPSD; i++){
				if(TitColfPSD[i] != NULL){
					for(int j=0; j<NTitColfPSD[i]; j++)
						MT.Free(TitColfPSD[i][j], NCharTitCol*sizeof(char));
				}
				MT.Free(TitColfPSD[i], NTitColfPSD[i]*sizeof(char*));
			}
			MT.Free(TitColfPSD, NfNPSD*sizeof(char**));
			MT.Free(NTitColfPSD, NfNPSD*sizeof(int));
		}
		
		if(fPSD != NULL){
			for(int i=0; i<NfNPSD; i++)
				if(fPSD[i] != NULL)
					delete fPSD[i];
			MT.Free(fPSD, NfNPSD*sizeof(LCDataFileRead*));
		}
		
		if(ConfNoise != NULL)
			delete ConfNoise;
		
		if(TDIName != NULL){
			for(int i=0; i<NTDI; i++)
				if(TDIName[i] != NULL)
					MT.Free(TDIName[i], 16*sizeof(char));
			MT.Free(TDIName, NTDI*sizeof(char*));
		}
		
		if(fTDIPSD != NULL)
			MT.Free(fTDIPSD, NTDI*sizeof(LCDataFileRead*));
		
		if(ifTDIPSD != NULL)
			MT.Free(ifTDIPSD, NTDI*sizeof(int));
		
		if(tSigDer != NULL){
			for(int i=0; i<NTDI; i++)
				if(tSigDer[i] != NULL)
					MT.Free(tSigDer[i], NtData*sizeof(double));
			MT.Free(tSigDer, NTDI*sizeof(double*));
		}
		
		if(tSigMinD != NULL){
			for(int i=0; i<NTDI; i++)
				if(tSigMinD[i] != NULL)
					MT.Free(tSigMinD[i], NtData*sizeof(double));
			MT.Free(tSigMinD, NTDI*sizeof(double*));
		}
		
		if(fSigDer != NULL){
			for(int iP=0; iP<NParStudy; iP++){
				if(fSigDer[iP] != NULL)
					for(int i=0; i<NTDI; i++)
						if(fSigDer[iP][i] != NULL)
							MT.Free(fSigDer[iP][i], NfData*sizeof(dcomplex));
				MT.Free(fSigDer[iP], NTDI*sizeof(dcomplex*));
			}
			MT.Free(fSigDer, NParStudy*sizeof(dcomplex*));
		}
		
		if(iParamRand != NULL)
			MT.Free(iParamRand,NParamRand*sizeof(int));
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "LC2FIM: error: " << e.what()<<Endl;
		std::cerr << "LC2FIM: abort!" << Endl;
		exit(1);
	}
	return(0);
};


void LoadConfusionNoise(LCTools * MT, LCSerie2 * & ConfNoise, int KeyModel, double fMin, double fMax, double df)
{
	double f(fMin);
	double Res;
	int NDat(MT->iceil((fMax-fMin)/df));
	
	double L_s, phiL, SxG;
	
	Cout << "Loading confusion noise ..."; MT->o->flush();
	
	
	ConfNoise = new LCSerie2(MT, fMin, df, NDat);
	ConfNoise->allocAll();
	
	for(int iF=0; iF<NDat; iF++){
		
		switch (KeyModel) {
			case 0:
				//! *** LISA
				L_s = 5.0e9/LC::c_SI;
				phiL = 2.0*M_PI*L_s*f;
				SxG = 0.;
				
				if((f>1.0e-4)&&(f<=1.0e-3))
					SxG = pow(10.0,-44.62) * pow(f,-2.3);
				if((f>1.0e-3)&&(f<=pow(10.0,-2.7)))
					SxG = pow(10.0,-50.92) * pow(f,-4.4);
				if((f>pow(10.0,-2.7))&&pow(10.0,-2.4))
					SxG = pow(10.0,-62.8) * pow(f,-8.8);
				if((f>pow(10.0,-2.4))&&(f<=1.0e-2))
					SxG = pow(10.0,-89.68) * pow(f,-20.0);
				
				Res = 16.0 * L_s*L_s * (2.0*M_PI*f)*(2.0*M_PI*f) * pow(sin(phiL),2.0) * SxG ;
				break;
				
			case 1:
				//! *** New LISA configurtion 2
				L_s = 1.0e9/LC::c_SI;
				phiL = 2.0*M_PI*L_s*f;
				SxG = 0.;
				
				if((f>4.5e-4)&&(f<=5.3e-4))
					SxG = 1e-13 * pow(f,7.);
				if((f>5.3e-4)&&(f<=2.2e-3))
					SxG = 2.9174e-47 * pow(f,-3.235);
				if((f>2.2e-3)&&(f<=4.e-3))
					SxG = 1.517e-51 * pow(f,-4.85);
				if((f>4.0e-3)&&(f<=5.88e-3))
					SxG = 6.706e-58 * pow(f,-7.5);
				
				Res = L_s*L_s * (2.0*M_PI*f)*(2.0*M_PI*f) * (3./5.) * pow(sin(phiL),2.0) * SxG ;
				break;
				
			default:
				
				throw std::invalid_argument("ERROR : Invalid key number of confusion noise model !");
				break;
		}
		
		ConfNoise->setBinValue(iF,Res);
		
		
		f += df;
	}
	
	Cout << "  --> OK." << Endl;
	
}




/** \}*/

// end of LISACODE-FIM.cpp
