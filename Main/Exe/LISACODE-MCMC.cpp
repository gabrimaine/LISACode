// $Id:  $
/*
 *  LISACODE- MCMC.cpp
 *  V 2.0
 *
 *  Created on 02/06/2005 by  Antoine Petiteau (AEI)
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
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
#include "LISACODE-ModSig.h"



/** \ingroup main Main 
 * \{
 */




/** \brief The executable compute the SNR curve and by product the GW and noise response
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
		
		int Nsrc(1);
		int iNStart(0);
		int iRealiStudy(-1);
		int NRand(1);
		
		LCModSig X(&MT), Xp(&MT);
		
		X.Detector = LISA;
		X.GalacticNoise = false;
		X.sSrcType(SPINBBHHHARM);
		
		double SNRth(6.);
		char CharIndPar[1064];
		strcpy(CharIndPar, "None");
		int MCMCNStepMax(100000);
		double MCMCscale(1.);
		double MCMCheat(2.);
		int MCMCFIMRecomp(-1);
		
		
		//! *** Name of TDI observables
		char TDINameAll[16384];
		strcpy(TDINameAll, "X");
		//strcpy(TDINameAll, "Am,Em");
		
		
		// *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution :" << Endl;
			Coutm << "\t\t(./)LC2 MCMC [Options] ConfigFile1.xml ConfigFile2.xml ..." << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << Endl << "\t\t * ConfigFileI.xml (required) : Configuration fileX. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -o %sOutFile  : Output file containing results [default: SNRResultX.txt]"  << Endl ;
			
			Coutm << "\t\t * -G key  : Type of GW source. [default : MBHBInsPrecHH ] "  << Endl ;
			Coutm << "\t\t\t GB or GalBin        : galactic binary,"  << Endl ;
			Coutm << "\t\t\t SH or MBHBInsPrecHH : Massive black hole binary inspiral phase with spin precession and higher harmonicX."  << Endl ;
			Coutm << "\t\t * -C %sCatalog : Catalogue containing the list of parameters : if not random source. [default : None] "  << Endl ;
			Coutm << "\t\t * -N %dNSrc    : Number of sources to study (-1 : number of source in catalogue). [default : " << Nsrc << " ]"  << Endl ;
			Coutm << "\t\t * -f %diFSrc   : Fisrt source to study is iFSrc*NSrc. [default : "  << iNStart << "]"  << Endl ;
			Coutm << "\t\t * -R %diReali  : Study only the source from the realization %diReali (-1:not used). [default : "<< iRealiStudy << "]"  << Endl ;
			//Coutm << "\t\t * -r %dNrand   : Make Nrand randomization of free paramters for each source. [default : "<< NRand << "]"  << Endl ;
			
			Coutm << "\t\t * -dt %step : Minimal time step. [default : " << X.dtMesMin << "]"  << Endl ;
			
			Coutm << "\t\t * -D %dConfig : For computing analytic noise : LISA,C1,C2,C3,C4,C5 [default : LISA]"  << Endl ;
			Coutm << "\t\t * -cn         : Include analytic confusion noise  [default : None]"  << Endl ;
			
			
			Coutm << "\t\t * -tdi %sObsName : Names of TDI generators [default : " << TDINameAll << " ]"  << Endl ;
			//Coutm << "\t\t * -S %fSNRth : SNR threshold for selected for FIM compute [default : " << SNRth << "]"  << Endl ;
			Coutm << "\t\t * -P %diP1,%diP2  : Index of parameters to study ('L' before the index mean that we have to use the log). [default : 0,1,14,13,12,4,L5,L15,17,6,7,8,9,10,11 if src MBHBInsPrecHH , 0,1,2,3,L4,6,7 if src GalBin ]"  << Endl ;
			//Coutm << "\t\t * -cv %diP %fDeltaMin %fDeltaMax %fDeltaLogStep : Test convergence for parameter iP with value of Delta between DeltaMin and DeltaMax with step log step of DeltaLogStep [default : None]"  << Endl ;
			Coutm << "\t\t * -M %fMax %sDir : Maximum memory in RAM and directory for swaping [default : Max = " << X.gMaxMemRAM() << ", Dir = " << X.gDirSwapMem() << " ]"  << Endl ;
			
			Coutm << "\t\t * -n %dN : Maximal number of MCMC steps [default : " << MCMCNStepMax << "]"  << Endl ;
			Coutm << "\t\t * -E %fscale : Scale used for MCMC steps [default : " << MCMCscale << "]"  << Endl ;
			Coutm << "\t\t * -F %dFreqFIM : Frequency of FIM recomputation : recompute each %dFreqFIM MCMC accepted steps [default : " << MCMCFIMRecomp << "]"  << Endl ;
			
			
			Coutm << "\t\t * -co       : Write the check output fileX. [default : false]"  << Endl ;
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
		if(((argc>1)&&(strcmp(argv[1],"--version")==0))||((argc>1)&&(strcmp(argv[1],"-V")==0))){
			Coutm << " ----- VERSION -----" << Endl;
			Coutm << " LC2MCMC : Compute the signal to noise ratio (SNR) and Fisher Matrix (FIM)   - LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		
		
		
		
		// ***** Main declaration *****
		
		
		//! *** Output (results) file
		char bNOut[16384];
		strcpy(bNOut, "ResMCMC");
		
		//! *** Catalogue and parameters of sources
		char fNInCatalogue[16384];
		int NRecCat(-1); 
		int NRowsCat(-1);
		std::ifstream fInCatalogue;
		int TypeCat(0);	// 0->standard, 1->result file 
		bool StopEndFile(true);
		std::string junk;
		bool CheckOutput(false);		
		std::vector<char*> Words(0);
		
		//! **** Selection parameters
		//! ** Threshold in mass ratio
		double qth(1./20.);
		qth = 1./50.;
		
		strcpy(fNInCatalogue, "None");
		
		//! *** For SNR, FIM, errors
		LCMatrix FIM;
		LCMatrix CovM;
		double SNR;
		double * Err(NULL);
		double * ptry(NULL);
		
		//! *** Data
		dcomplex ** fDat(NULL);
		int NfDat, NTDI, iFmin;
		double df;
		
		
		//! ********** Read options **********
		
		for(int iarg=1; iarg<argc; iarg++){
			
			if((argc>1)&&(strcmp(argv[iarg],"-o")==0)){
				// *** Read base name of output file in option
				strcpy(bNOut, argv[iarg+1]);
				nOptions +=2;
			}
			
			
			if((argc>1)&&(strcmp(argv[iarg],"-G")==0)){
				bool NotFound(true);
				if( (MT.wcmp(argv[iarg+1],"GB")) || (MT.wcmp(argv[iarg+1],"GalBin")) ){
					X.sSrcType(GALBIN);
					NotFound = false;
				}
				if( (MT.wcmp(argv[iarg+1],"SH")) || (MT.wcmp(argv[iarg+1],"MBHBInsPrecHH")) ){
					X.sSrcType(SPINBBHHHARM);
					NotFound = false;
				}
				if( (MT.wcmp(argv[iarg+1],"NR")) || (MT.wcmp(argv[iarg+1],"NRwave")) ){
					X.sSrcType(NRWAVE);
					NotFound = false;
				}
				if(NotFound){
					Coutm << "ERROR : The source type " << argv[iarg+1] << " is unknown !" << Endl;
					throw std::invalid_argument("ERROR : The source type is unknown !");
				}
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-C")==0)){
				strcpy(fNInCatalogue, argv[iarg+1]);
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-N")==0)){
				Nsrc = atoi(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-f")==0)){
				iNStart = atoi(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-R")==0)){
				iRealiStudy = atoi(argv[iarg+1]);
				nOptions +=2;
			}
			/*if((argc>1)&&(strcmp(argv[iarg],"-r")==0)){
			 NRand = atoi(argv[iarg+1]);
			 nOptions +=2;
			 }*/
			
			if((argc>1)&&(strcmp(argv[iarg],"-dt")==0)){
				X.dtMesMin = atof(argv[iarg+1]);
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-D")==0)){
				X.Detector = UNKNOWNDC;
				if(MT.wcmp(argv[iarg+1], "LISA")){
					X.Detector = LISA;
				}
				if(MT.wcmp(argv[iarg+1], "ELISA")){
					X.Detector = ELISA;
				}
				if(MT.wcmp(argv[iarg+1], "C1")){
					X.Detector = C1;
				}
				if(MT.wcmp(argv[iarg+1], "C2")){
					X.Detector = C2;
				}
				if(MT.wcmp(argv[iarg+1], "C3")){
					X.Detector = C3;
				}
				if(MT.wcmp(argv[iarg+1], "C4")){
					X.Detector = C4;
				}
				if(MT.wcmp(argv[iarg+1], "C5")){
					X.Detector = C5;
				}
				
				if(X.Detector == UNKNOWNDC)
					throw std::invalid_argument("Unknow type of detector.");
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-cn")==0)){
				X.GalacticNoise = true;
				nOptions ++;
			}
			
			
			if((argc>1)&&(strcmp(argv[iarg],"-tdi")==0)){
				strcpy(TDINameAll, argv[iarg+1]);
				nOptions +=2;
			}
			
			/*if((argc>1)&&(strcmp(argv[iarg],"-S")==0)){
			 SNRth = atof(argv[iarg+1]);
			 nOptions +=2;
			 }*/
			
			
			if((argc>1)&&(strcmp(argv[iarg],"-P")==0)){
				//! **** Read list of parameters to study
				strcpy(CharIndPar, argv[iarg+1]);
				nOptions +=2;
			}
			
			/*if((argc>1)&&(strcmp(argv[iarg],"-cv")==0)){
			 if(argc-iarg<4)
			 throw std::invalid_argument("ERROR : Need 4 parameters for testing the convergence !");
			 ConviP = atoi(argv[iarg+1]);
			 ConvDMin = atof(argv[iarg+2]);
			 ConvDMax = atof(argv[iarg+3]);
			 ConvDLogStep = atof(argv[iarg+4]);
			 NConvStep = MT.ifloor((log10(ConvDMax)-log10(ConvDMin))/ConvDLogStep);
			 nOptions += 5;
			 }*/
			
			if((argc>1)&&(strcmp(argv[iarg],"-M")==0)){
				if(argc-iarg<4)
					throw std::invalid_argument("ERROR : Need 2 parameters (path of directory and max RAM size) for configuring swap !");
				X.sSwap(argv[iarg+2], atof(argv[iarg+1]));
				nOptions +=3;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-n")==0)){
				MCMCNStepMax = atoi(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-E")==0)){
				MCMCscale = atof(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-F")==0)){
				MCMCFIMRecomp = atoi(argv[iarg+1]);
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-co")==0)){
				CheckOutput = true;
				nOptions++;
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
			if((argc>1)&&(strcmp(argv[iarg],"-v")==0)){
				MT.setDispDetails();
				nOptions++;
			}
		} 
		MT.setRandSeed(SeedRand);
		
		
		//! **********  Configuration  **********
		
		//! ***** Decode a char string to obtain the list of parameters to study 
		if(MT.wcmp(CharIndPar, "None")){
			if(X.gSrcType()==SPINBBHHHARM)
				strcpy(CharIndPar, "0,1,14,13,12,4,L5,L15,17,6,7,8,9,10,11");
			if(X.gSrcType()==GALBIN)
				strcpy(CharIndPar, "0,1,2,3,L4,6,7");
		}
		X.sParStudy(CharIndPar);
		
		//! ***** Read configuration files in arguments
		for(int ia=1+nOptions; ia<argc; ia++){
			X.addfCfg(argv[ia]);
		}
		
		// ***** Display for checking which files we are using
		Coutm << "Configuration files : " << Endl;
		for(int i=0; i<X.NfNCfg; i++){
			Coutm << "\t + " << X.fNCfg[i] ; MT.CheckFile("", X.fNCfg[i]);
		}
		
		//! *** Read the bloc name of TDI generators and obtain the individual names and the number of generators
		X.sTDI(TDINameAll);
		Coutm << "TDI generators : ";
		for(int i=0; i<X.gNTDI(); i++)
			Coutm << " " << X.gTDIName(i);
		Coutm << Endl;
		
		
		
		//! **********  Initialization  **********
		
		X.AllocBaseMem();
		Xp.CopyBase(X);
		
		
		FIM.init(&MT, X.gNParStudy(), X.gNParStudy());
		CovM.init(&MT, X.gNParStudy(), X.gNParStudy());
		Err = (double*) MT.AllocMemory(X.gNParStudy()*sizeof(double));
		for(int i=0; i<X.gNParStudy();  i++)
			Err[i] = 0.;
		
		
		//! **** Open and analyse the catalogue
		if(!MT.wcmp(fNInCatalogue,"None")){
			Coutm << "Open the catalogue file " << fNInCatalogue << " ... " << Endl;
			fInCatalogue.open(fNInCatalogue);
			if(!fInCatalogue)
				throw std::invalid_argument("Problem in openning parameters file.");
			
			//! ** Count number of record
			while(fInCatalogue.peek() == '#')
				fInCatalogue.ignore(16384,'\n');
			NRecCat = 0;
			while((fInCatalogue.peek() != '\n')&&(!fInCatalogue.eof())){
				fInCatalogue >> junk;
				if(fInCatalogue.peek() == ' '){
					int ipos(fInCatalogue.tellg());
					ipos++;
					fInCatalogue.seekg(ipos);
				}
				NRecCat++;
			}
			fInCatalogue.close();
			fInCatalogue.clear();
			
			
			//! ** Count the number of rows
			fInCatalogue.open(fNInCatalogue);
			char fCatHead[16384];
			while(fInCatalogue.peek() == '#')
				fInCatalogue.getline(fCatHead, 16384,'\n');
			MT.wextract(fCatHead, Words);
			if(Words[4][0]=='z')
				TypeCat = 1;
			NRowsCat = 0;
			while(!fInCatalogue.eof()){
				fInCatalogue.ignore(16384,'\n');
				if(!fInCatalogue.eof()){	
					NRowsCat++;
				}
			}
			fInCatalogue.close();
			fInCatalogue.clear();
			if(Nsrc == -1)
				Nsrc = NRowsCat;
		}
		
		if(NRecCat!=-1){
			fInCatalogue.open(fNInCatalogue);
			while(fInCatalogue.peek() == '#')
				fInCatalogue.ignore(16384,'\n');
			if(iNStart*Nsrc<NRowsCat){
				int iS(0);
				while((!fInCatalogue.eof())&&(iS < iNStart*Nsrc)){
					fInCatalogue.ignore(16384,'\n');
					iS++;
				}
			}
		}
		
		
		//! ***** Prepare the result output file
		
		//! ** Common part
		int iTitCol(0);
		int iTitColBase(0);
		char TitColBase[10000];
		if(X.gSrcType()==SPINBBHHHARM){
			strcpy(TitColBase,"#realID[1] srcID[2] Name[3] iStep[4] z[5] M1[6] q[7] lamS[8] betS[9] phL[10] thL[11] phi0[12] tc[13] a1[14] a2[15] thS1[16] thS2[17] phS1[18] phS2[19] per[20] e[21] DL0[22] Mcz[23] eta[24] thBS1[25] thBS2[26] phBS1[27] phBS2[28]");
			iTitColBase = 29;
		}
		if(X.gSrcType()==GALBIN){
			strcpy(TitColBase,"#realID[1] srcID[2] Name[3] iStep[4] M1[5] M2[6] P[7] Pdot[8] e[9] i[10] lam[11] bet[12] d[13] psi[14] Amp[15] f[16] fdot[17] phi0[18]");
			iTitColBase = 19;
		}
		if(X.gSrcType()==NRWAVE){
			strcpy(TitColBase,"#realID[1] srcID[2] Name[3] RandID[4] z[5] M1[6] q[7] lamS[8] betS[9] phL[10] thL[11] phi0[12] tc[13] a1[14] a2[15] thS1[16] thS2[17] phS1[18] phS2[19] per[20] e[21] DL0[22] Mcz[23] eta[24] thBS1[25] thBS2[26] phBS1[27] phBS2[28]");
			iTitColBase = 29;
		}
		
		//! *** Open file containing error FIM results and write header
		char fNOutErr[10000];
		sprintf(fNOutErr,"%s-Err.txt",bNOut);
		std::ofstream fOutErr(fNOutErr);
		iTitCol = iTitColBase;
		fOutErr << TitColBase;
		for(int i=0; i<X.gNTDI(); i++)
		for(int iP=0; iP<X.gNParStudy(); iP++){
			fOutErr << " Par" << X.giParStudy(iP) << "[" << iTitCol++ << "]";
			fOutErr << " Err" << X.giParStudy(iP) << "[" << iTitCol++ << "]";
		}
		fOutErr << " Likelihood[" << iTitCol++ << "]" << Endl;
		fOutErr.precision(12);
		
		
		//! *** Open FIM output file and write header
		char fNOutFIM[10000];
		sprintf(fNOutFIM, "%s-FIM.txt", bNOut);
		std::ofstream fOutFIM(fNOutFIM);
		iTitCol = iTitColBase;
		fOutFIM << TitColBase;
		for(int iP1=0; iP1<X.gNParStudy(); iP1++)
			for(int iP2=0; iP2<X.gNParStudy(); iP2++)
				fOutFIM << " FIM" << X.giParStudy(iP1) << X.giParStudy(iP2) << "[" << iTitCol++ << "]";
		for(int iP1=0; iP1<X.gNParStudy(); iP1++)
			for(int iP2=0; iP2<X.gNParStudy(); iP2++)
				fOutFIM << " CovM" << X.giParStudy(iP1) << X.giParStudy(iP2) << "[" << iTitCol++ << "]";
		fOutFIM << Endl;
		fOutFIM.precision(12);
		
		
		//! *** Open MCMC output file and write header
		char fNOutMCMC[10000];
		sprintf(fNOutMCMC, "%s-MCMC.txt", bNOut);
		std::ofstream fOutMCMC(fNOutMCMC);
		iTitCol = iTitColBase;
		fOutMCMC << TitColBase << " LogL[" << iTitCol << "] MCMCratio[" << iTitCol+1 << "] MCMCRand[" << iTitCol+2 << "] Accepted[" << iTitCol+3 << "]";
		iTitCol += 4;
		for(int iP=0; iP<X.gNParStudy(); iP++)
			fOutMCMC << " pTry" << X.giParStudy(iP) << "[" << iTitCol++ << "]"; 
		fOutMCMC << Endl;
		fOutMCMC.precision(12);
			
		//! ********** Main part : loop on sources **********
		
		int iSrc(0);
		bool CContinue(true);
		while(CContinue){
			
			int realID(-1), srcID(iSrc);
			double  z(0.), M1(0.), q(0.), betS(0.), lamS(0.), phL(0.), thL(0.), phi0(0.), tc(0.), a1(0.), a2(0.), thS1(0.), thS2(0.), phS1(0.), phS2(0.), per(0.), ecc(0.), DL0(0.), DL0d(0.), Mcz(0.), eta(0.), thBS1(0.), thBS2(0.), phBS1(0.), phBS2(0.);
			char srcName[100];
			double lamSd(0.), betSd(0.), M2(0.), Period(0.), Pdot(0.), inc(0.), incd(0.), psi(0.), Amp(0.), freq(0.), fdot(0.);
			char StrParam[10000];
			
			time_t tstart, tendSNR, tend;
			
			
			//! **** Read source in catalogue
			if((NRecCat!=-1)&&(!fInCatalogue.eof())){
				if(X.gSrcType()==SPINBBHHHARM){
					if(NRecCat>15){
						int iRStart(26);
						//! ** All parameters are in the catalogue
						fInCatalogue >> realID >> srcID;
						if(TypeCat==1){
							fInCatalogue >> junk >> junk;
							iRStart += 2;
						}
						fInCatalogue >> z >> M1 >> q >> lamS >> betS >> phL >> thL >> phi0 >> tc >> a1 >> a2 >> thS1 >> thS2 >> phS1 >> phS2 >> per >> ecc >> DL0d >> Mcz >> eta >> thBS1 >> thBS2 >> phBS1 >> phBS2;
						Coutm << phBS2 << Endl;
						for(int iR=iRStart; iR<NRecCat; iR++){
							fInCatalogue >> junk;
							//Coutm << iR << " --> " << junk << Endl;
						}
						if(TypeCat==1)
							DL0d *= 1.e-3; //! Converion from kpc to Mpc because we use Mpc as input and we write kpc so when we use output of previous run as catalogue, we have to convert
						
					}else{
						//! ** Some parameters are in the catalogue but not the 
						fInCatalogue >> realID >> srcID >> z >> M1 >> q  >> lamSd >> betSd >> phL >> thL >> phi0 >> tc >> DL0d >> Mcz >> eta;
					}
					sprintf(srcName,"Real%d-%d", realID, srcID);
				}
				if(X.gSrcType()==GALBIN){
					fInCatalogue >> M1 >> M2 >> Period >> Pdot >> ecc >> incd >> lamSd >> betSd >> DL0 >> srcID >> srcID >> srcName;
				}
				Coutm << "===================== SOURCE " << iSrc << " : src " << srcID << " of realization " << realID << " =====================" << Endl;
				
				//! *** Stop running if we reach the end of the file
				if((StopEndFile)&&(fInCatalogue.eof()))
					CContinue = false;
			}else{
				Coutm << "===================== SOURCE " << iSrc << " :  random =====================" << Endl;
				q = MT.RandUniform(1./20., 1.);
				M1 = pow(10,MT.RandUniform(5,8));
			}
			
			
			
			if(CContinue){
				
				if((iRealiStudy==-1)||(iRealiStudy==realID)){
					
					
					//! ******** Compute only the source which pass the threshold
					if( (X.gSrcType()==GALBIN) || ((X.gSrcType()==SPINBBHHHARM)&&(q>qth)) ){
						
						X.NotValidPar = false;
						
						time(&tstart);
						
						//! ****** Set the parameters
						
						if(X.gSrcType()==SPINBBHHHARM){
							if(NRecCat>0){
								X.sP(0, betS);
								X.sP(1, lamS);
								X.sP(14, phL);
								X.sP(13, thL);
								X.sP(12, phi0);
								X.sP(4, tc);
								X.sP(5, DL0d*1.e3);
								X.sP(15, Mcz);
								X.sP(17, eta);
								if(NRecCat>15){
									X.sP(6, a1);
									X.sP(7, a2);
									X.sP(8, thBS1);
									X.sP(9, thBS2);
									X.sP(10, phBS1);
									X.sP(11, phBS2);
								}
							}else{
								X.sP(2, M1);
								X.sP(3, M1*q);
							}
						}
						
						if(X.gSrcType()==GALBIN){
							X.sP(8, M1);
							X.sP(9, M2);
							X.sP(5, 2./Period);
							X.sP(12, Pdot);
							X.sP(1, lamSd*M_PI/180.);
							X.sP(0, betSd*M_PI/180.);
							X.sP(10, DL0);
							if(!MT.deq(incd,60.))
								X.sP(3, incd*M_PI/180.);
						}
						
						if(X.gSrcType()==NRWAVE){
							if(NRecCat>0){
								X.sP(0, betS);
								X.sP(1, lamS);
								X.sP(3, M1*(1+q)*(1+z));
								X.sP(4, tc);
								X.sP(5, DL0d*1.e3);
								//! TO DO : COMPUTE Thd, Phd AND Polarization from betS, lamS, thL AND phL
								X.sP(2, 0.);	//! Set polarization
								X.sP(6, 0.);	//! Set Thd
								X.sP(7, 0.);	//! Set Phd
							}
						}
						
						
                        X.DispDetails = true;
                        if(CheckOutput){
                            X.CheckOut = true;
                        }
                        sprintf(X.BaseNameOut, "%s-%d", bNOut, iSrc);
                        
						//! ****** Compute SNR
						
						X.DispDetails = true;
						
						for(int iP=0; iP<X.gNp(); iP++)
							Coutm << X.gP(iP) << " " ;
						Coutm << Endl << "\t ==>";
						for(int i=0; i<X.gNTDI(); i++)
							Coutm << " SNR_" << X.gTDIName(i) << " = " <<  X.gSNR(i);
						if(X.gNTDI()==1)
							SNR = X.gSNR4();
						if(X.gNTDI()==2)
							SNR = X.gSNR6();
						Coutm << " : SNR = " <<  SNR << Endl;
						
						
						//! ****** Recover the parameters
						if(X.gSrcType()==SPINBBHHHARM){
							betS = X.gP(0);
							lamS = X.gP(1);
							phL = X.gP(14);
							thL = X.gP(13);
							phi0 = X.gP(12);
							tc = X.gP(4);
							DL0 = X.gP(5);
							Mcz = X.gP(15);
							eta = X.gP(17);
							//M1 = X.gP(2);
							a1 = X.gP(6);
							a2 = X.gP(7);
							thBS1 = X.gP(8);
							thBS2 = X.gP(9);
							phBS1 = X.gP(10);
							phBS2 = X.gP(11);
							
							sprintf(StrParam,"%d %d %s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", realID, srcID, srcName, 0, z, M1, q, lamS, betS, phL, thL, phi0, tc, a1, a2, thS1, thS2, phS1, phS2, per, ecc, DL0, Mcz, eta, thBS1, thBS2, phBS1, phBS2);
						}
						
						if(X.gSrcType()==GALBIN){
							M1 = X.gP(8);
							M2 = X.gP(9);
							//Period = 1./(X.gP(5));
							//Pdot = X.gP(12);
							inc = X.gP(3);
							lamS = X.gP(1);
							betS = X.gP(0);
							DL0 = X.gP(10);
							
							psi = X.gP(2);
							Amp = X.gP(4);
							freq = X.gP(5);
							fdot = X.gP(6);
							phi0 = X.gP(7);
							ecc = X.gP(11);
							
							sprintf(StrParam,"%d %d %s %d %.8lf %.8lf %.8lf %.10e %lf %.8lf %.8lf %.8lf %.10e %.8lf %.10e %.10e %.10e %.8lf", realID, srcID, srcName, 0, M1, M2, Period, Pdot, ecc, inc, lamS, betS, DL0, psi, Amp, freq, fdot, phi0);
						}
						
						if(X.gSrcType()==NRWAVE){
							double Pol, Thd, Phd, Mz;
							betS = X.gP(0);
							lamS = X.gP(1);
							Pol = X.gP(2);
							Mz = X.gP(3);
							tc = X.gP(4);
							DL0 = X.gP(5);
							Thd = X.gP(6);
							Phd = X.gP(7);
							q = X.gP(8);
							a1 = X.gP(9);
							a2 = X.gP(10);
							thS1 = X.gP(11);
							thS2 = X.gP(12);
							phS1 = X.gP(13);
							phS2 = X.gP(14);
							phi0 = X.gP(15);
							thBS1 = X.gP(16);
							thBS2 = X.gP(17);
							phBS1 = X.gP(18);
							phBS2 = X.gP(19);
							thL	= X.gP(20);
							phL	= X.gP(21);
							
							
							z = MT.redshift(DL0, 70., 0.3, 0.7, 0.);
							eta = q/((1.+q)*(1.+q)); 
							Mcz = Mz*pow(eta,0.6);
							M1 = Mz/((1+q)*(1+z));
							
							sprintf(StrParam,"%d %d %s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", realID, srcID, srcName, 0, z, M1, q, lamS, betS, phL, thL, phi0, tc, a1, a2, thS1, thS2, phS1, phS2, per, ecc, DL0, Mcz, eta, thBS1, thBS2, phBS1, phBS2);
						}
						
						
						fOutMCMC << StrParam << " " << SNR*SNR/2. << " 1. 1. 1";
						X.gPStudy(ptry);
						for(int iP=0; iP<X.gNParStudy(); iP++)
							fOutMCMC << " " << ptry[iP];
						fOutMCMC << Endl;
						
						
						//! ****** Save the signal as data
						X.gfSig(fDat, NfDat, NTDI, iFmin, df);
						
						
						time(&tendSNR);
						MT.MemDisplay();
						
						
						
						//! ****** Compute the FIM
						X.DispDetails = false;
						X.CheckOut = false;
						Coutm << Endl << ">>>>>>>>>>>>>>>>>>>>>>>> Computation of FIM and errors ..." << Endl;
						
						
						
						
						//X.CheckOut = true;
						
						if(X.gNTDI()==1){
							FIM = X.gFIM4();
							CovM = X.gCovM4();
							for(int i=0; i<X.gNParStudy(); i++)
								Err[i] = X.gErr4(i);
						}
						if(X.gNTDI()==2){
							FIM = X.gFIM6();
							CovM = X.gCovM6();
							for(int i=0; i<X.gNParStudy(); i++)
								Err[i] = X.gErr6(i);
						}
												
						//! **** Write FIM and CovM
						fOutFIM << StrParam;
						for(int iP=0; iP<X.gNParStudy(); iP++)
							for(int jP=0; jP<X.gNParStudy(); jP++)
								fOutFIM << " " << FIM(iP,jP);
						for(int iP=0; iP<X.gNParStudy(); iP++)
							for(int jP=0; jP<X.gNParStudy(); jP++)
								fOutFIM << " " << CovM(iP,jP);
						fOutFIM << Endl; fOutFIM.flush();
						
						Coutm << "Covariance matrix error estimate : " << Endl;
						for(int iP=0; iP<X.gNParStudy(); iP++)
							Coutm << "\t - " << iP << " : " << X.giParStudy(iP) << " : " <<  X.gP(X.giParStudy(iP)) << " @ " << Err[iP] << Endl;
							
						
						fOutErr << StrParam ;
						for(int iP=0; iP<X.gNParStudy(); iP++)
							fOutErr << " " << X.gP(X.giParStudy(iP)) << " " << Err[iP];
						fOutErr << " " << SNR*SNR/2. << Endl; fOutErr.flush();
						
						
						time(&tend);
						
						Coutm << "Time usage : tot = " << (tend-tstart) << " s" << Endl;
						Coutm << "\tTU:\t SNR               : " << 100.*(tendSNR-tstart)/(tend-tstart) << " %  ( " << (tendSNR-tstart) << " s)" << Endl;
						Coutm << "\tTU:\t FIM               : " << 100.*(tend-tendSNR)/(tend-tstart) << " %  ( " << (tend-tendSNR) << " s)" << Endl;
						
						
						
						
						//! ****** Run MCMC
						
						double lLX, lLXp;
						double qX, qXp;
						bool MCMCContinue(true);
						double MCMCratio(0.), MCMCalpha(0.), MCMCbeta(0.);
						int iMCMC(1), NAccept(0), move(0);
						
						
						//Y = X;
						
						qX = 1.;
						
						X.Reset(false);
						lLX = X.LogLikelihood(fDat, NTDI, NfDat, iFmin, df);
						
						Xp.CopyParam(X);
						lLXp = lLX;
						qXp = qX;
						
						
						iMCMC = 0;
						while(MCMCContinue){
							
							//! **** Jump : new point X
							X.MCMCjump(MCMCheat, MCMCscale);
							
							
							//! **** Compute MCMC quantities
							
							lLX = X.LogLikelihood(fDat, NTDI, NfDat, iFmin, df);
							if(lLX<-1.e30){
								throw std::invalid_argument("ERROR found with success !!!");
							}
							
							MCMCratio = ( exp((lLX - lLXp)/MCMCheat) * qX / qXp ) ;
							
							MCMCalpha = MIN(1., MCMCratio);
							
							MCMCbeta = MT.RandUniform(0., 1.0);
							
							X.gPStudy(ptry);
							
							Coutm << iMCMC << " MCMC (scale = " << MCMCscale << " , heat = " << MCMCheat << ") >>> lLX = " << lLX << "  ==> MCMCratio = " << MCMCratio << "  ==> MCMCalpha = " << MCMCalpha << " vs " << MCMCbeta;
							
							if( MCMCalpha > MCMCbeta ){ //!< *** Accept the jump
								move = 1;
								Xp.CopyParam(X);
								lLXp = lLX;
								qXp = qX;
								NAccept++;
								Coutm << " ==> Accepted ( " << NAccept << " )";
							}else{
								move = 0;
								X.CopyParam(Xp);
								Coutm << " ==> Rejected ";
							}
	
							
							
							
							
							//! **** Record
							//! ** Recover the parameters
							if(X.gSrcType()==SPINBBHHHARM){
								betS = X.gP(0);
								lamS = X.gP(1);
								phL = X.gP(14);
								thL = X.gP(13);
								phi0 = X.gP(12);
								tc = X.gP(4);
								DL0 = X.gP(5);
								Mcz = X.gP(15);
								eta = X.gP(17);
								M1 = X.gP(2);
								a1 = X.gP(6);
								a2 = X.gP(7);
								thBS1 = X.gP(8);
								thBS2 = X.gP(9);
								phBS1 = X.gP(10);
								phBS2 = X.gP(11);
								M1 = X.gP(2);
								q = X.gP(3)/M1;
								sprintf(StrParam,"%d %d %s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", realID, srcID, srcName, iMCMC, z, M1, q, lamS, betS, phL, thL, phi0, tc, a1, a2, thS1, thS2, phS1, phS2, per, ecc, DL0, Mcz, eta, thBS1, thBS2, phBS1, phBS2);
							}
							
							if(X.gSrcType()==GALBIN){
								M1 = X.gP(8);
								M2 = X.gP(9);
								//Period = 1./(X.gP(5));
								//Pdot = X.gP(12);
								inc = X.gP(3);
								lamS = X.gP(1);
								betS = X.gP(0);
								DL0 = X.gP(10);
								psi = X.gP(2);
								Amp = X.gP(4);
								freq = X.gP(5);
								fdot = X.gP(6);
								phi0 = X.gP(7);
								ecc = X.gP(11);
								sprintf(StrParam,"%d %d %s %d %.8lf %.8lf %.8lf %.10e %lf %.8lf %.8lf %.8lf %.10e %.8lf %.10e %.10e %.10e %.8lf", realID, srcID, srcName, iMCMC, M1, M2, Period, Pdot, ecc, inc, lamS, betS, DL0, psi, Amp, freq, fdot, phi0);
							}
							fOutMCMC << StrParam << " " << lLX << " " << MCMCratio << " " << MCMCbeta << " " << move;
							for(int iP=0; iP<X.gNParStudy(); iP++)
								fOutMCMC << " " << ptry[iP];
							fOutMCMC << Endl;
							
							
							
							
							
							iMCMC++;
							
							Coutm << " (accpetance rate = " << 100.*NAccept/((double)iMCMC) << " %)." << Endl << Endl << Endl;
							
							//! **** Test stop conditions
							if(iMCMC > MCMCNStepMax)
								MCMCContinue = false;
							
							if((move)&&(MCMCFIMRecomp>0)&&(NAccept%MCMCFIMRecomp==0)){
								X.ComputeFIM();	
							}
							
						};
						

						
					}
					
					
					
					
				}else{
					X.NotValidPar = true;
					Coutm << ">>> Rejected because q = " << q << " < qth = " << qth << Endl; 
				}
				
				
			}else{
				//! ** Not increase the index of source
				iSrc--;
			}
			
			//! Increase the index of source
			iSrc++;
			if(iSrc>=Nsrc)
				CContinue = false;
			
		}
		
		fOutErr.close();
		fOutFIM.close();
		fOutMCMC.close();
		
		if(NRecCat!=-1){
			fInCatalogue.close();
			fInCatalogue.clear();
		}
		
		for(int iW=0; iW<Words.size(); iW++)
			MT.Free(Words[iW],256*sizeof(char));
		
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
		
		
		//! ******************** Free memory
		
		if(Err!=NULL)
			MT.Free(Err, X.gNParStudy()*sizeof(double));
		
		if(ptry!=NULL)
			MT.Free(ptry, X.gNParStudy()*sizeof(double));
		
		if(fDat!=NULL){
			for(int iS=0; iS<NTDI; iS++)
				if(fDat[iS]!=NULL)
					MT.Free(fDat[iS], NfDat*sizeof(dcomplex));
			MT.Free(fDat, NTDI*sizeof(dcomplex*));
		}
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "LC2MCMC: error: " << e.what()<<Endl;
		std::cerr << "LC2MCMC: abort!" << Endl;
		exit(1);
	}
	return(0);
};



/** \}*/

// end of LISACODE-SNR.cpp
