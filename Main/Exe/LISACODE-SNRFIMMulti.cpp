// $Id:  $
/*
 *  LISACODE- SNRFIMMulti.cpp
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




/** \brief The executable computes the SNR and the Fisher Matrix for multiple GW sources
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
		
		int Nsrc(-1);
		int iNStart(0);
		int iRealiStudy(-1);
		int NRand(1);
		
		LCModSig S(&MT);
		
		S.Detector = ELISA;
		S.GalacticNoise = false;
		S.sSrcType(SPINBBHHHARM);
		
		double SNRth(6.);
		char CharIndPar[1064];
		strcpy(CharIndPar, "None");
		
		/*
		int NSC(6);
		for(int i=0; i<(NSC*(NSC-1)); i++){
			int emC(1+floor(i/(NSC-1)));
			int reC((i+emC)%NSC+((i+emC)%NSC>=emC?1:0));
			Coutm << i << " : " << emC << " --> " << reC << Endl;
		}
		for(int ie=1; ie<=NSC; ie++){
			for(int ir=1; ir<=NSC; ir++){
				if(ie!=ir)
					Coutm << ie << " --> " << ir << " : " <<  (ie-1)*(NSC-1)+(ir<ie?ir-1:ir-2) << "   " <<  (ie-1)*(NSC-1) << "   " << (ir<ie?ir-1:ir-2) << Endl;
			}
		}
		return 0;
		*/
		
		// *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution :" << Endl;
			Coutm << "\t\t(./)LC2 SNRFIMMulti [Options] ConfigFile1.xml ConfigFile2.xml ..." << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << Endl << "\t\t * ConfigFileI.xml (required) : Configuration files. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -o %sOutFile  : Output file containing results [default: SNRResults.txt]"  << Endl ;
			
			Coutm << "\t\t * -G key  : Type of GW source. [default : MBHBInsPrecHH ] "  << Endl ;
			Coutm << "\t\t\t GB or GalBin           : galactic binary,"  << Endl ;
			Coutm << "\t\t\t SH or MBHBInsPrecHH    : Massive black hole binary inspiral phase with spin precession and higher harmonics starting integration at t = 0s. Standard catalog (no thJ phJ)."  << Endl ;
			//Coutm << "\t\t\t S9 or MBHBInsPrecHH9M  : Massive black hole binary inspiral phase with spin precession and higher harmonics starting integration at t = 9M."  << Endl ;
			Coutm << "\t\t\t SW or MBHBInsPrecHHOmegaM %fOmegaM : Massive black hole binary inspiral phase with spin precession and higher harmonics starting integration using a PN matching of omega M , i.e. orbital angular momentum times total mass in seconds. Catalog with extra parameters thJ and phJ."  << Endl ;
			Coutm << "\t\t\t NR or NRwave           : Waveform based on numerical relativity data."  << Endl ;
			Coutm << "\t\t\t NI or NRwaveInsp       : Waveform based on numerical relativity data but only the inspiral part (the same taper as the one use in SH is applied)."  << Endl ;
			Coutm << "\t\t * -C %sCatalog : Catalogue containing the list of parameters : if not random source. [default : None] "  << Endl ;
			Coutm << "\t\t * -N %dNSrc    : Number of sources to study (-1 : number of source in catalogue). [default : " << Nsrc << " ]"  << Endl ;
			Coutm << "\t\t * -f %diFSrc   : Fisrt source to study is iFSrc*NSrc. [default : "  << iNStart << "]"  << Endl ;
			Coutm << "\t\t * -R %diReali  : Study only the source from the realization %diReali (-1:not used). [default : "<< iRealiStudy << "]"  << Endl ;
			Coutm << "\t\t * -r %dNrand   : Make Nrand randomization of free paramters for each source. [default : "<< NRand << "]"  << Endl ;
			
			Coutm << "\t\t * -dt %step : Minimal time step. [default : " << S.dtMesMin << "]"  << Endl ;
			
			Coutm << "\t\t * -D %dConfig : For computing analytic noise : LISA,ELISA,C1,C2,C3,C4,C5,C6,SGOF1,SGOM1 [default : ELISA]"  << Endl ;
			Coutm << "\t\t * -cn         : Include analytic confusion noise  [default : None]"  << Endl ;
			
			//Coutm << "\t\t * -tdi %sObsName : Names of TDI generators [default : X,Am,Em]"  << Endl ;
			
			Coutm << "\t\t * -S %fSNRth : SNR threshold for selected for FIM compute [default : " << SNRth << "]"  << Endl ;
			Coutm << "\t\t * -P %diP1,%diP2  : Index of parameters to study ('L' before the index mean that we have to use the log). [default : 0,1,14,13,12,4,L5,L15,17,6,7,8,9,10,11 if src MBHBInsPrecHH , 0,1,2,3,L4,6,7 if src GalBin , 0,1,2,L3,4,L5,6,7 if src NRwave ]"  << Endl ;
			Coutm << "\t\t * -cv %diP %fDeltaMin %fDeltaMax %fDeltaLogStep : Test convergence for parameter iP with value of Delta between DeltaMin and DeltaMax with step log step of DeltaLogStep [default : None]"  << Endl ;
			Coutm << "\t\t * -M %fMax %sDir : Maximum memory in RAM and directory for swaping [default : Max = " << S.gMaxMemRAM() << ", Dir = " << S.gDirSwapMem() << " ]"  << Endl ;
			Coutm << "\t\t * -fd %FactStepD : Factor on step used in numerical derivative" << Endl;
            
			Coutm << "\t\t * -co       : Write the check output files. [default : false]"  << Endl ;
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
			Coutm << " LC2SNRFIMMulti : Compute the signal to noise ratio (SNR) and Fisher Matrix (FIM)   - LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		time_t jobstarttime;
		char hostname[256];
		time ( &jobstarttime );
		gethostname(hostname, 256);
		Coutm << "Job start on " << hostname << " at " << ctime (&jobstarttime);
		
		
		
		// ***** Main declaration *****
		
		//! *** Name of TDI observables
		char TDINameAll[16384];
		strcpy(TDINameAll, "X,Am,Em");
		
		
		//! *** Output (results) file
		char bNOut[16384];
		strcpy(bNOut, "ResSNRFIM");
		
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
		//qth = 1.;
		
		strcpy(fNInCatalogue, "None");
        bool ExtraParamJ(false);
		
		//! *** For testing convergence
		double ConvDMin(-1.), ConvDMax(-1.), ConvDLogStep(-1.);
		int ConviP(-1);
		int NConvStep(1); //!< Number of step for testing te convergence
		
		//! *** For FIM
		LCMatrix FIM4, FIM6;
		LCMatrix CovM4, CovM6;
		
		
		//! ********** Read options **********
		
		for(int iarg=1; iarg<argc; iarg++){
			
			//Coutm << argv[iarg] << " " << nOptions << Endl;
			
			if((argc>1)&&(strcmp(argv[iarg],"-o")==0)){
				// *** Read base name of output file in option
				strcpy(bNOut, argv[iarg+1]);
				nOptions +=2;
			}
			
			
			if((argc>1)&&(strcmp(argv[iarg],"-G")==0)){
				bool NotFound(true);
				if( (MT.wcmp(argv[iarg+1],"GB")) || (MT.wcmp(argv[iarg+1],"GalBin")) ){
					S.sSrcType(GALBIN);
					NotFound = false;
					nOptions +=2;
				}
				if( (MT.wcmp(argv[iarg+1],"SH")) || (MT.wcmp(argv[iarg+1],"MBHBInsPrecHH")) ){
					S.sSrcType(SPINBBHHHARM);
					NotFound = false;
					nOptions +=2;
				}
                /*
				if( (MT.wcmp(argv[iarg+1],"S9")) || (MT.wcmp(argv[iarg+1],"MBHBInsPrecHH9M")) ){
					S.sSrcType(SPINBBHHHARM);
					S.AddGWSpePar(1, 9.);
					NotFound = false;
					nOptions +=2;
				}
                 */
				if( (MT.wcmp(argv[iarg+1],"SW")) || (MT.wcmp(argv[iarg+1],"MBHBInsPrecHHOmega")) ){
					S.sSrcType(SPINBBHHHARM);
					S.AddGWSpePar(3, atof(argv[iarg+2]));
                    ExtraParamJ = true;
					NotFound = false;
					nOptions +=3;
				}
				if( (MT.wcmp(argv[iarg+1],"NR")) || (MT.wcmp(argv[iarg+1],"NRwave")) ){
					S.sSrcType(NRWAVE);
					NotFound = false;
					nOptions +=2;
				}
				if( (MT.wcmp(argv[iarg+1],"NI")) || (MT.wcmp(argv[iarg+1],"NRwaveInsp")) ){
					S.sSrcType(NRWAVE);
					S.AddGWSpePar(0, 1.);
					NotFound = false;
					nOptions +=2;
				}
				if(NotFound){
					Coutm << "ERROR : The source type " << argv[iarg+1] << " is unknown !" << Endl;
					throw std::invalid_argument("ERROR : The source type is unknown !");
				}
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
			if((argc>1)&&(strcmp(argv[iarg],"-r")==0)){
				NRand = atoi(argv[iarg+1]);
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-dt")==0)){
				S.dtMesMin = atof(argv[iarg+1]);
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-D")==0)){
				S.Detector = UNKNOWNDC;
				if(MT.wcmp(argv[iarg+1], "LISA")){
					S.Detector = LISA;
				}
				if(MT.wcmp(argv[iarg+1], "ELISA")){
					S.Detector = ELISA;
				}
				if(MT.wcmp(argv[iarg+1], "C1")){
					S.Detector = C1;
				}
				if(MT.wcmp(argv[iarg+1], "C2")){
					S.Detector = C2;
				}
				if(MT.wcmp(argv[iarg+1], "C3")){
					S.Detector = C3;
				}
				if(MT.wcmp(argv[iarg+1], "C4")){
					S.Detector = C4;
				}
				if(MT.wcmp(argv[iarg+1], "C5")){
					S.Detector = C5;
				}
				if(MT.wcmp(argv[iarg+1], "C6")){
					S.Detector = C6;
				}
				if(MT.wcmp(argv[iarg+1], "SGOF1")){
					S.Detector = SGOF1;
				}
				if(MT.wcmp(argv[iarg+1], "SGOM1")){
					S.Detector = SGOM1;
				}
				
				if(S.Detector == UNKNOWNDC)
					throw std::invalid_argument("Unknow type of detector.");
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-cn")==0)){
				S.GalacticNoise = true;
				nOptions ++;
			}
			
			/*
			 if((argc>1)&&(strcmp(argv[iarg],"-tdi")==0)){
			 strcpy(TDINameAll, argv[iarg+1]);
			 nOptions +=2;
			 }*/
			
			if((argc>1)&&(strcmp(argv[iarg],"-S")==0)){
				SNRth = atof(argv[iarg+1]);
				nOptions +=2;
			}
			
			
			if((argc>1)&&(strcmp(argv[iarg],"-P")==0)){
				//! **** Read list of parameters to study
				strcpy(CharIndPar, argv[iarg+1]);
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-cv")==0)){
				if(argc-iarg<4)
					throw std::invalid_argument("ERROR : Need 4 parameters for testing the convergence !");
				ConviP = atoi(argv[iarg+1]);
				ConvDMin = atof(argv[iarg+2]);
				ConvDMax = atof(argv[iarg+3]);
				ConvDLogStep = atof(argv[iarg+4]);
				NConvStep = MT.ifloor((log10(ConvDMax)-log10(ConvDMin))/ConvDLogStep);
				nOptions += 5;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-M")==0)){
				if(argc-iarg<4)
					throw std::invalid_argument("ERROR : Need 2 parameters (path of directory and max RAM size) for configuring swap !");
				S.sSwap(argv[iarg+2], atof(argv[iarg+1]));
				nOptions +=3;
			}
            
            if((argc>1)&&(strcmp(argv[iarg],"-fd")==0)){
				S.FactStepD = atof(argv[iarg+1]);
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
				//char cmd[10000];
				//sprintf(cmd, "hostname > Host_%s", argv[iarg+1]);
				//system(cmd);
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
			if(S.gSrcType()==SPINBBHHHARM)
				strcpy(CharIndPar, "0,1,12,4,L5,L15,17,6,7,8,9,10,11,21,22"); // changed by Sofiane
				//strcpy(CharIndPar, "0,1,14,13,12,4,L5,L15,17,6,7,8,9,10,11,21,22"); // changed by Sofiane
			if(S.gSrcType()==GALBIN)
				strcpy(CharIndPar, "0,1,2,3,L4,6,7");
			if(S.gSrcType()==NRWAVE)
				strcpy(CharIndPar, "0,1,L3,4,L5,6,7");  // changed by Sofiane polarization (p = 2) is not a free parameter
		}
		S.sParStudy(CharIndPar);
		
		//! ***** Read configuration files in arguments
		for(int ia=1+nOptions; ia<argc; ia++){
			S.addfCfg(argv[ia]);
		}
		
		// ***** Display for checking which files we are using
		Coutm << "Configuration files : " << Endl;
		for(int i=0; i<S.NfNCfg; i++){
			Coutm << "\t + " << S.fNCfg[i] ; MT.CheckFile("", S.fNCfg[i]);
		}
		
		//! *** Read the bloc name of TDI generators and obtain the individual names and the number of generators
		S.sTDI(TDINameAll);
		Coutm << "TDI generators : ";
		for(int i=0; i<S.gNTDI(); i++)
			Coutm << " " << S.gTDIName(i);
		Coutm << Endl;
		
		
		
		//! **********  Initialization  **********
		
		S.AllocBaseMem();
		
		FIM4.init(&MT, S.gNParStudy(), S.gNParStudy());
		FIM6.init(&MT, S.gNParStudy(), S.gNParStudy());
		CovM4.init(&MT, S.gNParStudy(), S.gNParStudy());
		CovM6.init(&MT, S.gNParStudy(), S.gNParStudy());
		
		
		
		
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
				//Coutm << NRecCat << " --> " << junk << Endl;
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
		
		Coutm << "\t\t ==> in catalogue : " << NRecCat << " records and " << Nsrc << " sources." << Endl;
		
		//! ***** Prepare the result output file
		
		char fNOutMain[10000];
		sprintf(fNOutMain,"%s.txt",bNOut);
		std::ofstream fOutMain(fNOutMain);
		
		//! *** Write the titles
		
		//! ** Common part
		int iTitCol(0);
		int iTitColBase(0);
		char TitColBase[10000];
		if(S.gSrcType()==SPINBBHHHARM){
			strcpy(TitColBase,"#realID[1] srcID[2] Name[3] RandID[4] z[5] M1[6] q[7] lamS[8] betS[9] phL[10] thL[11] phi0[12] tc[13] a1[14] a2[15] thS1[16] thS2[17] phS1[18] phS2[19] per[20] e[21] DL0[22] Mcz[23] eta[24] thBS1[25] thBS2[26] phBS1[27] phBS2[28] thBJ[29] phBJ[30]");
			iTitColBase = 31;
		}
		if(S.gSrcType()==GALBIN){
			strcpy(TitColBase,"#realID[1] srcID[2] Name[3] RandID[4] M1[5] M2[6] P[7] Pdot[8] e[9] i[10] lam[11] bet[12] d[13] psi[14] Amp[15] f[16] fdot[17] phi0[18]");
			iTitColBase = 19;
		}
		if(S.gSrcType()==NRWAVE){
			strcpy(TitColBase,"#realID[1] srcID[2] Name[3] RandID[4] z[5] M1[6] q[7] lamS[8] betS[9] phL[10] thL[11] phi0[12] tc[13] a1[14] a2[15] thS1[16] thS2[17] phS1[18] phS2[19] per[20] e[21] DL0[22] Mcz[23] eta[24] thBS1[25] thBS2[26] phBS1[27] phBS2[28] thBJ[29] phBJ[30]");
			iTitColBase = 31;
		}
		
		//! ** Write header of the main file containing resuls
		iTitCol = iTitColBase;
		fOutMain << TitColBase;
		for(int i=0; i<S.gNTDI(); i++)
			fOutMain << " SNR[" << iTitCol++ << "]";
		fOutMain << " SNR4[" << iTitCol++ << "]";
		fOutMain << " SNR6[" << iTitCol++ << "]";
		for(int iP=0; iP<S.gNParStudy(); iP++)
			fOutMain << " Err4P" << S.giParStudy(iP) << "[" << iTitCol++ << "]";
		fOutMain << " Err4Sky[" << iTitCol++ << "]";
		fOutMain << " Err4MajAx[" << iTitCol++ << "]";
		fOutMain << " Err4MinAx[" << iTitCol++ << "]";
		for(int iP=0; iP<S.gNParStudy(); iP++)
			fOutMain << " Err6P" << S.giParStudy(iP) << "[" << iTitCol++ << "]";
		fOutMain << " Err6Sky[" << iTitCol++ << "]";
		fOutMain << " Err6MajAx[" << iTitCol++ << "]";
		fOutMain << " Err6MinAx[" << iTitCol++ << "]";
		fOutMain << Endl;
		
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
		
		
		//! *** Open FIM output file
		char fNOutFIM[10000];
		if(ConviP>=0)
			sprintf(fNOutFIM, "%s-FIM-CP%d.txt", bNOut, ConviP);
		else
			sprintf(fNOutFIM, "%s-FIM.txt", bNOut);
		
		std::ofstream fOutFIM(fNOutFIM);
		
		//! *** Write header of FIM output file
		iTitCol = iTitColBase;
		fOutFIM << TitColBase;
		if(ConviP>=0){
			fOutFIM << " D_" << ConviP << "[" << iTitCol++ << "]";
			//fOutFIM << " P-D_" << ConviP << "[" << iTitCol++ << "]";
			//fOutFIM << " P+D_" << ConviP << "[" << iTitCol++ << "]";
		}
        for(int i=0; i<S.gNTDI(); i++)
			fOutFIM << " SNR[" << iTitCol++ << "]";
		fOutFIM << " SNR4[" << iTitCol++ << "]";
		fOutFIM << " SNR6[" << iTitCol++ << "]";
		for(int iP=0; iP<S.gNParStudy(); iP++)
			fOutFIM << " Err4P" << S.giParStudy(iP) << "[" << iTitCol++ << "]";
		fOutFIM << " Err4Sky[" << iTitCol++ << "]";
		fOutFIM << " Err4MajAx[" << iTitCol++ << "]";
		fOutFIM << " Err4MinAx[" << iTitCol++ << "]";
		for(int iP=0; iP<S.gNParStudy(); iP++)
			fOutFIM << " Err6P" << S.giParStudy(iP) << "[" << iTitCol++ << "]";
		fOutFIM << " Err6Sky[" << iTitCol++ << "]";
		fOutFIM << " Err6MajAx[" << iTitCol++ << "]";
		fOutFIM << " Err6MinAx[" << iTitCol++ << "]";
		for(int iP1=0; iP1<S.gNParStudy(); iP1++)
			for(int iP2=0; iP2<S.gNParStudy(); iP2++)
				fOutFIM << " FIM4P" << S.giParStudy(iP1) << S.giParStudy(iP2) << "[" << iTitCol++ << "]";
		for(int iP1=0; iP1<S.gNParStudy(); iP1++)
			for(int iP2=0; iP2<S.gNParStudy(); iP2++)
				fOutFIM << " Cov4P" << S.giParStudy(iP1) << S.giParStudy(iP2) << "[" << iTitCol++ << "]";
		for(int iP1=0; iP1<S.gNParStudy(); iP1++)
			for(int iP2=0; iP2<S.gNParStudy(); iP2++)
				fOutFIM << " FIM6P" << S.giParStudy(iP1) << S.giParStudy(iP2) << "[" << iTitCol++ << "]";
		for(int iP1=0; iP1<S.gNParStudy(); iP1++)
			for(int iP2=0; iP2<S.gNParStudy(); iP2++)
				fOutFIM << " Cov6P" << S.giParStudy(iP1) << S.giParStudy(iP2) << "[" << iTitCol++ << "]";
		fOutFIM << Endl;
		fOutFIM.precision(12);
		
		
		//! ********** Main part : loop on sources **********
		
		int iSrc(0);
		bool CContinue(true);
		while(CContinue){
			
			int realID(-1), srcID(iSrc);
			double  z(0.), M1(0.), q(0.), betS(0.), lamS(0.), phL(0.), thL(0.), phi0(0.), tc(0.), a1(0.), a2(0.), thS1(0.), thS2(0.), phS1(0.), phS2(0.), per(0.), ecc(0.), DL0(0.), DL0d(0.), Mcz(0.), eta(0.), thBS1(0.), thBS2(0.), phBS1(0.), phBS2(0.), thBJ(0.), phBJ(0.);
			char srcName[100];
			double lamSd(0.), betSd(0.), M2(0.), Period(0.), Pdot(0.), inc(0.), incd(0.), psi(0.), Amp(0.), freq(0.), fdot(0.);
			char StrParam[10000];
			
			time_t tstart, tendSNR, tend;
			
			strcpy(srcName,"None");
			
			//! **** Read source in catalogue
			if((NRecCat!=-1)&&(!fInCatalogue.eof())){
				if((S.gSrcType()==SPINBBHHHARM)||(S.gSrcType()==NRWAVE)){
					if(NRecCat>15){
						int iRStart(28);               // changed by Sofiane
						//! ** All parameters are in the catalogue
						fInCatalogue >> realID >> srcID;
						if(TypeCat==1){
							fInCatalogue >> junk >> junk;
							iRStart += 2;
						}
                        if(ExtraParamJ)
                            fInCatalogue >> z >> M1 >> q >> lamS >> betS >> phL >> thL >> phi0 >> tc >> a1 >> a2 >> thS1 >> thS2 >> phS1 >> phS2 >> per >> ecc >> DL0d >> Mcz >> eta >> thBS1 >> thBS2 >> phBS1 >> phBS2 >> thBJ >> phBJ;
                        else
                            fInCatalogue >> z >> M1 >> q >> lamS >> betS >> phL >> thL >> phi0 >> tc >> a1 >> a2 >> thS1 >> thS2 >> phS1 >> phS2 >> per >> ecc >> DL0d >> Mcz >> eta >> thBS1 >> thBS2 >> phBS1 >> phBS2;
						//Coutm << phBJ << Endl;
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
				if(S.gSrcType()==GALBIN){
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
				
				for(int iRand=0; iRand<NRand; iRand++){
					
					if(NRand>1){
						if(MT.wcmp(srcName,"None")){
							char tmpsrcName[100];
							strcpy(tmpsrcName,srcName);
							sprintf(srcName,"%s-Rnd%d", tmpsrcName, iRand);
						}else{
							sprintf(srcName,"Rnd%d", iRand);
						}
					}
					
					
					if((iRealiStudy==-1)||(iRealiStudy==realID)){
						
						Coutm << "================ >>>>> Random " << iRand << Endl;
                       
						
						if(S.gSrcType()==SPINBBHHHARM){
                            if(ExtraParamJ)
                                sprintf(StrParam,"%d %d %s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", realID, srcID, srcName, iRand, z, M1, q, lamS, betS, phL, thL, phi0, tc, a1, a2, thS1, thS2, phS1, phS2, per, ecc, DL0d*1.e3, Mcz, eta, thBS1, thBS2, phBS1, phBS2, thBJ, phBJ);  // changed by Sofiane
                            else
                                sprintf(StrParam,"%d %d %s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", realID, srcID, srcName, iRand, z, M1, q, lamS, betS, phL, thL, phi0, tc, a1, a2, thS1, thS2, phS1, phS2, per, ecc, DL0d*1.e3, Mcz, eta, thBS1, thBS2, phBS1, phBS2);
						}
						
						if(S.gSrcType()==GALBIN)
							sprintf(StrParam,"%d %d %s %d %.8lf %.8lf %.8lf %.10e %lf %.8lf %.8lf %.8lf %.10e %.8lf %.10e %.10e %.10e %.8lf", realID, srcID, srcName, iRand, M1, M2, Period, Pdot, ecc, inc, lamS, betS, DL0, psi, Amp, freq, fdot, phi0);
						
                      
						if(S.gSrcType()==NRWAVE)
							sprintf(StrParam,"%d %d %s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", realID, srcID, srcName, iRand, z, M1, q, lamS, betS, phL, thL, phi0, tc, a1, a2, thS1, thS2, phS1, phS2, per, ecc, DL0d*1.e3, Mcz, eta, thBS1, thBS2, phBS1, phBS2, thBJ, phBJ);  // changed by Sofiane
                  
						
						Coutm << StrParam << Endl;
						
						
						//! **** Compute only the source which pass the threshold
						if( (S.gSrcType()==GALBIN) || (S.gSrcType()==NRWAVE)  || ((S.gSrcType()==SPINBBHHHARM)&&(q>qth)) ){
							
							S.NotValidPar = false;
							
							time(&tstart);
							
							S.Reset(true);
							if(S.gSrcType()==SPINBBHHHARM){
								if(NRecCat>0){
									
									S.sP(0, betS);
									S.sP(1, lamS);
									S.sP(14, phL);
									S.sP(13, thL);
									S.sP(12, phi0);
									S.sP(4, tc);
									S.sP(5, DL0d*1.e3);
									S.sP(15, Mcz);  
									S.sP(17, eta);
									if(NRecCat>15){
										S.sP(6, a1);
										S.sP(7, a2);
										S.sP(8, thBS1);
										S.sP(9, thBS2);
										S.sP(10, phBS1);
										S.sP(11, phBS2);
                                        if(ExtraParamJ){
                                            S.sP(21, thBJ);    // added by Sofiane
                                            S.sP(22, phBJ);    // added by Sofiane
                                        }
									}
									
									/*
									//S.sP(0, betS);
									//S.sP(1, lamS);
									//S.sP(14, phL);
									//S.sP(13, thL);
									//S.sP(12, phi0);
									S.sP(4, tc);
									S.sP(5, DL0d*1.e3);
									S.sP(15, Mcz);
									S.sP(17, eta);
									if(NRecCat>15){
										S.sP(6, a1);
										S.sP(7, a2);
										//S.sP(8, thBS1);
										//S.sP(9, thBS2);
										//S.sP(10, phBS1);
										//S.sP(11, phBS2);
									}*/
									 
								}else{
									S.sP(2, M1);
									S.sP(3, M1*q);
								}
							}
							
							if(S.gSrcType()==GALBIN){
								S.sP(8, M1);
								S.sP(9, M2);
								S.sP(5, 2./Period);
								S.sP(12, Pdot);
								S.sP(1, lamSd*M_PI/180.);
								S.sP(0, betSd*M_PI/180.);
								S.sP(10, DL0);
								if(!MT.deq(incd,60.))
									S.sP(3, incd*M_PI/180.);
							}
							
							if(S.gSrcType()==NRWAVE){
								if(NRecCat>0){
									S.sP(0, betS);
									S.sP(1, lamS);
									S.sP(3, M1*(1+q)*(1+z));
									S.sP(4, tc);
									S.sP(5, DL0d*1.e3);
									S.sP(16, thBS1);
									S.sP(17, thBS2);
									S.sP(18, phBS1);
									S.sP(19, phBS2);
                                    S.sP(6, thBJ);
                                    S.sP(7, phBJ);
                                    S.sP(8, 1./q);
                                    S.sP(9, a1);
                                    S.sP(10, a2);
									S.sP(20, thL);   // added by Sofiane
									S.sP(21, phL);   // added by Sofiane
								}
							}
							
							
							S.DispDetails = true;
							if(CheckOutput){
								S.CheckOut = true;
							}
							sprintf(S.BaseNameOut, "%s-%d-%d", bNOut, iSrc, iRand);
							
							double SNR4(S.gSNR4());
							for(int iP=0; iP<S.gNp(); iP++)
								Coutm << S.gP(iP) << " " ;
							Coutm << Endl << "\t ==>";
							for(int i=0; i<S.gNTDI(); i++)
								Coutm << " SNR_" << S.gTDIName(i) << " = " <<  S.gSNR(i);
							Coutm << " ==> SNR_4links = " << S.gSNR4() << "   and  SNR_6links = " << S.gSNR6() << Endl;
							
							
							//! *** Recover the parameters
							if(S.gSrcType()==SPINBBHHHARM){
								betS = S.gP(0);
								lamS = S.gP(1);
								phL = S.gP(14);
								thL = S.gP(13);
								phi0 = S.gP(12);
								tc = S.gP(4);
								DL0 = S.gP(5);
								Mcz = S.gP(15);
								eta = S.gP(17);
								M1 = S.gP(2)/(1.+z);
								q = (1-sqrt(1-4*eta)-2.*eta)/(2.*eta);
								a1 = S.gP(6);
								a2 = S.gP(7);
								thBS1 = S.gP(8);
								thBS2 = S.gP(9);
								phBS1 = S.gP(10);
								phBS2 = S.gP(11);
                                thBJ = S.gP(21);    // added by Sofiane
                                phBJ = S.gP(22);    // added by Sofiane
								
								sprintf(StrParam,"%d %d %s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", realID, srcID, srcName, iRand, z, M1, q, lamS, betS, phL, thL, phi0, tc, a1, a2, thS1, thS2, phS1, phS2, per, ecc, DL0, Mcz, eta, thBS1, thBS2, phBS1, phBS2, thBJ, phBJ);  // changed by Sofiane
							}
							
							if(S.gSrcType()==GALBIN){
								M1 = S.gP(8);
								M2 = S.gP(9);
								//Period = 1./(S.gP(5));
								//Pdot = S.gP(12);
								inc = S.gP(3);
								lamS = S.gP(1);
								betS = S.gP(0);
								DL0 = S.gP(10);
								
								psi = S.gP(2);
								Amp = S.gP(4);
								freq = S.gP(5);
								fdot = S.gP(6);
								phi0 = S.gP(7);
								ecc = S.gP(11);
								
								sprintf(StrParam,"%d %d %s %d %.8lf %.8lf %.8lf %.10e %lf %.8lf %.8lf %.8lf %.10e %.8lf %.10e %.10e %.10e %.8lf", realID, srcID, srcName, iRand, M1, M2, Period, Pdot, ecc, inc, lamS, betS, DL0, psi, Amp, freq, fdot, phi0);
							}
							
							if(S.gSrcType()==NRWAVE){
								double Pol, Mz; // Pol, Thd, Phd, Mz;
								betS = S.gP(0);
								lamS = S.gP(1);
								Pol = S.gP(2);
								Mz = S.gP(3);
								tc = S.gP(4);
								DL0 = S.gP(5);
                                thBJ = S.gP(6);       //Thd = S.gP(6);  // changed by Sofiane
								phBJ = S.gP(7);       // Phd = S.gP(7); // changed by Sofiane
								q = 1./S.gP(8);
								a1 = S.gP(9);
								a2 = S.gP(10);
								thS1 = S.gP(11);
								thS2 = S.gP(12);
								phS1 = S.gP(13);
								phS2 = S.gP(14);
								phi0 = S.gP(15);
								thBS1 = S.gP(16);
								thBS2 = S.gP(17);
								phBS1 = S.gP(18);
								phBS2 = S.gP(19);
								thL	= S.gP(20);
								phL	= S.gP(21);
								
								
								z = MT.redshift(DL0, 70., 0.3, 0.7, 0.);
								eta = q/((1.+q)*(1.+q)); 
								Mcz = Mz*pow(eta,0.6);
								M1 = Mz/((1+q)*(1+z));
								
								sprintf(StrParam,"%d %d %s %d %lf %lf %lf %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f", realID, srcID, srcName, iRand, z, M1, q, lamS, betS, phL, thL, phi0, tc, a1, a2, thS1, thS2, phS1, phS2, per, ecc, DL0, Mcz, eta, thBS1, thBS2, phBS1, phBS2, thBJ, phBJ);  // changed by Sofiane
							}
							
							
							
							time(&tendSNR);
							
							
							MT.MemDisplay();
							
							//! ******* Compute the FIM if SNR threshold
							
							if(S.gSNR6() > SNRth){
								
								//S.DispDetails = true;
								//S.CheckOut = true;
								
								
								Coutm << Endl << ">>>>>>>>>>>>>>>>>>>>>>>> Source selected : Computation of FIM and errors ..." << Endl;
								
								
								for(int iConv=0; iConv<NConvStep; iConv++){
								
									fOutFIM << StrParam;
									
									if(ConviP>=0){
										S.Reset(false);
										S.DeltaParSpecInd = ConviP;
										S.DeltaParSpecVal = pow(10., log10(ConvDMin) + iConv*ConvDLogStep);
										fOutFIM << " " << S.DeltaParSpecVal;
									}
									
									
									//S.CheckOut = true;
									
									FIM4 = S.gFIM4();
									
									Coutm << "============================ Fisher Information Matrix for 4 links ============================" << Endl;
									Coutm << FIM4 << Endl;
									//! ** Mathematica format 
									Coutm << "FIM = ";
									for(int iP1=0; iP1<S.gNParStudy(); iP1++){
										if(iP1==0)
											Coutm << "{{";
										else 
											Coutm << ",{";
										for(int iP2=0; iP2<S.gNParStudy(); iP2++){
											if(iP2==0)
												Coutm <<  FIM4(iP1,iP2);
											else
												Coutm << "," <<  FIM4(iP1,iP2);
										}
										Coutm << "}";
									}
									Coutm << "}" << Endl;
									
									
									FIM6 = S.gFIM6();
									
									Coutm << "============================ Fisher Information Matrix for 6 links ============================" << Endl;
									Coutm << FIM6 << Endl;
									//! ** Mathematica format 
									Coutm << "FIM = ";
									for(int iP1=0; iP1<S.gNParStudy(); iP1++){
										if(iP1==0)
											Coutm << "{{";
										else 
											Coutm << ",{";
										for(int iP2=0; iP2<S.gNParStudy(); iP2++){
											if(iP2==0)
												Coutm <<  FIM6(iP1,iP2);
											else
												Coutm << "," <<  FIM6(iP1,iP2);
										}
										Coutm << "}";
									}
									Coutm << "}" << Endl;
									
									
									CovM4 = S.gCovM4();
									CovM6 = S.gCovM6();
									
									
									if(!S.gFIMnonan())
										std::cerr << "WARINING : PROBLEM OF NAN IN FIM OF SOURCE " << iSrc << " : src " << srcID << " of realization " << realID << " =====================" << Endl;
									
                                    
                                    for(int i=0; i<S.gNTDI(); i++)
                                        fOutFIM << " " << S.gSNR(i);
                                    fOutFIM << " " << S.gSNR4() << " " << S.gSNR6();
                                    
                                    double ErrMajAxis, ErrMinAxis;
                                    double ErrArea( S.gErrSky(4, ErrMajAxis, ErrMinAxis) );
                                    for(int iP=0; iP<S.gNParStudy(); iP++)
                                        fOutFIM << " " << S.gErr4(iP);
                                    fOutFIM << " " << ErrArea << " " << ErrMajAxis << " " << ErrMinAxis;
                                    for(int iP=0; iP<S.gNParStudy(); iP++)
                                        fOutFIM << " " << S.gErr6(iP);
                                    ErrArea = S.gErrSky(6, ErrMajAxis, ErrMinAxis) ;
                                    fOutFIM << " " << ErrArea << " " << ErrMajAxis << " " << ErrMinAxis;
                                    
                                    
									//! **** Write FIM and CovM
									for(int iP=0; iP<S.gNParStudy(); iP++)
										for(int jP=0; jP<S.gNParStudy(); jP++)
											fOutFIM << " " << FIM4(iP,jP);
									for(int iP=0; iP<S.gNParStudy(); iP++)
										for(int jP=0; jP<S.gNParStudy(); jP++)
											fOutFIM << " " << CovM4(iP,jP);
									
									//! **** Write FIM and CovM
									for(int iP=0; iP<S.gNParStudy(); iP++)
										for(int jP=0; jP<S.gNParStudy(); jP++)
											fOutFIM << " " << FIM6(iP,jP);
									for(int iP=0; iP<S.gNParStudy(); iP++)
										for(int jP=0; jP<S.gNParStudy(); jP++)
											fOutFIM << " " << CovM6(iP,jP);
									
									fOutFIM << Endl;
								}
								
							}
							//! ******* End of FIM part
							
							time(&tend);
							
							Coutm << "Time usage : tot = " << (tend-tstart) << " s" << Endl;
							Coutm << "\tTU:\t SNR               : " << 100.*(tendSNR-tstart)/(tend-tstart) << " %  ( " << (tendSNR-tstart) << " s)" << Endl;
							Coutm << "\tTU:\t FIM               : " << 100.*(tend-tendSNR)/(tend-tstart) << " %  ( " << (tend-tendSNR) << " s)" << Endl;
							
							
						}else{
							S.NotValidPar = true;
							Coutm << ">>> Rejected because q = " << q << " < qth = " << qth << Endl; 
						}
						
						fOutMain << StrParam;
						
						for(int i=0; i<S.gNTDI(); i++)
							fOutMain << " " << S.gSNR(i);
						fOutMain << " " << S.gSNR4() << " " << S.gSNR6();
						if(S.gSNR6() > SNRth){
							double ErrMajAxis, ErrMinAxis;
							double ErrArea( S.gErrSky(4, ErrMajAxis, ErrMinAxis) );
							for(int iP=0; iP<S.gNParStudy(); iP++)
								fOutMain << " " << S.gErr4(iP);
							fOutMain << " " << ErrArea << " " << ErrMajAxis << " " << ErrMinAxis;
							for(int iP=0; iP<S.gNParStudy(); iP++)
								fOutMain << " " << S.gErr6(iP);
							ErrArea = S.gErrSky(6, ErrMajAxis, ErrMinAxis) ;
							fOutMain << " " << ErrArea << " " << ErrMajAxis << " " << ErrMinAxis;
						}else {
							for(int iP=0; iP<2*(S.gNParStudy()+3); iP++)
								fOutMain << " " << 0.;
						}
						
						fOutMain << Endl;
					}else{
						//! ** Not increase the index of source
						iSrc--;
					}
				}
				
			}
			
			//! Increase the index of source
			iSrc++;
			if(iSrc>=Nsrc)
				CContinue = false;
			
		}
		
		fOutMain.close();
		fOutFIM.close();
		
		if(NRecCat!=-1){
			fInCatalogue.close();
			fInCatalogue.clear();
		}
		
		for(int iW=0; iW<Words.size(); iW++)
			MT.Free(Words[iW],256*sizeof(char));
		
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "LC2SNRFIMMulti: error: " << e.what()<<Endl;
		std::cerr << "LC2SNRFIMMulti: abort!" << Endl;
		exit(1);
	}
	return(0);
};



/** \}*/

// end of LISACODE-SNR.cpp
