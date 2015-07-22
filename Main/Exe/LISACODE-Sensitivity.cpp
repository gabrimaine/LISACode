// $Id:  $
/*
 *  LISACODE-Sensitivity.h
 *  V 2.0
 *
 *  Created on 25/04/2005 by  Antoine Petiteau (AEI)
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
#include "LISACODE-LISACode.h"


/** \ingroup main Main 
 * \{
 */


/** \brief The executable compute the sensitivity curve and by product the GW and noise response
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

		
		
		
		// *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution :" << Endl;
			Coutm << "\t\t(./)LC2Sensitivity [Options] ConfigFile1.xml ConfigFile2.xml ..." << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << Endl << "\t\t * ConfigFileI.xml (required) : Configuration files. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -nf %sPathNoiseFile : Noise file [required or option -nd]"  << Endl ;
			Coutm << "\t\t * -nd %sPathNoiseDir  : Directory containing noise files [required or option -nf]"  << Endl ;
			Coutm << "\t\t * -gf %sPathGWStoch192file : File used for describing an isotropic source distribution [default : $LISACODEDIR/ConfigFiles/Refs/Stoch-192.xml ]"  << Endl ;
			Coutm << "\t\t * -G \t\t: Read response to GW in file. [default: false]"  << Endl ;
			Coutm << "\t\t * -o %sOutBaseName  : Directory containing noise files [required or option -nf]"  << Endl ;
			Coutm << "\t\t * -tdi %sObsName : Names of TDI generators [default : X,Am,Em]"  << Endl ;
			Coutm << "\t\t * -t0 %ft0 : Initial time [default : 0.]"  << Endl ;
			Coutm << "\t\t * -it0 %dit0 : Initial time = it0 * duration [default : 0]"  << Endl ;
			Coutm << "\t\t * -T %fT   : Duration [default : 1048576.0 (=2^20)]"  << Endl ;
			Coutm << "\t\t * -D \t: Record delays. [default: false]"  << Endl ;
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
			Coutm << " LC2Sensitivity : Compute the sensitivity curve - LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		
		// ***** Main declaration *****
		LISACode * Sim(NULL);
		long SeedRand((long)time(NULL));
		double dtMes(0.3);
		double t0(0.);
		double tDur(1048576.0);
		int it0(0);
		bool Deft0(false), DeftDur(false);
		bool ReadGWResp(false);
		
		int NData(-1), tmpNData(-1);
		int NDatPSD(0);
		double dfPSD(1.);
		int PSDNMean(10);
		
		char LCDir[16384];
		char NoiseDir[16384];
		std::ofstream fOut;
		
		//! ** Base name of outpu
		char bNOut[16384];
		strcpy(bNOut,"Sens");
		
		//! File containing the gravitational wave sources
		char fNGWs[16384];
		strcpy(fNGWs, "None");
		
		//! List of names of files describing the noises 
		char ** fNNoises(NULL);
		int NfNNoises(0);
		
		//! List of names of files describing the configuration
		char ** fNConfig(NULL);
		int NfNConfig(0);
		
		//! Bloc list of TDI observables
		double ** tTDI(NULL);
		char TDINameAll[16384];
		char ** TDIName(NULL);
		int NTDI(0);
		strcpy(TDINameAll, "X,Am,Em");
		
		double ** tDelay(NULL);
		int NDataDelay;
		bool RecordDelay(false);
		
		double ** PSDGW(NULL);
		
		
		// ***** Read options *****
		for(int iarg=1; iarg<argc; iarg++){
			
			if((argc>1)&&(strcmp(argv[iarg],"-nf")==0)&&(fNNoises==NULL)){
				
				// ***** Read one noise file in option
				NfNNoises++;
				fNNoises = (char**)MT.AllocMemory(NfNNoises*sizeof(char*));
				fNNoises[0] = (char*)MT.AllocMemory(2048*sizeof(char));
				strcpy(fNNoises[0], argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-nd")==0)&&(fNNoises==NULL)){
				
				// ***** Read several noise files in directory specified in option
				struct dirent **listfiles;
				int n;
				strcpy(NoiseDir, argv[iarg+1]);
				n = scandir(NoiseDir,&listfiles,0,alphasort);
				if(n<0)
					throw std::invalid_argument("Problem to open the noise directory !");
				for(int i=0; i<n; i++){
					if(listfiles[i]->d_name[0] != '.'){
						NfNNoises++;
						fNNoises = (char**)MT.ReAllocMemory(fNNoises, (NfNNoises-1)*sizeof(char*), NfNNoises*sizeof(char*));
						fNNoises[NfNNoises-1] = (char*)MT.AllocMemory(2048*sizeof(char));
						strcpy(fNNoises[NfNNoises-1], listfiles[i]->d_name);
					}
				}
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-gf")==0)){
				// ***** Read gravitational waves file in option
				strcpy(fNGWs, argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-G")==0)){
				ReadGWResp = true;
				nOptions ++;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-o")==0)){
				// ***** Read base name of output file in option
				strcpy(bNOut, argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-tdi")==0)){
				// ***** Read gravitational waves file in option
				strcpy(TDINameAll, argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-t0")==0)){
				t0 = atof(argv[iarg+1]); 
				Deft0 = true;
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-it0")==0)){
				it0 = atof(argv[iarg+1]); 
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-T")==0)){
				tDur = atof(argv[iarg+1]); 
				DeftDur = true;
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-D")==0)){
				RecordDelay = true;
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
		if(it0!=0){
			t0 = it0*tDur;
			Deft0 = true;
		}
		MT.setRandSeed(SeedRand);
		
		
		// ***** If the gravitational waves configuration file is still None then try to load the default one from $LISACODEDIR/ConfigFiles/Refs/Stoch-192.xml
		if(strcmp(fNGWs,"None")==0){
			if (getenv("LISACODEDIR") != NULL){
				strcpy(LCDir,getenv("LISACODEDIR"));
				sprintf(fNGWs, "%s/ConfigFiles/Refs/Stoch-192.xml", LCDir);
			}else{
				throw std::invalid_argument("No configuration file for gravitational waves !");
			}
		}
		
		
		// ***** Read configuration files in arguments
		
		for(int ia=1+nOptions; ia<argc; ia++){
			NfNConfig++;
			fNConfig = (char**)MT.ReAllocMemory(fNConfig, (NfNConfig-1)*sizeof(char*), NfNConfig*sizeof(char*));
			fNConfig[NfNConfig-1] = (char*)MT.AllocMemory(2048*sizeof(char));
			strcpy(fNConfig[NfNConfig-1], argv[ia]);
		}
		
		
		
		// ***** Display for checking which files we are using
		Coutm << "Configuration files : " << Endl;
		Coutm << "\t - gravitational wave : " << fNGWs; MT.CheckFile("", fNGWs);
		Coutm << "\t - noises (in " << NoiseDir << " ) : " << Endl;
		for(int i=0; i<NfNNoises; i++){
			Coutm << "\t\t + " << fNNoises[i]; MT.CheckFile(NoiseDir, fNNoises[i]);
		}
		Coutm << "\t - detector and the rest : " << Endl;
		for(int i=0; i<NfNConfig; i++){
			Coutm << "\t\t + " << fNConfig[i] ; MT.CheckFile("", fNConfig[i]);
		}

		
		
		
		//! ********* Obtain the response to gravitational waves
		
		Coutm << Endl << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Compute the response to gravitational waves ..." << Endl;
		
		
		//! *** Read the bloc name of TDI generators and obtain the individual names and the number of generators
		int ipTDINameAll(0);
		MT.wextractcoma(TDINameAll, ipTDINameAll, TDIName, NTDI,16);
		Coutm << "TDI generators : ";
		for(int i=0; i<NTDI; i++)
			Coutm << " " << TDIName[i];
		Coutm << Endl;
		
		
		//! *** Allocate memory for the list of series
		tTDI = (double**) MT.AllocMemory(NTDI*sizeof(double*));
		PSDGW = (double**) MT.AllocMemory(NTDI*sizeof(double*));
		for(int i=0; i<NTDI; i++){
			tTDI[i] = NULL;
			PSDGW[i] = NULL;
		}
		
		//! *** Name of GW response file
		char fNRespGW[16384];
		sprintf(fNRespGW, "%s-GWs.txt", bNOut);
		
		if(ReadGWResp){
			
			//! *** Read in file
			LCDataFileRead fInRespGW(&MT, fNRespGW);
			fInRespGW.init();
			NDatPSD = fInRespGW.getNDat();
			dfPSD = fInRespGW.getdx();
			for(int iTDI=0; iTDI<NTDI; iTDI++){
				PSDGW[iTDI] = (double*) MT.AllocMemory(NDatPSD*sizeof(double));
				for(int i=0; i<NDatPSD; i++){
					PSDGW[iTDI][i] = fInRespGW.gDataBin(iTDI, i);
				}
			}
			
		}else{
			
			//! *** Compute the gravitational wave response
			Sim  = new LISACode(&MT);
			Sim->setTimeInfo(dtMes, t0, tDur);
			
			//! *** Configure the simulator using the GW file and the configuration files
			Sim->config(fNGWs);
			for(int i=0; i<NfNConfig; i++)
				Sim->config(fNConfig[i]);
			
			//! *** Set and/or update time informations
			if(Deft0)
				Sim->sett0(t0);
			if(DeftDur)
				Sim->settDur(tDur);
			
			t0 = Sim->gett0();
			dtMes = Sim->getdtMes();
			tDur = Sim->gettDur();
			
			std::cerr << "Time informations : t0 = " << t0 << "s  , dt = " << dtMes << "s  , duration = " << tDur << " s" << Endl;
			
			//! *** Initialization of simulation
			Sim->init();
			
			//! *** Add output series for TDI with allocation of the series
			for(int i=0; i<NTDI; i++){
				Sim->AddSerieOut(TDIName[i], 1, tmpNData, true);
				if(NData == -1)
					NData = tmpNData;
				if(NData != tmpNData){
					Coutm << "ERROR : The numbers of data in TDI time vector are not the same : previously it was " << NData << "but it's " << tmpNData << "for the last one (" << TDIName[i] << ") !" << Endl;
					throw std::invalid_argument("ERROR : The numbers of data in TDI time vector are not the same !");
				}
			}
			//! *** Link output series
			for(int i=0; i<NTDI; i++)
				Sim->LinkSerieOut(i, tTDI[i]);
			
			
			//! *** Add and link output series  for delay with allocation of the series
			if(RecordDelay){
				Sim->AddSerieOut("D1", 2, NDataDelay, true);
				Sim->AddSerieOut("D2", 2, NDataDelay, true);
				Sim->AddSerieOut("D3", 2, NDataDelay, true);
				tDelay = (double**)MT.AllocMemory(3*sizeof(double*));
				for(int i=0; i<3; i++)
					Sim->LinkSerieOut(NTDI+i, tDelay[i]);
			}
			
			Sim->DispInfo("");
			
			//! *** Compute the time series : Run the simulation
			Sim->Run();
			/*
			 if(tTDI != NULL){
			 if(tTDI[0] != NULL){
			 std::ofstream fTCheck("CheckTime_TDI.txt");
			 for(int iT=0; iT<NData; iT++){
			 fTCheck << iT*dtMes;
			 for(int i=0; i<NTDI; i++)
			 fTCheck << " " << tTDI[i][iT];
			 fTCheck << Endl;
			 }
			 }
			 }
			 */
			
			//! *** Delete the simulation
			delete Sim;
			Sim = NULL;
			
			//! *** Record delay
			if(RecordDelay){
				char fNOutDelay[16000];
				sprintf(fNOutDelay, "%s-Delay.txt", bNOut);
				std::ofstream fOutDelay(fNOutDelay);
				fOutDelay << "#t D1 D2 D3" << Endl;
				fOutDelay.precision(10);
				for(int i=0; i<NDataDelay; i++){
					fOutDelay << t0+i*dtMes;
					for(int iR=0; iR<3; iR++)
						fOutDelay << " " << tDelay[iR][i];
					fOutDelay << Endl;
				}
				fOutDelay.close();
				
			}
			
			//! ****** Compute the power spectral density
			for(int iTDI=0; iTDI<NTDI; iTDI++)
				MT.PSD(tTDI[iTDI], NData, dtMes, 4, PSDNMean, PSDGW[iTDI], NDatPSD, dfPSD, true);
			
			
			//! ***** Record the response to gravitational waves
			Coutm << "Record gravitational wave response in " << fNRespGW << " ..." << Endl;
			fOut.open(fNRespGW);
			fOut << "#f";
			for(int iTDI=0; iTDI<NTDI; iTDI++)
				fOut << " " << TDIName[iTDI];
			fOut << Endl;
			fOut.precision(10);
			for(int iF=0; iF<NDatPSD; iF++){
				fOut << iF*dfPSD;
				for(int iTDI=0; iTDI<NTDI; iTDI++)
					fOut << " " << PSDGW[iTDI][iF];
				fOut << Endl;
			}
			fOut.close();
			fOut.clear();
			
		}
		
		/*
		 std::ofstream fPCheck("CheckPSD_TDI.txt");
		 for(int iF=0; iF<NDatPSD; iF++){
		 fPCheck << iF*dfPSD;
		 for(int i=0; i<NTDI; i++)
		 fPCheck << " " << PSDGW[i][iF];
		 fPCheck << Endl;
		 }
		 */
		
		//! ********* Compute the noise response
		char fNTmpNoise[16384];
		double ** PSDNoise(NULL);
		char fNOutNoise[16384];
		char fNTmp[16384];
		int ipS(0), ipE(0);
		
		//! ********* Compute the response to each noise individually : Loop on noise configuration file
		
		for(int iN=0; iN<NfNNoises; iN++){
			
			Coutm << Endl << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Compute the response to individual noise " << fNNoises[iN] << " ..." << Endl;
			
			Sim  = new LISACode(&MT);
			
			//! *** Configure the simulator using the noise file and the configuration files
			sprintf(fNTmpNoise, "%s/%s", NoiseDir, fNNoises[iN]);
			Sim->config(fNTmpNoise);
			for(int i=0; i<NfNConfig; i++)
				Sim->config(fNConfig[i]);
			
			//! *** Initialization of simulation
			Sim->setTimeInfo(dtMes, t0, tDur);
			Sim->init();
			
			//! *** Add output series with allocation of the series
			for(int i=0; i<NTDI; i++){
				Sim->AddSerieOut(TDIName[i], 1, tmpNData, ((ReadGWResp)&&(iN==0)));
				if(NData == -1)
					NData = tmpNData;
				if(NData != tmpNData){
					Coutm << "ERROR : The numbers of data in TDI time vector are not the same : previously it was " << NData << " but it's " << tmpNData << " for the last one (" << TDIName[i] << ") !" << Endl;
					throw std::invalid_argument("ERROR : The numbers of data in TDI time vector are not the same !");
				}
			}
			//! *** Link output series
			for(int i=0; i<NTDI; i++)
				Sim->LinkSerieOut(i, tTDI[i]);
			
			Sim->DispInfo("");
			
			//! *** Compute the time series : Run the simulation
			Sim->Run();
			
			//! *** Delete the simulation
			delete Sim;
			Sim = NULL;
			
			//! ****** Compute the power spectral density with allocation if it's the first noise study
			if(iN==0)
				PSDNoise = (double**) MT.AllocMemory(NTDI*sizeof(double*));
			for(int iTDI=0; iTDI<NTDI; iTDI++){
				int tmpNDatPSD;
				MT.PSD(tTDI[iTDI], NData, dtMes, 4, PSDNMean, PSDNoise[iTDI], tmpNDatPSD, dfPSD, (iN==0));
				if(tmpNDatPSD != NDatPSD){
					Coutm << "ERROR : The numbers of data in PSD time vector are not the same : previously it was " << NDatPSD << " but it's " << tmpNDatPSD << " for the last one (" << TDIName[iTDI] << ") !" << Endl;
					throw std::invalid_argument("ERROR : The numbers of data in TDI time vector are not the same !");
				}
			}
			
			
			//! ***** Record sensitivity and response to noise
			
			//! ** Extract the inportant part of noise configuration file
			ipS = 0;
			ipE = strlen(fNNoises[iN])-1;
			if(strlen(fNNoises[iN])>7){
				if(strncmp(fNNoises[iN], "Config-", 7)==0)
					ipS = 7;
				if(strncmp(fNNoises[iN]+strlen(fNNoises[iN])-4, ".xml", 4)==0)
					ipE = strlen(fNNoises[iN])-4;
			}
			strncpy(fNTmp, fNNoises[iN]+ipS, ipE-ipS);
			fNTmp[ipE-ipS] = '\0';
			
			sprintf(fNOutNoise, "%s-Noise-%s.txt", bNOut, fNTmp);
			Coutm << "Record response to individual noise " << fNNoises[iN] << " in " << fNOutNoise << " ..." << Endl;
			fOut.open(fNOutNoise);
			fOut << "#f";
			for(int iTDI=0; iTDI<NTDI; iTDI++)
				fOut << " Resp_" << TDIName[iTDI];
			for(int iTDI=0; iTDI<NTDI; iTDI++)
				fOut << " Sens_" << TDIName[iTDI];
			fOut << Endl;
			fOut.precision(10);
			for(int iF=0; iF<NDatPSD; iF++){
				fOut << iF*dfPSD;
				for(int iTDI=0; iTDI<NTDI; iTDI++)
					fOut << " " << PSDNoise[iTDI][iF];
				//! *** Here computation of the sensitivity
				for(int iTDI=0; iTDI<NTDI; iTDI++)
					fOut << " " << 5.*sqrt(PSDNoise[iTDI][iF]/(LC::Yr_SI*PSDGW[iTDI][iF]));
				fOut << Endl;
			}
			fOut.close();
			fOut.clear();
			
			
		}
		
		
		//! ********* Compute the response to all the noises together 
		
		Coutm << Endl << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Compute the response to all the noises together ..." << Endl;
		
		Sim  = new LISACode(&MT);
		Sim->setTimeInfo(dtMes, t0, tDur);
		
		//! *** Configure the simulator using each noise files and the configuration files
		for(int iN=0; iN<NfNNoises; iN++){
			sprintf(fNTmpNoise, "%s/%s", NoiseDir, fNNoises[iN]);
			Sim->config(fNTmpNoise);
		}
		for(int i=0; i<NfNConfig; i++)
			Sim->config(fNConfig[i]);
		
		//! *** Initialization of simulation
		Sim->setTimeInfo(dtMes, t0, tDur);
		Sim->init();
		
		//! *** Add output series with allocation of the series
		for(int i=0; i<NTDI; i++){
			Sim->AddSerieOut(TDIName[i], 1, tmpNData, false);
			if(NData == -1)
				NData = tmpNData;
			if(NData != tmpNData){
				Coutm << "ERROR : The numbers of data in TDI time vector are not the same : previously it was " << NData << "but it's " << tmpNData << "for the last one (" << TDIName[i] << ") !" << Endl;
				throw std::invalid_argument("ERROR : The numbers of data in TDI time vector are not the same !");
			}
		}
		//! *** Link output series
		for(int i=0; i<NTDI; i++)
			Sim->LinkSerieOut(i, tTDI[i]);
		
		Sim->DispInfo("");
		
		//! *** Compute the time series : Run the simulation
		Sim->Run();
		
		//! *** Delete the simulation
		delete Sim;
		Sim = NULL;
		
		//! ****** Compute the power spectral density with allocation if it's the first noise study
		for(int iTDI=0; iTDI<NTDI; iTDI++)
			MT.PSD(tTDI[iTDI], NData, dtMes, 4, PSDNMean, PSDNoise[iTDI], NDatPSD, dfPSD, false);
		
		
		//! ***** Record sensitivity and response to noise
		
		sprintf(fNOutNoise, "%s-AllNoises.txt", bNOut, fNTmp);
		Coutm << "Record response to all noises together in " << fNOutNoise << " ..." << Endl;
		fOut.open(fNOutNoise);
		fOut << "#f";
		for(int iTDI=0; iTDI<NTDI; iTDI++)
			fOut << " Resp_" << TDIName[iTDI];
		for(int iTDI=0; iTDI<NTDI; iTDI++)
			fOut << " Sens_" << TDIName[iTDI];
		fOut << Endl;
		fOut.precision(10);
		for(int iF=0; iF<NDatPSD; iF++){
			fOut << iF*dfPSD;
			for(int iTDI=0; iTDI<NTDI; iTDI++)
				fOut << " " << PSDNoise[iTDI][iF];
			//! *** Here computation of the sensitivity
			for(int iTDI=0; iTDI<NTDI; iTDI++)
				fOut << " " << 5.*sqrt(PSDNoise[iTDI][iF]/(LC::Yr_SI*PSDGW[iTDI][iF]));
			fOut << Endl;
		}
		fOut.close();
		fOut.clear();
		
		
		
		// ***** Free memory *****
		
		if(Sim != NULL)
			delete Sim;
		Sim = NULL;
		
		if(fNNoises != NULL){
			for(int iF=0; iF<NfNNoises; iF++)
				if(fNNoises[iF] != NULL)
					MT.Free(fNNoises[iF], 2048*sizeof(char));
			MT.Free(fNNoises, NfNNoises*sizeof(char*));
		}
		
		if(fNConfig != NULL){
			for(int iF=0; iF<NfNConfig; iF++)
				if(fNConfig[iF] != NULL)
					MT.Free(fNConfig[iF], 2048*sizeof(char));
			MT.Free(fNConfig, NfNConfig*sizeof(char*));
		}
		
		if(TDIName != NULL){
			for(int i=0; i<NTDI; i++)
				if(TDIName[i] != NULL)
					MT.Free(TDIName[i], 16*sizeof(char));
			MT.Free(TDIName, NTDI*sizeof(char*));
		}
		
		if(tTDI != NULL){
			for(int i=0; i<NTDI; i++)
				if(tTDI[i] != NULL)
					MT.Free(tTDI[i], NData*sizeof(double));
			MT.Free(tTDI, NTDI*sizeof(double*));
		}
		
		if(tDelay != NULL){
			for(int iR=0; iR<3; iR++)
				if(tDelay[iR] != NULL)
					MT.Free(tDelay[iR], NDataDelay*sizeof(double));
			MT.Free(tDelay, 3*sizeof(double*));
		}
		
		if(PSDGW != NULL){
			for(int i=0; i<NTDI; i++)
				if(PSDGW[i] != NULL)
				MT.Free(PSDGW[i], NDatPSD*sizeof(double));
			MT.Free(PSDGW, NTDI*sizeof(double*));
		}
		
		if(PSDNoise != NULL){
			for(int i=0; i<NTDI; i++)
				if(PSDNoise[i] != NULL)
					MT.Free(PSDNoise[i], NDatPSD*sizeof(double));
			MT.Free(PSDNoise, NTDI*sizeof(double*));
		}
			
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "LC2Sensitivity: error: " << e.what()<<Endl;
		std::cerr << "LC2Sensitivity: abort!" << Endl;
		exit(1);
	}
	return(0);
};


/** \}*/

// end of LISACODE-Sensitivity.cpp
