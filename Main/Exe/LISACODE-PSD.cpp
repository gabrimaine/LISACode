// $Id:  $
/*
 *  LISACODE-PSD.h
 *  V 2.0
 *
 *  Created on 26/04/2005 by  Antoine Petiteau (AEI)
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


/** \brief The executable compute the power spectral density of TDI output 
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
		
		
		
		// *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution :" << Endl;
			Coutm << "\t\t(./)LC2PSD [Options] ConfigFile1.xml ConfigFile2.xml ..." << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << Endl << "\t\t * ConfigFileI.xml (required) : Configuration files. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -o %sOutBaseName  : Directory containing noise files [required or option -nf]"  << Endl ;
			Coutm << "\t\t * -tdi %sObsName : Names of TDI generators [default : X,Y,Z]"  << Endl ;
			Coutm << "\t\t * -r %dNbRun : Number of run [default : 10]"  << Endl ;
			Coutm << "\t\t * -t0 %ft0 : Initial time [default : 0.]"  << Endl ;
			Coutm << "\t\t * -T %fT   : Duration [default : 1048576.0 (=2^20)]"  << Endl ;
			Coutm << "\t\t * -s %dseed : Seed for random gennerator. [default : current time]"  << Endl ;
			Coutm << "\t\t\t  NB: If not specified the random seed is the time machine." << Endl;
			Coutm << "\t\t * -tdi %sObsName : Names of TDI generators [default : X,Am,Em]"  << Endl ;
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
			Coutm << " LC2PSD : Compute the power spectral density - LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		
		// ***** Main declaration *****
		LISACode * Sim(NULL);
		long SeedRand((long)time(NULL));
		double dtMes(0.3);
		double t0(0.);
		double tDur(1048576.0);
		int NData(-1), tmpNData(-1);
		int NDatPSD(0);
		double dfPSD(1.);
		
		int NRun(10);
		
		
		//! ** Base name of outpu
		char fNOut[16384];
		strcpy(fNOut,"PSD-Test.txt");
		std::ofstream fOut;
		
		//! List of names of files describing the configuration
		char ** fNConfig(NULL);
		int NfNConfig(0);
		
		//! Bloc list of TDI observables
		double ** tTDI(NULL);
		char TDINameAll[16384];
		char ** TDIName(NULL);
		int NTDI(0);
		strcpy(TDINameAll, "X,Am,Em");
		
		
		
		// ***** Read options *****
		for(int iarg=1; iarg<argc; iarg++){
			
			if((argc>1)&&(strcmp(argv[iarg],"-o")==0)){
				// ***** Read base name of output file in option
				strcpy(fNOut, argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-tdi")==0)){
				strcpy(TDINameAll, argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-r")==0)){
				NRun = atoi(argv[iarg+1]);
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
		
		
		// ***** Read configuration files in arguments
		
		for(int ia=1+nOptions; ia<argc; ia++){
			NfNConfig++;
			fNConfig = (char**)MT.ReAllocMemory(fNConfig, (NfNConfig-1)*sizeof(char*), NfNConfig*sizeof(char*));
			fNConfig[NfNConfig-1] = (char*)MT.AllocMemory(2048*sizeof(char));
			strcpy(fNConfig[NfNConfig-1], argv[ia]);
		}
		
		
		
		// ***** Display for checking which files we are using
		Coutm << "Configuration files : " << Endl;
		for(int i=0; i<NfNConfig; i++){
			Coutm << "\t + " << fNConfig[i] ; MT.CheckFile("", fNConfig[i]);
		}
		
		
		
		double ** PSDTmp(NULL);
		double ** PSDRes(NULL);
		char fNOutTmp[16384];
		sprintf(fNOutTmp, "%s.tmp", fNOut);
		
		//! ********* Compute the response to gravitational waves
		
		for(int iRun=0; iRun<NRun; iRun++){
			std::cerr << Endl << ">>>>>>> Run " << iRun << Endl;
			Coutm << Endl << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Compute the power spectral density : Run " << iRun << " ..." << Endl;
			
			Sim  = new LISACode(&MT);
			Sim->setTimeInfo(dtMes, t0, tDur);
			
			if(iRun==0){
				//! *** Read the bloc name of TDI generators and obtain the individual names and the number of generators
				int ipTDINameAll(0);
				MT.wextractcoma(TDINameAll, ipTDINameAll, TDIName, NTDI,16);
				Coutm << "TDI generators : ";
				for(int i=0; i<NTDI; i++)
					Coutm << " " << TDIName[i];
				Coutm << Endl;
				
				
				//! *** Allocate memory
				tTDI = (double**) MT.AllocMemory(NTDI*sizeof(double*));
				PSDRes = (double**) MT.AllocMemory(NTDI*sizeof(double*));
				PSDTmp = (double**) MT.AllocMemory(NTDI*sizeof(double*));
				for(int i=0; i<NTDI; i++){
					tTDI[i] = NULL;
					PSDRes[i] = NULL;
					PSDTmp[i] = NULL;
				}
			}
			
			//! *** Configure the simulator using the configuration files
			for(int i=0; i<NfNConfig; i++)
				Sim->config(fNConfig[i]);
			
			//! *** Update time informations
			if(iRun == 0){
				t0 = Sim->gett0();
				dtMes = Sim->getdtMes();
			}else{
				if(!MT.deq(t0,Sim->gett0()))
					throw std::invalid_argument("ERROR : The initial time changes !");
				if(!MT.deq(dtMes,Sim->getdtMes()))
					throw std::invalid_argument("ERROR : The time step changes !");
			}
			
			//! *** Initialization of simulation
			Sim->init();
			
			//! *** Add output series with allocation of the series if needed
			for(int i=0; i<NTDI; i++){
				Sim->AddSerieOut(TDIName[i], 1, tmpNData, (iRun == 0));
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
			
			
			//! ****** Compute the power spectral density
			
			for(int iTDI=0; iTDI<NTDI; iTDI++)
				MT.PSD(tTDI[iTDI], NData, dtMes, 4, 10, PSDTmp[iTDI], NDatPSD, dfPSD, (iRun == 0));
			
			//! *** For the first run initialize the result PSD 
			if(iRun == 0){
				for(int iTDI=0; iTDI<NTDI; iTDI++){
					PSDRes[iTDI] = (double*) MT.AllocMemory(NDatPSD*sizeof(double));
					for(int i=0; i<NDatPSD; i++)
						PSDRes[iTDI][i] = 0.; 
				}
			}
			
			//! *** Add the PSD of the run to the result PSD
			for(int iTDI=0; iTDI<NTDI; iTDI++)
				for(int i=0; i<NDatPSD; i++)
					PSDRes[iTDI][i] += PSDTmp[iTDI][i];
			
			//! *** Temporary record
			fOut.open(fNOutTmp);
			for(int iF=0; iF<NDatPSD; iF++){
				fOut << iF*dfPSD;
				for(int iTDI=0; iTDI<NTDI; iTDI++)
					fOut << " " << PSDRes[iTDI][iF]/(iRun+1);
				fOut << Endl;
			}
			fOut.close();
			fOut.clear();
			
		}
		
		//! *** Normalize the result
		for(int iTDI=0; iTDI<NTDI; iTDI++)
			for(int i=0; i<NDatPSD; i++)
				PSDRes[iTDI][i] /= NRun;
		
		sprintf(cmd, "rm %s", fNOutTmp);
		system(cmd);
		
		//! ***** Record the ressults
		Coutm << "Record the result PSD in " << fNOut << " ..." << Endl;
		fOut.open(fNOut);
		fOut << "#f";
		for(int iTDI=0; iTDI<NTDI; iTDI++)
			fOut << " " << TDIName[iTDI];
		fOut << Endl;
		for(int iF=0; iF<NDatPSD; iF++){
			fOut << iF*dfPSD;
			for(int iTDI=0; iTDI<NTDI; iTDI++)
				fOut << " " << PSDRes[iTDI][iF];
			fOut << Endl;
		}
		fOut.close();
		fOut.clear();
		
		
		
		
		// ***** Free memory *****
		
		if(Sim != NULL)
			delete Sim;
		Sim = NULL;
		
		
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
		
		if(PSDRes != NULL){
			for(int i=0; i<NTDI; i++)
				if(PSDRes[i] != NULL)
				MT.Free(PSDRes[i], NDatPSD*sizeof(double));
			MT.Free(PSDRes, NTDI*sizeof(double*));
		}
		
		if(PSDTmp != NULL){
			for(int i=0; i<NTDI; i++)
				if(PSDTmp[i] != NULL)
					MT.Free(PSDTmp[i], NDatPSD*sizeof(double));
			MT.Free(PSDTmp, NTDI*sizeof(double*));
		}
			
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "LC2PSD: error: " << e.what()<<Endl;
		std::cerr << "LC2PSD: abort!" << Endl;
		exit(1);
	}
	return(0);
};


void CheckFile(LCTools * MT, char * fDir, char * fName)
{
	std::ifstream fIn;
	
	char TmpFullName[16384];
	if(MT->wcmp(fDir,""))
		strcpy(TmpFullName, fName);
	else 
		sprintf(TmpFullName, "%s/%s", fDir, fName);
	
	fIn.open(TmpFullName);
	if(fIn == NULL){
		Cout << Endl << "ERROR: Can not open the file " << TmpFullName << " !" << Endl;
		throw std::invalid_argument("ERROR: Can not open a file.");
	}
	fIn.close();
	fIn.clear();
	Cout << " --> OK" << Endl;
}



/** \}*/

// end of LISACODE-PSD.cpp
