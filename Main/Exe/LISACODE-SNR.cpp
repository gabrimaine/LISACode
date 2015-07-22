// $Id:  $
/*
 *  LISACODE-SNR.h
 *  V 2.0
 *
 *  Created on 07/05/2005 by  Antoine Petiteau (AEI)
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
#include "LISACODE-ModSig.h"


/** \ingroup main Main 
 * \{
 */


/*! Load confusion noise model in a serie */ 
void LoadConfusionNoise(LCTools * MT, LCSerie2 * & ConfNoise, int KeyModel, double fMin, double fMax, double df);

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
		double dtMes(0.3);
		double t0(0.);
		double tDur(63115200.0);
		int NtData(-1), tmpNtData(-1);
		int NfData(-1);
        bool RecCumul(false);
        bool SnAnalytic(true);
        
        LCModSig S(&MT);
        S.Detector = ELISA;
        S.GalacticNoise = 0;
		
		int NRun(1);
		
		// *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution :" << Endl;
			Coutm << "\t\t(./)LC2SNR [Options] ConfigFile1.xml ConfigFile2.xml ..." << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << Endl << "\t\t * ConfigFileI.xml (required) : Configuration files. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -nf %sPSDFileName : File containing the PSD [default: None ==> analytic]"  << Endl ;
			Coutm << "\t\t * -nfc %sPSDFileName %sCodeColumn : File containing the PSD. Example : -nfc PSDFile.txt f,n,X,Y,Z [default: None ==> analytic]"  << Endl ;
			Coutm << "\t\t * -ca %dKeyModel : Analytic model of confusion noise  [default : None]"  << Endl ;
			Coutm << "\t\t * -o %sOutFile  : Output file containing results [default: SNRResults.txt]"  << Endl ;
			Coutm << "\t\t * -tdi %sObsName : Names of TDI generators [default : X,Am,Em]"  << Endl ;
            Coutm << "\t\t * -D %dConfig : For computing analytic noise : LISA,ELISA,C1,C2,C3,C4,C5,C6,SGOF1,SGOM1. Should not be used with tdi.[default : ELISA]"  << Endl ;
			Coutm << "\t\t * -t0 %ft0 : Initial time [default : " << t0 <<"]"  << Endl ;
			Coutm << "\t\t * -T %fT   : Duration [default : " << tDur  << " (= "  << tDur/LC::Yr_SI << " yr)]"  << Endl ;
            Coutm << "\t\t * -dt %fdt : Time step [default : " << dtMes  << " s ]"  << Endl ;
			Coutm << "\t\t * -pr %diP  : Index of parameters to randomize [default : None]"  << Endl ;
			Coutm << "\t\t * -s %dseed : Seed for random gennerator. [default : current time]"  << Endl ;
			Coutm << "\t\t\t  NB: If not specified the random seed is the time machine." << Endl;
			Coutm << "\t\t * -cu \t\t: Record cumulative SNR in file CumulSNR.txt (with other checking values). [default: false]"  << Endl ;
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
			Coutm << " LC2SNR : Compute the signal to noise ratio (SNR) - LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		
		// ***** Main declaration *****
		LISACode * Sim(NULL);
        double fLowMax(0.01); // WARNING: fixed by hard
		
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
		strcpy(TDINameAll, "X");
		
		// *** List of pointer on read file and index of signal in it. Size NTDI
		LCDataFileRead ** fTDIPSD(NULL);
		int * ifTDIPSD(NULL);
		
		//! *** Time and frequency data of the TDI signal of GW. Size NTDI
		double ** tSig(NULL);
		dcomplex ** fSig(NULL);
		
		//! *** SNR values
		double * SNR;
		double SNRAll;
		double df;
		int iFmin, iFmax;
		dcomplex tmpRes;
		double tmpSP;
		
		//! *** Output (results) file
		char fNOut[16384];
		std::ofstream fOut;
		strcpy(fNOut, "SNRResults.txt");
		
		//! *** Index of parameters to be chosen randomly
		int * iParamRand(NULL);
		int NParamRand(0);
		
		// ***** Read options *****
		for(int iarg=1; iarg<argc; iarg++){
			
			if((argc>1)&&(strcmp(argv[iarg],"-nf")==0)){
                SnAnalytic = false;
                
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
                SnAnalytic = false;
                
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
            
			if((argc>1)&&(strcmp(argv[iarg],"-o")==0)){
				// ***** Read base name of output file in option
				strcpy(fNOut, argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-t0")==0)){
				t0 = atof(argv[iarg+1]); 
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-T")==0)){
				tDur = atoi(argv[iarg+1]); 
				nOptions +=2;
			}
            if((argc>1)&&(strcmp(argv[iarg],"-dt")==0)){
				dtMes = atoi(argv[iarg+1]); 
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
            if((argc>1)&&(strcmp(argv[iarg],"-cu")==0)){
				RecCumul = true;
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
			NfNGWCfg++;
			fNGWCfg = (char**)MT.ReAllocMemory(fNGWCfg, (NfNGWCfg-1)*sizeof(char*), NfNGWCfg*sizeof(char*));
			fNGWCfg[NfNGWCfg-1] = (char*)MT.AllocMemory(2048*sizeof(char));
			strcpy(fNGWCfg[NfNGWCfg-1], argv[ia]);
		}
		
		
		
		// ***** Display for checking which files we are using
		Coutm << "Configuration files : " << Endl;
		for(int i=0; i<NfNGWCfg; i++){
			Coutm << "\t + " << fNGWCfg[i] ; MT.CheckFile("", fNGWCfg[i]);
		}
		
        
		//! *** Read the bloc name of TDI generators and obtain the individual names and the number of generators
		int ipTDINameAll(0);
		MT.wextractcoma(TDINameAll, ipTDINameAll, TDIName, NTDI, NCharTitCol);
		Coutm << "TDI generators : ";
		for(int i=0; i<NTDI; i++)
			Coutm << " " << TDIName[i];
		Coutm << Endl;
		
		//! *** Allocate memory related TDI size
		tSig = (double**) MT.AllocMemory(NTDI*sizeof(double*));
		fSig = (dcomplex**) MT.AllocMemory(NTDI*sizeof(dcomplex*));
		SNR = (double*) MT.AllocMemory(NTDI*sizeof(double));
		for(int i=0; i<NTDI; i++){
			tSig[i] = NULL;
			fSig[i] = NULL;
			SNR[i] = 0.;
		}
        
        S.sTDI(TDINameAll);
        S.AllocBaseMem();
		
        
        if(SnAnalytic){
            if(ModelConfNoise>=0)
                S.GalacticNoise = 1;
        }else{
            //! *** Read the PSD files
            fTDIPSD = (LCDataFileRead **) MT.AllocMemory(NTDI*sizeof(LCDataFileRead*));
            ifTDIPSD = (int *) MT.AllocMemory(NTDI*sizeof(int));
            for(int i=0; i<NTDI; i++){
                fTDIPSD[i] = NULL;
                ifTDIPSD[i] = 0;
            }
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
                    if(TDIName[iTDI]>=0){
                        Coutm << Endl << "\t - " << TDIName[iTDI] << " : signal " << ifTDIPSD[iTDI] << " of :";
                        fTDIPSD[iTDI]->ControlDisplay();
                    }else{
                        std::cerr << "ERROR :Index of column (" << ifTDIPSD[iTDI] << " ) for TDI " << TDIName << " have to be >= 0 ! Did you really include f as first column for example : example f,X,Am,Em ." << Endl;
                        throw std::invalid_argument("ERROR :Index of column have to be >= 0 ! !");
                    }
                }else{
                    std::cerr << "ERROR : No PSD of noises for TDI " << TDIName << " !" << Endl;
                    throw std::invalid_argument("ERROR : A PSD of noise file is missing for one TDI !");
                }
            }
            
            
            //! *** Load Confusion noise if needed
            if(ModelConfNoise>=0)
                LoadConfusionNoise(&MT, ConfNoise, ModelConfNoise, fTDIPSD[0]->getx0(), fTDIPSD[0]->getxend(), fTDIPSD[0]->getdx());
		}
		
		fOut.open(fNOut);
		
		//! ********* Compute the response to gravitational waves
		
		for(int iRun=0; iRun<NRun; iRun++){
			std::cerr << Endl << ">>>>>>> Run " << iRun << Endl;
			
			Sim  = new LISACode(&MT);
			Sim->setTimeInfo(dtMes, t0, tDur);
			
			
			//! *** Configure the simulator using the configuration files
			for(int i=0; i<NfNGWCfg; i++)
				Sim->config(fNGWCfg[i]);
			
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
			
			//! *** Choose randomly the parameters to randomize
			if(iRun!=0){
				for(int iP=0; iP<NParamRand; iP++){
					//! *** THIS CONDITION IS FOR VERIFICATION BINARIES (if inclination = 60 deg, then randomize it) SO IT SHOULD BE REMOVED FOR ANY OTHER SOURCE !!! 
					if((iParamRand[iP]!=3)||(((iParamRand[iP]==3))&&(MT.deq(Sim->GWgetParam(0,iParamRand[iP]),60.*M_PI/180.))))
						Sim->GWRandParam(iParamRand[iP]);
				}
			}
				
			//! *** Write header of output file
			if(iRun == 0){
				fOut << "#i";
				Sim->GWDispParamName(&fOut);
				for(int i=0; i<NTDI; i++)
					fOut << " SNR_" <<  TDIName[i];
				fOut << " SNRAll" << Endl;
			}
			
			//! *** Initialization of simulation
			Sim->init();
            //fLowMax=Sim->get... // WARNING: fLowMax fixed by hard
			
			//! *** Add output series with allocation of the series if needed
			for(int i=0; i<NTDI; i++){
				Sim->AddSerieOut(TDIName[i], 1, tmpNtData, (iRun == 0));
				if(NtData == -1)
					NtData = tmpNtData; 
				if(NtData != tmpNtData){
					Coutm << "ERROR : The numbers of data in TDI time vector are not the same : previously it was " << NtData << "but it's " << tmpNtData << "for the last one (" << TDIName[i] << ") !" << Endl;
					throw std::invalid_argument("ERROR : The numbers of data in TDI time vector are not the same !");
				}
			}
			//! *** Link output series
			for(int i=0; i<NTDI; i++)
				Sim->LinkSerieOut(i, tSig[i]);
			
			Sim->DispInfo("");
			
			
			//! *** Compute the time series : Run the simulation
			Sim->Run();
			
			/*
			std::ofstream fPtCheck("CheckSNR_tTDI.txt");
			for(int iT=0; iT<NtData; iT++){
				fPtCheck << iT*dtMes;
				for(int i=0; i<NTDI; i++)
					fPtCheck << " " << tSig[i][iT];
				fPtCheck << Endl;
			}
			*/
			
			//! ****** Compute the Fourrier transform
			//! *** Prepare computation
			if(iRun == 0){
				if(NRun<10)
					MT.unsetBestFTFwd();
				
				NfData = MT.getNfFTreal(NtData);
				df = 1./(dtMes*NtData);
                if(SnAnalytic){
                    iFmin = 2;
                    iFmax = NfData - 2;
                }else{
                    iFmin = MAX(MT.iceil(fTDIPSD[0]->getx0()/df), 0)+2;
                    iFmax = MIN(MT.iceil(fTDIPSD[0]->getxend()/df), NfData)-2;
                }
				for(int i=0; i<NTDI; i++)
					fSig[i] = (dcomplex*) MT.AllocMemory(NfData*sizeof(dcomplex));
				
				/*
				//! *** Check the PSD of noise and confusion noise
				std::ofstream fCheckSn("PSDCheck.txt");
				fCheckSn.precision(12);
				for(int iF=iFmin; iF<iFmax; iF++){
					fCheckSn << iF*df << " " << ConfNoise->gData(iF*df, LIN, 0);
					for(int iS=0; iS<NTDI; iS++)
						fCheckSn << " " << fTDIPSD[iS]->gData(ifTDIPSD[iS], iF*df, LIN, 0) << " " <<  (fTDIPSD[iS]->gData(ifTDIPSD[iS], iF*df, LIN, 0) + ConfNoise->gData(iF*df, LIN, 0));
					fCheckSn << Endl;
				}
				fCheckSn.close();
				 */
			}
			
			//! *** Compute Fourrier transform
			MT.FTMakeFwdMulti(tSig, fSig, NtData, NTDI);
			
			//! *** Normalize Fourrier transform
			for(int iF=0; iF<NfData; iF++)
				for(int iS=0; iS<NTDI; iS++)
					fSig[iS][iF] *= dtMes;
			
				
			
			Coutm << "Frequency frame : df = " << df << " Hz ,  Fmin = " << iFmin*df << " Hz ,  Fmax = " << iFmax*df << " Hz" << Endl; 
			
			/*
			std::ofstream fPfCheck("CheckSNR_fTDI.txt");
			for(int iF=0; iF<NfData; iF++){
				fPfCheck << iF*df;
				for(int i=0; i<NTDI; i++)
					fPfCheck << " " << norm(fSig[i][iF]);
				fPfCheck << Endl;
			}
			*/
			
			
			//! ***** Compute the SNR
			std::ofstream * fCheck(NULL);
            double Sn,SnMin(1.e30);
            if(RecCumul){
                fCheck = new std::ofstream("CumulSNR.txt");
                fCheck->precision(12);
            }
			SNRAll =0.;
			for(int iS=0; iS<NTDI; iS++){
				tmpRes = 0.;
				for(int iF=iFmin; iF<iFmax; iF++){
					//if(iS==0)
					//	Coutm << iF << " " << iF*df << " " << real(fSig[iS][iF]) << " " << imag(fSig[iS][iF]) << " " << real(fSig[iS][iF] * conj(fSig[iS][iF])) << " " << fTDIPSD[iS]->gData(ifTDIPSD[iS], iF*df, LIN, 0) << " " << real(tmpRes) << " " << imag(tmpRes) << Endl;
					if(SnAnalytic){
                        Sn = S.PSDNoise(iF*df,0);
                    }else{
                        if(ConfNoise != NULL){
                            Sn = fTDIPSD[iS]->gData(ifTDIPSD[iS], iF*df, LIN, 0) + ConfNoise->gData(iF*df, LIN, 0);
                        }else{
                            Sn = fTDIPSD[iS]->gData(ifTDIPSD[iS], iF*df, LIN, 0);
                        }
                        if (iF*df<fLowMax){
                            if(SnMin>Sn)
                                SnMin = Sn;
                        }else{
                            if(Sn<SnMin)
                                Sn = SnMin;
                        }
                    }
                    
                    tmpRes += fSig[iS][iF] * conj(fSig[iS][iF]) / Sn ;
                    if(RecCumul)
                        (*fCheck) << iF*df << " " <<  4.0 * df * tmpRes.real() << " " << Sn << " " << fSig[iS][iF].real() << " " << fSig[iS][iF].imag() << " " << (fSig[iS][iF] * conj(fSig[iS][iF])).real() << Endl;
                    
                }
                tmpSP = 4.0 * df * tmpRes.real();
                SNR[iS] = sqrt(tmpSP);
                SNRAll += SNR[iS]*SNR[iS];
            }
            if(RecCumul) fCheck->close();
            
			SNRAll = sqrt(SNRAll);
			
			//! **** Display results
			Coutm << Endl << "Run " << iRun << " : ";
			Sim->GWDispParam(MT.o);
			Coutm << Endl << "\t ==>";
			for(int i=0; i<NTDI; i++)
				Coutm << " SNR_" << TDIName[i] << " = " <<  SNR[i];
			Coutm << " ==> SNR_All = " << SNRAll << Endl;
			
			//! *** Record results
			fOut << iRun;
			Sim->GWDispParam(&fOut);
			for(int i=0; i<NTDI; i++)
				fOut << " " <<  SNR[i];
			fOut << " " << SNRAll << Endl;
			
			//! *** Delete the simulation
			delete Sim;
			Sim = NULL;
			
		}
		
		
		fOut.close();
		
		
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
			
		if(tSig != NULL){
			for(int i=0; i<NTDI; i++)
				if(tSig[i] != NULL)
					MT.Free(tSig[i], NtData*sizeof(double));
			MT.Free(tSig, NTDI*sizeof(double*));
		}
		
		if(fSig != NULL){
			for(int i=0; i<NTDI; i++)
				if(fSig[i] != NULL)
					MT.Free(fSig[i], NfData*sizeof(dcomplex));
			MT.Free(fSig, NTDI*sizeof(dcomplex*));
		}
		
		if(SNR != NULL)
			MT.Free(SNR, NTDI*sizeof(double));

		if(iParamRand != NULL)
			MT.Free(iParamRand,NParamRand*sizeof(int));
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "LC2SNR: error: " << e.what()<<Endl;
		std::cerr << "LC2SNR: abort!" << Endl;
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

// end of LISACODE-SNR.cpp
