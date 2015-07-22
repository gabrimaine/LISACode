// $Id:  $
/*
 *  LISACODE-Run.cpp
 *  V 2.0
 *
 *  Created on 14/04/2005 by  Antoine Petiteau (APC)
 *  Last modification on 14/02/2011 by Antoine Petiteau (AEI)
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

/**\mainpage Developers' Manual of LISACode 
 This document provides a description of the LISACode software for future developpers. Some information about the code execution (\ref use) are also provided to users. 
 This manual is divided in three sections:
 - \ref intro 
 - \ref codeDesc 
 - \ref use
 
 */

/**\page intro Introduction
 
 LISACode[<a href="Bibliography.html#Petiteau-2008">Petiteau &nbsp; 2008</a>] is a LISA mission simulator.
 More details on LISACode and news concerning this software are given on LISA website: http://www.aei.mpg.de/~petiteau/LISACode/Home.html
 It is highly structured and programmed in C++. 
 The simulator has the purpose to bridge the gap between the basic principles of LISA and a sophisticated end-to-end engineering level simulator.  
 This software package, which runs on most computer platforms, can be downloaded from the download page on LISACode web site ( http://www.aei.mpg.de/~petiteau/LISACode/Download.html ).
 
 \section techDesc LISACode technical description
 
 \todo Rewrite LISACode technical description 
 
 LISACode simulates the LISA gravitational wave (GW) detector (see http://www.esa.int/esaSC/SEMEJRR1VED_index_0.html and [<a href="Bibliography.html#LISAPPA">Bender &nbsp; 1998</a>]) and modelizes the waveform of GW emitted from large kind of sources.
 
 It does not aim at simulating the LISA detector in detail but rather it uses the response function of its main components, particularly because they will affect the noise level of the detector response. 
 It also includes an implementation of the TDI (Time Delay Interferometry, [<a href="Bibliography.html#TDIRevue">Tinto &nbsp; 2004</a>] [<a href="Bibliography.html#TDIVinet">Dhurandhar &nbsp; 2002</a>]) 
 technique which allows to suppress the noise introduced by lasers frequency instability. 
 
 The main inputs and outputs of LISACode are time-dependent sequences. 
 Input sequences describe the GW strain and output sequences describe the phasemeters response or their treatment via various TDI combination.
 
 A number of GW signals can be defined, but the main aim of the code is to be used in conjunction with more sophisticated GW simulators via intermediate data files.
 
 */


/**\page codeDesc A description of the Code
 \todo Rewrite description of the code
 
 \section org Code organisation 
 LISACode is written in C++ and has a very modular structure. 
 The main structure of LISACode is shown in the figure below. 
 This structure maps the main components of the LISA detector as well as its physical inputs (see details in [<a href="Bibliography.html#LISACode">Petiteau &nbsp; 2008</a>][<a href="Bibliography.html#LISACodeLISASymp6">Petiteau &nbsp; 2006</a>]).
 Its main components are: 
 \arg a variety of GW inputs, 
 \arg a detailed description of the orbits [<a href="Bibliography.html#OrbitLISA">Dhurandhar &nbsp; 2004</a>] of the three satellites (including the breathing and rotation
 modes of the LISA triangle), 
 \arg the different noise sources (lasers, DFS and the
 interferometric measurements),
 \arg the phasemetre measurements and 
 \arg the Ultra Stable Oscillator (USO) clock performances.
 
 \image html "Structure_EN.pdf"
 \image latex "Structure_EN.eps" "Structure of LISACode" width=12cm
 
 
 \todo ... continue ...
 
 */

/**\page use Use of LISACode and parameters
 
 \todo Write use of LISACode and parameters
 
 */

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "LISACODE-Constants.h"
#include "LISACODE-Tools.h"
#include "LISACODE-LISACode.h"


/** \defgroup main Main 
 * This group manage the main elements
 * \{
 */

/** \brief Main of LISACode.
 * \author A. Petiteau
 * \version 2.0
 * \date 14/02/2011
 * 
 * Running :
 * \arg Help.
 * \arg Version.
 * \arg Read options.
 *
 *
 */
int main (int argc, char * const argv[])
{
	try {
		LCTools MT;
		int nOptions(0);
		
		long SeedRand((long)time(NULL));
		double t0(0.);
		double tDur(1048576.0);
		int it0(0);
		bool Deft0(false), DeftDur(false);
		bool DefOut(false);
		char OutKey[100];
		char OutBName[10000];
        char XMLHeadName[1024];
		
        strcpy(XMLHeadName, "None");
        
		// *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution :" << Endl;
			Coutm << "\t\t(./)LISACode2 [Options] ConfigFile1.xml ConfigFile2.xml ..." << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << Endl << "\t\t * ConfigFileI.xml (required) : Configuration files. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -s %dseed : Seed for random gennerator. [default : current time]"  << Endl ;
			Coutm << "\t\t\t  NB: If not specified the random seed is the time machine." << Endl;
			Coutm << "\t\t * -dn \t\t: No screen display. [default: false]"  << Endl ;
			Coutm << "\t\t * -t0 %ft0 : Initial time [default : use the value read in config]"  << Endl ;
			Coutm << "\t\t * -it0 %dit0 : Initial time = it0 * duration [default : use the value read in config Sim or " << t0 << " (if no config Sim)]"  << Endl ;
			Coutm << "\t\t * -T %fT : Duration [default : use the value read in config Sim or " << tDur << " (if no config Sim)]"  << Endl ;
            Coutm << "\t\t * -x %s  : Name of global xml header [ default : no xml header]"  << Endl ;
			Coutm << "\t\t * -o %sKey %sBaseName : Shrotcut for defining common base name for output file and the type of output from"  << Endl ;
			Coutm << "\t\t\t key coding : P->Phasemeter(all), T->TDI(X,Y,Z,X2,Y2,Z2), D->Delay(D1,D2,D3,D1p,D2p,D3p) " << Endl;
			Coutm << "\t\t\t For example : '-o PTD Truc' [default : use the value read in config Sim or nothing (if no config Sim)]"  << Endl ;
			Coutm << "\t\t * -dl %sfile \t: Write standard output in a file. [default: no file]"  << Endl ;
			Coutm << "\t\t * -v \t\t: Verbose mode : display full details. [default: false]"  << Endl ;
			Coutm << "\t\t * -h \t\t: This help."  << Endl ;
			Coutm << "\t\t * -V \t\t: Version."  << Endl ;
			Coutm << Endl << "\tNOTE :" << Endl;
			Coutm <<  "\t\t There are more details in the user's manual (which is not yet ready, sorry !)." << Endl;
			Coutm <<  "\t\t Technical details in Documentation/latex/refman.pdf." << Endl;
			Coutm << Endl << " ----------------" << Endl;
			return 0;
			
		}
		
		// *********** Version *************
		if(((argc>1)&&(strcmp(argv[1],"--version")==0))||((argc>1)&&(strcmp(argv[1],"-V")==0))){
			Coutm << " ----- VERSION -----" << Endl;
			Coutm << " LISACode2 : main executable of LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		
		// ***** Main declaration *****
		LISACode Sim(&MT);
		
		
		// ***** Read options *****
		for(int iarg=1; iarg<argc; iarg++){
			if((argc>1)&&(strcmp(argv[iarg],"-s")==0)){
				if(iarg+1>=argc)
					throw std::invalid_argument("ERROR : We need a value after option -s !");
				SeedRand = atoi(argv[iarg+1]); 
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-t0")==0)){
				if(iarg+1>=argc)
					throw std::invalid_argument("ERROR : We need a value after option -t0 !");
				t0 = atof(argv[iarg+1]); 
				Deft0 = true;
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-it0")==0)){
				if(iarg+1>=argc)
					throw std::invalid_argument("ERROR : We need a value after option -it0 !");
				it0 = atof(argv[iarg+1]); 
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-T")==0)){
				if(iarg+1>=argc)
					throw std::invalid_argument("ERROR : We need a value after option -T !");
				tDur = atof(argv[iarg+1]); 
				DeftDur = true;
				nOptions +=2;
			}
            if((argc>1)&&(strcmp(argv[iarg],"-x")==0)){
				if(iarg+1>=argc)
					throw std::invalid_argument("ERROR : We need a value after option -T !");
				strcpy(XMLHeadName,argv[iarg+1]); 
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-o")==0)){
				if(iarg+2>=argc)
					throw std::invalid_argument("ERROR : We need 2 words after option -o !");
				strcpy(OutKey, argv[iarg+1]);
				strcpy(OutBName, argv[iarg+2]);
				DefOut = true;
				nOptions +=3;
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
		
		// ***** Configuration using xml file in arguments *****
		if (!MT.wcmp(XMLHeadName, "None")) {
            Coutm << "Global XML output file : " << XMLHeadName << Endl;
            Sim.setGlobalXMLHeader(XMLHeadName);
        }
        
		for(int ia=1+nOptions; ia<argc; ia++){
			Sim.config(argv[ia]);
		}
		
		//! *** If initial time and/or duration are define in options replace the current value by the one read in options 
		if(Deft0)
			Sim.sett0(t0);
		if(DeftDur)
			Sim.settDur(tDur);
		
		if(DefOut){
			//! *** Configure output file from shortcut
			for(int iC=0; iC<strlen(OutKey); iC++){
				if(OutKey[iC]=='P'){
					char * TmpN;
					char ** TmpO;
					TmpN = (char*) MT.AllocMemory(10000*sizeof(char));
					TmpO = (char**) MT.AllocMemory(4*sizeof(char*));
					for(int iO=0; iO<4; iO++)
						TmpO[iO] = (char*) MT.AllocMemory(16*sizeof(char));
					sprintf(TmpN,"%s-SC1.txt",OutBName);
					strcpy(TmpO[0],"sci1");
					strcpy(TmpO[1],"sci1s");
					strcpy(TmpO[2],"tau1");
					strcpy(TmpO[3],"tau1s");
					Sim.AddOutputFile(0, TmpN, ASCII, TmpO, 4);
					
					sprintf(TmpN,"%s-SC2.txt",OutBName);
					strcpy(TmpO[0],"sci2");
					strcpy(TmpO[1],"sci2s");
					strcpy(TmpO[2],"tau2");
					strcpy(TmpO[3],"tau2s");
					Sim.AddOutputFile(0, TmpN, ASCII, TmpO, 4);
					
					sprintf(TmpN,"%s-SC3.txt",OutBName);
					strcpy(TmpO[0],"sci3");
					strcpy(TmpO[1],"sci3s");
					strcpy(TmpO[2],"tau3");
					strcpy(TmpO[3],"tau3s");
					Sim.AddOutputFile(0, TmpN, ASCII, TmpO, 4);
				}
				if(OutKey[iC]=='T'){
					char * TmpN;
					char ** TmpO;
					TmpN = (char*) MT.AllocMemory(10000*sizeof(char));
					TmpO = (char**) MT.AllocMemory(6*sizeof(char*));
					for(int iO=0; iO<6; iO++)
						TmpO[iO] = (char*) MT.AllocMemory(16*sizeof(char));
					sprintf(TmpN,"%s-TDI.txt",OutBName);
					strcpy(TmpO[0],"X");
					strcpy(TmpO[1],"Y");
					strcpy(TmpO[2],"Z");
					strcpy(TmpO[3],"X2");
					strcpy(TmpO[4],"Y2");
					strcpy(TmpO[5],"Z2");
					Sim.AddOutputFile(1, TmpN, ASCII, TmpO, 6);
				}
				if(OutKey[iC]=='D'){
					char * TmpN;
					char ** TmpO;
					TmpN = (char*) MT.AllocMemory(10000*sizeof(char));
					TmpO = (char**) MT.AllocMemory(6*sizeof(char*));
					for(int iO=0; iO<6; iO++)
						TmpO[iO] = (char*) MT.AllocMemory(16*sizeof(char));
					sprintf(TmpN,"%s-Delay.txt",OutBName);
					strcpy(TmpO[0],"D1");
					strcpy(TmpO[1],"D2");
					strcpy(TmpO[2],"D3");
					strcpy(TmpO[3],"D1p");
					strcpy(TmpO[4],"D2p");
					strcpy(TmpO[5],"D3p");
					Sim.AddOutputFile(2, TmpN, ASCII, TmpO, 6);
				}
			}
			
		}
		
		//Sim.DispInfo("");
		
		
		
		
		// ***** Initialization *****
		
		Sim.init();
		Sim.DispInfo("");
		
		
		// ***** Run *****
		
		Sim.Run();
		
		
		
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "LISACode2:Run: error: " << e.what()<<Endl;
		std::cerr << "LISACode2:Run: abort!" << Endl;
		exit(1);
	}
	return(0);
};


/** \}*/

// end of LISACODE-Run.cpp
