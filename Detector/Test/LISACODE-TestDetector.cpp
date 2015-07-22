// $Id:  Exp $
/*
 *  LISACODE-TestDetector.cpp
 *
 *  Created by Antoine Petiteau on 25/02/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "LISACODE-Constants.h"
#include "LISACODE-Tools.h"
#include "LISACODE-GW.h"
#include "LISACODE-GWFile.h"
#include "LISACODE-GWStochastic.h"
#include "LISACODE-GWGalBin.h"
#include "LISACODE-DataFileWrite.h"
#include "LISACODE-OrbitsAnaLISA.h"
#include "LISACODE-ArmResp.h"

/**\ingroup OBPM
 * \{
 */


/** \brief Main of Code for testing measurement system : optical bench path + photodiode + phasemeter 
 * \author A. Petiteau
 * \version 2.0
 * \date 25/04/2011
 *
 *
 */
int main (int argc, char * const argv[])
{
	try {
		LCTools MT;
		int nOptions(0);
		
		throw std::invalid_argument("Sorry but this test module is not yet ready");
		
		//! ********* Declaration of varaiables 
		long SeedRand((long)time(NULL));
		char fOut[1024];
		char fInXML[1024];
		double t0(0.0), t;
		double dt(5.);
		int NDat(MT.ifloor(31557600./dt)); //! One year
		
		
		
		//! *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution : Compute gravitational wave signal for the 6 arms" << Endl;
			Coutm << "\t\t(./)LC2TestDetector [Options] " << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << "\t\t * %sfOut       : Output file name : t hBp1 hBc1 hBp2 hBc2 ... [default: TestDetector]. " << Endl;
			Coutm << "\t\t * %sfInXMLFile : XML input file [default: TestDetector.xml]. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -s %dseed : Seed for random gennerator. [default : current time]"  << Endl ;
			Coutm << "\t\t * -t %fdtmes : Time step. [default : 1 ]"  << Endl ;
			Coutm << "\t\t * -dn \t\t: No screen display. [default: false]"  << Endl ;
			Coutm << "\t\t * -dl %sfile \t: Write standard output in a file. [default: no file]"  << Endl ;
			Coutm << "\t\t * -v \t\t: Verobse : display full details. [default: false]"  << Endl ;
			Coutm << Endl << " ----------------" << Endl;
			return 0;
			
		}
		
		//! *********** Version *************
		if(((argc>1)&&(strcmp(argv[1],"--version")==0))&&((argc>1)&&(strcmp(argv[1],"-v")==0))){
			Coutm << " ----- VERSION -----" << Endl;
			Coutm << " LC2TestDetector : executable testing measurement system of LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		for(int iarg=1; iarg<argc; iarg++){
			if((argc>1)&&(strcmp(argv[iarg],"-s")==0)){
				SeedRand = atoi(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-t")==0)){
				dt = atof(argv[iarg+1]);
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
		
		
		
		//! ***** Initialization of variable
		strcpy(fOut,"TestDetector.txt");
		strcpy(fInXML,"/Users/petiteau/Applications/src/LISACode/LISACode_2_0/LISACode/Main/Test/Config-.xml");
		
		
		if(argc-nOptions>1){
			strcpy(fOut, argv[1+nOptions]);
		}
		
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "TestDetector: error: " << e.what()<<Endl;
		std::cerr << "TestDetector: abort!" << Endl;
		exit(1);
	}
	return(0);
};


/** \}*/

// end of LISACODE-TestDetector.cpp