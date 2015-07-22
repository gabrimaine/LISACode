// $Id:  Exp $
/*
 *  LISACODE-zDL.cpp
 *
 *  Created by Antoine Petiteau on 14/02/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "LISACODE-Constants.h"
#include "LISACODE-Tools.h"

/**\ingroup ToolBox 
 * \{
 */


/** \brief Main of Code for testing toolbox.
 * \author A. Petiteau
 * \version 2.0
 * \date 07/02/2012
 *
 *
 */
int main (int argc, char * const argv[])
{
	try {
		LCTools MT;
		int nOptions(0);
		
		long SeedRand((long)time(NULL));
		
		//! *** Cosmological parameters
		double H0(70.);
		double Omm(0.3);
		double Oml(0.7);
		double w(0.);
		
		
		// *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution :" << Endl;
			Coutm << "\t\t(./)LC2zDL [Options] keyW value  " << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << "\t\t * keyW describe type of computation :" << Endl;
			Coutm << "\t\t 'z' : redshift z as input --> compute luminosity distance DL " << Endl;
			Coutm << "\t\t 'D' : luminosity distance DL as input --> compute redhsift z " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -H0  %fValue : Current value of Hubble expansion parameter. [default : " << H0 << " ]"  << Endl ;
			Coutm << "\t\t * -Omm %fValue : Omega matter. [default : " << Omm << " ]"  << Endl ;
			Coutm << "\t\t * -Oml %fValue : Omega lambda. [default : " << Oml << " ]"  << Endl ;
			Coutm << "\t\t * -w   %fValue : (- 1 - omDE) : parameter of dark energy equation of state . [default : " << w << " ]"  << Endl ;
			Coutm << "\t\t * -v \t\t: Verobse : display full details. [default: false]"  << Endl ;
			Coutm << Endl << " ----------------" << Endl;
			return 0;
			
		}
		
		// *********** Version *************
		if(((argc>1)&&(strcmp(argv[1],"--version")==0))&&((argc>1)&&(strcmp(argv[1],"-v")==0))){
			Coutm << " ----- VERSION -----" << Endl;
			Coutm << " LC2zDL : Compute relation between redshift z and luminosity distance DL - LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		
		
		
		for(int iarg=1; iarg<argc; iarg++){
			if((argc>1)&&(strcmp(argv[iarg],"-H0")==0)){
				H0 = atof(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-Omm")==0)){
				Omm = atof(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-Oml")==0)){
				Oml = atof(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-w")==0)){
				w = atof(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-v")==0)){
				MT.setDispDetails();
				nOptions++;
			}
		} 
		MT.setRandSeed(SeedRand);
		
		
		//! ***** Some declarations
		double z(1.), DL(1.);
		char keyW;
		double Val(0.);

		
		
		if(argc-nOptions<3){
			throw std::invalid_argument("ERROR : we need two arguments : keyWord (z or D) and value");
		}
		
		keyW = argv[1+nOptions][0];
		Val  = atof(argv[2+nOptions]);
		

		Coutm << "Cosmological parameters :" << Endl;
		Coutm << "\t - Current value of Hubble expansion parameter : H0 = " << H0 << Endl;
		Coutm << "\t - Omega matter : Omm = " << Omm << Endl;
		Coutm << "\t - Omega lambda : Oml = " << Oml << Endl;
		Coutm << "\t - Parameter of dark energy equation of state : w = " << w << Endl;
		
		
		if(keyW == 'D'){
			DL = Val;
			z = MT.redshift(DL, H0, Omm, Oml, w);
			Coutm << "For the luminosity distance DL = " << DL << " kpc , the redshift is z = " << z << " ." << Endl;
		}
		
		
		if(keyW == 'z'){
			z = Val;
			DL = MT.DL(z, H0, Omm, Oml, w);
			Coutm << "For the redshift is z = " << z << " , the luminosity distance DL = " << DL << " kpc." << Endl;
		}
		
	}
	
	
	
	
	
	catch(std::exception & e ) {
		std::cerr << "zDL: error: " << e.what()<<Endl;
		std::cerr << "zDL: abort!" << Endl;
		exit(1);
	}
	return(0);
};


/** \}*/

// end of