// $Id:  $
/*
 *  LISACODE-MeanFiles.h
 *
 *  Created on 03/05/2005 by  Antoine Petiteau (AEI)
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
#include "LISACODE-DataFileRead.h"

/** \ingroup main Main 
 * \{
 */

/** \brief The executable compute mean between several input files with exactly the same size
 * \author A. Petiteau
 * \version 2.0
 * \date 03/05/2011
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
			Coutm << "\t\t(./)LC2MeanFiles [Options] File1 File2 ..." << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << Endl << "\t\t * FileI (at least one required) : input files. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -o %sOutputFile  : output file [default: MeanRes.txt]"  << Endl ;
			Coutm << "\t\t * -v \t\t: Verbose mode : display full details. [default: false]"  << Endl ;
			Coutm << "\t\t * -h \t\t: This help."  << Endl ;
			Coutm << "\t\t * -V \t\t: Version."  << Endl ;
			Coutm << Endl << " ----------------" << Endl;
			return 0;
			
		}
		
		// *********** Version *************
		if(((argc>1)&&(strcmp(argv[1],"--version")==0))||((argc>1)&&(strcmp(argv[1],"-V")==0))){
			Coutm << " ----- VERSION -----" << Endl;
			Coutm << " LC2MeanFiles : Just mean several files - LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		
		// ***** Main declaration *****
		char fNOut[16384];
		strcpy(fNOut,"MeanRes.txt");
		char ** fNIn(NULL);
		int NfNIn(0);
		
		double ** Dat(NULL);
		int NDat, NRec;
		double x0, dx;
		
		
		// ***** Read options *****
		for(int iarg=1; iarg<argc; iarg++){
			
			if((argc>1)&&(strcmp(argv[iarg],"-o")==0)){
				// ***** Read base name of output file in option
				strcpy(fNOut, argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-v")==0)){
				MT.setDispDetails();
				nOptions++;
			}
		}
		
		
		// ***** Read configuration files in arguments
		if (argc<=1+nOptions)
			throw std::invalid_argument("ERROR : At least one input file is needed");
		
		for(int ia=1+nOptions; ia<argc; ia++){
			NfNIn++;
			fNIn = (char**)MT.ReAllocMemory(fNIn, (NfNIn-1)*sizeof(char*), NfNIn*sizeof(char*));
			fNIn[NfNIn-1] = (char*)MT.AllocMemory(2048*sizeof(char));
			strcpy(fNIn[NfNIn-1], argv[ia]);
		}
		
		
		//! **** Loop into input files
		
		LCDataFileRead * fIn;
		for(int iF=0; iF<NfNIn; iF++){
			fIn = new LCDataFileRead(&MT, fNIn[iF]);
			fIn->init();
			
			if(iF==0){
				//! *** In the first input file
				//! ** Read size of data and time info 
				x0 = fIn->getx0();
				dx = fIn->getdx();
				NDat = fIn->getNDat();
				NRec = fIn->getNRec();
				
				//! ** Allocate memory for output
				Dat = (double**) MT.AllocMemory(NRec*sizeof(double*));
				for(int iR=0; iR<NRec; iR++){
					Dat[iR] = (double*) MT.AllocMemory(NDat*sizeof(double));
					for(int i=0; i<NDat; i++)
						Dat[iR][i] = 0.;
				}
			}
			
			if(!MT.deq(x0,fIn->getx0()))
				throw std::invalid_argument("ERROR : The last file don't have the same initial time as the first one !");
			
			if(!MT.deq(dx,fIn->getdx()))
				throw std::invalid_argument("ERROR : The last file don't have the same time step as the first one !");
			
			if(NDat != fIn->getNDat())
				throw std::invalid_argument("ERROR : The last file don't have the same number of data as the first one !");
			
			if(NRec != fIn->getNRec())
				throw std::invalid_argument("ERROR : The last file don't have the same number of records as the first one !");
			
			//! *** Read data
			for(int iR=0; iR<NRec; iR++){
				for(int i=0; i<NDat; i++){
					Dat[iR][i] += fIn->gDataBin(iR, i);
				}
			}
			
			delete fIn;
		}
		
		
		//! **** Compute the mean
		for(int iR=0; iR<NRec; iR++)
			for(int i=0; i<NDat; i++)
				Dat[iR][i] /= (double)NfNIn;
			
		
		//! **** Record in output
		std::ofstream fOut(fNOut);
		fOut.precision(12);
		Coutm << "Record the results in " << fNOut << " ... " << Endl;
		for(int i=0; i<NDat; i++){
			fOut << x0 + i*dx;
			for(int iR=0; iR<NRec; iR++){
				fOut << " " << Dat[iR][i];
			}
			fOut << Endl;
		}
		
		fOut.close();
		
		
		
		
		//! ** Free memory
		if(Dat!=NULL){
			for(int iR=0; iR<NRec; iR++)
				MT.Free(Dat[iR],NDat*sizeof(double));
			MT.Free(Dat,NRec*sizeof(double*));
		}
			
		if(fNIn != NULL){
			for(int iF=0; iF<NfNIn; iF++)
				if(fNIn[iF] != NULL)
					MT.Free(fNIn[iF], 2048*sizeof(char));
			MT.Free(fNIn, NfNIn*sizeof(char*));
		}
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "LC2MeanFiles: error: " << e.what()<<Endl;
		std::cerr << "LC2MeanFiles: abort!" << Endl;
		exit(1);
	}
	return(0);
};


// end of LISACODE-MeanFiles.cpp