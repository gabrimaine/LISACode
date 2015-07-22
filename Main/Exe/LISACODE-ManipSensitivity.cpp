// $Id:  $
/*
 *  LISACODE-ManipSensitivity.h
 *
 *  Created on 10/09/2011 by  Antoine Petiteau (AEI)
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

/** \brief The executable manipulate sensitivity data
 * \author A. Petiteau
 * \version 2.0
 * \date 10/09/2011
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
		
		int Nsm(1);
		int Nds(1);
		double fStartSmooth(2.e-4);
		double fMin(1.e-6), fMax(1.5);
		int SensType(0); //! Type of sensitivity (see option -O in help for details) : 0->S5, 1->S1, 2->Oh, 3->Nf, 4->Np
		char fNGalBin[16384];
		double L0m(1.e9);
		double Fact2SensStd((sqrt(LC::Yr_SI)/5.));
		
		strcpy(fNGalBin,"None");
		
		// *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution :" << Endl;
			Coutm << "\t\t(./)LC2ManipSensitivity [Options] InFile iCol OutFile" << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << Endl << "\t\t * InFile  : Input file containing sensitivity data . ";
			Coutm << Endl << "\t\t * iCol    : Index of column to read :";
			Coutm << "\t\t\t 2,...,NTDI+1 for noise PSD -> options '-O Nf' and '-O Np')" << Endl;
			Coutm << "\t\t\t NTDI+2,...,2*NTDI+1 for sensitivity -> options '-O S5', '-O S1' and '-O Oh')" << Endl;
			Coutm << Endl << "\t\t * OutFile : Output file. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -S %dNsm   : Number of data used to smooth (typical 11). [ default : " << Nsm << "]"  << Endl ;
			Coutm << "\t\t * -D %dNds   : Down-sampling : keep one data each Nds (typical 5). [ default : " << Nds << "]"  << Endl ;
			Coutm << "\t\t * -i %fx     : Value of reference after which one we start smoothing. [ default : " << fStartSmooth << " Hz ]"  << Endl ;
			Coutm << "\t\t * -f %ffMin  : Minimal frequency. [ default : " << fMin << " Hz ]"  << Endl ;
			Coutm << "\t\t * -F %ffMax  : Maximal frequency. [ default : " << fMax << " Hz ]"  << Endl ;
			Coutm << "\t\t * -O Key : Type of output. [ default : " << fMax << " Hz ]"  << Endl ;
			Coutm << "\t\t\t S5 : Sensitivity for SNR=5 and one year integration. [ default ]"  << Endl ;
			Coutm << "\t\t\t S1 : Standard sensitivity, i.e. corresponding to sci. req. doc. (SNR=1 and one sec integration.)."  << Endl ;
			Coutm << "\t\t\t Oh : Omega h2 : energy density."  << Endl ;
			Coutm << "\t\t\t Nf : Noise power spectral density in relative frequency unit. "  << Endl ;
			Coutm << "\t\t\t Np : Noise power spectral density in phase unit. "  << Endl ;
			Coutm << "\t\t * -C %fL %sFile : Armlength and file containing the confusion noise from Galactic binaries (format Tyson). [ default : L = " << L0m << " m , file = " << fNGalBin << " ]"  << Endl ;
			Coutm << "\t\t * -v \t\t: Verbose mode : display full details. [default: false]"  << Endl ;
			Coutm << "\t\t * -h \t\t: This help."  << Endl ;
			Coutm << "\t\t * -V \t\t: Version."  << Endl ;
			Coutm << Endl << " ----------------" << Endl;
			return 0;
			
		}
		
		// *********** Version *************
		if(((argc>1)&&(strcmp(argv[1],"--version")==0))||((argc>1)&&(strcmp(argv[1],"-V")==0))){
			Coutm << " ----- VERSION -----" << Endl;
			Coutm << " LC2ManipSensitivity : Manipulate sensitivity data - LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		
		// ***** Main declaration *****
		char fNIn[16384];
		char fNOut[16384];
		
		int iCol;
		
		
		
		// ***** Read options *****
		for(int iarg=1; iarg<argc; iarg++){
			
			if((argc>1)&&(strcmp(argv[iarg],"-S")==0)){
				Nsm = atoi(argv[iarg+1]);
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-i")==0)){
				fStartSmooth = atof(argv[iarg+1]);
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-D")==0)){
				Nds = atof(argv[iarg+1]);
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-f")==0)){
				fMin = atof(argv[iarg+1]);
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-F")==0)){
				fMax = atof(argv[iarg+1]);
				nOptions +=2;
			}
			
			if((argc>1)&&(strcmp(argv[iarg],"-O")==0)){
				bool NotFound(true);
				if( (MT.wcmp(argv[iarg+1],"S5")) ){
					SensType = 0;
					NotFound = false;
				}
				if( (MT.wcmp(argv[iarg+1],"S1")) ){
					SensType = 1;
					NotFound = false;
				}
				if( (MT.wcmp(argv[iarg+1],"Oh")) ){
					SensType = 2;
					NotFound = false;
				}
				if( (MT.wcmp(argv[iarg+1],"Nf")) ){
					SensType = 3;
					NotFound = false;
				}
				if( (MT.wcmp(argv[iarg+1],"Np")) ){
					SensType = 4;
					NotFound = false;
				}
				if(NotFound){
					Coutm << "ERROR : The sensitivity type " << argv[iarg+1] << " is unknown  (only S5, S1 or Oh) !" << Endl;
					throw std::invalid_argument("ERROR : The source type is unknown !");
				}
				nOptions +=2;
			}
			
			if((argc>2)&&(strcmp(argv[iarg],"-C")==0)){
				L0m = atof(argv[iarg+1]);
				strcpy(fNGalBin, argv[iarg+2]);
				nOptions +=3;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-v")==0)){
				MT.setDispDetails();
				nOptions++;
			}
		}
		
		
		//! ***** Read configuration files in arguments 
		if(argc > nOptions+3){
			strcpy(fNIn, argv[nOptions+1]);
			iCol = atoi(argv[nOptions+2]);
			strcpy(fNOut, argv[nOptions+3]);
		}else{
			throw std::invalid_argument("ERROR : We need 3 arguments : input file name, column index and output file name.");
		}
		
		
		//! ******* Load data
		Coutm << ">>> Load data ..." << Endl;
		LCDataFileRead fIn(&MT, fNIn);
		fIn.init();
		
		
		
		//! ******* Smooth data
		Coutm << ">>> Smooth data over " << Nsm << " points ..." << Endl;
		double * DatSm;
		int NDatSm(fIn.getNDat()-1);
		int NSmo2(Nsm/2), iSmin, iSmax, DiS;
		DatSm = (double*) MT.AllocMemory(NDatSm*sizeof(double));
		
		for(int iF=1; iF<fIn.getNDat(); iF++){
			iSmin = MAX(iF-NSmo2,0);
			iSmax = MIN(iF+NSmo2+1,fIn.getNDat());
			DiS = MAX(iSmax - iSmin, 1);
			if(fIn.gRefBin(iF)<fStartSmooth){
				DatSm[iF-1] = fIn.gDataBin(iCol-2, iF);
			}else{
				DatSm[iF-1] = 0.;
				for(int iS=iSmin; iS<iSmax; iS++){
					DatSm[iF-1] += fIn.gDataBin(iCol-2, iS);
					//Coutm << fIn.gDataBin(iCol-1, iF) << " "; 
				}
				DatSm[iF-1] /= ((double)(DiS));
				//if(iF > fIn.getNDat()-3)
				//	Coutm << fIn.gRefBin(iF) << " = " << iF << " => " << fIn.gDataBin(iCol-2, iF) << " [ " << iSmin << " : " << iSmax << " ] = " << DiS << " " << DatSm[iF] << Endl;
			}
		}
		

		
		//! ******* Down-sampling data
		Coutm << ">>> Down-sampling data at 1/" << Nds << " ..." << Endl;
		double * DatDs;
		double * RefDs;
		int NDatDs(MT.ifloor(NDatSm/Nds)), iDDs(0), iDownS(Nds);
		DatDs = (double*) MT.AllocMemory(NDatDs*sizeof(double));
		RefDs = (double*) MT.AllocMemory(NDatDs*sizeof(double));
		for(int iF=0; iF<NDatDs; iF++){
			RefDs[iF] = 0.;
			DatDs[iF] = 0.;
		}
		for(int iF=0; iF<NDatSm; iF++){
			if(iDownS==Nds){
				RefDs[iDDs] = fIn.gRefBin(iF+1);
				DatDs[iDDs] = DatSm[iF];
				iDDs++;
				iDownS = 0;
			}
			iDownS++;
		}
		//Coutm << iDDs << " / " << NDatDs << Endl;
		
		
		
		//! ******* Add confusion noise
		if(!MT.wcmp(fNGalBin, "None")){
			
			Coutm << ">>> Add confusion noise from Galactic binaries read in " << fNGalBin << " (using armlength = "  << L0m << " m)  ..." << Endl; 
			double L0s(L0m/LC::c_SI), sphiL;
			LCDataFileRead fGalBin(&MT, fNGalBin);
			fGalBin.init();
			if((SensType==0)||(SensType==1)||(SensType==2)){
				for(int iF=0; iF<NDatDs; iF++){
					if( (RefDs[iF]>fGalBin.getx0()) && (RefDs[iF]<fGalBin.getxend()) ){
						sphiL = sin( 2.0*M_PI*RefDs[iF]*L0s);
						//DatDs[iF] = sqrt( DatDs[iF]*DatDs[iF] + fGalBin.gData(3, RefDs[iF], LIN, 0)  * (5./3.) / pow( sin(2.*M_PI*RefDs[iF]*L0m/LC::c_SI), 2.) / (Fact2SensStd*Fact2SensStd) ) ;
						DatDs[iF] = sqrt( DatDs[iF]*DatDs[iF] + fGalBin.gData(3, RefDs[iF], LIN, 0)  * (5./3.) / ( sphiL*sphiL) / (Fact2SensStd*Fact2SensStd) ) ;
						//Coutm << RefDs[iF] << " " << DatDs[iF] << " " << fGalBin.gData(3, RefDs[iF], LIN, 0) << " " << sqrt(fGalBin.gData(3, RefDs[iF], LIN, 0) * (5./3.) / pow( sin(2.*M_PI*RefDs[iF]*L0/LC::c_SI), 2.)) << " " << sqrt(fGalBin.gData(3, RefDs[iF], LIN, 0) * (5./3.) / pow( sin(2.*M_PI*RefDs[iF]*L0/LC::c_SI), 2.)) / Fact2SensStd << Endl;
					}
				}
			}
			
			if((SensType==3)||(SensType==4)){
				for(int iF=0; iF<NDatDs; iF++){
					if( (RefDs[iF]>fGalBin.getx0()) && (RefDs[iF]<fGalBin.getxend()) ){
						//Coutm << RefDs[iF] << " " << DatDs[iF] << " " << 4.0 * L0s*L0s * (2.0*M_PI*RefDs[iF])*(2.0*M_PI*RefDs[iF]) * fGalBin.gData(3, RefDs[iF], LIN, 0) << Endl;
						DatDs[iF] =  DatDs[iF] + 4.0 * L0s*L0s * (2.0*M_PI*RefDs[iF])*(2.0*M_PI*RefDs[iF]) * fGalBin.gData(3, RefDs[iF], LIN, 0) ;
					}
				}
			}
		}

		
		//! ******* Conversion
		if(SensType==1){
			Coutm << ">>> Convertion of sensitvity to standard sensitivity ..." << Endl;
			for(int iF=0; iF<NDatDs; iF++)
				DatDs[iF] *= Fact2SensStd;
		}
		
		if(SensType==2){
			Coutm << ">>> Convertion of sensitvity to Omega h2 ..." << Endl;
			double H0std (100.);
			double H0Hz (H0std / LC::kpc_m);
			double FactConv( Fact2SensStd*Fact2SensStd* 4. * M_PI*M_PI / (3. * H0Hz*H0Hz ) );
			
			//! Conversion : \f[ S_{\Omega h^2} = C S_{std}^2 * f^3 \f] 
			//! with \f[ C = { 4 \pi^2 \over 3 H_{0,Hz}^2  \f]
			
			for(int iF=0; iF<NDatDs; iF++)
				DatDs[iF] = FactConv * DatDs[iF]*DatDs[iF] * RefDs[iF]*RefDs[iF]*RefDs[iF]  ;
		}
		
		if(SensType==4){
			Coutm << ">>> Convertion of noise PSD in relative frequency unit to phase unit ..." << Endl;
			for(int iF=0; iF<NDatDs; iF++)
				DatDs[iF] *= (LC::c_SI*LC::c_SI) / (4.*M_PI*M_PI * RefDs[iF]*RefDs[iF]) ;
		}
		
		
		//! ******* Record results
		Coutm << ">>> Record results ..." << Endl;
		std::ofstream fOut(fNOut);
		fOut.precision(12);
		for(int iF=0; iF<NDatDs; iF++)
			if((RefDs[iF]>fMin)&&(RefDs[iF]<fMax))
				fOut << RefDs[iF] << " " << DatDs[iF] << Endl;
		fOut.close();

		
		
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
		
		//! ******* Free
		
		MT.Free(DatSm, NDatSm*sizeof(double));
		MT.Free(DatDs, NDatDs*sizeof(double));
		MT.Free(RefDs, NDatDs*sizeof(double));
		
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "LC2ManipSensitivity: error: " << e.what()<<Endl;
		std::cerr << "LC2ManipSensitivity: abort!" << Endl;
		exit(1);
	}
	return(0);
};


// end of LISACODE-ManipSensitivity.cpp