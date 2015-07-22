// $Id:  Exp $
/*
 *  LISACODE-TestTools.cpp
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
#include "LISACODE-Serie2.h"
#include "LISACODE-DataFileWrite.h"
#include "LISACODE-DataFileRead.h"

/**\ingroup ToolBox 
 * \{
 */


/** \brief Main of Code for testing toolbox.
 * \author A. Petiteau
 * \version 2.0
 * \date 14/02/2011
 *
 *
 */
int main (int argc, char * const argv[])
{
	try {
		LCTools MT;
		int nOptions(0);
		
		long SeedRand((long)time(NULL));
		
		// *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution :" << Endl;
			Coutm << "\t\t(./)LC2TestTools [Options] %sfOutBaseName " << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << "\t\t * %sfOutBaseName : Base for output file name [default: TestToolBox]. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -s %dseed : Seed for random gennerator. [default : current time]"  << Endl ;
			Coutm << "\t\t * -dn \t\t: No screen display. [default: false]"  << Endl ;
			Coutm << "\t\t * -dl %sfile \t: Write standard output in a file. [default: no file]"  << Endl ;
			Coutm << "\t\t * -v \t\t: Verobse : display full details. [default: false]"  << Endl ;
			Coutm << Endl << " ----------------" << Endl;
			return 0;
			
		}
		
		// *********** Version *************
		if(((argc>1)&&(strcmp(argv[1],"--version")==0))&&((argc>1)&&(strcmp(argv[1],"-v")==0))){
			Coutm << " ----- VERSION -----" << Endl;
			Coutm << " TestTools : executable testing toolbox of LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		for(int iarg=1; iarg<argc; iarg++){
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
		
		
		//! ***** Declaration of variable
		char fBOut[1024];
		
		//! ***** Initialization of variable
		strcpy(fBOut,"TestToolBox");
		
		
		if(argc-nOptions>1){
			strcpy(fBOut, argv[1+nOptions]);
		}
		
		
		
		double t0(-10.), dt(2.);
		int NDatInSerie(20);
		int NDat(10000);
		double Delay(5.e9/LC::c_SI);
		
		Coutm << "Serie of " << NDat << " data and record current and delay by " << Delay << " ..." << Endl;
		
		
		LCSerie2 SerieTest(&MT, t0, dt, NDatInSerie);
		
		char fNOutSerie[1024];
		sprintf(fNOutSerie, "%s-Serie.bin", fBOut);
		double * Rec0;
		double * RecD;
		LCDataFileWrite fOutSerie(&MT, fNOutSerie, XML);
		Rec0 = fOutSerie.AddRecord("Current");
		RecD = fOutSerie.AddRecord("Delay");
		fOutSerie.sett0(t0);
		fOutSerie.setdt(dt);
		fOutSerie.setNDatExpect(NDat);
		fOutSerie.init(NULL,0);
		
		
		char fNOutSerietxt[1024];
		sprintf(fNOutSerietxt, "%s-Serie.txt", fBOut);
		double * RecT0;
		double * RecTD;
		LCDataFileWrite fOutSerietxt(&MT, fNOutSerietxt, ASCII);
		RecT0 = fOutSerietxt.AddRecord("Current");
		RecTD = fOutSerietxt.AddRecord("Delay");
		fOutSerietxt.sett0(t0);
		fOutSerietxt.setdt(dt);
		fOutSerietxt.setNDatExpect(NDat);
		fOutSerietxt.init(NULL,0);
		
		//! *** Initialization of the serie
		for(int iD=0; iD<NDatInSerie; iD++) {
			double tmp = MT.RandGaussian(0.0, 10.0);
			SerieTest.addData(tmp);
		}
		
		//! *** Running
		for(int iD=0; iD<NDat; iD++) {
			double tmp = MT.RandGaussian(0.0, 10.0);
			SerieTest.addData(tmp);
			(*Rec0) = SerieTest.gData(0., LAG, 7);
			(*RecD) = SerieTest.gData(Delay, LAG, 7);
			fOutSerie.RecordData();
			(*RecT0) = (*Rec0);
			(*RecTD) = (*RecD);
			fOutSerietxt.RecordData();

		}
		
		//! *** Write the xml header
		char fNOutSerieXML[1024];
		sprintf(fNOutSerieXML, "%s-Serie.xml", fBOut);
		std::ofstream fOutSerieXML(fNOutSerieXML);
		
		
		fOutSerieXML << "<?xml version=\"1.0\"?>" << Endl;
		fOutSerieXML << "<!DOCTYPE XSIL SYSTEM \"http://www.vallis.org/lisa-xml.dtd\">" << Endl;
		fOutSerieXML << "<?xml-stylesheet type=\"text/xsl\" href=\"lisa-xml.xsl\"?>" << Endl;
		fOutSerieXML << "<XSIL>" << Endl;
		fOutSerieXML << "\t<Param Name=\"Author\">" << Endl;
        fOutSerieXML << "\t\tAntoine Petiteau" << Endl;
		fOutSerieXML << "\t</Param>" << Endl;
		fOutSerieXML << "\t<Param Name=\"GenerationDate\" Type=\"ISO-8601\">" << Endl;
		fOutSerieXML << "\t\t" << MT.TimeISO8601() << Endl;
		fOutSerieXML << "\t</Param>" << Endl;
		fOutSerieXML << "\t<XSIL Type=\"TestData\">" << Endl;
		fOutSerie.WriteXMLHeader(&fOutSerieXML,2);
		fOutSerieXML << "\t</XSIL>" << Endl;
		fOutSerieXML << "</XSIL>" << Endl;
		
		fOutSerieXML.close();
		
		
		
		//! ***** Read the data again for testing the I/O system
		
		//! ***** Read the data
		LCDataFileRead fInSerie(&MT, fNOutSerieXML);
		fInSerie.init();
		fInSerie.ControlDisplay();
		
		//! ***** Write
		char fNOutSerieASCII[1024];
		double t02(100.), dt2(2.);
		int NDat2(9000);
		sprintf(fNOutSerieASCII, "%s-Serie2.txt", fBOut);
		double * RecB0;
		double * RecBD;
		LCDataFileWrite fOutSerieASCII(&MT, fNOutSerieASCII, ASCII);
		RecB0 = fOutSerieASCII.AddRecord("Current");
		RecBD = fOutSerieASCII.AddRecord("Delay");
		fOutSerieASCII.sett0(t02);
		fOutSerieASCII.setdt(dt2);
		fOutSerieASCII.setNDatExpect(NDat2);
		fOutSerieASCII.init(NULL,0);
		
		//! *** Read
		double t;
		for(int iD=0; iD<NDat2; iD++) {
			t = t02 + iD*dt2;
			(*RecB0) = fInSerie.gData(0, t, TRU, 7);
			(*RecBD) = fInSerie.gData(0, t-Delay, TRU, 7);
			fOutSerieASCII.RecordData();
		}
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
	}
	
	
	
	
	
	catch(std::exception & e ) {
		std::cerr << "TestTools: error: " << e.what()<<Endl;
		std::cerr << "TestTools: abort!" << Endl;
		exit(1);
	}
	return(0);
};


/** \}*/

// end of