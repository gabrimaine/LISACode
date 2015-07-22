// $Id:  Exp $
/*
 *  LISACODE-TestUSO.cpp
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
#include "LISACODE-USO.h"
#include "LISACODE-USONoiseDeriv.h"
#include "LISACODE-DataFileWrite.h"

/**\ingroup Noise 
 * \defgroup TestUSO TestUSO
 * \{
 */


/** \brief Main of Code for testing USO clocks.
 * \author A. Petiteau
 * \version 2.0
 * \date 12/04/2011
 *
 *
 */
int main (int argc, char * const argv[])
{
	try {
		LCTools MT;
		int nOptions(0);
		
		//! ********* Declaration of varaiables 
		long SeedRand((long)time(NULL));
		char fBOut[1024];
		char fInXML[1024];
		double t0(0.0), t;
		double dtPhy(0.3);
		double dtMes(3.0);
		int NDat(10000);
		double Delay(5.e9/LC::c_SI);
		
		//! *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution :" << Endl;
			Coutm << "\t\t(./)LC2TestUSO [Options] " << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << "\t\t * %sfInXMLFile : XML input file [default: TestUSO.xml]. " << Endl;
			Coutm << "\t\t * %sfOutBaseName : Base for output file name [default: TestUSO]. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -s %dseed : Seed for random gennerator. [default : current time]"  << Endl ;
			Coutm << "\t\t * -tp %fdtphy : Physical time step. [default : 0.1 ]"  << Endl ;
			Coutm << "\t\t * -tm %fdtmes : Measurement time step. [default : 3 ]"  << Endl ;
			Coutm << "\t\t * -dn \t\t: No screen display. [default: false]"  << Endl ;
			Coutm << "\t\t * -dl %sfile \t: Write standard output in a file. [default: no file]"  << Endl ;
			Coutm << "\t\t * -v \t\t: Verobse : display full details. [default: false]"  << Endl ;
			Coutm << Endl << " ----------------" << Endl;
			return 0;
			
		}
		
		//! *********** Version *************
		if(((argc>1)&&(strcmp(argv[1],"--version")==0))&&((argc>1)&&(strcmp(argv[1],"-v")==0))){
			Coutm << " ----- VERSION -----" << Endl;
			Coutm << " TestUSO : executable testing USO clock of LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		for(int iarg=1; iarg<argc; iarg++){
			if((argc>1)&&(strcmp(argv[iarg],"-s")==0)){
				SeedRand = atoi(argv[iarg+1]); 
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-tp")==0)){
				dtPhy = atof(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-tm")==0)){
				dtMes = atof(argv[iarg+1]);
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
		strcpy(fBOut,"TestUSO");
		strcpy(fInXML,"/Users/petiteau/Data/LISACode/TestLC20/Test3Noise.xml");
		
		
		if(argc-nOptions>1){
			strcpy(fBOut, argv[1+nOptions]);
		}
		
		if(argc-nOptions>2){
			strcpy(fInXML, argv[2+nOptions]);
		}
		
		
		
		//! ********************* Configuration
		
		//! ***** Output file
		char fNOutUSO[1024];
		double ** RecT;
		double ** RecN0;
		double ** RecND;
		sprintf(fNOutUSO, "%s.txt", fBOut);
		Coutm << "Write output in " << fNOutUSO << " ..." << Endl;
		LCDataFileWrite fOutTest(&MT, fNOutUSO, ASCII);
		fOutTest.sett0(t0);
		fOutTest.setdt(dtMes);
		fOutTest.setNDatExpect(NDat);
		
		
		
		//! ***** Noise
		
		//! *** Declaration
		LCUSO ** USOs;
		int NbUSOs(0);
		
		
		//! *** Configuration with XML file
		//! ** Read XML file
		ezxml_t tree, section;
		tree = ezxml_parse_file(fInXML);
		const char *type;
		for (section = ezxml_child(tree, "XSIL"); section; section = section->next) {
			type = ezxml_attr(section, "Type");
			if(MT.wcmp(type,"LISACode")){
				ezxml_t usodata; 
				for (usodata = ezxml_child(section, "XSIL"); usodata; usodata = usodata->next) {
					type = ezxml_attr(usodata, "Type");
					//Coutm << "type  =  " << type << Endl;
					if(MT.wcmp(type,"USO")){
						ezxml_t param;
						char * USOType(NULL);
						for(param = ezxml_child(usodata,"Param"); param; param = param->next){
							//Coutm << "ezxml_attr(param,\"Name\")  =  " << ezxml_attr(param,"Name") << Endl;
							if(MT.wcmp(ezxml_attr(param,"Name"),"USOType")){
								MT.stripcopy((*param).txt, USOType);
								if(MT.wcmp(USOType, "NoiseDeriv")){
									NbUSOs++;
									USOs = (LCUSO**) MT.ReAllocMemory(USOs, (NbUSOs-1)*sizeof(LCUSO*), NbUSOs*sizeof(LCUSO*));
									USOs[NbUSOs-1] = new LCUSONoiseDeriv(&MT);
									USOs[NbUSOs-1]->config(usodata);
									
									//USOs[NbUSOs-1]->setPhaseOutput();
									USOs[NbUSOs-1]->setFreqOutput();
								}
							}
						}
						if(USOType!=NULL)
							MT.Free(USOType, (strlen(USOType)+1) * sizeof(char));
					}
				}
			}
		}
		ezxml_free(tree);
		
		//! ** Set the time informations
		for (int iN=0; iN<NbUSOs; iN++)
			USOs[iN]->setTimeInfo(dtMes, dtPhy, 0., 30.);
		
		char NameRecT[128], NameRecN0[128], NameRecND[128];
		RecT = (double**) MT.AllocMemory(NbUSOs*sizeof(double*));
		RecN0 = (double**) MT.AllocMemory(NbUSOs*sizeof(double*));
		RecND = (double**) MT.AllocMemory(NbUSOs*sizeof(double*));
		for (int iN=0; iN<NbUSOs; iN++){
			sprintf(NameRecT, "TimeShift%d", iN);
			sprintf(NameRecN0, "CurrentNoise%d", iN);
			sprintf(NameRecND, "DelayedNoise%d", iN);
			RecT[iN] = fOutTest.AddRecord(NameRecT);
			RecN0[iN] = fOutTest.AddRecord(NameRecN0);
			RecND[iN] = fOutTest.AddRecord(NameRecND);
		}
		
		
		
		//! ********************* Initialization
		
		
		//! ** Initialization of noises
		for (int iN=0; iN<NbUSOs; iN++){
			USOs[iN]->init();
			Coutm << Endl;
			USOs[iN]->DispInfo("");
		}
		
		//! ** Initialization of output file
		fOutTest.init(NULL,0);
		
		
		//! ********************* Running
		for(int iT=0; iT<NDat; iT++){
			t = t0 + iT*dtMes;
			for (int iN=0; iN<NbUSOs; iN++){
				//! * Run one step of noise
				USOs[iN]->RunStep(t);
				//! * Copy the result in output
				(*RecT[iN]) = USOs[iN]->gT(t);
				(*RecN0[iN]) = USOs[iN]->gN(0.);
				(*RecND[iN]) = USOs[iN]->gN(0.);
			}
			//! * Write the output
			fOutTest.RecordData();
		}
		
		for (int iN=0; iN<NbUSOs; iN++)
			delete USOs[iN];
		MT.Free(USOs, NbUSOs*sizeof(LCUSO*));
		MT.Free(RecT, NbUSOs*sizeof(double*));
		MT.Free(RecN0, NbUSOs*sizeof(double*));
		MT.Free(RecND, NbUSOs*sizeof(double*));
		ezxml_free(tree);
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "TestUSO: error: " << e.what()<<Endl;
		std::cerr << "TestUSO: abort!" << Endl;
		exit(1);
	}
	return(0);
};


/** \}*/

// end of