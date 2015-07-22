// $Id:  Exp $
/*
 *  LISACODE-TestNoises.cpp
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
#include "LISACODE-Noise.h"
#include "LISACODE-NoiseWhite.h"
#include "LISACODE-NoiseFilter.h"
#include "LISACODE-NoiseFile.h"
#include "LISACODE-DataFileWrite.h"

/**\ingroup Noise 
 * \defgroup TestNoises TestNoises
 * \{
 */


/** \brief Main of Code for testing noises.
 * \author A. Petiteau
 * \version 2.0
 * \date 08/04/2011
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
		double dtMes(0.3);
		int NDat(10000);
		double tDur(NDat*dtMes);
		double Delay(16.66);
		
		//! *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution :" << Endl;
			Coutm << "\t\t(./)LC2TestNoise [Options] %sfOutBaseName %sfInXMLFile " << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << "\t\t * %sfOutBaseName : Base for output file name [default: TestNoises]. " << Endl;
			Coutm << "\t\t * %sfInXMLFile : XML input file [default: TestNoises.xml]. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -s %dseed : Seed for random gennerator. [default : current time]"  << Endl ;
			Coutm << "\t\t * -tp %fdtphy : Physical time step. [default : 0.3 ]"  << Endl ;
			Coutm << "\t\t * -tm %fdtmes : Measurement time step. [default : 0.3 ]"  << Endl ;
			Coutm << "\t\t * -T %fdur    : Duration. [default : 1000*dtMes ]"  << Endl ;
			Coutm << "\t\t * -dn \t\t: No screen display. [default: false]"  << Endl ;
			Coutm << "\t\t * -dl %sfile \t: Write standard output in a file. [default: no file]"  << Endl ;
			Coutm << "\t\t * -v \t\t: Verobse : display full details. [default: false]"  << Endl ;
			Coutm << Endl << " ----------------" << Endl;
			return 0;
			
		}
		
		//! *********** Version *************
		if(((argc>1)&&(strcmp(argv[1],"--version")==0))&&((argc>1)&&(strcmp(argv[1],"-v")==0))){
			Coutm << " ----- VERSION -----" << Endl;
			Coutm << " LC2TestNoise : executable testing noises of LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
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
			if((argc>1)&&(strcmp(argv[iarg],"-T")==0)){
				tDur = atof(argv[iarg+1]);
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
		NDat = MT.ifloor(tDur/dtMes);
		Coutm << "Number of data = " << NDat << Endl;
		MT.setRandSeed(SeedRand);
		
		
		//! ***** Initialization of variable
		strcpy(fBOut,"TestNoises");
		strcpy(fInXML,"Config-Noises_std.xml");
		
		
		if(argc-nOptions>1){
			strcpy(fBOut, argv[1+nOptions]);
		}
		
		if(argc-nOptions>2){
			strcpy(fInXML, argv[2+nOptions]);
		}
		

		
		//! ********************* Configuration
		
		//! ***** Output file
		char fNOutNoise[1024];
		double ** Rec0;
		double ** RecD;
		sprintf(fNOutNoise, "%s.txt", fBOut);
		Coutm << "Write output in " << fNOutNoise << " ..." << Endl;
		LCDataFileWrite fOutTest(&MT, fNOutNoise, ASCII);
		fOutTest.sett0(t0);
		fOutTest.setdt(dtMes);
		fOutTest.setNDatExpect(NDat);
		
		
		
		//! ***** Noise
		
		//! *** Declaration
		LCNoise ** Noises;
		int NNoises(0);
		/*
		Noises = (LCNoise**)MT.AllocMemory(sizeof(LCNoise*));
		Noises[0] = new LCNoiseWhite(&MT);
		
		//! *** Configuration with given data
		//! ** Set square root of PSD
		Noises[0]->config(0, 1.e-13);
		//! ** Set the name
		Noises[0]->setName("x1s");
		//! ** Set the time informations
		Noises[0]->setTimeInfo(dtMes, dtPhy, 0., 30.);
		NNoises++;
		*/
		
		//! *** Configuration with XML file
		//! ** Read XML file
		ezxml_t tree, section;
		tree = ezxml_parse_file(fInXML);
		const char *type;
		for (section = ezxml_child(tree, "XSIL"); section; section = section->next) {
			type = ezxml_attr(section, "Type");
			if(MT.wcmp(type,"NoiseData")){
				ezxml_t noisedata; 
				for (noisedata = ezxml_child(section, "XSIL"); noisedata; noisedata = noisedata->next) {
					ezxml_t noisedata;
                    for (noisedata = ezxml_child(section, "XSIL"); noisedata; noisedata = noisedata->next) {
                        ezxml_t param;
                        char * SourceType(NULL);
                        char * SpectralType(NULL);
                        
                        
                        for(param = ezxml_child(noisedata,"Param"); param; param = param->next){
                            if(MT.wcmp(ezxml_attr(param,"Name"),"SourceType"))
                                MT.stripcopy((*param).txt, SourceType);
                            if(MT.wcmp(ezxml_attr(param,"Name"),"SpectralType"))
                                MT.stripcopy((*param).txt, SpectralType);
                        }
                        
                        if((MT.wcmp(SourceType, "PseudoRandomNoise"))&&(MT.wcmp(SpectralType, "White"))){
                            NNoises++;
                            Noises = (LCNoise**) MT.ReAllocMemory(Noises, (NNoises-1)*sizeof(LCNoise*), NNoises*sizeof(LCNoise*));
                            Noises[NNoises-1] = new LCNoiseWhite(&MT);
                            Noises[NNoises-1]->config(noisedata);
                        }
                        
                        if((MT.wcmp(SourceType, "PseudoRandomNoise"))&&((MT.wcmp(SpectralType, "Filter_f"))
                                                                         ||(MT.wcmp(SpectralType, "Filter_1of"))
                                                                         ||(MT.wcmp(SpectralType, "WhitePhase"))
                                                                         ||(MT.wcmp(SpectralType, "BlueFrequency"))
                                                                         ||(MT.wcmp(SpectralType, "WhiteFrequency"))
                                                                         ||(MT.wcmp(SpectralType, "RedFrequency"))
                                                                         ||(MT.wcmp(SpectralType, "PinkFrequency"))
                                                                         ||(MT.wcmp(SpectralType, "PinkAcceleration"))
                                                                         ||(MT.wcmp(SpectralType, "Filter_1of_1of2"))
                                                                         ||(MT.wcmp(SpectralType, "Filter_1of_1of32"))
                                                                         ||(MT.wcmp(SpectralType, "PreStabLaserNoiseFreq")))){
                            NNoises++;
                            Noises = (LCNoise**) MT.ReAllocMemory(Noises, (NNoises-1)*sizeof(LCNoise*), NNoises*sizeof(LCNoise*));
                            Noises[NNoises-1] = new LCNoiseFilter(&MT);
                            Noises[NNoises-1]->config(noisedata);
                        }
                        
                        if(MT.wcmp(SourceType, "DataFile")){
                            NNoises++;
                            Noises = (LCNoise**) MT.ReAllocMemory(Noises, (NNoises-1)*sizeof(LCNoise*), NNoises*sizeof(LCNoise*));
                            Noises[NNoises-1] = new LCNoiseFile(&MT);
                            Noises[NNoises-1]->config(noisedata);
                        }
                        
                        //! ** HERE, Add new type of noise
                        
                        if(SpectralType != NULL)
                            MT.Free(SpectralType, (strlen(SpectralType)+1) * sizeof(char));
                        if(SourceType != NULL)
                            MT.Free(SourceType, (strlen(SourceType)+1) * sizeof(char));
                    }
				}
			}
        }
        ezxml_free(tree);
		
		
		//! ** Set the time informations
		for (int iN=0; iN<NNoises; iN++)
			Noises[iN]->setTimeInfo(dtMes, dtPhy, -10., 100.);
		
        //! ** Check display
        for (int iN=0; iN<NNoises; iN++)
            Noises[iN]->DispInfo("\t");
        
        
		char NameRec0[128], NameRecD[128];
		Rec0 = (double**) MT.AllocMemory(NNoises*sizeof(double*));
		RecD = (double**) MT.AllocMemory(NNoises*sizeof(double*));
		for (int iN=0; iN<NNoises; iN++){
			sprintf(NameRec0, "Noise%d-current", iN);
			sprintf(NameRecD, "Noise%d-delay", iN);
			Rec0[iN] = fOutTest.AddRecord(NameRec0);
			RecD[iN] = fOutTest.AddRecord(NameRecD);
		}
		
		
		
		//! ********************* Initialization
		
		
		//! ** Initialization of noises
		for (int iN=0; iN<NNoises; iN++){
			Noises[iN]->init();
			Coutm << Endl;
			Noises[iN]->DispInfo("");
		}
		
		//! ** Initialization of output file
		fOutTest.init(NULL,0);
		
		
		//! ********************* Running
		for(int iT=0; iT<NDat; iT++){
			t = t0 + iT*dtMes;
			//Coutm << t << Endl;
			for (int iN=0; iN<NNoises; iN++){
				//! * Run one step of noise
				Noises[iN]->RunStep();
				//Noises[iN]->DispData();
				//! * Copy the result in output
				(*Rec0[iN]) = Noises[iN]->gN(0.);
				(*RecD[iN]) = Noises[iN]->gN(Delay);
			}
			//! * Write the output
			fOutTest.RecordData();
		}
		
		for (int iN=0; iN<NNoises; iN++)
			delete Noises[iN];
		MT.Free(Noises, NNoises*sizeof(LCNoise*));
		MT.Free(Rec0, NNoises*sizeof(double*));
		MT.Free(RecD, NNoises*sizeof(double*));
		ezxml_free(tree);
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "TestNoises: error: " << e.what()<<Endl;
		std::cerr << "TestNoises: abort!" << Endl;
		exit(1);
	}
	return(0);
};


/** \}*/

// end of