// $Id:  Exp $
/*
 *  LISACODE-TestArmResp.cpp
 *
 *  Created by Antoine Petiteau on 23/02/11.
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
#include "LISACODE-GWSpinBBHHHarm1.h"
#include "LISACODE-DataFileWrite.h"
#include "LISACODE-OrbitsAnaLISA.h"
#include "LISACODE-OrbitsAnaMLDC.h"
#include "LISACODE-OrbitsData.h"
#include "LISACODE-OrbitsAnaOctahedron.h"
#include "LISACODE-ArmResp.h"

/**\ingroup GW
 * \{
 */


/** \brief Main of Code for testing gravitaional waves.
 * \author A. Petiteau
 * \version 2.0
 * \date 23/04/2011
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
		char fOut[1024];
		char fInXML[1024];
		double t0(0.0), t;
		double dt(5.);
		int NDat(MT.ifloor(31557600./dt)); //! One year
		
		//! *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution : Compute gravitational wave signal for all arms" << Endl;
			Coutm << "\t\t(./)LC2GWArmResp [Options] fOut fInXMLFile " << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << "\t\t * %sfOut       : Output file name : t hBp1 hBc1 hBp2 hBc2 ... [default: TestArmResp]. " << Endl;
			Coutm << "\t\t * %sfInXMLFile : XML input file [default: TestArmResp.xml]. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -s %dseed : Seed for random gennerator. [default : current time]"  << Endl ;
			Coutm << "\t\t * -t %fdtmes : Time step. [default : 1 ]"  << Endl ;
			Coutm << "\t\t * -T %fdur   : Duration. [default : 1000*dtMes ]"  << Endl ;
			Coutm << "\t\t * -dn \t\t: No screen display. [default: false]"  << Endl ;
			Coutm << "\t\t * -dl %sfile \t: Write standard output in a file. [default: no file]"  << Endl ;
			Coutm << "\t\t * -v \t\t: Verobse : display full details. [default: false]"  << Endl ;
			Coutm << Endl << " ----------------" << Endl;
			return 0;
			
		}
		
		//! *********** Version *************
		if(((argc>1)&&(strcmp(argv[1],"--version")==0))&&((argc>1)&&(strcmp(argv[1],"-v")==0))){
			Coutm << " ----- VERSION -----" << Endl;
			Coutm << " LC2GWArmResp : executable computing gravitational wave signal for the 6 arms - LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
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
			if((argc>1)&&(strcmp(argv[iarg],"-T")==0)){
				NDat = MT.ifloor(atof(argv[iarg+1])/dt);
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
		strcpy(fOut,"TestArmResp.txt");
		strcpy(fInXML,"/Users/petiteau/Applications/src/LISACode/LISACode_2_0/LISACode/Main/Test/Config-GWs_Test.xml");
		
		
		if(argc-nOptions>1){
			strcpy(fOut, argv[1+nOptions]);
		}
		
		
		
		
		
		//! ********************* Configuration
		
		//! ***** Output file
		double ** Rec;
		Coutm << "Write output in " << fOut << " ..." << Endl;
		LCDataFileWrite fOutTest(&MT, fOut, ASCII);
		fOutTest.sett0(t0);
		fOutTest.setdt(dt);
		fOutTest.setNDatExpect(NDat);
		

		
		//! *** Declaration
		LCGW ** GWs(NULL);
		int NGWs(0);
		LCOrbits * Orb(NULL);
		LCArmResp sGW(&MT);
		bool GWOptimizeTime(false);
		
		int iarg(1+nOptions);
		while (iarg<argc-1){
			iarg++;
			strcpy(fInXML, argv[iarg]);
			//! *** Configuration with XML file
			//! ** Read XML file
			ezxml_t tree, section, orbitdata;
			tree = ezxml_parse_file(fInXML);
			const char *type;
			for (section = ezxml_child(tree, "XSIL"); section; section = section->next) {
				type = ezxml_attr(section, "Type");
				if(MT.wcmp(ezxml_attr(section, "Type"),"LISAData")){
					for(orbitdata = ezxml_child(section, "XSIL"); orbitdata; orbitdata = orbitdata->next) {
						if(MT.wcmp(ezxml_attr(orbitdata,"Type"),"LISACode_Orbits")){
                            if (Orb != NULL)
                                throw std::invalid_argument("ERROR in LISACode::config : The orbit have already been created !");
                            Orb = new LCOrbitsAnaLISA(&MT);
                            Orb->config(orbitdata);
                        }
                        if(MT.wcmp(ezxml_attr(orbitdata,"Type"),"MLDC_Orbits")){
                            if (Orb != NULL)
                                throw std::invalid_argument("ERROR in LISACode::config : The orbit have already been created !");
                            Orb = new LCOrbitsAnaMLDC(&MT);
                            Orb->config(orbitdata);
                        }
                        if(MT.wcmp(ezxml_attr(orbitdata,"Type"),"OrbitsFile")){
                            Orb = new LCOrbitsData(&MT);
                            Orb->config(orbitdata);
                        }
                        if(MT.wcmp(ezxml_attr(orbitdata,"Type"),"Octahedron_FirstOrbits")){
                            if (Orb != NULL)
                                throw std::invalid_argument("ERROR in LISACode::config : The orbit have already been created !");
                            Orb = new LCOrbitsAnaOctahedron(&MT);
                            Orb->config(orbitdata);
                        }
					}
				}
				if(MT.wcmp(type,"SourceData")){
					ezxml_t gwdata; 
					for (gwdata = ezxml_child(section, "XSIL"); gwdata; gwdata = gwdata->next) {
						ezxml_t param;
						
						//! * Configure GW read in file
						if(MT.wcmp(ezxml_attr(gwdata,"Type"),"SampledPlaneWave")){
							NGWs++;
							GWs = (LCGW**) MT.ReAllocMemory(GWs, (NGWs-1)*sizeof(LCGW*), NGWs*sizeof(LCGW*));
							GWs[NGWs-1] = new LCGWFile(&MT);
							GWs[NGWs-1]->config(gwdata);
						}else{
							if(MT.wcmp(ezxml_attr(gwdata,"Type"),"PlaneWave")){
								
								for(param = ezxml_child(gwdata,"Param"); param; param = param->next){
									//! * Configure GW modeled in the code
									
									char * SourceType(NULL);
									if(MT.wcmp(ezxml_attr(param,"Name"),"SourceType")){
										MT.stripcopy((*param).txt, SourceType);
										if(MT.wcmp(SourceType, "Stochastic")){
											NGWs++;
											GWs = (LCGW**) MT.ReAllocMemory(GWs, (NGWs-1)*sizeof(LCGW*), NGWs*sizeof(LCGW*));
											GWs[NGWs-1] = new LCGWStochastic(&MT);
											GWs[NGWs-1]->config(gwdata);
											GWOptimizeTime = true;
										}
										if((MT.wcmp(SourceType, "GalacticBinaryV"))||(MT.wcmp(SourceType, "GalacticBinary"))){
											NGWs++;
											GWs = (LCGW**) MT.ReAllocMemory(GWs, (NGWs-1)*sizeof(LCGW*), NGWs*sizeof(LCGW*));
											GWs[NGWs-1] = new LCGWGalBin(&MT);
											GWs[NGWs-1]->config(gwdata);
										}
										if((MT.wcmp(SourceType, "SpinBBHHighHarm"))||(MT.wcmp(SourceType, "SpinBBHHH"))){
											NGWs++;
											GWs = (LCGW**) MT.ReAllocMemory(GWs, (NGWs-1)*sizeof(LCGW*), NGWs*sizeof(LCGW*));
											GWs[NGWs-1] = new LCGWSpinBBHHHarm1(&MT);
											GWs[NGWs-1]->config(gwdata);
										}
										
									}
									if(SourceType != NULL)
										MT.Free(SourceType, (strlen(SourceType)+1) * sizeof(char));
								}
							}
						}
					}
				}
			}
			ezxml_free(tree);
		}
		
		if(Orb == NULL)
			throw std::invalid_argument("ERROR in main : The orbits has not been defined !");
		
		
		
		//! ** Link the orbits 
		sGW.LinkOrb(Orb);
		sGW.LinkGWs(GWs, NGWs);
		
		
		
		//for (int iN=0; iN<NGWs; iN++)
		//	GWs[iN]->setTimeInfo(0., dt, 900., NDat*dt, 2.*17. );
		
		char NameRec[128];
		/*
		Rec = (double**) MT.AllocMemory(6*sizeof(double*));
		for (int iSC=0; iSC<3; iSC++){
			sprintf(NameRec, "Arm_%d->%d", iSC+1, 1+(iSC+1)%3);
			Rec[iSC] = fOutTest.AddRecord(NameRec);
		}
		for (int iSC=0; iSC<3; iSC++){
			sprintf(NameRec, "Arm-_%d->%d", iSC+1, 1+(iSC+2)%3);
			Rec[iSC+3] = fOutTest.AddRecord(NameRec);
		}*/
		Rec = (double**) MT.AllocMemory(Orb->getNArm()*sizeof(double*));
		int iArm(0);
		for (int iSCe=1; iSCe<=Orb->getNSC(); iSCe++){
			for (int iSCr=1; iSCr<=Orb->getNSC(); iSCr++){
				if(iSCe!=iSCr){
					sprintf(NameRec, "Arm_%d->%d", iSCe, iSCr);
					Rec[iArm++] = fOutTest.AddRecord(NameRec);
				}
			}
		}
		
		
		//! ********************* Initialization
		
		//! ** Initialization of orbits
		Orb->init();
		
		
		//! ** Initialization and setting time informations of gravitational wave 
		for (int iN=0; iN<NGWs; iN++){
			double tDOrbMax(900.);
			double tGWmin(t0-tDOrbMax);
			if(GWOptimizeTime){
				double tOrbMin(0.), tOrbMax(0.);
				Orb->tGWMinMax(GWs[iN]->getParam(0), GWs[iN]->getParam(1), t0, NDat*dt, dt, tOrbMin, tOrbMax);
				tGWmin = -tOrbMax;
				tDOrbMax = tOrbMax - tOrbMin;
				if(MT.Disp())
					Coutm << "Optimize GW time : tsmallest = " << tGWmin << " s ,  maximal time difference in orbits = " <<  tDOrbMax << " s." << Endl;
			}
			GWs[iN]->setTimeInfo(t0, dt, NDat*dt, tGWmin, tDOrbMax, 2.*Orb->getNominalArm() );
			GWs[iN]->init();
			Coutm << Endl;
			GWs[iN]->DispInfo("");
		}
		
		
		
		//! ** Initalization of arm response
		sGW.init();
		
		//! ** Initialization of output file
		fOutTest.init(NULL,0);
		
		
		//! ********************* Running
		std::cerr << "[........10........20........30........40........50........60........70........80........90.......100]" << Endl;
		std::cerr << "["; MT.o->flush();
		int iDisp = 0;
		iArm = 0;
		for(int iT=0; iT<NDat; iT++){
			t = t0 + iT*dt;
			//Coutm << t << Endl;
			/*
			for(int iSC=0; iSC<3; iSC++){
				(*Rec[iSC]) = sGW.gS(1+(iSC+1)%3, iSC+1, t);
				(*Rec[iSC+3]) = sGW.gS(1+(iSC+2)%3, iSC+1, t);
			}
			 */
			iArm=0;
			for (int iSCe=1; iSCe<=Orb->getNSC(); iSCe++){
				for (int iSCr=1; iSCr<=Orb->getNSC(); iSCr++){
					if(iSCe!=iSCr)
						(*Rec[iArm++]) = sGW.gS(iSCe, iSCr, t);
				}
			}
			
			iDisp++;
			if(iDisp>NDat/100){
				iDisp = 0;
				std::cerr << "=";
				fflush(stderr);
			}
			
			//! * Write the output
			fOutTest.RecordData();
		}
		std::cerr << "=]" << Endl;
		
		
		MT.Free(Rec, Orb->getNArm()*sizeof(double*));
		for (int iN=0; iN<NGWs; iN++)
			delete GWs[iN];
		MT.Free(GWs, NGWs*sizeof(LCGW*));
		delete Orb;
		
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "TestArmResp: error: " << e.what()<<Endl;
		std::cerr << "TestArmResp: abort!" << Endl;
		exit(1);
	}
	return(0);
};


/** \}*/

// end of