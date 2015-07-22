// $Id:  Exp $
/*
 *  LISACODE-TestOrbits.cpp
 *
 *  Created by Antoine Petiteau on 10/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "LISACODE-OrbitsAnaLISA.h"
#include "LISACODE-OrbitsAnaMLDC.h"
#include "LISACODE-OrbitsData.h"
#include "LISACODE-OrbitsAnaOctahedron.h"
#include "LISACODE-DataFileWrite.h"

/**\ingroup Orbits
 * \defgroup TestOrbits TestOrbits
 * \{
 */


/** \brief Main of Code for testing orbits.
 * \author A. Petiteau
 * \version 2.0
 * \date 10/04/2011
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
		double dt(3600.0);
        double TObs(63115200.);
		int NDat(TObs/dt);
		double Delay(5.e9/LC::c_SI);
		bool ExeRun2(false);
		
		//! *********** Help *************
		if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
			Coutm << " ----- HELP -----" << Endl;
			Coutm << Endl << "\tExecution : generation of orbits" << Endl;
			Coutm << "\t\t(./)LC2Orbits [Options] %sfOutBaseName %sfInXMLFile" << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << "\t\t * %sfOutBaseName : Base for output file name [default: TestOrbits]. " << Endl;
			Coutm << "\t\t * %sfInXMLFile : XML input file [default: TestOrbits.xml]. " << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -s %dseed : Seed for random gennerator. [default : current time]"  << Endl ;
			Coutm << "\t\t * -dt %fdt : Time step (in sec). [default : 3600.0 ]"  << Endl ;
			Coutm << "\t\t * -T %fT : Duration (in sec). [default : 63115200.0 = 2 years ]"  << Endl ;
			Coutm << "\t\t * -r2 : Test. [default : 3600.0 ]"  << Endl ;
			Coutm << "\t\t * -dn \t\t: No screen display. [default: false]"  << Endl ;
			Coutm << "\t\t * -dl %sfile \t: Write standard output in a file. [default: no file]"  << Endl ;
			Coutm << "\t\t * -v \t\t: Verobse : display full details. [default: false]"  << Endl ;
			Coutm << Endl << " ----------------" << Endl;
			return 0;
			
		}
		
		//! *********** Version *************
		if(((argc>1)&&(strcmp(argv[1],"--version")==0))&&((argc>1)&&(strcmp(argv[1],"-v")==0))){
			Coutm << " ----- VERSION -----" << Endl;
			Coutm << " LC2Orbits : executable for generating orbits - LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		for(int iarg=1; iarg<argc; iarg++){
			if((argc>1)&&(strcmp(argv[iarg],"-s")==0)){
				SeedRand = atoi(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-dt")==0)){
				dt = atof(argv[iarg+1]);
				nOptions +=2;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-r2")==0)){
				ExeRun2 = true;
				nOptions ++;
			}
			if((argc>1)&&(strcmp(argv[iarg],"-T")==0)){
				TObs = atof(argv[iarg+1]);
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
		
        
		NDat = MT.ifloor(TObs/dt);
		
		//! ***** Initialization of variable
		strcpy(fBOut,"TestOrb");
		strcpy(fInXML,"/Users/petiteau/Applications/src/LISACode/LISACode_2_0/LISACode/ConfigFiles/ConfigNewLISA/ToyLISALike5/Config-Orbits_std.xml");
		
		
		if(argc-nOptions>1){
			strcpy(fBOut, argv[1+nOptions]);
		}
		
		if(argc-nOptions>2){
			strcpy(fInXML, argv[2+nOptions]);
		}
		
		
		
		//! ********************* Configuration
		
		//! ***** Output file
		char fNOutOrbPos[1024];
		char fNOutOrbArm[1024];
		double ** RecPos;
		double ** RecArm;
		
		sprintf(fNOutOrbPos, "%s-Pos.txt", fBOut);
		Coutm << "Write position in " << fNOutOrbPos << " ..." << Endl;
		LCDataFileWrite fOutPos(&MT, fNOutOrbPos, ASCII);
		fOutPos.sett0(t0);
		fOutPos.setdt(dt);
		fOutPos.setNDatExpect(NDat);
		
		sprintf(fNOutOrbArm, "%s-Arm.txt", fBOut);
		Coutm << "Write position in " << fNOutOrbArm << " ..." << Endl;
		LCDataFileWrite fOutArm(&MT, fNOutOrbArm, ASCII);
		fOutArm.sett0(t0);
		fOutArm.setdt(dt);
		fOutArm.setNDatExpect(NDat);
		
		
		//! ***** Orbits
		
		//! *** Declaration
		LCOrbits ** Orb;
		int NOrb(0);
		
		
		//! *** Configuration with XML file
		//! ** Read XML file
		Coutm << "Configuration using " << fInXML << " ... " << Endl; 
		ezxml_t tree, section;
		tree = ezxml_parse_file(fInXML);
		const char *type;
		for (section = ezxml_child(tree, "XSIL"); section; section = section->next) {
			type = ezxml_attr(section, "Type");
			if(MT.wcmp(type,"LISAData")){
				ezxml_t orbitdata; 
				for (orbitdata = ezxml_child(section, "XSIL"); orbitdata; orbitdata = orbitdata->next) {
					if(MT.wcmp(ezxml_attr(orbitdata,"Type"),"LISACode_Orbits")){
						NOrb++;
						Orb = (LCOrbits**) MT.ReAllocMemory(Orb, (NOrb-1)*sizeof(LCOrbits*), NOrb*sizeof(LCOrbits*));
						Orb[NOrb-1] = new LCOrbitsAnaLISA(&MT);
						Orb[NOrb-1]->config(orbitdata);
					}
                    if(MT.wcmp(ezxml_attr(orbitdata,"Type"),"MLDC_Orbits")){
						NOrb++;
						Orb = (LCOrbits**) MT.ReAllocMemory(Orb, (NOrb-1)*sizeof(LCOrbits*), NOrb*sizeof(LCOrbits*));
						Orb[NOrb-1] = new LCOrbitsAnaMLDC(&MT);
						Orb[NOrb-1]->config(orbitdata);
					}
					if(MT.wcmp(ezxml_attr(orbitdata,"Type"),"OrbitsFile")){
						NOrb++;
						Orb = (LCOrbits**) MT.ReAllocMemory(Orb, (NOrb-1)*sizeof(LCOrbits*), NOrb*sizeof(LCOrbits*));
						Orb[NOrb-1] = new LCOrbitsData(&MT);
						Orb[NOrb-1]->config(orbitdata);
					}
					if(MT.wcmp(ezxml_attr(orbitdata,"Type"),"Octahedron_FirstOrbits")){
						NOrb++;
						Orb = (LCOrbits**) MT.ReAllocMemory(Orb, (NOrb-1)*sizeof(LCOrbits*), NOrb*sizeof(LCOrbits*));
						Orb[NOrb-1] = new LCOrbitsAnaOctahedron(&MT);
						Orb[NOrb-1]->config(orbitdata);
					}
				}
			}
		}
		ezxml_free(tree);
		
		if(NOrb == 0)
			throw std::invalid_argument("No orbit description found in te xml file !");
		
		int NSCTotCur(0), NArmTotCur(0);
		
		for(int iO=0; iO<NOrb; iO++){
			NSCTotCur += Orb[iO]->getNSC();
			NArmTotCur += Orb[iO]->getNArm();
		}
		
		//! ** Set the time informations
		char TmpName[128];
		RecPos = (double**) MT.AllocMemory(NOrb*3*NSCTotCur*sizeof(double*));
		NSCTotCur = 0;
		for(int iO=0; iO<NOrb; iO++){
			for(int iSC=0; iSC<Orb[iO]->getNSC(); iSC++){
				sprintf(TmpName, "x%d-%d", iSC+1, iO);
				RecPos[3*NSCTotCur + 3*iSC + 0] = fOutPos.AddRecord(TmpName);
				sprintf(TmpName, "y%d-%d", iSC+1, iO);
				RecPos[3*NSCTotCur + 3*iSC + 1] = fOutPos.AddRecord(TmpName);
				sprintf(TmpName, "z%d-%d", iSC+1, iO);
				RecPos[3*NSCTotCur + 3*iSC + 2] = fOutPos.AddRecord(TmpName);
			}
			NSCTotCur += Orb[iO]->getNSC();
		}
		RecArm = (double**) MT.AllocMemory(NOrb*NArmTotCur*sizeof(double*));
		NArmTotCur = 0;
		for (int iO=0; iO<NOrb; iO++){
			for(int iArm=0; iArm<Orb[iO]->getNArm(); iArm++){
				sprintf(TmpName, "L%d-%d", iArm, iO);
				RecArm[NArmTotCur+iArm] = fOutArm.AddRecord(TmpName);
			}
			NArmTotCur += Orb[iO]->getNArm();
		}
		
		
		
		//! ********************* Initialization
		
		
		//! ** Initialization of noises
		for (int iO=0; iO<NOrb; iO++){
			Orb[iO]->init();
			Coutm << Endl;
			Orb[iO]->DispInfo("");
		}
		
		//! ** Initialization of output file
		fOutPos.init(NULL,0);
		fOutArm.init(NULL,0);
		
		
		//! ********************* Running
		Coutm << "Running (" << NDat << " steps) ... " << Endl;
		LCVector tmpP(&MT);
		for(int iT=0; iT<NDat; iT++){
			t = t0 + iT*dt;
			NSCTotCur = 0;
			NArmTotCur = 0;
			for (int iO=0; iO<NOrb; iO++){
				for (int iSC=1; iSC<=Orb[iO]->getNSC(); iSC++) {
					tmpP = Orb[iO]->Pos(iSC, t);
					(*RecPos[NSCTotCur+3*(iSC-1)+0]) = tmpP(0);
					(*RecPos[NSCTotCur+3*(iSC-1)+1]) = tmpP(1);
					(*RecPos[NSCTotCur+3*(iSC-1)+2]) = tmpP(2);
				}
				if(Orb[iO]->getNSC()==3){
					(*RecArm[NArmTotCur+0]) = Orb[iO]->Arm(3, 2, t);
					(*RecArm[NArmTotCur+1]) = Orb[iO]->Arm(1, 3, t);
					(*RecArm[NArmTotCur+2]) = Orb[iO]->Arm(2, 1, t);
					(*RecArm[NArmTotCur+3]) = Orb[iO]->Arm(2, 3, t);
					(*RecArm[NArmTotCur+4]) = Orb[iO]->Arm(3, 1, t);
					(*RecArm[NArmTotCur+5]) = Orb[iO]->Arm(1, 2, t);
				}else{
					int iArm(0);
					for(int iem=1; iem<=Orb[iO]->getNSC(); iem++)
						for(int ire=1; ire<=Orb[iO]->getNSC(); ire++)
							if(iem!=ire)
								(*RecArm[NArmTotCur+iArm++]) = Orb[iO]->Arm(iem, ire, t);
				}
				NSCTotCur += Orb[iO]->getNSC();
				NArmTotCur += Orb[iO]->getNArm();
			}
			//! * Write the output
			fOutPos.RecordData();
			fOutArm.RecordData();
		}
		
		
		
		
		
		
		//! *** Write the xml header
		char fNOutOrbXML[1024];
		sprintf(fNOutOrbXML, "%s-Arm.xml", fBOut);
		std::ofstream fOutOrbXML(fNOutOrbXML);
		
		
		fOutOrbXML << "<?xml version=\"1.0\"?>" << Endl;
		fOutOrbXML << "<!DOCTYPE XSIL SYSTEM \"http://www.vallis.org/lisa-xml.dtd\">" << Endl;
		fOutOrbXML << "<?xml-stylesheet type=\"text/xsl\" href=\"lisa-xml.xsl\"?>" << Endl;
		fOutOrbXML << "<XSIL>" << Endl;
		fOutOrbXML << "\t<Param Name=\"Author\">" << Endl;
        fOutOrbXML << "\t\tAntoine Petiteau" << Endl;
		fOutOrbXML << "\t</Param>" << Endl;
		fOutOrbXML << "\t<Param Name=\"GenerationDate\" Type=\"ISO-8601\">" << Endl;
		fOutOrbXML << "\t\t" << MT.TimeISO8601() << Endl;
		fOutOrbXML << "\t</Param>" << Endl;
		fOutOrbXML << "\t<XSIL Type=\"LISAData\">" << Endl;
		fOutOrbXML << "\t\t<XSIL Name=\"LISACode orbits\" Type=\"OrbitsFile\">" << Endl;
		fOutPos.WriteXMLHeader(&fOutOrbXML,3);
		fOutOrbXML << "\t\t</XSIL>" << Endl;
		fOutOrbXML << "\t</XSIL>" << Endl;
		fOutOrbXML << "</XSIL>" << Endl;
		
		fOutOrbXML.close();
		
		
		
		//! ************* Free memory of the first run
		
		for (int iO=0; iO<NOrb; iO++)
			delete Orb[iO];
		MT.Free(Orb, NOrb*sizeof(LCOrbits*));
		MT.Free(RecPos, NOrb*9*sizeof(double*));
		MT.Free(RecArm, NOrb*6*sizeof(double*));
		
		
		//! ***** If we are in a LISA case we also test the reading of data file 
		
		if(NSCTotCur==3){
			
			if(ExeRun2){
				
				//! ********************* Read the orbits in the file for testing orbit file reader
				
				double dt2(2000.0);
				double t02(t0+2*dt2);
				int NDat2(MT.iceil((NDat*dt-5*dt2)/dt2));
				LCOrbits * Orb2;
				
				
				// ** Prepare output files
				
				sprintf(fNOutOrbPos, "%s-Pos-2.txt", fBOut);
				Coutm << "Write position in " << fNOutOrbPos << " ..." << Endl;
				LCDataFileWrite fOutPos2(&MT, fNOutOrbPos, ASCII);
				fOutPos2.sett0(t02);
				fOutPos2.setdt(dt2);
				fOutPos2.setNDatExpect(NDat2);
				
				sprintf(fNOutOrbArm, "%s-Arm-2.txt", fBOut);
				Coutm << "Write position in " << fNOutOrbArm << " ..." << Endl;
				LCDataFileWrite fOutArm2(&MT, fNOutOrbArm, ASCII);
				fOutArm2.sett0(t02);
				fOutArm2.setdt(dt2);
				fOutArm2.setNDatExpect(NDat2);
				
				char TmpName2[128];
				double ** RecPos2;
				double ** RecArm2;
				RecPos2 = (double**) MT.AllocMemory(9*sizeof(double*));
				for(int iSC=0; iSC<3; iSC++){
					sprintf(TmpName2, "x%d", iSC+1);
					RecPos2[3*iSC+0] = fOutPos2.AddRecord(TmpName2);
					sprintf(TmpName2, "y%d", iSC+1);
					RecPos2[3*iSC+1] = fOutPos2.AddRecord(TmpName2);
					sprintf(TmpName2, "z%d", iSC+1);
					RecPos2[3*iSC+2] = fOutPos2.AddRecord(TmpName2);
				}
				RecArm2 = (double**) MT.AllocMemory(6*sizeof(double*));
				RecArm2[0] = fOutArm2.AddRecord("L1");
				RecArm2[1] = fOutArm2.AddRecord("L2");
				RecArm2[2] = fOutArm2.AddRecord("L3");
				RecArm2[3] = fOutArm2.AddRecord("L1'");
				RecArm2[4] = fOutArm2.AddRecord("L2'");
				RecArm2[5] = fOutArm2.AddRecord("L3'");
				
				fOutPos2.init(NULL,0);
				fOutArm2.init(NULL,0);
				
				
				
				//! *** Read XML file
				ezxml_t tree2, section2;
				tree2 = ezxml_parse_file(fNOutOrbXML);
				const char *type2;
				for (section2 = ezxml_child(tree2, "XSIL"); section2; section2 = section2->next) {
					type2 = ezxml_attr(section2, "Type");
					if(MT.wcmp(type2,"LISAData")){
						ezxml_t orbitdata2; 
						for (orbitdata2 = ezxml_child(section2, "XSIL"); orbitdata2; orbitdata2 = orbitdata2->next) {
							if(MT.wcmp(ezxml_attr(orbitdata2,"Type"),"OrbitsFile")){
								Orb2 = new LCOrbitsData(&MT);
								Orb2->config(orbitdata2);
							}
						}
					}
				}
				ezxml_free(tree2);
                
				Orb2->init();
				
				
				
				//! ********************* Running
				Coutm << "Running ..." << Endl;
				LCVector tmpP2(&MT);
				for(int iT=0; iT<NDat2; iT++){
					t = t02 + iT*dt2;
					for (int iSC=1; iSC<4; iSC++) {
						tmpP2 = Orb2->Pos(iSC, t);
						(*RecPos2[3*(iSC-1)+0]) = tmpP2(0);
						(*RecPos2[3*(iSC-1)+1]) = tmpP2(1);
						(*RecPos2[3*(iSC-1)+2]) = tmpP2(2);
					}
					(*RecArm2[0]) = Orb2->Arm(3, 2, t);
					(*RecArm2[1]) = Orb2->Arm(1, 3, t);
					(*RecArm2[2]) = Orb2->Arm(2, 1, t);
					(*RecArm2[3]) = Orb2->Arm(2, 3, t);
					(*RecArm2[4]) = Orb2->Arm(3, 1, t);
					(*RecArm2[5]) = Orb2->Arm(1, 2, t);
					//! * Write the output
					fOutPos2.RecordData();
					fOutArm2.RecordData();
				}
				
				
				//! ********************* Free memory of the second run
				
				delete Orb2;
				MT.Free(RecPos2, 9*sizeof(double*));
				MT.Free(RecArm2, 6*sizeof(double*));
				
			}
			
		}
		
		if (MT.Disp())
			Coutm << Endl << "End." << Endl;
	}
	
	
	catch(std::exception & e ) {
		std::cerr << "TestOrbits: error: " << e.what()<<Endl;
		std::cerr << "TestOrbits: abort!" << Endl;
		exit(1);
	}
	return(0);
};


/** \}*/

// end of