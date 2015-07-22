// $Id:  Exp $
/*
 *  LISACODE-TestGWs.cpp
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
#include "LISACODE-GWSpinBBHNR1.h"
#include "LISACODE-DataFileWrite.h"
#include "LISACODE-GWCosmicString.h"

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
    double tDur(NDat*dt);
    bool FFTCompute(false);
    
    //! *********** Help *************
    if(((argc>1)&&(strcmp(argv[1],"--help")==0))||((argc>1)&&(strcmp(argv[1],"-h")==0))){
      Coutm << " ----- HELP -----" << Endl;
      Coutm << Endl << "\tExecution : Compute hB+ and hBx of each source" << Endl;
      Coutm << "\t\t(./)LC2GWbary [Options] %sfOut %sfInXMLFile " << Endl;
      Coutm << Endl << "\tArguments :" << Endl;
      Coutm << "\t\t * %sfOut       : Output file name : t hBp1 hBc1 hBp2 hBc2 ... [default: TestGWs]. " << Endl;
      Coutm << "\t\t * %sfInXMLFile : XML input file [default: TestGWs.xml]. " << Endl;
      Coutm << Endl << "\tOptions :" << Endl;
      Coutm << "\t\t * -s %dseed : Seed for random gennerator. [default : current time]"  << Endl ;
      Coutm << "\t\t * -t %fdtmes : Time step. [default : 1 ]"  << Endl ;
      Coutm << "\t\t * -T %fdur   : Duration. [default : 1000*dtMes ]"  << Endl ;
      Coutm << "\t\t * -t0 %ft0 : Initial time. [default : 0 ]"  << Endl ;
      Coutm << "\t\t * -fft : Also compute fft of h+ and hx. [default : false ]"  << Endl ;
      Coutm << "\t\t * -dn \t\t: No screen display. [default: false]"  << Endl ;
      Coutm << "\t\t * -dl %sfile \t: Write standard output in a file. [default: no file]"  << Endl ;
      Coutm << "\t\t * -v \t\t: Verobse : display full details. [default: false]"  << Endl ;
      Coutm << Endl << " ----------------" << Endl;
      return 0;
      
    }
    
    //! *********** Version *************
    if(((argc>1)&&(strcmp(argv[1],"--version")==0))&&((argc>1)&&(strcmp(argv[1],"-v")==0))){
      Coutm << " ----- VERSION -----" << Endl;
      Coutm << " LC2GWbary : executable computing hB+ and hBx of each source - LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
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
	tDur = atof(argv[iarg+1]);
	nOptions +=2;
      }
      if((argc>1)&&(strcmp(argv[iarg],"-t0")==0)){
	t0 = atof(argv[iarg+1]);
	nOptions +=2;
      }
      if((argc>1)&&(strcmp(argv[iarg],"-fft")==0)){
	FFTCompute = true;
	nOptions++;
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
    NDat = MT.ifloor(tDur/dt);
    
    Coutm << "NData = " << NDat  << " :  from " <<  t0 << " s  to " <<  t0 + NDat*dt << " s  with step of " << dt  << " s." << Endl;  
    
    //! ***** Initialization of variable
    strcpy(fOut,"TestGWs.txt");
    strcpy(fInXML,"/Users/petiteau/Applications/src/LISACode/LISACode_2_0/LISACode/Main/Test/Config-GWs_Test.xml");
    
    
    if(argc-nOptions>1){
      strcpy(fOut, argv[1+nOptions]);
    }else{
      throw std::invalid_argument("ERROR : Need output filename (and xml input filename as second argument).");
    }
    
    if(argc-nOptions>2){
      strcpy(fInXML, argv[2+nOptions]);
    }else{
      throw std::invalid_argument("ERROR : Need xml input filename.");
    }
    
    
    
    //! ********************* Configuration
    
    //! ***** Output file
    double ** RecP;
    double ** RecC;
    char tmpN[10000];
    sprintf(tmpN, "%s.txt", fOut);
    Coutm << "Write output in " << tmpN << " ..." << Endl;
    LCDataFileWrite fOutTest(&MT, tmpN, ASCII);
    fOutTest.sett0(t0);
    fOutTest.setdt(dt);
    fOutTest.setNDatExpect(NDat);
    
    //! ***** Noise
    
    //! *** Declaration
    LCGW ** GWs(NULL);
    int NGWs(0);
    
    //! *** Configuration with XML file
    //! ** Read XML file
    ezxml_t tree, section;
    tree = ezxml_parse_file(fInXML);
    const char *type;
    for (section = ezxml_child(tree, "XSIL"); section; section = section->next) {
      type = ezxml_attr(section, "Type");
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
		  if((MT.wcmp(SourceType, "NRwave"))){
		    NGWs++;
		    GWs = (LCGW**) MT.ReAllocMemory(GWs, (NGWs-1)*sizeof(LCGW*), NGWs*sizeof(LCGW*));
		    GWs[NGWs-1] = new LCGWSpinBBHNR1(&MT);
		    GWs[NGWs-1]->config(gwdata);
		  }
		  
		  
		  if((MT.wcmp(SourceType, "CosmicString"))||(MT.wcmp(SourceType, "CosmicString"))){
		    NGWs++;
		    GWs = (LCGW**) MT.ReAllocMemory(GWs, (NGWs-1)*sizeof(LCGW*), NGWs*sizeof(LCGW*));
		    GWs[NGWs-1] = new LCGWCosmicString(&MT);
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
    
    //! ** Set the time informations
    for (int iN=0; iN<NGWs; iN++)
      GWs[iN]->setTimeInfo(t0, dt, NDat*dt, 0., 20.*dt, 2.*20.*dt);
    
    char NameRecP[128], NameRecC[128];
    RecP = (double**) MT.AllocMemory(NGWs*sizeof(double*));
    RecC = (double**) MT.AllocMemory(NGWs*sizeof(double*));
    for (int iN=0; iN<NGWs; iN++){
      sprintf(NameRecP, "%s_hBc", GWs[iN]->getName());
      sprintf(NameRecC, "%s_hBx", GWs[iN]->getName());
      RecP[iN] = fOutTest.AddRecord(NameRecP);
      RecC[iN] = fOutTest.AddRecord(NameRecC);
    }
    
    
    
    //! ********************* Initialization
    
    
    //! ** Initialization of noises
    for (int iN=0; iN<NGWs; iN++){
      GWs[iN]->init();
      Coutm << Endl;
      GWs[iN]->DispInfo("");
    }
    
    //! ** Initialization of output file
    fOutTest.init(NULL,0);
    
    
    //! ************ if FFT compute create an array to store tDat
    double ** tHpc;
    dcomplex ** fHpc;
    double tmpHp, tmpHc;
    int NfDat(MT.getNfFTreal(NDat));
    double df(1./(dt*NDat));
    tHpc = (double**) MT.AllocMemory(2*sizeof(double*));
    fHpc = (dcomplex**) MT.AllocMemory(2*sizeof(dcomplex*));
    if(FFTCompute){
      tHpc[0] = (double*) MT.AllocMemory(NDat*sizeof(double));
      tHpc[1] = (double*) MT.AllocMemory(NDat*sizeof(double));
      fHpc[0] = (dcomplex*) MT.AllocMemory(NfDat*sizeof(dcomplex));
      fHpc[1] = (dcomplex*) MT.AllocMemory(NfDat*sizeof(dcomplex));
    }
    
    
    //! ********************* Running
    std::cerr << "[........10........20........30........40........50........60........70........80........90.......100]" << Endl;
    std::cerr << "["; MT.o->flush();
    int iDisp = 0;
    for(int iT=0; iT<NDat; iT++){
      t = t0 + iT*dt;
      //Coutm << t << Endl;
      for (int iN=0; iN<NGWs; iN++){
	tmpHp = GWs[iN]->hBp(t);
	tmpHc = GWs[iN]->hBc(t);
	(*RecP[iN]) = tmpHp ;
	(*RecC[iN]) = tmpHc ;
	if(FFTCompute){
	  tHpc[0][iT] = tmpHp ;
	  tHpc[1][iT] = tmpHc ;
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
    
    
    //! ***** Compute Fourrier transform
    if(FFTCompute){
      
      //for(int i=0; i<NDat; i++)
      //	Coutm << i*dt << " " << tHpc[0][i] << " " << tHpc[1][i] << Endl;
      
      //for(int i=0; i<NfDat; i++)
      //	Coutm << i*df << " " <<  norm(fHpc[0][i]) << " " <<  arg(fHpc[1][i])<< " " <<  norm(fHpc[1][i]) << " " <<  arg(fHpc[1][i]) << Endl;
      
      //! *** Compute Fourrier transform
      MT.unsetBestFTFwd();
      MT.FTMakeFwdMulti(tHpc, fHpc, NDat, 2);
      
      //! *** Normalize Fourrier transform
      for(int iF=0; iF<NfDat; iF++){
	fHpc[0][iF] *= dt;
	fHpc[1][iF] *= dt;
      }
      
      //! *** The results are recorded for checking 
      
      sprintf(tmpN, "%s_fSig.txt", fOut);
      Coutm << "Write frequency output in " << tmpN << " ..." << Endl;
      std::ofstream ffout(tmpN);
      ffout.precision(12);
      ffout << "#f norm[h+] arg[h+] norm[hx] arg[hx]" << Endl;
      for(int i=0; i<NfDat; i++){
	ffout << i*df << " " <<  norm(fHpc[0][i]) << " " <<  arg(fHpc[1][i])<< " " <<  norm(fHpc[1][i]) << " " <<  arg(fHpc[1][i]) << Endl;
      }
      ffout.close();
      
      
      MT.Free(tHpc[0], NDat*sizeof(double));
      MT.Free(tHpc[1], NDat*sizeof(double));
      MT.Free(fHpc[0], NfDat*sizeof(dcomplex));
      MT.Free(fHpc[1], NfDat*sizeof(dcomplex));
      MT.Free(tHpc, 2*sizeof(double*));
      MT.Free(fHpc, 2*sizeof(dcomplex*));
    }
    MT.Free(tHpc,2*sizeof(double*));
    MT.Free(fHpc,2*sizeof(dcomplex*));
    
    for (int iN=0; iN<NGWs; iN++)
      delete GWs[iN];
    
    

    MT.Free(GWs, NGWs*sizeof(LCGW*));
    MT.Free(RecP, NGWs*sizeof(double*));
    MT.Free(RecC, NGWs*sizeof(double*));


    
    if (MT.Disp())
      Coutm << Endl << "End." << Endl;
  }
  
  
  catch(std::exception & e ) {
    std::cerr << "TestGWs: error: " << e.what()<<Endl;
    std::cerr << "TestGWs: abort!" << Endl;
    exit(1);
  }
  return(0);
};


/** \}*/

// end of
