// $Id:  Exp $
/*
 *  LISACODE-xml2ascii.c
 *
 *  Created by Antoine Petiteau on 30/08/13.
 *  Copyright 2013 APC - CNRS - University Paris Diderot (Paris,France). All rights reserved.
 *  Inspred by lisatools/lisaXML/C-examples/xml2ascii.c --- Copyright (c) 2006 Michele Vallisneri
 */


#include "LISACODE-Tools.h"
#include "LISACODE-DataFileRead.h"
#include "LISACODE-DataFileWrite.h"



/**\ingroup ToolBox 
 * \{
 */


/** \brief Convert xml data in ascii
 * \author A. Petiteau
 * \version 3.0
 * \date 30/08/2013
 *
 *
 */


void LoopConvert(LCTools * MT, ezxml_t secRoot){
    ezxml_t sec;
    for (sec = ezxml_child(secRoot, "XSIL"); sec; sec = sec->next) {
        if(MT->wcmp(ezxml_attr(sec, "Type"),"TimeSeries")){
            char Names[512];
            char ** Titles;
            int NTitles(0), NCharTitles(64), iTi(0);
            strcpy(Names,ezxml_attr(sec, "Name"));
            
            Cout << "Converting time serie named " << Names << Endl;
            MT->wextractcoma(Names, iTi, Titles, NTitles, NCharTitles);
            
            
            //! **** Read and load the file
            LCDataFileRead fIn(MT);
            int NDat, NRec;
            double t(0.) , t0(0.), dt(0.);
            fIn.config(sec);
            fIn.init();
            fIn.ControlDisplay();
            NDat = fIn.getNDat();
            NRec = fIn.getNRec();
            t0 = fIn.getx0();
            dt = fIn.getdx();
            
            //! **** Defined ascii file name
            char fNIn[1024], fNBase[1024], fNOut[1024];
            fIn.getFileName(fNIn);
            strncpy(fNBase,fNIn,strlen(fNIn)-4);
            sprintf(fNOut,"%s.txt",fNBase);
            
            //! **** Prepare output file
            LCDataFileWrite fOut(MT, fNOut, ASCII);
            double ** RecPtr;
            RecPtr = (double**)MT->AllocMemory(NRec*sizeof(double*));
            fOut.sett0(t0);
            fOut.setdt(dt);
            fOut.setNDatExpect(NDat);
            for (int iR=0; iR<NRec; iR++)
                RecPtr[iR] = fOut.AddRecord(Titles[MAX(0,NTitles-NRec+iR)]);
            fOut.init(NULL,0);
            
            
            //! ********************* Running
            for(int iT=0; iT<NDat; iT++){
                t = t0 + iT*dt;
                for (int iR=0; iR<NRec; iR++)
                    (*RecPtr[iR]) = fIn.gDataBin(iR,iT);
                //! * Write the output
                fOut.RecordData();
            }
            
            for(int iW=0; iW<NTitles; iW++)
                MT->Free(Titles[iW],64*sizeof(char));
            
            
        }
        LoopConvert(MT,sec);
    }
}




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
			Coutm << "\t\t(./)LC2zDL [Options] InputFile.xml  " << Endl;
			Coutm << Endl << "\tArguments :" << Endl;
			Coutm << "\t\t * InputFile.xml : XML input file from which one the data will be extracted." << Endl;
			Coutm << Endl << "\tOptions :" << Endl;
			Coutm << "\t\t * -v \t\t: Verobse : display full details. [default: false]"  << Endl ;
			Coutm << Endl << " ----------------" << Endl;
			return 0;
			
		}
		
		// *********** Version *************
		if(((argc>1)&&(strcmp(argv[1],"--version")==0))&&((argc>1)&&(strcmp(argv[1],"-v")==0))){
			Coutm << " ----- VERSION -----" << Endl;
			Coutm << " LC2-xml2ascii : Convert xml data in ascii - LISACode package - version " << LC::LCVersion << " at " << LC::DateOfLastUpdate << Endl;
			Coutm << " ----------------" << Endl;
			return 0;
		}
		
		
		//! ***** Options 
		for(int iarg=1; iarg<argc; iarg++){
			if((argc>1)&&(strcmp(argv[iarg],"-v")==0)){
				MT.setDispDetails();
				nOptions++;
			}
		} 
		MT.setRandSeed(SeedRand);
		
		
		//! ***** Some declarations
        char fNIn[1024];
        
		
		if(argc-nOptions<1){
			throw std::invalid_argument("ERROR : we need to know the input file");
		}
        strcpy(fNIn,argv[1+nOptions]);
		
        
		ezxml_t tree;
		tree = ezxml_parse_file(fNIn);

        LoopConvert(&MT,tree);

        
        ezxml_free(tree);
        
		
	}
	
	
	
	
	
	catch(std::exception & e ) {
		std::cerr << "xml2ascii: error: " << e.what()<<Endl;
		std::cerr << "xml2ascii: abort!" << Endl;
		exit(1);
	}
	return(0);
};


/** \}*/

/*
int main(int argc,char **argv) {
    TimeSeries *timeseries;
    LISASource *lisasource;
	
    FILE *outfile;
	
    int i,j;
    int nOpt;
    int Osrcf;
	int iTimeS;
	int ReadContinue;
	
	char outfilename[512];
	char tmpstring[512];
	
	
    if((argc < 2)||(strcmp(argv[1],"--help")==0)||(strcmp(argv[1],"-h")==0)) {
		printf("Version 2.1 (modified by Antoine Petiteau\n)");
		printf("Usage: %s [options] source.xml\n",argv[0]);
		printf("     options : -s : Read source file (h+ hx)\n");
		exit(0);
    }
	
    nOpt = 0;
    Osrcf = 0;
	
    for(i=1; i<argc; i++ ){
		if(strcmp(argv[i],"-s")==0){
			Osrcf = 1;
			fprintf(stderr, " Option : Read source file data\n");
			nOpt++;
		}
    };
	
	iTimeS = 0;
	ReadContinue = 1;
	do{
		if(Osrcf){ // Read source datafile
			if(iTimeS < 1){
				lisasource = getLISASources(argv[1+nOpt]);
				if(lisasource == 0)
					ReadContinue = 0;
				else
					timeseries = lisasource->TimeSeries;
			}else{
				ReadContinue = 0;
			}
		}else{ // Read TDI datafile
			timeseries = getmultipleTDIdata(argv[1+nOpt],iTimeS);
			if(timeseries == 0)
				ReadContinue = 0;
		}
		
		if(ReadContinue){
			strncpy(tmpstring, timeseries->FileName, strlen(timeseries->FileName)-4);
			//printf("%s\n", tmpstring);
			tmpstring[strlen(timeseries->FileName)-4] = '\0';
			//printf("%s\n", tmpstring);
			sprintf(outfilename, "%s.txt", tmpstring); 
			outfile = fopen(outfilename,"w");
			
			fprintf(stderr,"-> Read TimeSeries '%s' from file '%s':\n",timeseries->Name,timeseries->FileName);
			fprintf(stderr,"-> %d x %d values (%s)\n",timeseries->Length,timeseries->Records,timeseries->Name);
			
			
			for(i=0;i<timeseries->Records;i++) {
				fprintf(stderr,"-> column %d: %s\n",i,timeseries->Data[i]->Name);
			}
			
			fprintf(stderr,"-> Write data in '%s' ...\n",outfilename);
			
			if(argc > 2) {
				fprintf(outfile,"#\n#----------\n");
			}
			
			for(j=0;j<timeseries->Length;j++) {
				if(Osrcf){
					fprintf(outfile,"%.16g\t",timeseries->TimeOffset + j*timeseries->Cadence);
				}
				for(i=0;i<timeseries->Records;i++) {
					fprintf(outfile,"%.16g",timeseries->Data[i]->data[j]);
					
					if(i+1 < timeseries->Records) {
						fprintf(outfile,"\t");
					} else {
						fprintf(outfile,"\n");
					}
				}
			}
			
			freeTimeSeries(timeseries);
			
			fclose(outfile);
			
			iTimeS++;
		}
		
	}while( ReadContinue );
	
	if(iTimeS == 0){
		if(Osrcf)
			fprintf(stderr,"No source data found. Maybe you are looking for TDI data ? In this case, remove -s option.\n");
		else
			fprintf(stderr,"No TDI data found. Maybe you are looking for source data ? In this case, use -s option.\n");
	}
	return 0;
}
 */

// end of LISACODE-xml2ascii.c