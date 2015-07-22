/*
 *  LISACODE-DataFileWrite.cpp
 *  LC20_xcode
 *
 *  Created by Antoine Petiteau on 15/02/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-DataFileWrite.h"


// *****************
// *  Constructor  *
// *****************

/**
 *
 */
LCDataFileWrite::LCDataFileWrite()
{
	MT = new LCTools;
	MT->LocTools = true;
	initNULL(false);
}


/**
 *
 */
LCDataFileWrite::LCDataFileWrite(LCTools * MT_n, const char * fNOut_n, TypeFile typeOutFile_n)
{
	MT = MT_n;
	initNULL(false);
	
	strcpy(fNOut, fNOut_n);
	typeOutFile = typeOutFile_n;

}


/**
 *
 */
LCDataFileWrite::~LCDataFileWrite()
{
	initNULL(true);
}


// ********************
// *  Access methods  *
// ********************





// ***************************************
// * Linking and initialization methods  *
// ***************************************


/**
 *
 */
void LCDataFileWrite::initNULL(bool CleanMem)
{
	/** *** Clean memory if required */ 
	if(CleanMem){
		if(iDatWrite != -1){
			fOut.close();
			fOut.clear();
			if(NDatExpect != iDatWrite)
				std::cerr << "WARNING in LCDataFileWrite::initNULL : For the file " << fNOut << " , the number of data expected (" << NDatExpect << ") does not correspond to the number of data written (" << iDatWrite << ") !" << Endl;
		}
		
		if(RecData != NULL){
			for (int i=0; i<NRec; i++)
				if(RecData[i] != NULL)
					MT->Free(RecData[i], sizeof(double));
			MT->Free(RecData, NRec*sizeof(double*));
		}
		
		if(RecName != NULL){
			for (int i=0; i<NRec; i++)
				MT->Free(RecName[i], 128*sizeof(char));
			MT->Free(RecName, NRec*sizeof(char*));
		}
	}
	
	/** *** Set all values at 0 or NULL */
	
	NRec = 0;
	RecData = NULL;
	RecName = NULL;
	strcpy(fNOut,"");
	typeOutFile = UNDEFINED;
	
	WriteRef = true;
	t0 = 0.;
	dt = 1.;
	NDatExpect = 0;
	
	iDatWrite = -1;
}




/**
 *
 */
double * LCDataFileWrite::AddRecord(char * RecNameAdd)
{
	NRec++;
	
	//! *** Add memory for the new data 
	RecData = (double**) MT->ReAllocMemory(RecData, (NRec-1)*sizeof(double*), NRec*sizeof(double*));
	RecData[NRec-1] = (double*) MT->AllocMemory(sizeof(double));
	(*RecData[NRec-1]) = 0.;
	
	// *** Add memory for the name of new data
	RecName = (char**) MT->ReAllocMemory(RecName, (NRec-1)*sizeof(char*), NRec*sizeof(char*));
	RecName[NRec-1] = (char*) MT->AllocMemory(128*sizeof(char));
	strcpy(RecName[NRec-1], RecNameAdd);
	
	return(RecData[NRec-1]);
}


/**
 *
 */
void LCDataFileWrite::init(std::ostream * GlobalHeader, int indent)
{
	//! *** If the file is already openned, close it
	if(iDatWrite != -1){
		fOut.close();
		fOut.clear();
		if(NDatExpect != iDatWrite)
			std::cerr << "WARNING in LCDataFileWrite::initNULL : For the file " << fNOut << " , the number of data expected (" << NDatExpect << ") does not correspond to the number of data written (" << iDatWrite << ") !" << Endl;
	}
	
	//! *** Open the file
	if(typeOutFile == ASCII )
		fOut.open(fNOut);
	else
		fOut.open(fNOut, std::ios::binary);
	
	iDatWrite = 0;
	
	//! *** Write the header 
	//! ** For ASCII file
	if(typeOutFile == ASCII){
		//! ** Just the titles
		fOut << "#";
		if(WriteRef)
			fOut << "Ref ";
		for(int iR=0; iR<NRec; iR++){
			if(iR!=0)
				fOut << " ";
			fOut << RecName[iR];
		}
		fOut << Endl;
		fOut.precision(12);
	}
	
	//! ** For binary file which requires an header
	if(typeOutFile == BINARY){
		fOut << "#TITLE Ref";
		for(int iR=0; iR<NRec; iR++)
			fOut << " " << RecName[iR];
		fOut << Endl; 
		fOut << "#RECORD " << NRec << " " << NDatExpect <<  Endl;
		fOut << "#TIME " << t0 << " " << dt << " " << t0+NDatExpect*dt << Endl;
		fOut.precision(10);
	}
	if(typeOutFile == XML){
        if(GlobalHeader != NULL)
            WriteXMLHeader(GlobalHeader, indent);
    }
	
}



// *********************
// *  Running methods  *
// *********************

void LCDataFileWrite::RecordData()
{
	//! *** If reference required, compute the value and pass to the next function 
	if(WriteRef)
		RecordData(t0+iDatWrite*dt);
	else
		RecordData(0.);
}


void LCDataFileWrite::RecordData(double vRef)
{
	if(iDatWrite==-1)
		throw std::invalid_argument("ERROR in LCDataFileWrite::RecordData : The file is not open ! You maybe forget the init() !");
	
	//! **** Write in ASCII 
	if(typeOutFile==ASCII){
		if(WriteRef)
			fOut << vRef << " ";
		for(int iR=0; iR<NRec; iR++){
			if(iR!=0)
				fOut << " ";
			fOut << (*RecData[iR]);
		}
		fOut << Endl;
	}
	
	//! **** Write in binary
	if((typeOutFile==BINARY)||(typeOutFile==XML)){
		if(WriteRef)
			fOut.write((char*) &vRef, sizeof(double));
		for(int iR=0; iR<NRec; iR++)
			fOut.write((char*) RecData[iR], sizeof(double));
	}
	
	iDatWrite++;
	
}


void LCDataFileWrite::WriteXMLHeader(std::ostream * fHeader, int NIndent)
{
	char BTab[2048];
	char ObsDescr[2048];
	int checksum(0); // TODO : Give the good checksum
	
	if(fHeader == NULL)
		throw std::invalid_argument("LCDataFileWrite::WriteXMLHeader : The output file is NULL !");
	
	strcpy(ObsDescr,"");
	if(WriteRef)
		strcat(ObsDescr, "t");
	for(int iR=0; iR<NRec; iR++){
		if((WriteRef)||(iR!=0))
			strcat(ObsDescr, ",");
		strcat(ObsDescr, RecName[iR]);
	}
	
	
	//! *** Prepare the indentation
	strcpy(BTab,"");
	for(int i=0; i<NIndent; i++)
		strcat(BTab,"\t");
	
	//! *** 
    //! * Note : all data are type of TDIObservable to simplify and match lisatools format, but later it should be better to change that !
    fHeader->precision(15);
	(*fHeader) << BTab << "<XSIL Name=\"" << ObsDescr << "\" Type=\"TDIObservable\">" << Endl;
    (*fHeader) << BTab <<LC::xmlind<<"<XSIL Name=\"" << ObsDescr << "\" Type=\"TimeSeries\">" << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<"<Param Name=\"TimeOffset\" Unit=\"Second\">" << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<"" << t0 << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<"</Param>" << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<"<Param Name=\"Cadence\" Unit=\"Second\">" << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<"" << dt << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<"</Param>" << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<"<Param Name=\"Duration\" Unit=\"Second\">" << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<"" << dt*NDatExpect << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<"</Param>" << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<"<Array Name=\"" << ObsDescr << "\" Type=\"double\" Unit=\"\">" << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<"<Dim Name=\"Length\">" << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<LC::xmlind<<"" << NDatExpect << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<"</Dim>" << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<"<Dim Name=\"Records\">" << Endl;
	if(WriteRef)
		(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<LC::xmlind<<"" << NRec+1 << Endl;
	else
		(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<LC::xmlind<<"" << NRec << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<"</Dim>" << Endl;
	
	if(typeOutFile==ASCII)
		(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<"<Stream Checksum=\"" <<  checksum << "\" Type=\"Remote\" Encoding=\"ASCII\">" << Endl;
	
	if((typeOutFile==BINARY)||(typeOutFile==XML)){
		if(MT->IsLocalBigEndian())
			(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<"<Stream Checksum=\"" <<  checksum << "\" Type=\"Remote\" Encoding=\"Binary,BigEndian\">" << Endl;
		else
			(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<"<Stream Type=\"Remote\" Checksum=\"" <<  checksum << "\" Encoding=\"Binary,LittleEndian\">" << Endl;
	}
	
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<LC::xmlind<<"" << fNOut << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<LC::xmlind<<"</Stream>" << Endl;
	(*fHeader) << BTab <<LC::xmlind<<LC::xmlind<<"</Array>" << Endl;
	(*fHeader) << BTab <<LC::xmlind<<"</XSIL>" << Endl;
    (*fHeader) << BTab << "</XSIL>" << Endl;
	
}



void LCDataFileWrite::DispInfo(char * BTab)
{
	if(MT->Disp()){
		Cout << BTab <<LC::xmlind<<"- File = " << fNOut << Endl;
		Cout << BTab <<LC::xmlind<<"- Records =";
		for(int iR=0; iR<NRec; iR++)
			Cout << " " << RecName[iR];
		Cout << Endl;
	}
}



// end of LISACODE-DataFileWrite.cpp

