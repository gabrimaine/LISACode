/*
 *  LISACODE-DataFileRead.cpp
 *  LC20_xcode
 *
 *  Created by Antoine Petiteau on 15/02/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-DataFileRead.h"


// *****************
// *  Constructor  *
// *****************


LCDataFileRead::LCDataFileRead()
{
	MT = new LCTools;
	initNULL(false);
}


LCDataFileRead::LCDataFileRead(LCTools * MT_n)
{
	MT = MT_n;
	initNULL(false);
}


LCDataFileRead::LCDataFileRead(LCTools * MT_n, char * NewFileName)
{
	MT = MT_n;
	initNULL(false);
	setFile(NewFileName);
}


LCDataFileRead::~LCDataFileRead()
{
	initNULL(true);
}



void LCDataFileRead::initNULL(bool CleanMem)
{
	if(CleanMem){
		if(Data!=NULL){
			for(int iR=0; iR<NRec; iR++)
				if(Data[iR]!= NULL)
					delete Data[iR];
			MT->Free(Data, NRec*sizeof(LCSerie2*));
		}
		for(int iW=0; iW<Header.size(); iW++){
			if(Header[iW] != NULL)
				MT->Free(Header[iW], 16384*sizeof(char));
		}
	}
	
	FileType = UNDEFINED;
	strcpy(fNIn, "None");
	
	NRec = 0;
	NData = 0;
	dx = 1.;
	x0 = 0.;
	xend = 0.;
	SkipFirstRec = false;
	Header.resize(0);
	PosStartBin = 0;
	FileEncoding = MT->testbyteorder();
	
	Data = NULL;
	
}


// ********************
// *  Access methods  *
// ********************

void LCDataFileRead::getHeader(char ** & Words, int & NWords, int NCharWord)
{
	NWords = Header.size();
	Words = (char**) MT->AllocMemory(NWords*sizeof(char*));
	for(int i=0; i<NWords; i++){
		Words[i] = (char*) MT->AllocMemory(NCharWord*sizeof(char));
		strncpy(Words[i], Header[i], NCharWord);
		Cout << Header[i] << " -> " << Words[i] << Endl;
	}
}

void LCDataFileRead::setFile(char * NewFileName)
{
	if(!MT->wcmp(NewFileName,fNIn)){
		//! *** If it's a new file
		//! ** Clean the previous data 
		initNULL(true);
		strcpy(fNIn, NewFileName);
		FileType = UNDEFINED;
	}
}




// ***************************************
// * Linking and initialization methods  *
// ***************************************

void LCDataFileRead::config(ezxml_t seriesxmlbloc)
{
	char DataName[2048];
	ezxml_t param, array, dim, stream;
	double tDur;
	
	if(seriesxmlbloc == NULL)
		throw std::invalid_argument("ERROR in LCDataFileRead::config : The xml bloc is NULL !");
	
	if(!MT->wcmp(ezxml_attr(seriesxmlbloc,"Type"),"TimeSeries"))
		throw std::invalid_argument(" ERROR in LCDataFileRead::config : The type is not TimeSeries ");
	
	
	strcpy(DataName,ezxml_attr(seriesxmlbloc,"Name"));
	
	for(param = ezxml_child(seriesxmlbloc,"Param"); param; param = param->next) {
		if(strcmp(ezxml_attr(param,"Name"),"TimeOffset")==0)
			x0 = MT->gXMLTime(param);
		if(strcmp(ezxml_attr(param,"Name"),"Cadence")==0)
			dx = MT->gXMLTime(param);
		if(strcmp(ezxml_attr(param,"Name"),"Duration")==0)
			tDur = MT->gXMLTime(param);
	}
	array = ezxml_child(seriesxmlbloc,"Array");
	for(dim = ezxml_child(array,"Dim"); dim; dim = dim->next) {
		if(strcmp(ezxml_attr(dim,"Name"),"Length")==0)
			NData = atoi(ezxml_txt(dim));
		if(strcmp(ezxml_attr(dim,"Name"),"Records")==0)
			NRec = atoi(ezxml_txt(dim));
	}
	stream = ezxml_child(array,"Stream");
	
	FileEncoding = MT->ReadEncodingConst(ezxml_attr(stream,"Encoding"));
	/*! *** If it's an "endian" file, it's a binary with or without header inside. 
	 *  So at the moment we suppose without (XML) but it can still be change during the file analysis if the header is detected. */
	if(FileEncoding != NOENDIAN)
		FileType = XML;
	
	
	char * TmpFileName(NULL);
	MT->stripcopy(ezxml_txt(stream), TmpFileName);
	strcpy(fNIn, TmpFileName);
	MT->Free(TmpFileName, (strlen(TmpFileName)+1) * sizeof(char) );
	
	
	//! *** Copy the data name in the header
	Header.push_back(NULL);
	Header[Header.size()-1] = (char*) MT->AllocMemory(16384*sizeof(char));
	strcpy(Header[Header.size()-1],DataName);
	
	//! *** Detect if the first record is the time
	if(Header[0]!=NULL){
		if(Header[0][0]=='t')
			SkipFirstRec = true;
	}
	if(SkipFirstRec)
		NRec--;
	if(!MT->deq(tDur, dx * NData))
		throw std::invalid_argument("LCDataFileRead::config : Duration is not coherent with the other parameters !");
	xend = x0 + dx * NData;
	   
}



void LCDataFileRead::init()
{
	//! *** Analysis the file
	FileAnalysis();
	ControlDisplay();
	
	//! *** Allocate memory
	Data = (LCSerie2**)MT->AllocMemory(sizeof(LCSerie2*) * NRec);
	for(int iR=0; iR< NRec; iR++){
		Data[iR] = new LCSerie2(MT, x0, dx, NData);
		Data[iR]->allocAll();
	}
	
	//! *** Read data
	ReadDataSeries(NULL);
	
}


void LCDataFileRead::init(double ** & DataExt, int & NRecExt, int & NDataExt)
{
	//! *** Analysis the file
	FileAnalysis();
	ControlDisplay();
	
	NRecExt = NRec;
	NDataExt = NData;
	
	//! *** Allocate memory
	DataExt = (double**)MT->AllocMemory(sizeof(double*) * NRecExt);
	for(int iR=0; iR< NRec; iR++){
		DataExt[iR] = (double*)MT->AllocMemory(sizeof(double) * NDataExt);
	}
	
	//! *** Read data
	ReadDataSeries(DataExt);
	
}



// *******************
// *  Other methods  *
// *******************

void LCDataFileRead::FileAnalysis()
{
	char * ptr;
	std::ifstream fIn;
	char Buf[2048];
	bool NeedAnalysis(true);
	
	//! ***** Control if file exist 
	fIn.open(fNIn);
	if(fIn == NULL){
		Cout << " ERROR: LCDataFileRead::FileAnalysis : Can not open the file " << fNIn << " ." << Endl;
		throw std::invalid_argument(" ERROR: LCDataFileRead::FileAnalysis : Can not open the file.");
	}
	fIn.close();
	fIn.clear();
	
	//! ***** Test XML : Read extention of file
	ptr = fNIn + strlen(fNIn);
	if(strncmp(ptr-4,".xml",4)==0) {
		//! *** XML file
		FileType = XML;
		
	}else{
		
		//! *** Test BINARY : Read the 3 first line of file (for binary file with header)
		fIn.open(fNIn);
		if(fIn.peek() == '#'){
			fIn >> Buf;
			if(MT->wcmp(Buf,"#TITLE")){
				fIn.ignore(16384,'\n');
				fIn >> Buf;
				if(MT->wcmp(Buf,"#RECORD")){
					fIn.ignore(16384,'\n');
					fIn >> Buf;
					if(MT->wcmp(Buf,"#TIME"))
						FileType = BINARY;
				}
			}
		}
		fIn.close();
		fIn.clear();
		
		/*! *** If at this point the file type is XML, it means that the xml bloc is included in another file and have already been read.
		 *  So we don't have to analyse the file again.
		 */
		if(FileType == XML)
			NeedAnalysis = false;
		
		//! ***** Test ASCII : Read the first data for testing if it's a good value
		if(FileType == UNDEFINED){
			fIn.open(fNIn);
			//! *** Junk header
			while(fIn.peek() == '#'){
				fIn.ignore(16384,'\n');
			};
			fIn >> Buf;
			if(MT->isValue(Buf))
				FileType = ASCII;
			fIn.close();
			fIn.clear();
		}
	}
	
	if(NeedAnalysis){
		switch (FileType) {
			case ASCII :
				FileAnalysisASCII() ;
				break;
			case BINARY :
				FileAnalysisBinary() ;
				break;
			case XML :
				FileAnalysisXML() ;
				break;
			default :
				throw std::invalid_argument("ERROR: LCDataFileRead::FileAnalysis : The type of the file can't be identified.") ;
		}
	}
}


/*! Analyse ASCII file */
void LCDataFileRead::FileAnalysisASCII()
{
	std::ifstream fIn;
	double tr0(0.0), tr1(0.0), trn(0.0);
	char junk[512];
	
	if(MT->Disp())
		Cout << "  Analysing of ASCII file " << fNIn << " ...";  MT->o->flush();
	
	//! *** Count the number of columns 
	fIn.open(fNIn);
	if (!fIn)
		throw std::invalid_argument("ERROR: LCDataFileRead::FileAnalysisASCII : Problem in openning data file.");
	while(fIn.peek() == '#'){
		Header.push_back(NULL);
		Header[Header.size()-1] = (char*) MT->AllocMemory(16384*sizeof(char));
		fIn.getline(Header[Header.size()-1],16384,'\n');
	};
	NRec = 0;
	while((fIn.peek() != '\n')&&(!fIn.eof())){
		fIn >> junk;
		//Cout << "- Record " << NRec << " : " << junk << " :: " << fIn.peek()  << Endl;
		// ** Considering the case " \n"
		if((fIn.peek() == ' ')||fIn.peek() == '\t'){
			int ipos(fIn.tellg());
			ipos++;
			fIn.seekg(ipos);
		}
		NRec++;
	}
	fIn.close();
	fIn.clear();
	
	//! *** Count the number of rows 
	fIn.open(fNIn);
	while(fIn.peek() == '#'){
		fIn.ignore(16384,'\n');
	};
	NData = 0;
	while(!fIn.eof()){
		//! ** Read reference
		switch (NData){
			case 0 :
				fIn >> tr0;
				//Cout << "x0 = " << x0 << Endl;
				//fIn >> junk;
				//Cout << "junk = " << junk << Endl;
				break;
			case 1 :
				fIn >> tr1;
				break;
			default :
				fIn >> trn;
		}
		
		/*! Read end of row */
		fIn.ignore(16384,'\n');
		if(!fIn.eof()){	
			NData++;
		}
	};
	fIn.close();
	fIn.clear();
	
	//! *** The first column is the reference so number of records = number of columns - 1
	SkipFirstRec = true;
	NRec--;
	
	x0 = tr0;
	dx = fabs(tr0-tr1);
	xend = trn;
	
	if(MT->Disp())
		Cout << "  --> OK." << Endl;
}


void LCDataFileRead::FileAnalysisBinary()
{
	if(MT->Disp()){
		Cout << "  Analysing of Binary file " << fNIn << " ..."; MT->o->flush();
	}
	
	/*! *** Read the title */
	std::ifstream BinFileTitle(fNIn);
	char tmpTitle[16384];
	std::vector<char*> Words(0);
	bool AddHeadLine(false);
	
	
	BinFileTitle.getline(tmpTitle,16384,'\n');
	BinFileTitle.close();
	BinFileTitle.clear();
	
	MT->wextract(tmpTitle, Words);
	if(Words.size()<2){
		throw std::invalid_argument("ERROR: LCDataFileRead::FileAnalysisBinary :  There no title (first line of the binary file).");
	}else{
		Header.push_back(NULL);
		Header[Header.size()-1] = (char*) MT->AllocMemory(16384*sizeof(char));
		strcpy(Header[Header.size()-1],Words[1]);
	}
	
	for(int iW=2; iW<Words.size(); iW++){
		Cout << "\t" << Words[iW] << Endl ;
		if(MT->wcmp(Words[iW],"]")){
			AddHeadLine = false;
		}
		if(AddHeadLine){
			strcat(Header[Header.size()-1],Words[iW]);
			strcat(Header[Header.size()-1]," ");
		}
		if(MT->wcmp(Words[iW],"[")){
			Header.push_back(NULL);
			Header[Header.size()-1] = (char*) MT->AllocMemory(16384*sizeof(char));
			strcpy(Header[Header.size()-1],"");
			AddHeadLine = true;
		}
	}
	//setTypeData(Words[1]);
	
	/*! *** Read the binary file */
	FILE * fIn;
	fIn = fopen(fNIn,"rb");
	if(fIn == NULL)
		throw std::invalid_argument("ERROR: LCDataFileRead::FileAnalysisBinary :  Problem in opening the file.");
	
	fseek(fIn , strlen(tmpTitle)+1, SEEK_SET);
	fscanf(fIn, "#RECORD %d %ld \n", &NRec, &NData);	
	fscanf(fIn, "#TIME %lf %lf %lf \n", &x0, &dx, &xend);
	Cout << "tend : read = " << xend << "  / computed = " << x0+dx*NData << Endl;
	
	PosStartBin = ftell(fIn);
	
	fclose(fIn);
	
		
	//! *** Detect if the first record is the time
	if(Header[0]!=NULL){
		if(Header[0][0]=='t')
			SkipFirstRec = true;
	}
	if(SkipFirstRec)
		NRec--;
	
	
	for(int iW=0; iW<Words.size(); iW++)
		MT->Free(Words[iW],256*sizeof(char));
	
	if(MT->Disp())
		Cout << "  --> OK." << Endl;
}



void LCDataFileRead::FileAnalysisXML()
{
	ezxml_t tree, section, sub1, sub2, sub3;
	bool FileNotFound(true);
	
	if(MT->Disp())
		Cout << "  Analysing of XML file " << fNIn << " ..." ;  MT->o->flush();
	
	//! **** Go through the xml file and configure using the FIRST XSIL bloc of type "TimeSeries"  
	tree = ezxml_parse_file(fNIn);
	
	//! *** Go through the section level (base level)
	for(section = ezxml_child(tree, "XSIL"); section; section = section->next) {
		if(FileNotFound){
			if(MT->wcmp(ezxml_attr(section, "Type"),"TimeSeries")){
				config(section);
				FileNotFound = false;
			}else{

				//! *** Go through the sub-level 1
				for(sub1 = ezxml_child(section, "XSIL"); sub1; sub1 = sub1->next) {
					if(FileNotFound){
						if(MT->wcmp(ezxml_attr(sub1, "Type"),"TimeSeries")){
							config(sub1);
							FileNotFound = false;
						}else{
							
							//! *** Go through the sub-level 2
							for(sub2 = ezxml_child(section, "XSIL"); sub2; sub2 = sub2->next) {
								if(FileNotFound){
									if(MT->wcmp(ezxml_attr(sub2, "Type"),"TimeSeries")){
										config(sub2);
										FileNotFound = false;
									}else{
										
										//! *** Go through the sub-level 3
										for(sub3 = ezxml_child(section, "XSIL"); sub3; sub3 = sub3->next) {
											if(FileNotFound){
												if(MT->wcmp(ezxml_attr(sub2, "Type"),"TimeSeries")){
													config(sub3);
													FileNotFound = false;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
    ezxml_free(tree);									

	if(MT->Disp())
		Cout << "  --> OK." << Endl;
}



void LCDataFileRead::ReadDataSeries(double ** ExtData)
{
	switch (FileType) {
		case ASCII :
			ReadDataSeriesASCII(ExtData);
			break;
		case BINARY :
			ReadDataSeriesBinary(ExtData);
			break;
		case XML :
			ReadDataSeriesBinary(ExtData);
			break;
		default :
			throw std::invalid_argument("ERROR in LCDataFileRead::ReadDataSeries : The type of the file can't be identified.") ;
	}
}



void LCDataFileRead::ReadDataSeriesASCII(double ** ExtData)
{
	std::ifstream fIn;
	char junk[512];
	int NCol, iRec;
	double tmpV;
	
	if(MT->Disp()){
		Cout << "  Loading of ASCII file " << fNIn << " ..." ;
		MT->o->flush();
	}
	
	//! *** Compute the number of columns
	if(SkipFirstRec)
		NCol = NRec+1;
	else 
		NCol = NRec;
	
	
	//! *** Open file and read header
	fIn.open(fNIn);
	if (!fIn)
		throw std::invalid_argument("DataFileManager::LCDataFileRead : Problem in openning data file.");
	while(fIn.peek() == '#'){
		fIn.ignore(16384,'\n');
	};
	
	
	//! *** Read and store data
	
	for(int nl=0; nl<NData; nl++){
		iRec = 0;
		for(int col=0; col<NCol; col++){
			if((col!=0)||((col==0)&&(!SkipFirstRec))){
				fIn >> tmpV;
				if(ExtData == NULL)
					Data[iRec]->setBinValue(nl, tmpV);
				else
					ExtData[iRec][nl] = tmpV;
				iRec++;
			}else{
				fIn >> junk;
			}
		}
	}
	
	fIn.close();
	fIn.clear();
	if(MT->Disp())
		Cout << "  --> OK." << Endl;
}



/*! Read data in Binary file and store only records specified by the vector of index iReadRec */
void LCDataFileRead::ReadDataSeriesBinary(double ** ExtData)
{
	FILE * fIn;
	int NCol;
	double * TmpData;
	
	if(MT->Disp()){
		Cout << "  Loading of Binary file " << fNIn << " ..." ;
		MT->o->flush();
	}
	
	//! *** Compute the number of columns
	if(SkipFirstRec)
		NCol = NRec+1;
	else 
		NCol = NRec;
	
	// ** Allocate memory for temporary buffer of data
	TmpData = (double*)MT->AllocMemory(sizeof(double) * NCol * NData);
	
	//! *** Open file
	fIn = fopen(fNIn,"rb");
	if (!fIn)
		throw std::invalid_argument("DataFileManager::ReadDataSeries : Problem in openning data file.");
	
	//! *** Read header if there is one
	if(FileType == BINARY){
		fseek(fIn , PosStartBin, SEEK_SET);
	}
	
	//! *** Read data
	fread(TmpData, sizeof(double), NCol * NData, fIn);
	fclose(fIn);
	
	//for(int i=0; i<10; i++)
	//	Cout << i << " " << TmpData[i] << Endl;
	
	
	//! Do the encoding switch if necessary
	if( (FileEncoding !=  MT->testbyteorder()) && (FileEncoding != NOENDIAN) ){
		Cout << " --> Converting endianness of binary data! ";
		for(int i=0; i<NCol*NData; i++)
			MT->convertendianness(&TmpData[i]);
	}
	
	//! *** Store data
	int pS;
	int pR(0);
	if(SkipFirstRec)
		pR = 1;
	for(int nl=0; nl<NData; nl++){
		pS = NCol*nl + pR;
		for(int iR=0; iR<NRec; iR++){
			if(ExtData == NULL)
				Data[iR]->setBinValue(nl, TmpData[ pS + iR ]);
			else
				ExtData[iR][nl] = TmpData[ pS + iR ];
		}
	}
	MT->Free(TmpData,sizeof(double) * NCol * NData);
	
	if(MT->Disp())
		Cout << "  --> OK." << Endl;
}



void LCDataFileRead::ControlDisplay()
{
	if(MT->Disp()){	
		Cout << "  File  " << fNIn << " : ";
		MT->DispFileType(FileType, MT->o);
		Cout << " (";
		MT->DispEncoding(FileEncoding, &Cout );
		Cout << ") : " << Endl;
		Cout << "    --> Data : "<< NRec << " records x " << NData << " rows" << Endl;
		Cout << "    --> Time : from "<< x0 << "s to "<< xend << "s  with step of " << dx << "s" <<  Endl;
		//Cout << Endl;
	}
}



// *************************
// *  Access data methods  *
// *************************


double LCDataFileRead::gData(int iS, double t, INTERP InterpType, double InterpUtilValue)
{
	//return( Data[iS]->gData(xend + x0 - dx - t, InterpType, InterpUtilValue) );
	return( Data[iS]->gData(t, InterpType, InterpUtilValue) );
}

double LCDataFileRead::gDataBin(int iS, int iT)
{
	return( Data[iS]->getBinValue(iT) );
}


double LCDataFileRead::DerivRawCur(int iS, int NDbin)
{
	return( Data[iS]->DerivRawCur(NDbin) );
}


double LCDataFileRead::DerivBackOrder2Cur(int iS)
{
	return( Data[iS]->DerivBackOrder2Cur() );
}


double LCDataFileRead::DerivRawSpe(int iS, double x, double Dx, INTERP InterpType, double InterpUtilValue)
{
	return( Data[iS]->DerivRawSpe(x, Dx, InterpType, InterpUtilValue) );
}


double LCDataFileRead::DerivBackOrder2Spe(int iS, double x, double Dx, INTERP InterpType, double InterpUtilValue)
{
	return( Data[iS]->DerivBackOrder2Spe(x, Dx, InterpType, InterpUtilValue) );
}



// end of LISACODE-DataFileRead.cpp