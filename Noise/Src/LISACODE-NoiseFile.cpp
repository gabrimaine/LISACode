/*
 *  LISACODE-NoiseFile.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 14/03/13.
 *  Copyright 2013 Laboratory AstroParicles and Cosmology - University Paris-Diderot (Paris, France). All rights reserved.
 *
 */


#include "LISACODE-NoiseFile.h"


// *********************
// ***  Constructor  ***
// *********************

/**
 *
 */
LCNoiseFile::LCNoiseFile()
: LCNoise()
{
	initNULL(false);
}


LCNoiseFile::LCNoiseFile(LCTools * MT_n)
: LCNoise(MT_n)
{
	initNULL(false);
}


LCNoiseFile::~LCNoiseFile()
{
	initNULL(true);
}



// ********************************************
// ***        Required methods              ***
// ********************************************

void LCNoiseFile::initNULL(bool CleanMem)
{
	initNULLBase(CleanMem);
	if(CleanMem){
		if(fDat != NULL)
			delete fDat;

		if(Filter != NULL){
			for(int iF=0; iF<NFil; iF++)
				delete Filter[iF];
			MT->Free(Filter, NFil*sizeof(LCFilter*));
		}
		if(TmpFilData != NULL){
			for(int iF=0; iF<NFil; iF++)
				delete TmpFilData[iF];
			MT->Free(TmpFilData, NFil*sizeof(LCSerie2*));
		}
    }
	
	fDat = NULL;
	InterpType = LAG;
	InterpUtilValue = 3;
	
	RawData = NULL;
	Filter = NULL;
	TmpFilData = NULL;
	NFil = 0;
    Factor = 1;
    IndexRecToRead = 1;
    tStartReadData=0.;
    tMaxReadData=0.;
    tReadData = 0.;
	
	
	NbDataStab = 1000;
	fSpec.resize(0);
	
}


void LCNoiseFile::config(ezxml_t xmlbloc)
{
	char * SpectralType(NULL);
    
    ezxml_t param;
	ezxml_t section;
    
	setName(ezxml_attr(xmlbloc,"Name"));

    
	bool NoDatFound(true);
	for (section = ezxml_child(xmlbloc, "XSIL"); section; section = section->next) {
		if((NoDatFound)&&(MT->wcmp(ezxml_attr(section,"Type"),"TimeSeries"))){
			NoDatFound = false;
			fDat = new LCDataFileRead(MT);
			fDat->config(section);
		}
	}
    
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		if(MT->wcmp(ezxml_attr(param,"Name"),"DataType")){
			MT->stripcopy((*param).txt, SpectralType);
			TypeFilter = -1;
            
			if(MT->wcmp(SpectralType, "Frequency")){
				TypeFilter = 0;
                Factor *= 1. ;
            }
            if(MT->wcmp(SpectralType, "Length")){
				TypeFilter = 1;
                Factor *= 2.*M_PI/LC::c_SI ;
            }
			
			if(MT->wcmp(SpectralType, "Phase")){
				TypeFilter = 1;
                Factor *= 2.*M_PI/LC::c_SI ; // TO BE CORRECTED 
            }
            
			if(TypeFilter == -1)
				throw std::invalid_argument("ERROR in LCNoiseFile::config : The type of input data of the xml bloc is unknow (use only Frequency, Phase or Length )");
		}
        
		if(MT->wcmp(ezxml_attr(param,"Name"),"Factor")){
			Factor = atof((*param).txt);
        }
        
        if(MT->wcmp(ezxml_attr(param,"Name"),"IndexRecord")){
			IndexRecToRead = atof((*param).txt);
        }
        
        if(MT->wcmp(ezxml_attr(param,"Name"),"StartShift")){
			tReadData += atof((*param).txt);
        }
	}
	
	if(SpectralType != NULL)
		MT->Free(SpectralType, (strlen(SpectralType)+1) * sizeof(char) );
}

void LCNoiseFile::config(int iParam, double ParamVal)
{
	switch (iParam) {
		default:
			Cout << "ERROR in LCNoiseFile::config : Unknow type of parameter : " << iParam << Endl;
			throw std::invalid_argument("ERROR in LCNoiseFile::config : The type of parameters is unknown !");
			break;
	}
}


void LCNoiseFile::init()
{
	initBase();
	
	
    //! ***** Prepare reading of input data
    if(fDat==NULL)
		throw std::invalid_argument("ERROR in LCNoiseFile::init : There is no file to read !");
	fDat->init();
    
    tStartReadData = fDat->getx0() + fDat->getdx() * InterpUtilValue;
    tMaxReadData = fDat->getxend() - fDat->getdx() * InterpUtilValue;
    tReadData += tStartReadData;
    
    if(IndexRecToRead>fDat->getNRec()){
        Cout << "ERROR in LCNoiseFile::init : Index to read ("<< IndexRecToRead << ") is higher than the number of records ("<<fDat->getNRec()<<") !" << Endl;
		throw std::invalid_argument("ERROR in LCNoiseFile::init : Index to read is higher than the number of records !");
    }
    
    
    
	//! ***** Computation of the filter coefficient
	//! Coefficient alpha and beta contained in 3 D structure : Number of filters x Number of cells in the filter x Index in the cell 
	std::vector< std::vector< std::vector<double> > > Alpha(0);
	std::vector< std::vector< std::vector<double> > > Beta(0);
	
    
	/*! ****** Filter with \f$ PSD(f) = PSD_0 f^2 \f$
	 * This filter is construct if it's :
	 * a direct noise in f
	 * or a white phase noise for producing frequency data = blue frequency noise
	 */
	if(TypeFilter==1){
		//! ** One filter
		Alpha.resize(1);
		Beta.resize(1);
		
		//! * One cell in the filter
		(Alpha[0]).resize(1);
		(Beta[0]).resize(1);
		
		//! * Coefficients : \f$ \alpha_0 = -1 , \beta_0 = - \beta_1 = PSD_0 / (\pi \Delta t_{phy}) \f$
		(Alpha[0][0]).push_back(-1.0);
		(Beta[0][0]).push_back(1/(M_PI*dtPhy));
		(Beta[0][0]).push_back(-1.0*Beta[0][0][0]);
		
	}
	
	
	/*! ****** Filter with \f$ PSD(f) = PSD_0 f^{-2} \f$
	 * This filter is construct if it's :
	 * a direct noise in 1/f
	 * or a red frequency noise
	 */
	if(TypeFilter==2){
		//! ** One filter
		Alpha.resize(1);
		Beta.resize(1);
		
		//! * One cell in the filter
		(Alpha[0]).resize(1);
		(Beta[0]).resize(1);
		
		//! * Coefficients : \f$ \alpha_0 = 1 , \beta_0 = \beta_1 = PSD_0 \pi \Delta t_{phy} \f$
		(Alpha[0][0]).push_back(1.0);
		(Beta[0][0]).push_back(M_PI*dtPhy);
		(Beta[0][0]).push_back(Beta[0][0][0]);	
	}
	
	
    
	
	
	//! *** Delete previous filter
	if(Filter != NULL){
		for(int iF=0; iF<NFil; iF++)
			delete Filter[iF];
		MT->Free(Filter, NFil*sizeof(LCFilter*));
	}
	if(TmpFilData != NULL){
		for(int iF=0; iF<NFil; iF++)
			delete TmpFilData[iF];
		MT->Free(TmpFilData, NFil*sizeof(LCSerie2*));
	}
	
	
	/*!
	 * Allocation and initialization. Number of #Filter , #TmpFilData and #RawData :
	 * - #NFil = 0 ->  0 #Filter , 0 #TmpFilData , (1) #RawData (pointer on #NoiseData)
	 * - #NFil = 1 ->  1 #Filter , 0 #TmpFilData , 1 #RawData 
	 * - #NFil = 2 ->  2 #Filter , 2 #TmpFilData , 2 #RawData 
	 * - #NFil = 3 ->  3 #Filter , 3 #TmpFilData , 3 #RawData 
	 */
	
	
	NFil = Alpha.size();
	
	if(NFil == 0){
		//! *** If no filter is needed, it's a white noise and the raw data are the noise data
		Filter = NULL;
		RawData = (LCSerie2**) MT->AllocMemory(sizeof(LCSerie2*));
		RawData[0] = NoiseData;
	}else{
		//! *** If filter is needed :
		//! ** create raw data (white noise)
		RawData = (LCSerie2**) MT->AllocMemory(NFil*sizeof(LCSerie2*));
		for(int iF=0; iF<NFil; iF++){
			RawData[iF] = new LCSerie2(MT, NoiseData->getRef(0), NoiseData->getRefStep(), NoiseData->getNmax());
			RawData[iF]->allocAll();
		}
		
		//! ** Create temporary filtered data
		if(NFil>1){
			TmpFilData = (LCSerie2**) MT->AllocMemory(NFil*sizeof(LCSerie2*));
			for(int iF=0; iF<NFil; iF++){
				TmpFilData[iF] = new LCSerie2(MT, NoiseData->getRef(0), NoiseData->getRefStep(), NoiseData->getNmax());
				TmpFilData[iF]->allocAll();
			}
		}
		
		//! ** create, link and initialize the filter
		Filter = (LCFilter**) MT->AllocMemory(NFil*sizeof(LCFilter*));
		for(int iF=0; iF<NFil; iF++){
			Filter[iF] = new LCFilter(MT);
			Filter[iF]->setRawdata(RawData[iF]);
			if(NFil==1)
				Filter[iF]->setFilteredData(NoiseData);
			else
				Filter[iF]->setFilteredData(TmpFilData[iF]);
			Filter[iF]->init(Alpha[iF], Beta[iF], NbDataStab);
		}
	}
	
	//! *** Delete memory for coefficients
	for(int iF=0; iF<NFil; iF++){
		for(int i=0; i<Alpha[iF].size(); i++)
			Alpha[iF][i].clear();
		Alpha[iF].clear();
		for(int i=0; i<Beta[iF].size(); i++)
			Beta[iF][i].clear();
		Beta[iF].clear();
	}
	
	
	
	//! *** Fill the noise serie
	if(Filter != NULL){
		int MaxDepth(1);
		for(int iF=0; iF<NFil; iF++){
			int tmpDepth(Filter[iF]->getDepth());
			if(MaxDepth < tmpDepth){
				MaxDepth = tmpDepth;
			}
		}
		generNoise(NoiseData->getNmax()-MaxDepth);
		
	}else{
		
		generNoise(NoiseData->getNmax()-1);
	}
	
	//! *** Stabilization of filter
	for(int i=0; i<NbDataStab; i++)
		generNoise(0);
	
}


void LCNoiseFile::generNoise(int iStartBin)
{
    //Cout << "generNoise " << tReadData << "  " << iStartBin << Endl;
	//! *** Read data from file and applying the factor 
	double tmpN, tmpR;
	int NRaw(MAX(NFil,1));
	for(int iF=0; iF<NRaw; iF++){
		for(int i=iStartBin; i>=0; i--){
			tmpN = Factor * fDat->gData(0, tReadData, InterpType, InterpUtilValue) ;
            tReadData += dtPhy;
            if(tReadData>tMaxReadData)
                tReadData = tStartReadData;
			RawData[iF]->addData(tmpN);
            //Cout << tReadData << " " << tmpN << Endl; 
		}
	}
	
    //! *** Filtering
	if(Filter != NULL){
		//! *** Application of the filter
		for(int iF=0; iF<NFil; iF++)
			Filter[iF]->App(iStartBin);
		if(NFil>1){
			for(int i=iStartBin; i>=0; i--){
				tmpR = 0.;
				for(int iF=0; iF<NFil; iF++)
					tmpR += TmpFilData[iF]->getBinValue(i);
				NoiseData->addData(tmpR);
			}
		}		
	}
	
}


void LCNoiseFile::DispInfo(char * BTab)
{
	if(MT->Disp()){
		DispInfoBase(BTab);
		Cout << BTab << "\t- type filter = " << TypeFilter << Endl;
		Cout << BTab << "\t- factor = " << Factor << Endl;
        Cout << BTab << "\t- index of record to read = " << IndexRecToRead << Endl;
        if(fDat != NULL)
			fDat->ControlDisplay();
		char BTab2[128];
		sprintf(BTab2, "%s\t" , BTab);
		if(Filter != NULL){
			for(int iF=0; iF<NFil; iF++)
				Filter[iF]->DispInfo(BTab2);
		}
        Cout << BTab << "\t- Use data from " << tStartReadData << " s to " << tMaxReadData << " s : tCurrent " << tReadData << " s " << Endl;
        
	}
}		

// ***********************
// ***  Local mehtods  ***
// ***********************


// ***************************************
// * Linking and initialization methods  *
// ***************************************





// end of LISACODE-Noise.cpp



