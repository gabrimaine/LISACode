// $Id:  Exp $
/*
 *  LISACode-Filter.cpp
 *  LISACode V 2.0
 *
 *  Created on 08/04/11 by Antoine PETITEAU (AEI)
 *  Last modification on 08/04/11 by Antoine PETITEAU (AEI)
 *
 */


#include "LISACODE-Filter.h"

// *****************
// *  Constructor  *
// *****************

LCFilter::LCFilter()
{
	MT = new LCTools;
	MT->LocTools = true;
	initNULL(false);
}


LCFilter::LCFilter(LCTools * MT_n)
{
	MT = MT_n;
	initNULL(false);
}


LCFilter::~LCFilter()
{
	initNULL(true);
}




// ***************************************
// * Linking and initialization methods  *
// ***************************************


void LCFilter::initNULL(bool CleanMem)
{
	if(CleanMem){
		
		if(Alpha != NULL){
			for(int iI=0; iI<NFil; iI++)
				if(Alpha[iI] != NULL)
					MT->Free(Alpha[iI], NAlpha[iI]*sizeof(double));
			MT->Free(Alpha, NFil*sizeof(double*));
		}
		
		if(NAlpha != NULL)
			MT->Free(NAlpha, NFil*sizeof(int));
		
		if(Beta != NULL){
			for(int iI=0; iI<NFil; iI++)
				if(Beta[iI] != NULL)
					MT->Free(Beta[iI], NBeta[iI]*sizeof(double));
			MT->Free(Beta, NFil*sizeof(double*));
		}
		
		if(NBeta != NULL)
			MT->Free(NBeta, NFil*sizeof(int));
		
		for(int iI=0; iI<NFil-1; iI++)
			delete TmpData[iI];
		MT->Free(TmpData, NFil*sizeof(LCSerie2*));
		
	}
	
	
	RawData = NULL;
	FilData = NULL;
	TmpData = NULL;
	NFil = 0;
	Alpha = NULL;
	NAlpha = NULL;
	Beta = NULL;
	NBeta = NULL;
	NbDataStab = 0;
}


void LCFilter::init(std::vector< std::vector<double> > Alpha_n,
					std::vector< std::vector<double> > Beta_n,
					int NbDataStabilization_n)
{
	//! ***** Allocation and copy of alpha coefficients
	NFil = Alpha_n.size();
	Alpha = (double**) MT->AllocMemory(NFil*sizeof(double*));
	NAlpha = (int*) MT->AllocMemory(NFil*sizeof(int));
	for(int iI=0; iI<NFil; iI++){
		NAlpha[iI] = Alpha_n[iI].size();
		Alpha[iI] = (double*) MT->AllocMemory(NAlpha[iI]*sizeof(double));
		for(int iC=0; iC<NAlpha[iI]; iC++)
			Alpha[iI][iC] = Alpha_n[iI][iC];
	}
	
	//! ***** Allocation and copy of beta coefficients
	if(NFil!=Beta_n.size()){
		std::cerr << "ERROR in LCFilter::init : The number of intermediate lists of coefficients in alpha (" << NFil << ") is not the same than in beta (" << Beta_n.size() << ") !" << Endl;
		throw std::invalid_argument("ERROR in LCFilter::init : The number of intermediate lists of coefficients in alpha and in beta !");
	}
	Beta = (double**) MT->AllocMemory(NFil*sizeof(double*));
	NBeta = (int*) MT->AllocMemory(NFil*sizeof(int));
	for(int iI=0; iI<NFil; iI++){
		NBeta[iI] = Beta_n[iI].size();
		Beta[iI] = (double*) MT->AllocMemory(NBeta[iI]*sizeof(double));
		for(int iC=0; iC<NBeta[iI]; iC++)
			Beta[iI][iC] = Beta_n[iI][iC];
	}
	
	//! ***** Allocation for intermediate file
	TmpData = (LCSerie2**) MT->AllocMemory(NFil*sizeof(LCSerie2*));
	for(int iI=0; iI<NFil-1; iI++){
		TmpData[iI] = new LCSerie2(MT, RawData->getRef(0), RawData->getRefStep(), RawData->getNmax());
		TmpData[iI]->allocAll();
	}
	
	//! ***** For stabilization of the filter
	NbDataStab = NbDataStabilization_n;
	
	
	
	
}



void LCFilter::init(double at, 
					double bp, 
					double fb,
					double fa)
{
	std::vector< std::vector<double> > alpha_n;
	std::vector< std::vector<double> > beta_n;
	double CellsCoef[30][5];
	int NCellsOut(0);
	std::vector<double> TmpVect;
	
	//! *** Sampling frequency in hertz
	double fe(1./RawData->getRefStep());
	
	if(MT->Disp())
		Cout << "\t >>> Make Elliptic filter : fe=" << fe << "Hz att=" << at << "dB osc=" << bp << "dB fb=" << fb << "Hz fa=" << fa << "Hz";
	
	
	//! *** Computed the cells of ecliptic filter
	CalcEllipticFilter4LISACode(fe, at, bp, fb, fa, 30, CellsCoef, &NCellsOut);
	
	if(MT->Disp())
		Cout << " " << NCellsOut << " Cells " << Endl; 
	for(int iCell=0; iCell<NCellsOut; iCell++){
		//cout << "   Cellule " << iCell << " : " << CellsCoef[iCell][0] << " , " << CellsCoef[iCell][1] << " , " << CellsCoef[iCell][2] << " , " << CellsCoef[iCell][3] << " , " << CellsCoef[iCell][4] << Endl;
		TmpVect.resize(0);
		TmpVect.push_back(-1.*CellsCoef[iCell][0]);
		TmpVect.push_back(-1.*CellsCoef[iCell][1]);
		alpha_n.push_back(TmpVect);
		TmpVect.resize(0);
		TmpVect.push_back(CellsCoef[iCell][2]);
		TmpVect.push_back(CellsCoef[iCell][3]);
		TmpVect.push_back(CellsCoef[iCell][4]);
		beta_n.push_back(TmpVect);
	}
	
	init(alpha_n, beta_n, 1000./fe);
	
}


// ********************
// *  Access methods  *
// ********************

int LCFilter::getDepth()
{
	int iRes(0);
	for(int i=0; i<NFil; i++){
		if(iRes<NAlpha[i])
			iRes = NAlpha[i];
		if(iRes<NBeta[i])
			iRes = NBeta[i];
	}
	return(iRes);
}


// *********************
// *  Running methods  *
// *********************

void LCFilter::App(int StartBin)
{
	double Num_tmp (0.0);
	double Den_tmp (0.0);
	LCSerie2 * InpDat;
	LCSerie2 * OutDat;
	
	
	for(int iFil=0; iFil<NFil; iFil++){
		
		//! *** Define the input data
		if(iFil == 0)
			InpDat = RawData;
		else 
			InpDat = TmpData[iFil-1];

		//! *** Define the output data
		if(iFil == NFil-1)
			OutDat = FilData;
		else 
			OutDat = TmpData[iFil];
		
		//! *** Apply filter : \f[ y_{n} = \sum_{k=1}^{N_{\alpha} } \alpha_{k} \;y_{n-k} +  \sum_{k=0}^{N_{\beta} } \; \beta_{k} x_{n-k}  \f]
		//if(MT->DispDet())
		//	Cout << "++++++++++++++++++++++++++++++++ Apply filter " << iFil << Endl;
		for(int i=StartBin; i>=0; i--){
			Num_tmp = 0.0;
			Den_tmp = 0.0;
			
			//! ** If there is enough data, compute the new filtered data  
			if(NAlpha[iFil]+1 < OutDat->getNbVal()){ 
				//! ** Computation of recursive part (from denominator): \f$ \sum_{k=1}^{N_{\alpha} } \; \alpha_{k} \; y_{n-k} \f$
				for(int k=0; k<NAlpha[iFil]; k++){
					Den_tmp += Alpha[iFil][k]*OutDat->getBinValue(k);
					//if(MT->DispDet())
					//	Cout << "\t Den: " << k << " : " << Alpha[iFil][k] << " x " << OutDat->getBinValue(k) << " = " << Alpha[iFil][k]*OutDat->getBinValue(k) << " ==> Sum = " << Den_tmp << Endl;
				}
				
				//! ** Computation of direct part (from numerator): \f$ \sum_{k=0}^{N_{\beta} } \; \beta_{k} \; x_{n-k} \f$
				for(int k=0; k<(int)(NBeta[iFil]); k++){
					Num_tmp += Beta[iFil][k]*InpDat->getBinValue(k+i);
					//if(MT->DispDet())
					//	Cout << "\t Num: " << k << " : " << Beta[iFil][k] << " x " << InpDat->getBinValue(k+i) << " = " << Beta[iFil][k]*InpDat->getBinValue(k+i) << " ==> Sum = " << Num_tmp << Endl;
				}
			}
			
			//! ** Computation of final data and add it to the data
			OutDat->addData(Num_tmp + Den_tmp);
			//if(MT->DispDet())
			//	Cout << i << " In = " << InpDat->getBinValue(i) << "  , Out = " << OutDat->getBinValue(0) << Endl;
		}
		
	}
	
}



void LCFilter::DispInfo(char * BTab)
{
	if(MT->Disp()){
		Cout << BTab << "Filter :" << Endl;
		for(int iC=0; iC<NFil; iC++){
			Cout << BTab << "\t- cell " << iC << " :" << Endl;
			Cout << BTab << "\t\t + alpha =";
			for(int iA=0; iA<NAlpha[iC]; iA++)
				Cout << " " << Alpha[iC][iA];
			Cout << Endl;
			Cout << BTab << "\t\t + beta =";
			for(int iB=0; iB<NBeta[iC]; iB++)
				Cout << " " << Beta[iC][iB];
			Cout << Endl;
		}
	}
}


// end of LISACODE-Filter.cpp
