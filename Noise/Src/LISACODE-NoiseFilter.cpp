/*
 *  LISACODE-NoiseFilter.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 09/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


#include "LISACODE-NoiseFilter.h"


// *********************
// ***  Constructor  ***
// *********************

/**
 *
 */
LCNoiseFilter::LCNoiseFilter()
: LCNoise()
{
	initNULL(false);
}


LCNoiseFilter::LCNoiseFilter(LCTools * MT_n)
: LCNoise(MT_n)
{
	initNULL(false);
}


LCNoiseFilter::~LCNoiseFilter()
{
	initNULL(true);
}



// ********************************************
// ***        Required methods              ***
// ********************************************

void LCNoiseFilter::initNULL(bool CleanMem)
{
	initNULLBase(CleanMem);
	
	if(CleanMem){
		if(RawData != NULL){
			if(NFil!=0){
				for(int iF=0; iF<NFil; iF++)
					delete RawData[iF];
				MT->Free(RawData, NFil*sizeof(LCSerie2*));
			}else{
				// ** Case of white noise : just a pointer of raw data (to noise data) as been allocated
				MT->Free(RawData, sizeof(LCSerie2*));
			}
		}
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
	
	RawData = NULL;
	Filter = NULL;
	TmpFilData = NULL;
	NFil = 0;
	
	SqPSD = 1.0;
	
	NbDataStab = 1000;
	TypeNoise = -1;
	fLow = 0.;
	fknee = 0.;
	slope = 0.;
	fSpec.resize(0);
	
}


void LCNoiseFilter::config(ezxml_t noisexmlbloc)
{
	ezxml_t param;
	
	setName(ezxml_attr(noisexmlbloc,"Name"));
	
	char * SpectralType(NULL);
	
	for(param = ezxml_child(noisexmlbloc,"Param"); param; param = param->next){
		if(MT->wcmp(ezxml_attr(param,"Name"),"SpectralType")){
			MT->stripcopy((*param).txt, SpectralType);
			
			if(MT->wcmp(SpectralType, "White"))
				TypeNoise = 0;
			
			if(MT->wcmp(SpectralType, "Filter_f"))
				TypeNoise = 1;
			
			if(MT->wcmp(SpectralType, "Filter_1of"))
				TypeNoise = 2;
			
			if(MT->wcmp(SpectralType, "Filter_1of_1of2"))
				TypeNoise = 3;
			
			if(MT->wcmp(SpectralType, "Filter_1of_1of32"))
				TypeNoise = 5;
			
			if((MT->wcmp(SpectralType, "WhitePhase"))||(MT->wcmp(SpectralType, "BlueFrequency")))
				TypeNoise = 10;
			
			if(MT->wcmp(SpectralType, "WhiteFrequency"))
				TypeNoise = 11;
			
			if(MT->wcmp(SpectralType, "RedFrequency"))
				TypeNoise = 12;
			
			if((MT->wcmp(SpectralType, "PinkFrequency"))||(MT->wcmp(SpectralType, "PinkAcceleration"))||(MT->wcmp(SpectralType, "PinkPhase")))
				TypeNoise = 13;
			
			if(MT->wcmp(SpectralType, "PreStabLaserNoiseFreq"))
				TypeNoise = 14;
			
			if(TypeNoise == -1)
				throw std::invalid_argument("ERROR in LCNoiseFilter::config : The spectral type of the xml bloc is unknow as potential filtered noise.");
		}
		if(MT->wcmp(ezxml_attr(param,"Name"),"PowerSpectralDensity")){
			double PSD;
			PSD = atof((*param).txt);
			SqPSD = sqrt(PSD);
		}
		if(MT->wcmp(ezxml_attr(param,"Name"),"SpectralSlope")){
			TypeNoise = 4;
			slope = atof((*param).txt);
		}
		if(MT->wcmp(ezxml_attr(param,"Name"),"Fknee")){
			fknee = atof((*param).txt);
		}
		if(MT->wcmp(ezxml_attr(param,"Name"),"Flow")){
			fLow = atof((*param).txt);
		}
		if(MT->wcmp(ezxml_attr(param,"Name"),"SpecificFrequency")){
			fSpec.push_back(atof((*param).txt));
		}
	}
	
	if(SpectralType != NULL)
		MT->Free(SpectralType, (strlen(SpectralType)+1) * sizeof(char) );
}

void LCNoiseFilter::config(int iParam, double ParamVal)
{
	switch (iParam) {
		default:
			Cout << "ERROR in LCNoiseFilter::config : Unknow type of parameter : " << iParam << Endl;
			throw std::invalid_argument("ERROR in LCNoiseFilter::config : The type of parameters is unknown !");
			break;
	}
}


void LCNoiseFilter::init()
{
	initBase();
	
	
	//! ***** Computation of the filter coefficient
	//! Coefficient alpha and beta contained in 3 D structure : Number of filters x Number of cells in the filter x Index in the cell 
	std::vector< std::vector< std::vector<double> > > Alpha(0);
	std::vector< std::vector< std::vector<double> > > Beta(0);
	
	// \todo TODO Include the full computation of this noise in case of phase data
	
	
	/*! ****** Filter with \f$ PSD(f) = PSD_0 f^2 \f$
	 * This filter is construct if it's :
	 * a direct noise in f
	 * or a white phase noise for producing frequency data = blue frequency noise
	 */
	if((TypeNoise==1)||((TypeNoise==10)&&(IsFreqOutput))){
		//! ** One filter
		Alpha.resize(1);
		Beta.resize(1);
		
		//! * One cell in the filter
		(Alpha[0]).resize(1);
		(Beta[0]).resize(1);
		
		//! * Coefficients : \f$ \alpha_0 = -1 , \beta_0 = - \beta_1 = PSD_0 / (\pi \Delta t_{phy}) \f$
		(Alpha[0][0]).push_back(-1.0);
		(Beta[0][0]).push_back(SqPSD/(M_PI*dtPhy));
		(Beta[0][0]).push_back(-1.0*Beta[0][0][0]);
		
	}
	
	
	/*! ****** Filter with \f$ PSD(f) = PSD_0 f^{-2} \f$
	 * This filter is construct if it's :
	 * a direct noise in 1/f
	 * or a red frequency noise
	 */
	if((TypeNoise==2)||((TypeNoise==12)&&(IsFreqOutput))){
		//! ** One filter
		Alpha.resize(1);
		Beta.resize(1);
		
		//! * One cell in the filter
		(Alpha[0]).resize(1);
		(Beta[0]).resize(1);
		
		//! * Coefficients : \f$ \alpha_0 = 1 , \beta_0 = \beta_1 = PSD_0 \pi \Delta t_{phy} \f$
		(Alpha[0][0]).push_back(1.0);
		(Beta[0][0]).push_back(SqPSD*M_PI*dtPhy);
		(Beta[0][0]).push_back(Beta[0][0][0]);	
	}
	
	
	
	/*! ****** Filter with \f$ PSD(f) = PSD_0 \left({1 \over f^{2}} + {f_{knee}^2 \over f^{4}} \right)  \f$
	 * This filter is construct if it's :
	 * a pink frequency noise
	 */
	if((TypeNoise==13)&&(IsFreqOutput)){
		
		//! ** Two filter
		Alpha.resize(2);
		Beta.resize(2);
		
		//! Correction to saturate at low frequencies : \f$ C_{sat} = 1 - 10^{-14} \f$
		double coef_alpha_ep(1-1.e-14); 
		
		//! ** First filter :  \f$ \sqrt{PSD_{0}} f^{-1}  \f$
		//! * One cell of the filter
		(Alpha[0]).resize(1);
		(Beta[0]).resize(1);
		
		//! * Coefficients : \f[ \alpha_0 = C_{sat} , \beta_0 = \beta_1 = \sqrt{PSD_0} \pi \Delta t \f]
		(Alpha[0][0]).push_back(coef_alpha_ep);
		(Beta[0][0]).push_back(SqPSD * M_PI * dtPhy);
		(Beta[0][0]).push_back(Beta[0][0][0]);
		
		//! ** Second filter : \f$  \sqrt{PSD_0} f_{knee} f^{-2} \f$
		//! * One cell of the filter
		(Alpha[1]).resize(1);
		(Beta[1]).resize(1);
		
		//! * Coefficients : \f[ \alpha_0 = 2 C_{sat} , \alpha_1 = - C_{sat} , \beta_0 = \beta_2 = C , \beta_1 = 2 C \f] with \f[ C = \sqrt{PSD_0} f_{knee} \pi^2 \Delta t_{phy}^2 \f]
		(Alpha[1][0]).push_back(2.0*coef_alpha_ep);
		(Alpha[1][0]).push_back(-1.0*coef_alpha_ep);
		(Beta[1][0]).push_back( SqPSD * fknee * M_PI*M_PI * dtPhy*dtPhy);
		(Beta[1][0]).push_back(2.0*Beta[1][0][0]);
		(Beta[1][0]).push_back(Beta[1][0][0]);
		
		
	}
	
	
	/*! ****** Filter with \f$ PSD(f) = PSD_0 {\left({1 \over f} + {f_{knee} \over f^{2}} \right)}^2  \f$
	 * This filter is construct if it's :
	 * a direct noise in 1/f + fk/f^2
	 */
	if((TypeNoise==3)){
		//! *** Tree filter
		Alpha.resize(3);
		Beta.resize(3);
		
		//! Correction to saturate at low frequencies : \f$ C_{sat} = 1 - 10^{-14} \f$
		double coef_alpha_ep(1-1.e-14); 
		
		int iF;
		
		//! ** First filter : \f$ \sqrt{PSD_{0}} f^{-1}  \f$
		iF = 0;
		//! * One cell of the filter
		(Alpha[iF]).resize(1);
		(Beta[iF]).resize(1);
		
		//! * Coefficients : \f[ \alpha_0 = C_{sat} , \beta_0 = \beta_1 = \sqrt{PSD_0} \pi \Delta t \f]
		(Alpha[iF][0]).push_back(coef_alpha_ep);
		(Beta[iF][0]).push_back(SqPSD * M_PI * dtPhy);
		(Beta[iF][0]).push_back(Beta[iF][0][0]);
		
		
		//! ** Second filter : \f$  2 \sqrt{PSD_0} f_{knee} f^{-3/2}\f$
		iF = 1;
		CoefOofFilter(Alpha[iF], Beta[iF], 1.e-6, 1./(2.*dtPhy), -3.0, sqrt(2.*fknee)*SqPSD * pow(2.*dtPhy, 3./2.) );
		
		
		//! ** Third filter : \f$  \sqrt{PSD_0} f^2_{knee} f^{-2} \f$
		iF = 2;
		//! * Two cells of the filter
		(Alpha[iF]).resize(1);
		(Beta[iF]).resize(1);
		
		
		//! * Coefficients : \f[ \alpha_0 = 2 C_{sat} , \alpha_1 = - C_{sat} , \beta_0 = \beta_2 = C , \beta_1 = 2 C \f] with \f[ C = \sqrt{PSD_0} f_{knee} \pi \Delta t_{phy}^2 \f]
		double C ( SqPSD * M_PI * M_PI * fknee * dtPhy * dtPhy  );
		(Alpha[iF][0]).push_back(2.*coef_alpha_ep);
		(Alpha[iF][0]).push_back(-1.*coef_alpha_ep);
		(Beta[iF][0]).push_back(C);
		(Beta[iF][0]).push_back(2.*C);
		(Beta[iF][0]).push_back(C);
		
	}
	
	/*! ****** Filter with \f$ PSD(f) = PSD_0 f^{a} \f$
	 * This filter is construct if it's :
	 * a direct noise in f^a
	 */
	if((TypeNoise==4)){
		//! *** One filter
		Alpha.resize(1);
		Beta.resize(1);
		
		//! ** First filter : \f$  2 \sqrt{PSD_0} f_{knee} f^{-3/2}\f$
		CoefOofFilter(Alpha[0], Beta[0], fLow, fknee, slope, SqPSD * pow(2.*dtPhy, -slope/2.) );
		
	}
	
	/*! ****** Filter with \f$ PSD(f) = PSD_0 \left({1 \over f^2} + {f_{knee} \over f^{3}} \right)  \f$
	 * This filter is construct if it's :
	 * a direct noise in 1/f + fk/f^2
	 */
	if((TypeNoise==5)){
		//! *** Two filter
		Alpha.resize(2);
		Beta.resize(2);
		
		//! Correction to saturate at low frequencies : \f$ C_{sat} = 1 - 10^{-14} \f$
		double coef_alpha_ep(1-1.e-14); 
		
		int iF;
		
		
		//! ** First filter : \f$ \sqrt{PSD_{0}} f^{-1}  \f$
		iF = 0;
		//! * One cell of the filter
		(Alpha[iF]).resize(1);
		(Beta[iF]).resize(1);
		
		//! * Coefficients : \f[ \alpha_0 = C_{sat} , \beta_0 = \beta_1 = \sqrt{PSD_0} \pi \Delta t \f]
		(Alpha[iF][0]).push_back(coef_alpha_ep);
		(Beta[iF][0]).push_back(SqPSD * M_PI * dtPhy);
		(Beta[iF][0]).push_back(Beta[iF][0][0]);
		
		
		//! ** Second filter : \f$  2 \sqrt{PSD_0} f_{knee} f^{-3/2}\f$
		iF = 1;
		CoefOofFilter(Alpha[iF], Beta[iF], 1.e-7, (1./(2.*dtPhy)), -3.0, sqrt(fknee)*SqPSD * pow(2.*dtPhy, 3./2.) );
		
		
	}
	
	
	/*! ****** Filter with \f$ \tilde{x}(f) = A { (s + \omega_2)^2 \over (s + \omega_1)^2 (s + \omega_3) } \f$
	 * This filter is construct if it's :
	 * a pre-stabilized laser noise in frequency
	 */
	if((TypeNoise==14)&&(IsFreqOutput)){
		if(fSpec.size()<3)
			throw std::invalid_argument("LCNoiseFilter::init : For a pre-stabilized laser noise in frequency we need the 3 specific frequencies !");
		
		double w1(2.0*M_PI*fSpec[0]), w2(2.0*M_PI*fSpec[1]), w3(2.0*M_PI*fSpec[2]);
		double cden(1.0/((1.+dtPhy*w1)*(1.+dtPhy*w1)*(1.+dtPhy*w3)));
		
		//! ** One filter
		Alpha.resize(1);
		Beta.resize(1);
		
		//! * Two cell in the filter
		(Alpha[0]).resize(5);
		(Beta[0]).resize(5);
		
		
		//! ** Filter 0 : \f$ ( s + \omega_3 ) x^{(1)} = A x^{(0)} \f$  
		(Alpha[0][0]).push_back(1./(1.+dtPhy*w3));
		(Beta[0][0]).push_back(dtPhy*pow(SqPSD,1./5.)/(1.+dtPhy*w3));
		
		//! ** Filter 1 : \f$ ( s + \omega_1 ) x^{(2)} = A x^{(1)} \f$  
		(Alpha[0][1]).push_back(1./(1.+dtPhy*w1));
		(Beta[0][1]).push_back(dtPhy*pow(SqPSD,1./5.)/(1.+dtPhy*w1));
		
		//! ** Filter 2 : \f$ ( s + \omega_1 ) x^{(3)} = A x^{(2)} \f$  
		(Alpha[0][2]).push_back(1./(1.+dtPhy*w1));
		(Beta[0][2]).push_back(dtPhy*pow(SqPSD,1./5.)/(1.+dtPhy*w1));
		
		//! ** Filter 3 : \f$  x^{(4)} = ( s + \omega_2 ) A x^{(3)}  \f$  
		(Beta[0][3]).push_back(pow(SqPSD,1./5.)*(1./dtPhy+w2));
		(Beta[0][3]).push_back(-1.*pow(SqPSD,1./5.)/dtPhy);
		
		//! ** Filter 4 : \f$  x^{(5)} = ( s + \omega_2 ) A x^{(4)}  \f$  
		(Beta[0][4]).push_back(pow(SqPSD,1./5.)*(1./dtPhy+w2));
		(Beta[0][4]).push_back(-1.*pow(SqPSD,1./5.)/dtPhy);
		
		
		
		// *** One filter and central derivative
		
		/*! * Recursive coefficients :
		 *	\f$ \alpha_0 = 1  \f$
		 *	\f$ \alpha_1 = -2 \Delta t_{phy} (2 \omega_1 + \omega_3) \f$
		 *	\f$ \alpha_2 = 3 - 4 \Delta t_{phy}^2 \omega_1 (\omega_1 + 2 \omega_3) \f$
		 *	\f$ \alpha_3 = 4 \Delta t_{phy} (2 \omega_1 + \omega_3 - 2 \Delta t_{phy}^2 \omega_1^2 \omega_3 ) \f$
		 *	\f$ \alpha_4 = - \alpha_2 \f$
		 *	\f$ \alpha_5 = - \alpha_1 \f$
		 *	\f$ \alpha_6 = 1 \f$
		 /
		(Alpha[0][0]).push_back(-2.0 * dtPhy * ( 2.0*w1 + w3 ) );
		(Alpha[0][0]).push_back( 3.0 - 4.0 * dtPhy*dtPhy * w1 * ( w1 + 2.0*w3 ));
		(Alpha[0][0]).push_back( 4.0 * dtPhy * (2.0*w1 + w3 - 2.0*dtPhy*dtPhy * w1*w1*w3 ) );
		(Alpha[0][0]).push_back( -1.0 * Alpha[0][0][1] );
		(Alpha[0][0]).push_back( -1.0 * Alpha[0][0][0] );
		(Alpha[0][0]).push_back( 1.0 );
		
		/*! * Direct coefficients : 
		 *	\f$ \beta_0 = 0 \f$
		 *	\f$ \beta_1 = 2 A \Delta t_{phy} \f$
		 *	\f$ \beta_2 = 8 A \Delta t_{phy}^2 \omega_2 \f$
		 *	\f$ \beta_3 = 2 A \Delta t_{phy} ( 4 \omega_2^2 \Delta t_{phy}^2 - 2 ) \f$
		 *	\f$ \beta_4 = - \beta_2 \f$
		 *	\f$ \beta_5 = \beta_1 \f$
		 /
		(Beta[0][0]).push_back(0.0);
		(Beta[0][0]).push_back(2.0*dtPhy*SqPSD);
		(Beta[0][0]).push_back(8.0*dtPhy*dtPhy*SqPSD*w2);
		(Beta[0][0]).push_back(2.0*dtPhy*SqPSD*(4.0*w2*w2*dtPhy*dtPhy-2.0));
		(Beta[0][0]).push_back(-1.0*Beta[0][0][2]);
		(Beta[0][0]).push_back(Beta[0][0][1]);
		 
		 
		// *** One filter and backward derivative 
		/*
		(Alpha[0][0]).push_back(cden*-1.0*(3.0+(w1+2.0*w3)*dtPhy));
		(Alpha[0][0]).push_back(cden*(3.0+(2.0*w1+w3)*dtPhy));
		(Alpha[0][0]).push_back(1.0);
		
		(Beta[0][0]).push_back(cden*SqPSD*dtPhy*(1.0+w2*dtPhy)*(1.0+w2*dtPhy));
		(Beta[0][0]).push_back(-2.0*cden*SqPSD*dtPhy*(1.0+w2*dtPhy));
		(Beta[0][0]).push_back(cden*SqPSD*dtPhy);
		 */
		
		
		
	}
	
	
	
	
	
	//! *** Delete previous raw data and filter
	if((RawData != NULL)&&(NFil!=0)){
		for(int iF=0; iF<NFil; iF++)
			delete RawData[iF];
		MT->Free(RawData, NFil*sizeof(LCSerie2*));
	}
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
	
	
	
	//! *** Compute #Sigma which multiply the output : 
	//! Take into account the time step and the square root of power spectral density for the white noise (No filter : Filter == NULL)
	if(Filter == NULL)
		Sigma = SqPSD/sqrt(2.0*dtPhy);
	else
		Sigma = 1./sqrt(2.0*dtPhy);
	
	
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


void LCNoiseFilter::generNoise(int iStartBin)
{
	//! *** Generation of white noise with sigma = 1
	double tmpN, tmpR;
	int NRaw(MAX(NFil,1));
	for(int iF=0; iF<NRaw; iF++){
		for(int i=iStartBin; i>=0; i--){
			tmpN = Sigma * MT->RandGaussian(0., 1.0) ;
			//tct += dtPhy;
			//tmpN = sin(tct/(20.*iSC+5*IndirectDir));
			//tmpN = 1.0;
			RawData[iF]->addData(tmpN);
		}
	}
	
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
	//! *** Multiply the result by the sigma for white noise
	//for(int i=iStartBin; i>=0; i--)
	//	NoiseData->setBinValue(i, Sigma * NoiseData->getBinValue(i));

	//Cout << NoiseData->getBinValue(0) << Endl;
}


void LCNoiseFilter::DispInfo(char * BTab)
{
	if(MT->Disp()){
		DispInfoBase(BTab);
		Cout << BTab << "\t- type = " << TypeNoise << Endl;
		Cout << BTab << "\t- Sigma = " << Sigma << Endl;
		Cout << BTab << "\t- sqrt(PSD) = " << SqPSD << Endl;
		Cout << BTab << "\t- flow  = " << fLow << " Hz" << Endl;
		Cout << BTab << "\t- fknee = " << fknee << " Hz" << Endl;
		Cout << BTab << "\t- slope = " << slope << Endl;
		for(int iFS=0; iFS<fSpec.size(); iFS++) {
			Cout << BTab << "\t- specific frequency " << iFS << " = " << fSpec[iFS] << " Hz" << Endl;
		}
		char BTab2[128];
		sprintf(BTab2, "%s\t" , BTab);
		if(Filter != NULL){
			for(int iF=0; iF<NFil; iF++)
				Filter[iF]->DispInfo(BTab2);
		}
	}
}		

// ***********************
// ***  Local mehtods  ***
// ***********************


// ***************************************
// * Linking and initialization methods  *
// ***************************************


void LCNoiseFilter::CoefOofFilter(std::vector< std::vector<double> > & alpha_n, std::vector< std::vector<double> > & beta_n, double fmin, double fmax, double fpow, double SqrPSDf0 )
{
	double w0,w1,wmax,wmin;
	double p0,z0,dp;
	double p,z,pcourant,zcourant,dpcourant,den,tStep_inv;
	int nproc;
	
	
	//! *** Calculation of boundaries : \f[ \ w_{0} = 2 \pi f_{min} , \quad w_1 = 2 \pi f_{max} , \quad w_{max} = log_{10} w_1  , \quad w_{min} = log_{10} w_0
	w0   = 2*M_PI*fmin;
	w1   = 2*M_PI*fmax;
	wmin = log10(w0);
	wmax = log10(w1);
	
	tStep_inv = 1./dtPhy;
	
	//! *** Number of poles : \f$ N_p = 2 ( w_{max} - w_{min}) + log_{10} \Delta t \f$
	nproc = (int)((wmax-wmin)*2.+log10(dtPhy));
	
	if(MT->Disp())
		Cout << "Filter 1/f^a : w0 = " << w0 << " , w1 = " << w1 << " ,  wmax = " << wmax << ", wmin = " << wmin << ", tStep = " << dtPhy << ", tStep_inv = " << tStep_inv << ", nproc = " << nproc << Endl;	
	
	alpha_n.resize(nproc);
	beta_n.resize(nproc);
	
	
	
	//! **** Calculation of poles and coefficients :
	/*! *** Initailization :
	 *	Step :\f[ \Delta p = { w_{min} - w_{max} \over N_p } \f]
	 *	Initial $p$ :\f[ p_0 = w_{min} + { 1 + { f_{pow} \over 2}  \over 2 \Delta p } \f]
	 *	Initial $z$ :\f[ z_0 = p_0 + { f_{pow} \over 2}  \Delta p \f]
	 *	Initial state : \f[ p_{cur} = p_0 , \quad z_{cur} = z_0 , \quad \Delta p_{cur} = 0 \f]
	 */
	
	dp = (wmax-wmin)/nproc;
	p0 = wmin+(1+fpow/2)/2*dp;
	z0 = p0-fpow/2*dp;
	pcourant = p0;
	zcourant = z0;
	dpcourant = 0.;
	
	
	//! *** Loop over the number of pole and zero
	for(int i=0; i<nproc; i++){
		/*!	** Calculation of pole and zero :
		 *	\f[ p_i = p_{cur} + \Delta p_{cur} \f]
		 *	\f[	z_i = p_i -	{ f_{pow} \over 2 \Delta p } \f]
		 *	\f[ w_0 = \Delta t * 10^{p_i}  \f]
		 *	\f[ w_1 = \Delta t * 10^{z_i}  \f]
		 */
		p = pcourant+dpcourant;
		dpcourant =	(wmax-wmin)/nproc;
		z = p-fpow/2*dp;
		w0 = (pow(10,p)/2)*dtPhy;
		w1 = (pow(10,z)/2)*dtPhy;
		pcourant = p;
		zcourant = z;
		
		/*!	** Calculation of coefficients :
		 *	\f[ \alpha_1 = { 1 - w_0 \over 1 + w_0 } , \quad
		 *		\beta_0 = { 1 + w_1 \over 1 + w_0 } , \quad
		 *		\beta_1 = - { 1 - w_1 \over 1 + w_0 } \f]
		 */
		den = 1+w0;
		alpha_n[i].resize(1);
		beta_n[i].resize(2);
		alpha_n[i][0]=(1-w0)/den;
		beta_n[i][0]=(1+w1)/den;
		beta_n[i][1]=-(1-w1)/den;
		
		if(MT->Disp())
			Cout << "\t + Filter " << i << " : alpha1 = " << alpha_n[i][0] << " , beta0 = " << beta_n[i][0] << " , beta1 = " << beta_n[i][1] << "  (pole = " << pow(10., pcourant) << " , zero = " << pow(10.,zcourant) << " )" << Endl;
	} 
	
	//! *** Apply the amplitude factor on the beta coefficient of first cell :
	if(beta_n.size()>0){
		beta_n[0][0] *= SqrPSDf0;
		beta_n[0][1] *= SqrPSDf0;
	}
	
}





// end of LISACODE-Noise.cpp



