/*
 *  LISACODE-OBPM.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 12/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


#include "LISACODE-OBPM.h"


// *********************
// ***  Constructor  ***
// *********************

LCOBPM::LCOBPM()
{
	MT = new LCTools;
	MT->LocTools = true;
	initNULLBase(false);
}


LCOBPM::LCOBPM(LCTools * MT_n)
{
	MT = MT_n;
	initNULLBase(false);
}



LCOBPM::~LCOBPM()
{
	initNULLBase(true);
}


void LCOBPM::initNULLBase(bool CleanMem)
{
	if(CleanMem){
		
		if(pn != NULL)
			MT->Free(pn, Npn*sizeof(LCNoise*));
		
		if(MesPho != NULL)
			delete  MesPho;
		MesPho = NULL;
		
		if(MesPha != NULL)
			delete MesPha;
		MesPha = NULL;
		
		if(PMFilter != NULL)
			delete PMFilter;
		
		if(MesFilter != NULL)
			delete MesFilter;
		MesFilter = NULL;
		
		if((LastMes != NULL)&&(LastMesLocalAlloc))
			MT->Free(LastMes, sizeof(double));
	}
	
	strcpy(Name,"");
	ArmGW = NULL;
	pn = NULL;
	Npn = 0;
	clockl = NULL;
	orb = NULL;
	MesPho = NULL;
	MesPha = NULL;
	PMFilter = NULL;
	MesFilter = NULL;
	LastMes = NULL;
	GWSignal = 0.;
	FactShotNoise = NULL;
	
	LastMesLocalAlloc = true;
	
	iSC = -1;
	IndirectDir = -1;
	
	PMFAttenuation  = 180; 
	PMFOscillation  = 0.1;
	PMFLowFreqFact  = 0.1;
	PMFHighFreqFact = 0.3;
	PMFalpha.resize(0);
	PMFbeta.resize(0);
	
	dtMes = 1.;
	dtPhy = 1.;
	tFirst = 0.;
	tLast = 30.;
	N2Add = 1;
	NDataPha = 30;
	
	
}



// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCOBPM::initNULL(bool CleanMem)
{
	initNULLBase(CleanMem);
}


void LCOBPM::config(ezxml_t xmlbloc)
{
	
}


void LCOBPM::config(int iParam, double ParamVal)
{
	
}

void LCOBPM::LinkNoise(LCNoise ** AllNoises, int NAllNoises)
{
	
}

void LCOBPM::LinkFactShotNoise(double * FactShotNoise_n)
{
	FactShotNoise = NULL;
}

void LCOBPM::init()
{
	
}


double LCOBPM::MeasurePho(double TimeLocal, double TShiftEff)
{
	return(0.);
}


double LCOBPM::MeasurePha(double tGlob)
{
	double t(tGlob), TimeLocal(tGlob), TShiftEff(0.);
	
	
	//! ** The GWs signal is computed here to avoid computation at each physical time step but it's less correct than computing in photodiode measurement 
	//GWSignal = sGW(clockl->gT(tGlob));
	
	
	//! **** Compute the photodiode measurement
	for(int i=N2Add-1; i>=0; i--){
		TimeLocal = tGlob - i*dtPhy;
		if(clockl != NULL){
			//! ** Compute the effective time of the measurement
			TimeLocal = clockl->gT(TimeLocal);
		}
		//! ** Shift between effective time #TimeLocal (local time) and time tag of the measurement #tGlob 
		TShiftEff = tGlob - TimeLocal;
		
		(*LastMes) = MeasurePho(TimeLocal, TShiftEff);
		
		//Cout << Name << " : t = " << TimeLocal << " : (*LastMes) = MeasurePho(TimeLocal, TShiftEff) = " << (*LastMes) << Endl;
		
		//! ** Measure by the photodiode
		MesPho->addData((*LastMes));
	}
	
	//! *** Measure by the phasemeter and application of its filter  
	if(PMFilter != NULL){
		
		//if((iSC==1)&&(IndirectDir==0)&&(tGlob>1000.0)){
		//	MT->setDispDetails();
		//}
		
		PMFilter->App(N2Add-1);
		(*LastMes) = MesFilter->getBinValue(0);
		
		//MT->unsetDispDetails();
	}
	
	//Cout << Name << " : t = " << t << " : if(PMFilter != NULL) = " << (*LastMes) << Endl;
	
	//! *** Add data to the buffer of phasemeter measurements 
	//Cout << (*LastMes) << Endl;
	MesPha->addData((*LastMes));
	
	return((*LastMes));
}


void LCOBPM::DispInfo(char * BTab, bool LinkedModule)
{
	
}



// ***********************
// ***  Local mehtods  ***
// ***********************


// ***************************
// *  Configuration methods  *
// ***************************

void LCOBPM::configFilter(ezxml_t xmlbloc)
{
	
	ezxml_t param;
	char * PMFilterType(NULL);
	bool ReqPMFilter(false);
	
	//! **** Delete previous filter
	if(PMFilter != NULL)
		delete PMFilter;
	
	
	//! **** Configuration of filter parameters
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		
		//! *** Pre-registered configurtation 
		if(MT->wcmp(ezxml_attr(param,"Name"),"FilterType")){
			MT->stripcopy((*param).txt, PMFilterType);
			if(MT->wcmp(PMFilterType,"Standard")){
				ReqPMFilter = true;
				PMFAttenuation  = 240; 
				PMFOscillation  = 0.1;
				PMFLowFreqFact  = 0.1;
				PMFHighFreqFact = 0.3;
			}
			if(MT->wcmp(PMFilterType,"Weak")){
				ReqPMFilter = true;
				PMFAttenuation  = 100; 
				PMFOscillation  = 0.1;
				PMFLowFreqFact  = 0.4;
				PMFHighFreqFact = 0.48;
			}
			if(MT->wcmp(PMFilterType,"CombTest")){
				double dtPhytmp(0.03), dtMestmp(0.3);
				int NFilter(10);
				ReqPMFilter = true;
				PMFalpha.resize(NFilter);
				PMFbeta.resize(NFilter);
				for(int i=0; i<NFilter; i++){
					PMFalpha[i].resize(1);
					PMFbeta[i].resize(1);
					PMFalpha[i][0] = 1. / (1. + dtPhytmp * 2.*M_PI*0.5/dtMestmp);
					PMFbeta[i][0] = dtPhytmp / (1. + dtPhytmp * 2.*M_PI*0.5/dtMestmp);
				}
			}
			
			
			if(MT->wcmp(PMFilterType,"Test")){
				ReqPMFilter = true;
				int iFil(0), NCells(21);
				
				
				iFil=0;
				PMFalpha.resize(NCells);
				PMFbeta.resize(NCells);
				for(int i=0; i<NCells; i++){
					PMFalpha[i].resize(2);
					PMFbeta[i].resize(3);
				}
				/*
				//! Filter type butterworth : dtphy=0.03, dtmes=0.3, HighFrF=0.5, LowFrF=2.0, Attenuation=100.0, OscilPass=0.1
				PMFalpha[iFil+0][0]=1.72647865016355295253e+00;
				PMFalpha[iFil+0][1]=-8.73616692264348349006e-01;
				PMFbeta[iFil+0][0]=3.67845105251988352402e-02;
				PMFbeta[iFil+0][1]=7.35690210503976704803e-02;
				PMFbeta[iFil+0][2]=3.67845105251988352402e-02;
				PMFalpha[iFil+1][0]=1.54320532956762135335e+00;
				PMFalpha[iFil+1][1]=-6.74724019781706485510e-01;
				PMFbeta[iFil+1][0]=3.28796725535213177349e-02;
				PMFbeta[iFil+1][1]=6.57593451070426354699e-02;
				PMFbeta[iFil+1][2]=3.28796725535213177349e-02;
				PMFalpha[iFil+2][0]=1.42029594771891654048e+00;
				PMFalpha[iFil+2][1]=-5.41339764235997611408e-01;
				PMFbeta[iFil+2][0]=3.02609541292702746706e-02;
				PMFbeta[iFil+2][1]=6.05219082585405493413e-02;
				PMFbeta[iFil+2][2]=3.02609541292702746706e-02;
				PMFalpha[iFil+3][0]=1.35011056176559107733e+00;
				PMFalpha[iFil+3][1]=-4.65172873517303409052e-01;
				PMFbeta[iFil+3][0]=2.87655779379280517050e-02;
				PMFbeta[iFil+3][1]=5.75311558758561034099e-02;
				PMFbeta[iFil+3][2]=2.87655779379280517050e-02;
				PMFalpha[iFil+4][0]=6.63665449906185700435e-01;
				PMFalpha[iFil+4][1]=-0.00000000000000000000e+00;
				PMFbeta[iFil+4][0]=1.68167275046907094271e-01;
				PMFbeta[iFil+4][1]=1.68167275046907094271e-01;
				PMFbeta[iFil+4][2]=0.00000000000000000000e+00;
				
				iFil+=5;
				
				
				
				//! Filter type butterworth : dtphy=0.03, dtmes=0.3, HighFrF=0.3, LowFrF=1.2, Attenuation=100.0, OscilPass=0.1
				PMFalpha[iFil+0][0]=1.86628585368501709318e+00;
				PMFalpha[iFil+0][1]=-9.25729096347488678020e-01;
				PMFbeta[iFil+0][0]=1.48608106656179066174e-02;
				PMFbeta[iFil+0][1]=2.97216213312358132348e-02;
				PMFbeta[iFil+0][2]=1.48608106656179066174e-02;
				PMFalpha[iFil+1][0]=1.74315654941839048853e+00;
				PMFalpha[iFil+1][1]=-7.98677989267037813370e-01;
				PMFbeta[iFil+1][0]=1.38803599621618312110e-02;
				PMFbeta[iFil+1][1]=2.77607199243236624220e-02;
				PMFbeta[iFil+1][2]=1.38803599621618312110e-02;
				PMFalpha[iFil+2][0]=1.65052545261811678401e+00;
				PMFalpha[iFil+2][1]=-7.03096490868682422182e-01;
				PMFbeta[iFil+2][0]=1.31427595626414026042e-02;
				PMFbeta[iFil+2][1]=2.62855191252828052084e-02;
				PMFbeta[iFil+2][2]=1.31427595626414026042e-02;
				PMFalpha[iFil+3][0]=1.58917018379793462124e+00;
				PMFalpha[iFil+3][1]=-6.39786989728784671883e-01;
				PMFbeta[iFil+3][0]=1.26542014827125039872e-02;
				PMFbeta[iFil+3][1]=2.53084029654250079744e-02;
				PMFbeta[iFil+3][2]=1.26542014827125039872e-02;
				PMFalpha[iFil+4][0]=1.55870820320291980643e+00;
				PMFalpha[iFil+4][1]=-6.08354761783443831469e-01;
				PMFbeta[iFil+4][0]=1.24116396451309802390e-02;
				PMFbeta[iFil+4][1]=2.48232792902619604780e-02;
				PMFbeta[iFil+4][2]=1.24116396451309802390e-02;
				
				iFil+=5;
				
				 
				//! Filter type butterworth : dtphy=0.03, dtmes=0.3, HighFrF=0.1, LowFrF=0.3, Attenuation=100.0, OscilPass=0.1
				PMFalpha[iFil+0][0]=1.97538958658917085387e+00;
				PMFalpha[iFil+0][1]=-9.81404537290183620613e-01;
				PMFbeta[iFil+0][0]=1.50373767525314701657e-03;
				PMFbeta[iFil+0][1]=3.00747535050629403314e-03;
				PMFbeta[iFil+0][2]=1.50373767525314701657e-03;
				PMFalpha[iFil+1][0]=1.94035623671806556345e+00;
				PMFalpha[iFil+1][1]=-9.46264512830026660595e-01;
				PMFbeta[iFil+1][0]=1.47706902799027059957e-03;
				PMFbeta[iFil+1][1]=2.95413805598054119914e-03;
				PMFbeta[iFil+1][2]=1.47706902799027059957e-03;
				PMFalpha[iFil+2][0]=1.90947346119909333595e+00;
				PMFalpha[iFil+2][1]=-9.15287700988539398850e-01;
				PMFbeta[iFil+2][0]=1.45355994736154781749e-03;
				PMFbeta[iFil+2][1]=2.90711989472309563498e-03;
				PMFbeta[iFil+2][2]=1.45355994736154781749e-03;
				PMFalpha[iFil+3][0]=1.88412361723844279382e+00;
				PMFalpha[iFil+3][1]=-8.89860668172213875593e-01;
				PMFbeta[iFil+3][0]=1.43426273344272880994e-03;
				PMFbeta[iFil+3][1]=2.86852546688545761988e-03;
				PMFbeta[iFil+3][2]=1.43426273344272880994e-03;
				PMFalpha[iFil+4][0]=1.86532971044030948526e+00;
				PMFalpha[iFil+4][1]=-8.71009534980037747331e-01;
				PMFbeta[iFil+4][0]=1.41995613493202865393e-03;
				PMFbeta[iFil+4][1]=2.83991226986405730787e-03;
				PMFbeta[iFil+4][2]=1.41995613493202865393e-03;
				PMFalpha[iFil+5][0]=1.85378682607320866005e+00;
				PMFalpha[iFil+5][1]=-8.59431503176255384702e-01;
				PMFbeta[iFil+5][0]=1.41116927576169634215e-03;
				PMFbeta[iFil+5][1]=2.82233855152339268429e-03;
				PMFbeta[iFil+5][2]=1.41116927576169634215e-03;
				PMFalpha[iFil+6][0]=9.24947826808538486887e-01;
				PMFalpha[iFil+6][1]=-0.00000000000000000000e+00;
				PMFbeta[iFil+6][0]=3.75260865957307912510e-02;
				PMFbeta[iFil+6][1]=3.75260865957307912510e-02;
				PMFbeta[iFil+6][2]=0.00000000000000000000e+00;
				
				iFil += 7;
				*/
				
				/*
				//! Filter type butterworth : dtphy=0.03, dtmes=0.3, HighFrF=0.1, LowFrF=0.4, Attenuation=200.0, OscilPass=0.1
				PMFalpha[iFil+0][0]=1.98292831533441438197e+00;
				PMFalpha[iFil+0][1]=-9.87834812904273884548e-01;
				PMFbeta[iFil+0][0]=1.22662439246492443491e-03;
				PMFbeta[iFil+0][1]=2.45324878492984886982e-03;
				PMFbeta[iFil+0][2]=1.22662439246492443491e-03;
				PMFalpha[iFil+1][0]=1.95945333823030898301e+00;
				PMFalpha[iFil+1][1]=-9.64301750030136961556e-01;
				PMFbeta[iFil+1][0]=1.21210294995700634466e-03;
				PMFbeta[iFil+1][1]=2.42420589991401268931e-03;
				PMFbeta[iFil+1][2]=1.21210294995700634466e-03;
				PMFalpha[iFil+2][0]=1.93756619553906950237e+00;
				PMFalpha[iFil+2][1]=-9.42360450458088916292e-01;
				PMFbeta[iFil+2][0]=1.19856372975483439769e-03;
				PMFbeta[iFil+2][1]=2.39712745950966879538e-03;
				PMFbeta[iFil+2][2]=1.19856372975483439769e-03;
				PMFalpha[iFil+3][0]=1.91782344840985374823e+00;
				PMFalpha[iFil+3][1]=-9.22568852475283085468e-01;
				PMFbeta[iFil+3][0]=1.18635101635733626081e-03;
				PMFbeta[iFil+3][1]=2.37270203271467252162e-03;
				PMFbeta[iFil+3][2]=1.18635101635733626081e-03;
				PMFalpha[iFil+4][0]=1.90069231145756201151e+00;
				PMFalpha[iFil+4][1]=-9.05395326758266061340e-01;
				PMFbeta[iFil+4][0]=1.17575382517607338918e-03;
				PMFbeta[iFil+4][1]=2.35150765035214677837e-03;
				PMFbeta[iFil+4][2]=1.17575382517607338918e-03;
				PMFalpha[iFil+5][0]=1.88655197098266680378e+00;
				PMFalpha[iFil+5][1]=-8.91219997854570733509e-01;
				PMFbeta[iFil+5][0]=1.16700671797601300723e-03;
				PMFbeta[iFil+5][1]=2.33401343595202601447e-03;
				PMFbeta[iFil+5][2]=1.16700671797601300723e-03;
				PMFalpha[iFil+6][0]=1.87569738335866587065e+00;
				PMFalpha[iFil+6][1]=-8.80338551968782412338e-01;
				PMFbeta[iFil+6][0]=1.16029215252915168625e-03;
				PMFbeta[iFil+6][1]=2.32058430505830337251e-03;
				PMFbeta[iFil+6][2]=1.16029215252915168625e-03;
				PMFalpha[iFil+7][0]=1.86834406966712540665e+00;
				PMFalpha[iFil+7][1]=-8.72967043461282621308e-01;
				PMFbeta[iFil+7][0]=1.15574344853933098573e-03;
				PMFbeta[iFil+7][1]=2.31148689707866197146e-03;
				PMFbeta[iFil+7][2]=1.15574344853933098573e-03;
				PMFalpha[iFil+8][0]=1.86463271749818049194e+00;
				PMFalpha[iFil+8][1]=-8.69246508035304565887e-01;
				PMFbeta[iFil+8][0]=1.15344763428106944318e-03;
				PMFbeta[iFil+8][1]=2.30689526856213888636e-03;
				PMFbeta[iFil+8][2]=1.15344763428106944318e-03;
				
				iFil += 9;
				 */
				/*
				//! Filter type butterworth : dtphy=0.03, dtmes=0.3, HighFrF=0.2, LowFrF=1.0, Attenuation=100.0, OscilPass=0.1
				PMFalpha[iFil+0][0]=1.90814003377194474353e+00;
				PMFalpha[iFil+0][1]=-9.39591852819607664671e-01;
				PMFbeta[iFil+0][0]=7.86295476191570946733e-03;
				PMFbeta[iFil+0][1]=1.57259095238314189347e-02;
				PMFbeta[iFil+0][2]=7.86295476191570946733e-03;
				PMFalpha[iFil+1][0]=1.80564257049575327407e+00;
				PMFalpha[iFil+1][1]=-8.35404926710106798815e-01;
				PMFbeta[iFil+1][0]=7.44058905358838552219e-03;
				PMFbeta[iFil+1][1]=1.48811781071767710444e-02;
				PMFbeta[iFil+1][2]=7.44058905358838552219e-03;
				PMFalpha[iFil+2][0]=1.72989131383837868583e+00;
				PMFalpha[iFil+2][1]=-7.58405064198417133703e-01;
				PMFbeta[iFil+2][0]=7.12843759000959115157e-03;
				PMFbeta[iFil+2][1]=1.42568751800191823031e-02;
				PMFbeta[iFil+2][2]=7.12843759000959115157e-03;
				PMFalpha[iFil+3][0]=1.68378496605175875800e+00;
				PMFalpha[iFil+3][1]=-7.11538746764985186033e-01;
				PMFbeta[iFil+3][0]=6.93844517830660787588e-03;
				PMFbeta[iFil+3][1]=1.38768903566132157518e-02;
				PMFbeta[iFil+3][2]=6.93844517830660787588e-03;
				PMFalpha[iFil+4][0]=8.34171061727226215154e-01;
				PMFalpha[iFil+4][1]=-0.00000000000000000000e+00;
				PMFbeta[iFil+4][0]=8.29144691363868924228e-02;
				PMFbeta[iFil+4][1]=8.29144691363868924228e-02;
				PMFbeta[iFil+4][2]=0.00000000000000000000e+00;
				
				iFil+=5;
				*/
				
				
				
				/*
				//! Filter type Chebyshev II : dtphy=0.03, dtmes=0.3, HighFrF=0.1, LowFrF=0.4, Attenuation=260.0, OscilPass=0.1
				PMFalpha[iFil+0][0]=1.98113205029021743897e+00;
				PMFalpha[iFil+0][1]=-9.86411929007989352058e-01;
				PMFbeta[iFil+0][0]=8.32348197910269677058e-02;
				PMFbeta[iFil+0][1]=-1.61189760864282133346e-01;
				PMFbeta[iFil+0][2]=8.32348197910269677058e-02;
				PMFalpha[iFil+1][0]=1.95480135414007660799e+00;
				PMFalpha[iFil+1][1]=-9.60043738953899716826e-01;
				PMFbeta[iFil+1][0]=7.65126858278178806350e-02;
				PMFbeta[iFil+1][1]=-1.47782986841812707945e-01;
				PMFbeta[iFil+1][2]=7.65126858278178806350e-02;
				PMFalpha[iFil+2][0]=1.92995173388151552984e+00;
				PMFalpha[iFil+2][1]=-9.35188170777667693301e-01;
				PMFbeta[iFil+2][0]=6.51100073796006678650e-02;
				PMFbeta[iFil+2][1]=-1.24983577863049102885e-01;
				PMFbeta[iFil+2][2]=6.51100073796006678650e-02;
				PMFalpha[iFil+3][0]=1.90720934369374184314e+00;
				PMFalpha[iFil+3][1]=-9.12464593252902300691e-01;
				PMFbeta[iFil+3][0]=5.05059129871344505980e-02;
				PMFbeta[iFil+3][1]=-9.57565764151083465006e-02;
				PMFbeta[iFil+3][2]=5.05059129871344505980e-02;
				PMFalpha[iFil+4][0]=1.88735064891626524997e+00;
				PMFalpha[iFil+4][1]=-8.92640265913224095762e-01;
				PMFbeta[iFil+4][0]=3.46706252483858085034e-02;
				PMFbeta[iFil+4][1]=-6.40516334998127712108e-02;
				PMFbeta[iFil+4][2]=3.46706252483858085034e-02;
				PMFalpha[iFil+5][0]=1.87127559020546030411e+00;
				PMFalpha[iFil+5][1]=-8.76604491840515587242e-01;
				PMFbeta[iFil+5][0]=1.98821495282552745953e-02;
				PMFbeta[iFil+5][1]=-3.44353974214552591149e-02;
				PMFbeta[iFil+5][2]=1.98821495282552745953e-02;
				PMFalpha[iFil+6][0]=1.85990110361867433397e+00;
				PMFalpha[iFil+6][1]=-8.65263614648361878423e-01;
				PMFbeta[iFil+6][0]=8.41923331261805810155e-03;
				PMFbeta[iFil+6][1]=-1.14759555955485023637e-02;
				PMFbeta[iFil+6][2]=8.41923331261805810155e-03;
				PMFalpha[iFil+7][0]=1.85399374602266298062e+00;
				PMFalpha[iFil+7][1]=-8.59375519284960143196e-01;
				PMFbeta[iFil+7][0]=2.15539803834323938359e-03;
				PMFbeta[iFil+7][1]=1.07097718561080935788e-03;
				PMFbeta[iFil+7][2]=2.15539803834323938359e-03;
				
				iFil += 8;
				
				*/
				
				//! Filter type Chebyshev II : dtphy=0.03, dtmes=0.3, HighFrF=0.1, LowFrF=0.4, Attenuation=300.0, OscilPass=0.1
				PMFalpha[iFil+0][0]=1.98349115190152391897e+00;
				PMFalpha[iFil+0][1]=-9.88355720609021548562e-01;
				PMFbeta[iFil+0][0]=7.68409151221110670038e-02;
				PMFbeta[iFil+0][1]=-1.48817261536724448900e-01;
				PMFbeta[iFil+0][2]=7.68409151221110670038e-02;
				PMFalpha[iFil+1][0]=1.96082831104362087515e+00;
				PMFalpha[iFil+1][1]=-9.65659344869815416956e-01;
				PMFbeta[iFil+1][0]=7.18165737066151221857e-02;
				PMFbeta[iFil+1][1]=-1.38802113587035702569e-01;
				PMFbeta[iFil+1][2]=7.18165737066151221857e-02;
				PMFalpha[iFil+2][0]=1.93924208933892749940e+00;
				PMFalpha[iFil+2][1]=-9.44061462768908565835e-01;
				PMFbeta[iFil+2][0]=6.32164940322232360437e-02;
				PMFbeta[iFil+2][1]=-1.21613614634465419528e-01;
				PMFbeta[iFil+2][2]=6.32164940322232360437e-02;
				PMFalpha[iFil+3][0]=1.91915080688949069021e+00;
				PMFalpha[iFil+3][1]=-9.23976756950907418542e-01;
				PMFbeta[iFil+3][0]=5.19339516319045554904e-02;
				PMFbeta[iFil+3][1]=-9.90419532023923410158e-02;
				PMFbeta[iFil+3][2]=5.19339516319045554904e-02;
				PMFalpha[iFil+4][0]=1.90104875053455901757e+00;
				PMFalpha[iFil+4][1]=-9.05894465838598672569e-01;
				PMFbeta[iFil+4][0]=3.91655717866221317336e-02;
				PMFbeta[iFil+4][1]=-7.34854282692046778536e-02;
				PMFbeta[iFil+4][2]=3.91655717866221317336e-02;
				PMFalpha[iFil+5][0]=1.88550454980879966804e+00;
				PMFalpha[iFil+5][1]=-8.90377187858027774325e-01;
				PMFbeta[iFil+5][0]=2.63299711349363023405e-02;
				PMFbeta[iFil+5][1]=-4.77873042206444428870e-02;
				PMFbeta[iFil+5][2]=2.63299711349363023405e-02;
				PMFalpha[iFil+6][0]=1.87312448918745699800e+00;
				PMFalpha[iFil+6][1]=-8.78024799276237000178e-01;
				PMFbeta[iFil+6][0]=1.49355168593089168011e-02;
				PMFbeta[iFil+6][1]=-2.49707236298379806061e-02;
				PMFbeta[iFil+6][2]=1.49355168593089168011e-02;
				PMFalpha[iFil+7][0]=1.86448333991856052272e+00;
				PMFalpha[iFil+7][1]=-8.69406080077212206625e-01;
				PMFbeta[iFil+7][0]=6.39640267854032357198e-03;
				PMFbeta[iFil+7][1]=-7.87006519842904997408e-03;
				PMFbeta[iFil+7][2]=6.39640267854032357198e-03;
				PMFalpha[iFil+8][0]=1.86003633121763423119e+00;
				PMFalpha[iFil+8][1]=-8.64971576082619875159e-01;
				PMFbeta[iFil+8][0]=1.82107336778795426448e-03;
				PMFbeta[iFil+8][1]=1.29309812940975777003e-03;
				PMFbeta[iFil+8][2]=1.82107336778795426448e-03;
				
				iFil += 9;
				
				
				//! Filter type butterworth : dtphy=0.03, dtmes=0.3, HighFrF=0.1, LowFrF=0.4, Attenuation=260.0, OscilPass=0.1
				//! Best
				PMFalpha[iFil+0][0]=1.98597059487818050627e+00;
				PMFalpha[iFil+0][1]=-9.90670513451891054935e-01;
				PMFbeta[iFil+0][0]=1.17497964342764085223e-03;
				PMFbeta[iFil+0][1]=2.34995928685528170446e-03;
				PMFbeta[iFil+0][2]=1.17497964342764085223e-03;
				PMFalpha[iFil+1][0]=1.96778317936738189609e+00;
				PMFalpha[iFil+1][1]=-9.72440056331016711155e-01;
				PMFbeta[iFil+1][0]=1.16421924090866039936e-03;
				PMFbeta[iFil+1][1]=2.32843848181732079872e-03;
				PMFbeta[iFil+1][2]=1.16421924090866039936e-03;
				PMFalpha[iFil+2][0]=1.95042197941471462386e+00;
				PMFalpha[iFil+2][1]=-9.55037770056965995558e-01;
				PMFbeta[iFil+2][0]=1.15394766056277310248e-03;
				PMFbeta[iFil+2][1]=2.30789532112554620497e-03;
				PMFbeta[iFil+2][2]=1.15394766056277310248e-03;
				PMFalpha[iFil+3][0]=1.93416735430789699102e+00;
				PMFalpha[iFil+3][1]=-9.38744677404533489629e-01;
				PMFbeta[iFil+3][0]=1.14433077415905960078e-03;
				PMFbeta[iFil+3][1]=2.28866154831811920156e-03;
				PMFbeta[iFil+3][2]=1.14433077415905960078e-03;
				PMFalpha[iFil+4][0]=1.91926737852683437779e+00;
				PMFalpha[iFil+4][1]=-9.23809439936766008117e-01;
				PMFbeta[iFil+4][0]=1.13551535248289543756e-03;
				PMFbeta[iFil+4][1]=2.27103070496579087512e-03;
				PMFbeta[iFil+4][2]=1.13551535248289543756e-03;
				PMFalpha[iFil+5][0]=1.90593740896341690139e+00;
				PMFalpha[iFil+5][1]=-9.10447924200560532171e-01;
				PMFbeta[iFil+5][0]=1.12762880928587191746e-03;
				PMFbeta[iFil+5][1]=2.25525761857174383493e-03;
				PMFbeta[iFil+5][2]=1.12762880928587191746e-03;
				PMFalpha[iFil+6][0]=1.89436049879259571505e+00;
				PMFalpha[iFil+6][1]=-8.98843616577189918893e-01;
				PMFbeta[iFil+6][0]=1.12077944614847875397e-03;
				PMFbeta[iFil+6][1]=2.24155889229695750794e-03;
				PMFbeta[iFil+6][2]=1.12077944614847875397e-03;
				PMFalpha[iFil+7][0]=1.88468836622633584277e+00;
				PMFalpha[iFil+7][1]=-8.89148594328870611747e-01;
				PMFbeta[iFil+7][0]=1.11505702563369463003e-03;
				PMFbeta[iFil+7][1]=2.23011405126738926005e-03;
				PMFbeta[iFil+7][2]=1.11505702563369463003e-03;
				PMFalpha[iFil+8][0]=1.87704265569559192173e+00;
				PMFalpha[iFil+8][1]=-8.81484789765402787509e-01;
				PMFbeta[iFil+8][0]=1.11053351745269085829e-03;
				PMFbeta[iFil+8][1]=2.22106703490538171658e-03;
				PMFbeta[iFil+8][2]=1.11053351745269085829e-03;
				PMFalpha[iFil+9][0]=1.87151626753484778831e+00;
				PMFalpha[iFil+9][1]=-8.75945323075485560160e-01;
				PMFbeta[iFil+9][0]=1.10726388515939287195e-03;
				PMFbeta[iFil+9][1]=2.21452777031878574390e-03;
				PMFbeta[iFil+9][2]=1.10726388515939287195e-03;
				PMFalpha[iFil+10][0]=1.86817457473941450630e+00;
				PMFalpha[iFil+10][1]=-8.72595721963545734035e-01;
				PMFbeta[iFil+10][0]=1.10528680603276877048e-03;
				PMFbeta[iFil+10][1]=2.21057361206553754096e-03;
				PMFbeta[iFil+10][2]=1.10528680603276877048e-03;
				PMFalpha[iFil+11][0]=9.33528194180002102165e-01;
				PMFalpha[iFil+11][1]=-0.00000000000000000000e+00;
				PMFbeta[iFil+11][0]=3.32359029099989281009e-02;
				PMFbeta[iFil+11][1]=3.32359029099989281009e-02;
				PMFbeta[iFil+11][2]=0.00000000000000000000e+00;
				
				iFil += 12;
				
			}
			
			
			
			if(MT->wcmp(PMFilterType,"Butter007503")){
				ReqPMFilter = true;
				//! Filter type butterworth : dtphy=0.075000, dtmes=0.300000, HighFrF=0.100000, LowFrF=0.500000, Attenuation=220.000000, OscilPass=0.100000
				PMFalpha.resize(9);
				PMFbeta.resize(9);
				for(int i=0; i<9; i++){
					PMFalpha[i].resize(2);
					PMFbeta[i].resize(3);
				}
				PMFalpha[0][0]=1.9381305653;
				PMFalpha[0][1]=-0.9683158828;
				PMFbeta[0][0]=0.0075463294;
				PMFbeta[0][1]=0.0150926588;
				PMFbeta[0][2]=0.0075463294;
				PMFalpha[1][0]=1.8795913781;
				PMFalpha[1][1]=-0.9088649800;
				PMFbeta[1][0]=0.0073184005;
				PMFbeta[1][1]=0.0146368010;
				PMFbeta[1][2]=0.0073184005;
				PMFalpha[2][0]=1.8272370226;
				PMFalpha[2][1]=-0.8556952342;
				PMFbeta[2][0]=0.0071145529;
				PMFbeta[2][1]=0.0142291058;
				PMFbeta[2][2]=0.0071145529;
				PMFalpha[3][0]=1.7819799831;
				PMFalpha[3][1]=-0.8097333412;
				PMFbeta[3][0]=0.0069383395;
				PMFbeta[3][1]=0.0138766791;
				PMFbeta[3][2]=0.0069383395;
				PMFalpha[4][0]=1.7444253545;
				PMFalpha[4][1]=-0.7715938200;
				PMFbeta[4][0]=0.0067921164;
				PMFbeta[4][1]=0.0135842327;
				PMFbeta[4][2]=0.0067921164;
				PMFalpha[5][0]=1.7149528366;
				PMFalpha[5][1]=-0.7416622838;
				PMFbeta[5][0]=0.0066773618;
				PMFbeta[5][1]=0.0133547236;
				PMFbeta[5][2]=0.0066773618;
				PMFalpha[6][0]=1.6937863201;
				PMFalpha[6][1]=-0.7201661104;
				PMFbeta[6][0]=0.0065949476;
				PMFbeta[6][1]=0.0131898952;
				PMFbeta[6][2]=0.0065949476;
				PMFalpha[7][0]=1.6810479643;
				PMFalpha[7][1]=-0.7072293618;
				PMFbeta[7][0]=0.0065453494;
				PMFbeta[7][1]=0.0130906987;
				PMFbeta[7][2]=0.0065453494;
				PMFalpha[8][0]=0.8383980874;
				PMFalpha[8][1]=-0.0000000000;
				PMFbeta[8][0]=0.0808009563;
				PMFbeta[8][1]=0.0808009563;
				PMFbeta[8][2]=0.0000000000;
			}
			if(MT->wcmp(PMFilterType,"Butter01503")){
				ReqPMFilter = true;
				//! Filter type butterworth : dtphy=0.150000, dtmes=0.300000, HighFrF=0.100000, LowFrF=0.400000, Attenuation=240.000000, OscilPass=0.100000
				PMFalpha.resize(10);
				PMFbeta.resize(10);
				for(int i=0; i<10; i++){
					PMFalpha[i].resize(2);
					PMFbeta[i].resize(3);
				}
				PMFalpha[0][0]=1.8338657265;
				PMFalpha[0][1]=-0.9483674610;
				PMFbeta[0][0]=0.0286254336;
				PMFbeta[0][1]=0.0572508672;
				PMFbeta[0][2]=0.0286254336;
				PMFalpha[1][0]=1.7448822180;
				PMFalpha[1][1]=-0.8538280571;
				PMFbeta[1][0]=0.0272364598;
				PMFbeta[1][1]=0.0544729196;
				PMFbeta[1][2]=0.0272364598;
				PMFalpha[2][0]=1.6669956004;
				PMFalpha[2][1]=-0.7710784048;
				PMFbeta[2][0]=0.0260207011;
				PMFbeta[2][1]=0.0520414022;
				PMFbeta[2][2]=0.0260207011;
				PMFalpha[3][0]=1.6000821292;
				PMFalpha[3][1]=-0.6999870331;
				PMFbeta[3][0]=0.0249762260;
				PMFbeta[3][1]=0.0499524519;
				PMFbeta[3][2]=0.0249762260;
				PMFalpha[4][0]=1.5438154543;
				PMFalpha[4][1]=-0.6402072155;
				PMFbeta[4][0]=0.0240979403;
				PMFbeta[4][1]=0.0481958806;
				PMFbeta[4][2]=0.0240979403;
				PMFalpha[5][0]=1.4977808301;
				PMFalpha[5][1]=-0.5912983110;
				PMFbeta[5][0]=0.0233793702;
				PMFbeta[5][1]=0.0467587404;
				PMFbeta[5][2]=0.0233793702;
				PMFalpha[6][0]=1.4615533950;
				PMFalpha[6][1]=-0.5528089305;
				PMFbeta[6][0]=0.0228138839;
				PMFbeta[6][1]=0.0456277677;
				PMFbeta[6][2]=0.0228138839;
				PMFalpha[7][0]=1.4347497683;
				PMFalpha[7][1]=-0.5243317560;
				PMFbeta[7][0]=0.0223954969;
				PMFbeta[7][1]=0.0447909939;
				PMFbeta[7][2]=0.0223954969;
				PMFalpha[8][0]=1.4170607816;
				PMFalpha[8][1]=-0.5055383158;
				PMFbeta[8][0]=0.0221193835;
				PMFbeta[8][1]=0.0442387671;
				PMFbeta[8][2]=0.0221193835;
				PMFalpha[9][0]=1.4082713013;
				PMFalpha[9][1]=-0.4962000435;
				PMFbeta[9][0]=0.0219821855;
				PMFbeta[9][1]=0.0439643711;
				PMFbeta[9][2]=0.0219821855;
				
			}
			
			
			
		}
		
		//! *** Individual parameters
		if(MT->wcmp(ezxml_attr(param,"Name"),"PMFAttenuation"))
			PMFAttenuation = atof((*param).txt);
		if(MT->wcmp(ezxml_attr(param,"Name"),"PMFOscillation"))
			PMFOscillation = atof((*param).txt);
		if(MT->wcmp(ezxml_attr(param,"Name"),"PMFLowFreqFact"))
			PMFLowFreqFact = atof((*param).txt);
		if(MT->wcmp(ezxml_attr(param,"Name"),"PMFHighFreqFact"))
			PMFHighFreqFact = atof((*param).txt);
	}
	
	//! *** Create the filter 
	if(ReqPMFilter)
		PMFilter = new LCFilter(MT);
	
	
	
	if(PMFilterType != NULL)
		MT->Free(PMFilterType, (strlen(PMFilterType)+1) * sizeof(char) );
	
}


// ***************************************
// * Linking and initialization methods  *
// ***************************************


void LCOBPM::LinkClock(LCUSO ** Clocks, int NClocks)
{
	for (int i=0; i<NClocks; i++) {
		if(Clocks[i] != NULL)
			if(Clocks[i]->getiSC()==iSC)
				clockl = Clocks[i];
	}
}



void LCOBPM::initBase()
{
	int NDataPho;
	
	//! *** Allocation of memory for the last measurement 
	if(LastMes == NULL)
		LastMes = (double*)MT->AllocMemory(sizeof(double));
	
	
	//! *** Compute some number of data
	N2Add = MT->ifloor(dtMes/dtPhy);
	NDataPho = MT->iceil(N2Add * NDataPha);
	
	//! *** Delete previous raw data and filter
	if(MesPho != NULL)
		delete MesPho;
	if(MesPha != NULL)
		delete MesPha;
	if(MesFilter != NULL)
		delete MesFilter;
		
	
	
	//! *** Create photodiode measurements (raw) data
	MesPho = new LCSerie2(MT, 0., dtPhy, NDataPho );
	
	//! *** Create photodiode measurements (raw) data
	MesPha = new LCSerie2(MT, 0., dtMes, NDataPha );
	
	//! *** Create the filter
	if(PMFilter != NULL){
		
		//! *** Create filtered photodiode data which are equivalent to the phasemeter measurements except the time step
		MesFilter = new LCSerie2(MT, 0., dtPhy, NDataPho);
		
		//! *** Linked the filter with the local series
		PMFilter->setInOutData(MesPho, MesFilter);
		
		//! *** Initialize the filter
		if(PMFalpha.size()==0)
			PMFilter->init(PMFAttenuation, PMFOscillation, PMFLowFreqFact/dtMes, PMFHighFreqFact/dtMes);
		else
			PMFilter->init(PMFalpha, PMFbeta, 10000);
	}
	
	
	//! *** Fill the input and output buffer by temporarly switching off the filter
	LCFilter * SaveFilter;
	SaveFilter = PMFilter;
	PMFilter = NULL;
	for(int i=0; i<MesPho->getNmax(); i++)
		MeasurePha(t0);
	PMFilter = SaveFilter;
	
	//! *** Stabilization of the filter
	if(PMFilter != NULL)
		for(int i=0; i<PMFilter->getNbDataStab(); i++)
			MeasurePha(t0);
	
		
}


void LCOBPM::LinkLastMesOutput(double * RecLoc)
{
	if((LastMes != NULL)&&(LastMesLocalAlloc))
		MT->Free(LastMes, sizeof(double));
	LastMesLocalAlloc = false;
	LastMes = RecLoc;
}


// ********************
// *  Access methods  *
// ********************

void LCOBPM::setNamePos(const char * NamePos)
{
	int ipos(strlen(NamePos)-1);
	char TmpiSC[2];
	
	//! *** Start from the end : detect if the direction of optical bench
	if(NamePos[ipos] == 's'){
		IndirectDir = 1;
		ipos--;
	}else{
		if(isdigit(NamePos[ipos])&&(!(isdigit(NamePos[ipos-1])))){
			IndirectDir = 0;
		}else{
			//! ** There are 2 digits at the end of the name , so it's not a LISA like detector (not 3 SCs)
			TmpiSC[0] = NamePos[ipos];
			TmpiSC[1] = '\0';
			IndirectDir = -1*atoi(TmpiSC);
			ipos--;
		}
	}
	
	//! *** Then : Spacecraft index
	TmpiSC[0] = NamePos[ipos];
	TmpiSC[1] = '\0';
	iSC = atoi(TmpiSC);
	//if((iSC<1)||(iSC>3))
	//	throw std::invalid_argument("ERROR in LCNoise::setNamePos : Index of spacecraft should be in 1,2 or 3 !" );
	
	
	int NChar(MIN(strlen(NamePos), 8));
	for(int i=0; i<NChar; i++)
		Name[i] = NamePos[i];
	Name[NChar] = '\0';
	
}


void LCOBPM::setTimeInfo( double t0_n, double dtMes_n, double dtPhy_n, double tFirst_n, double tLast_n, int NDataPha_n)
{
	t0 = t0_n;
	dtMes = dtMes_n;
	dtPhy = dtPhy_n;
	tFirst = tFirst_n;
	tLast = tLast_n;
	NDataPha = NDataPha_n;
}


bool LCOBPM::MatchName(const char * Cmpname)
{
	return(MT->wcmp(Name, Cmpname));
}


bool LCOBPM::CmpPosName(int iSC_cmp, int IndirectDir_cmp, const char * Name_cmp, int nchar)
{
	bool res(true);
	
	if(iSC != iSC_cmp){
		res = false;
	}else{
		if(IndirectDir != IndirectDir_cmp){
			res = false;
		}else{
			for(int i=0; i<nchar; i++)
				if(Name[i] != Name_cmp[i])
					res = false;
		}
	}
	return res;
}


// *******************
// *  Other methods  *
// *******************

double LCOBPM::sGW(double t)
{
	if(ArmGW == NULL)
		return(0.);
	else
		return(ArmGW->gS(iSCLink, iSC, t));
}


double LCOBPM::gN(int iN, double tDelay)
{
	if(pn[iN] == NULL)
		return( 0. );
	else
		return( pn[iN]->gN(tDelay) );
}


double LCOBPM::gFSN(double trec)
{
	if(FactShotNoise != NULL){
		return((*FactShotNoise));
	}else{
		return(1.);
	}
}


void LCOBPM::DispInfoBase(char * BTab, bool LinkedModule)
{
	if(MT->Disp()){
		char BTab2[1024];
		Cout << BTab << "\t- position : spacecraft = " << iSC << "  , direction = " << IndirectDir << Endl;
		
		Cout << BTab << "\t- position of the linked optical bench : spacecraft = " << iSCLink << "  , direction = " << IndirectDirLink << Endl;
		
		if(PMFilter != NULL){
			Cout << BTab << "\t- Filter :" << Endl;
			if(PMFalpha.size()==0){
				Cout << BTab << "\t\t+ Attenuation = " << PMFAttenuation << " dB" << Endl;
				Cout << BTab << "\t\t+ Oscilation  = " << PMFOscillation << " dB" << Endl;
				Cout << BTab << "\t\t+ Low transition frequency  = " << PMFLowFreqFact/dtMes << " Hz" << Endl;
				Cout << BTab << "\t\t+ High transition frequency = " << PMFHighFreqFact/dtMes << " Hz" << Endl;
			}
			strcpy(BTab2,BTab);
			strcat(BTab2, "\t\t");
			PMFilter->DispInfo(BTab2);
		}
		
		
		if(LinkedModule){
			strcpy(BTab2,BTab);
			strcat(BTab2, ">\t");
			if(clockl != NULL)
				clockl->DispInfo(BTab2);
			for(int i=0; i<Npn; i++)
				if(pn[i] != NULL)
					pn[i]->DispInfo(BTab2);
			
		}
	}
}



// end of LISACODE-OBPM.cpp

