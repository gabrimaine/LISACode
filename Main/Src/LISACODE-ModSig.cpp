/*
 *  LISACODE-ModSig.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 28/05/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */




#include "LISACODE-ModSig.h"

// *********************
// ***  Constructor  ***
// *********************	

LCModSig::LCModSig()
{
	MT = new LCTools;
	MT->LocTools = true;
	initNULL(false);
}


LCModSig::LCModSig(LCTools * MT_n)
{
	MT = MT_n;
	initNULL(false);
	
}


LCModSig::~LCModSig()
{
	initNULL(true);
}


void LCModSig::initNULL(bool CleanMem)
{
	if(CleanMem){
		
		//CleanNDatMem();
		CleanBaseMem();
		
		if(Err4 != NULL)
			MT->Free(Err4, Np*sizeof(double));
		if(Err6 != NULL)
			MT->Free(Err6, Np*sizeof(double));
		if(p != NULL)
			MT->Free(p, Np*sizeof(double));
		if(pS != NULL)
			MT->Free(pS, Np*sizeof(bool));
		
		if(Sim != NULL)
			delete Sim;
		
		if(TDIName != NULL){
			for(int i=0; i<NTDI; i++){
				if(TDIName[i] != NULL)
					MT->Free(TDIName[i], NCharTDIName*sizeof(char));
				TDIName[i] = NULL;
			}
			MT->Free(TDIName, NTDI*sizeof(char*));
			TDIName = NULL;
		}
		
		if(iParStudy != NULL)
			MT->Free(iParStudy, NParStudy*sizeof(int));
		if(lnParStudy != NULL)
			MT->Free(lnParStudy, NParStudy*sizeof(bool));
		if(ParStuDelta != NULL)
			MT->Free(ParStuDelta, NParStudy*sizeof(double));
		
		if(fNCfg != NULL){
			for(int iF=0; iF<NfNCfg; iF++)
				if(fNCfg[iF] != NULL)
					MT->Free(fNCfg[iF], NCharfNCfg*sizeof(char));
			MT->Free(fNCfg, NfNCfg*sizeof(char*));
		}
		
		if(GWmod != NULL)
			delete GWmod;
	}
	
	SNR = NULL;
	Err4 = NULL;
	Err6 = NULL;
	
	p = NULL;
	pS = NULL;
	Sim = NULL;
	SrcType = UNKNOWNGW;
	TDIName = NULL;
	NTDI = 0;
	
	tSig = NULL;
	tSig2 = NULL;
    tSig3 = NULL; // added by Sofiane
    tSig4 = NULL; // added by Sofiane
	NtData = 0;
	tEndSecur = 1.e5;
	
	fSig = NULL;
	fDSig = NULL;
	NfData = 0;
	NfDataLoc = 0;
	df = 1.e-8;
	iFmin = 0;
	iFmax = 1;
	
	iParStudy = NULL;
	NParStudy = 0;
	lnParStudy = NULL;
	ParStuDelta = NULL;
	
	NoBaseMem = true;
	NoIni = true;
	NotDat = true;
	NofDat = true;
	NoSNR = true;
	NoFIM = true;
	NoErr = true;
	
	strcpy(DirSwapMem,"./");
	MaxMemRAM = 16.;
	SwapMem = false;
	
	NoiseMinLowf = NULL;
	NCharfNCfg = 2048;
	NCharTDIName = 16;
	
	Detector = UNKNOWNDC;
	GalacticNoise = 1;
    dtMesMin = 4.; 
	tDurMax = 67108864.;
	t0 = 0.;
    FactStepD = 1.;
	fNCfg = NULL;
	NfNCfg = 0;
	CheckOut = false;
	DispDetails = false;
	strcpy(BaseNameOut,"DefOut");
	strcpy(CurrentNameOut,BaseNameOut);
	DeltaParSpecInd = -1;
	DeltaParSpecVal = 0.;
	NotValidPar = false;
	
	FIMnonan = true;
	
	GWmod = NULL;
	GWSpeParInd = NULL;
	GWSpeParVal = NULL;
	NGWSpePar = 0;
	
	//! ** Values with dependance
	dtMes = dtMesMin;
	tDur = tDurMax;

	
	
	
}


// ***************************
// *  Configuration methods  *
// ***************************
void LCModSig::addfCfg(char * fNCfgNew)
{
	NfNCfg++;
	fNCfg = (char**)MT->ReAllocMemory(fNCfg, (NfNCfg-1)*sizeof(char*), NfNCfg*sizeof(char*));
	fNCfg[NfNCfg-1] = (char*)MT->AllocMemory(NCharfNCfg*sizeof(char));
	strcpy(fNCfg[NfNCfg-1], fNCfgNew);
}



// ********************
// *  Access methods  *
// ********************
void LCModSig::sP(int iP, double Val)
{
	AllocBaseMem();
	
	if(iP>=Np){
		Cout << "ERROR in LCModSig::sP : index " << iP << " is not in the range [0," << Np << "]" << Endl; 
		throw std::invalid_argument("ERROR in LCModSig::sP : The index of parameters should be in the range [0,Nb param] !");
	}
	p[iP] = Val;
	pS[iP] = true;
	
	NoIni = true;
	NotDat = true;
	NofDat = true;
	NoSNR = true;
	NoFIM = true;
	NoErr = true;
}

void LCModSig::AddGWSpePar(int iSP, double Value)
{
	NGWSpePar++;
	GWSpeParInd = (int*) MT->ReAllocMemory(GWSpeParInd, (NGWSpePar-1)*sizeof(int), NGWSpePar*sizeof(int));
	GWSpeParVal = (double*) MT->ReAllocMemory(GWSpeParInd, (NGWSpePar-1)*sizeof(double), NGWSpePar*sizeof(double));
	GWSpeParInd[NGWSpePar-1] = iSP;
	GWSpeParVal[NGWSpePar-1] = Value; 
}


void LCModSig::sSrcType(GWSrcType SrcType_n)
{
	//! *** Clean 
	if(Err4 != NULL)
		MT->Free(Err4, Np*sizeof(double));
	if(Err6 != NULL)
		MT->Free(Err6, Np*sizeof(double));
	if(p != NULL)
		MT->Free(p, Np*sizeof(double));
	if(pS != NULL)
		MT->Free(pS, Np*sizeof(bool));
	
	SrcType = SrcType_n;
	switch (SrcType) {
		case SPINBBHHHARM:
			Np = 23;   // changed by Sofiane  19
			break;
		case GALBIN:
			Np = 13;
			break;
		case NRWAVE:
			Np = 22;
			break;
		default:
			Cout << "ERROR in LCModSig::sSrcType : Unknown type of source." << Endl;
			throw std::invalid_argument("ERROR in LCModSig::sSrcType : Unknown type of source.");
			break;
	}
	
	Err4 = (double*) MT->AllocMemory(Np*sizeof(double));
	Err6 = (double*) MT->AllocMemory(Np*sizeof(double));
	p  = (double*) MT->AllocMemory(Np*sizeof(double));
	pS = (bool*) MT->AllocMemory(Np*sizeof(bool));
    for (int i=0; i<Np; i++) {
        Err4[i] = 0.;
        Err6[i] = 0.;
        p[i]  = 0.;
        pS[i] = 0.;
    }
	
	Reset(true);
}


void LCModSig::sTDI(char * TDIAllName)
{
	//! *** Clean previous allocation
	if(TDIName != NULL){
		for(int i=0; i<NTDI; i++){
			if(TDIName[i] != NULL)
				MT->Free(TDIName[i], NCharTDIName*sizeof(char));
			TDIName[i] = NULL;
		}
		MT->Free(TDIName, NTDI*sizeof(char*));
	}
	
	//! *** Decode the char string and set the new TDI
	int ipTDINameAll(0);
	MT->wextractcoma(TDIAllName, ipTDINameAll, TDIName, NTDI, NCharTDIName);
	
	NoBaseMem = true;
}


void LCModSig::sParStudy(char * CharIndPar)
{
	//! *** Clean memory previously allocated
	if(iParStudy != NULL)
		MT->Free(iParStudy, NParStudy*sizeof(int));
	if(lnParStudy != NULL)
		MT->Free(lnParStudy, NParStudy*sizeof(bool));
	if(ParStuDelta != NULL)
		MT->Free(ParStuDelta, NParStudy*sizeof(double));
	
	//! *** Decode the char string
	char ** tmpWordInd(NULL);
	int iPtmp(0);
	MT->wextractcoma(CharIndPar, iPtmp, tmpWordInd, NParStudy, 10);
	
	//! *** Allocate and set the varaibles linked to the parameters to study
	iParStudy = (int*) MT->AllocMemory(NParStudy*sizeof(int));
	lnParStudy = (bool*) MT->AllocMemory(NParStudy*sizeof(bool));
	ParStuDelta = (double*) MT->AllocMemory(NParStudy*sizeof(double));
	for(int i=0; i<NParStudy; i++){
		if(tmpWordInd[i][0] == 'L'){
			iParStudy[i] = atoi(tmpWordInd[i]+1);
			lnParStudy[i] = true;
		}else{
			iParStudy[i] = atoi(tmpWordInd[i]);
			lnParStudy[i] = false;
		}
		ParStuDelta[i] = 0.;
	}
	for(int i=0; i<NParStudy; i++)
		MT->Free(tmpWordInd[i], 10*sizeof(char));
	MT->Free(tmpWordInd, NParStudy*sizeof(char*));
	NoBaseMem = true;
}


void LCModSig::sSwap(const char * DirSwapMem_n, double MaxMemRAM_n)
{
	MaxMemRAM = MaxMemRAM_n;
	strcpy(DirSwapMem, DirSwapMem_n);
}



void LCModSig::gPStudy(double * & pStudy)
{
	if(pStudy == NULL)
		pStudy = (double*) MT->AllocMemory(NParStudy*sizeof(double));
	for(int iP=0; iP<NParStudy; iP++)
		pStudy[iP] = p[iParStudy[iP]];
}


void LCModSig::gfSig(dcomplex ** &fSigCpy, int & NfDataLocCpy, int & NTDICpy, int & iFminCpy, double & dfCpy)
{
	if(fSigCpy!= NULL){
		for(int i=0; i<NTDICpy; i++)
			if(fSigCpy[i] != NULL)
				MT->Free(fSigCpy[i],NfDataLoc*sizeof(dcomplex));
		MT->Free(fSigCpy,NTDI*sizeof(dcomplex*));
	}
	NTDICpy = NTDI;
	NfDataLocCpy = NfDataLoc;
	iFminCpy = iFmin;
	dfCpy = df;
	fSigCpy = (dcomplex**) MT->AllocMemory(NTDICpy*sizeof(dcomplex*));
	for(int iS=0; iS<NTDI; iS++){
		fSigCpy[iS] = (dcomplex*) MT->AllocMemory(NfDataLoc*sizeof(dcomplex));
		memcpy(fSigCpy[iS], fSig[iS]+iFmin, NfDataLoc*sizeof(dcomplex));
	}
	
}


double LCModSig::gSNR(int iTDI)
{
	if(NotValidPar)
		return(0.);
	
	if(NoSNR)
		ComputeSNR();
	return(SNR[iTDI]);
}


double LCModSig::gSNR4()
{	
	if(NotValidPar)
		return(0.);
	
	double res(-1.);
	if(NoSNR)
		ComputeSNR();
	for(int i=0; i<NTDI; i++)
		if((MT->wcmp(TDIName[i],"X")))
			res = SNR[i];
	if(res<0){
		Cout << "ERROR in LCModSig::gSNR4 : Don't find the TDI X on which one the 4 links are based !" << Endl;
		throw std::invalid_argument("ERROR in LCModSig::gSNR4 : Don't find the TDI X on which one the 4 links are based !");
	}
	return(res);
}


double LCModSig::gSNR6()
{
	if(NotValidPar)
		return(0.);
	
	double resA(-1.), resE(-1.);
	if(NoSNR)
		ComputeSNR();
	for(int i=0; i<NTDI; i++){
		if((MT->wcmp(TDIName[i],"Am")))
			resA = SNR[i];
		if((MT->wcmp(TDIName[i],"Em")))
			resE = SNR[i];
	}
	if(resA<0){
		Cout << "ERROR in LCModSig::gSNR6 : Don't find the TDI Am need for the 6 links results !" << Endl;
		throw std::invalid_argument("ERROR in LCModSig::gSNR6 : Don't find the TDI A need for the 6 links results !");
	}
	if(resE<0){
		Cout << "ERROR in LCModSig::gSNR6 : Don't find the TDI Em need for the 6 links results !" << Endl;
		throw std::invalid_argument("ERROR in LCModSig::gSNR6 : Don't find the TDI E need for the 6 links results !");
	}
	
	return(sqrt( resA*resA + resE*resE ));
}


LCMatrix LCModSig::gFIM4()
{
	if(NoFIM)
		ComputeFIM();
	return(FIM4);
}


LCMatrix LCModSig::gFIM6()
{
	if(NoFIM)
		ComputeFIM();
	return(FIM6);
}


LCMatrix LCModSig::gCovM4()
{
	if(NoErr)
		ComputeErr();
	return(CovM4);
}


LCMatrix LCModSig::gCovM6()
{
	if(NoErr)
		ComputeErr();
	return(CovM6);	
}


double LCModSig::gErr4(int iP)
{
	if(NotValidPar)
		return(0.);
	if(NoErr)
		ComputeErr();
	return(Err4[iP]);
}


double LCModSig::gErr6(int iP)
{
	if(NotValidPar)
		return(0.);
	if(NoErr)
		ComputeErr();
	return(Err6[iP]);
}

double LCModSig::gErrSky(int Nlinks, double & ErrMajAxis, double & ErrMinAxis)
{
	double cB(0.);
	int iErrBet(-1), iErrLam(-1);
	double ErrArea(0.);
	double Cbb, Cll, Cbl;
	
	if(NotValidPar)
		return(0.);
	
	for(int iP=0; iP<NParStudy; iP++){
		if((SrcType==SPINBBHHHARM)||(SrcType==GALBIN)||(SrcType==NRWAVE)){
			cB = cos(p[0]);
			if(iParStudy[iP]==0)
				iErrBet = iParStudy[iP];
			if(iParStudy[iP]==1)
				iErrLam = iParStudy[iP];
		}
	}
	if((iErrBet<0)||(iErrLam<0))
		return(0.);
	
	if(NoErr)
		ComputeErr();
	
	if(Nlinks==4){
		Cbb = CovM4(iErrBet,iErrBet);
		Cll = CovM4(iErrLam,iErrLam);
		Cbl = CovM4(iErrBet,iErrLam);
	}
	if(Nlinks==6){
		Cbb = CovM6(iErrBet,iErrBet);
		Cll = CovM6(iErrLam,iErrLam);
		Cbl = CovM6(iErrBet,iErrLam);
	}
	
	//! * Denormalize
	Cbl *= sqrt(Cbb*Cll);
	
	ErrArea = 2.*M_PI*sqrt( (cB*cB*Cbb*Cll) - (cB*Cbl)*(cB*Cbl) ) ;
	ErrMajAxis =  2.*sqrt( cB*cB*Cbb + Cll + sqrt( (cB*cB*Cbb-Cll)*(cB*cB*Cbb-Cll) + 4*cB*cB*Cbl*Cbl ) ) ;
	ErrMinAxis =  2.*sqrt( cB*cB*Cbb + Cll - sqrt( (cB*cB*Cbb-Cll)*(cB*cB*Cbb-Cll) + 4*cB*cB*Cbl*Cbl ) ) ;
	
	return(ErrArea);
}


double LCModSig::gErr4Sky() 
{
	double ErrArea, ErrMajAxis, ErrMinAxis;
	ErrArea = gErrSky(4, ErrMajAxis, ErrMinAxis);
	return(ErrArea);
}

double LCModSig::gErr6Sky() 
{
	double ErrArea, ErrMajAxis, ErrMinAxis;
	ErrArea = gErrSky(6, ErrMajAxis, ErrMinAxis);
	return(ErrArea);
}


void LCModSig::CopyBase(const LCModSig &a)
{
	initNULL(true);
	
	Detector = a.Detector;
	GalacticNoise = a.GalacticNoise;
	sSwap(a.DirSwapMem, a.MaxMemRAM);
	dtMesMin = a.dtMesMin;
	tDurMax = a.tDurMax;
	tEndSecur = a.tEndSecur;
	t0 = a.t0;
	NCharfNCfg = a.NCharfNCfg;
	NCharTDIName = a.NCharTDIName;
	
	NfNCfg = 0;
	for(int i=0; i<a.NfNCfg; i++)
		addfCfg(a.fNCfg[i]);
	
	NTDI = a.NTDI;
	TDIName = (char**) MT->AllocMemory(NTDI*sizeof(char*));
	for(int i=0; i<NTDI; i++){
		TDIName[i] =  (char*) MT->AllocMemory(NCharTDIName*sizeof(char));
		for(int ic=0; ic<NCharTDIName; ic++)
			TDIName[i][ic] = a.TDIName[i][ic];
		
	}
	
	sSrcType(a.SrcType); //! < NOTE : One #AllocBaseMem() is done inside
	
	NParStudy = a.NParStudy;
	iParStudy = (int*) MT->AllocMemory(NParStudy*sizeof(int));
	lnParStudy = (bool*) MT->AllocMemory(NParStudy*sizeof(bool));
	ParStuDelta = (double*) MT->AllocMemory(NParStudy*sizeof(double));
	for(int i=0; i<NParStudy; i++){
		iParStudy[i] = a.iParStudy[i];
		lnParStudy[i] = a.lnParStudy[i];
		ParStuDelta[i] = a.ParStuDelta[i];
	}
	FIM4.init(MT, NParStudy, NParStudy);
	FIM6.init(MT, NParStudy, NParStudy);
	CovM4.init(MT, NParStudy, NParStudy);
	CovM6.init(MT, NParStudy, NParStudy);
	CheckInv.init(MT, NParStudy, NParStudy);
	
	
	AllocBaseMem();
}


void LCModSig::CopyParam(const LCModSig &a)
{
	for(int iP=0; iP<Np; iP++)
		p[iP] = a.p[iP];
}



// ***************
// *  Operators  *
// ***************

void LCModSig::operator=(const LCModSig &a)
{
	//! ***** Reallocate all the basis if needed 
	if((NTDI!=a.NTDI)||(NtData!=a.NtData)){
		
		initNULL(true);
		
		Detector = a.Detector;
		GalacticNoise = a.GalacticNoise;
		sSwap(a.DirSwapMem, a.MaxMemRAM);
		dtMesMin = a.dtMesMin;
		tDurMax = a.tDurMax;
		tEndSecur = a.tEndSecur;
		t0 = a.t0;
		NCharfNCfg = a.NCharfNCfg;
		NCharTDIName = a.NCharTDIName;
		
		NfNCfg = 0;
		for(int i=0; i<a.NfNCfg; i++)
			addfCfg(a.fNCfg[i]);
		
		NTDI = a.NTDI;
		TDIName = (char**) MT->AllocMemory(NTDI*sizeof(char*));
		for(int i=0; i<NTDI; i++){
			TDIName[i] =  (char*) MT->AllocMemory(NCharTDIName*sizeof(char));
			for(int ic=0; ic<NCharTDIName; ic++)
				TDIName[i] = a.TDIName[i];
			
		}
		
		sSrcType(a.SrcType); //! < NOTE : One #AllocBaseMem() is done inside
		
		NParStudy = a.NParStudy;
		iParStudy = (int*) MT->AllocMemory(NParStudy*sizeof(int));
		lnParStudy = (bool*) MT->AllocMemory(NParStudy*sizeof(bool));
		ParStuDelta = (double*) MT->AllocMemory(NParStudy*sizeof(double));
		for(int i=0; i<NParStudy; i++){
			iParStudy[i] = a.iParStudy[i];
			lnParStudy[i] = a.lnParStudy[i];
			ParStuDelta[i] = a.ParStuDelta[i];
		}
		FIM4.init(MT, NParStudy, NParStudy);
		FIM6.init(MT, NParStudy, NParStudy);
		CovM4.init(MT, NParStudy, NParStudy);
		CovM6.init(MT, NParStudy, NParStudy);
		CheckInv.init(MT, NParStudy, NParStudy);
		
		
		AllocBaseMem();
		
	}
	for(int iS=0; iS<NTDI; iS++)
		NoiseMinLowf[iS] = a.NoiseMinLowf[iS];
	
	
	
	
	for(int iP=0; iP<Np; iP++){
		p[iP] = a.p[iP];
		pS[iP] = a.pS[iP];
	}
	
	NtData = a.NtData;
	dtMes = a.dtMes;
	tDur = a.tDur;
	NfData = a.NfData;
	df = a.df;
	iFmin = a.iFmin;
	iFmax = a.iFmax;
	NfDataLoc = a.NfDataLoc ;
	
	
	NoIni = a.NoIni;
	NotDat = a.NotDat;
	if((!a.NoIni)&&(!a.NotDat)){ //!< If a has been initialized and a time data has been computed, copy the initialization parameters and the data
		
		for(int iS=0; iS<NTDI; iS++)
			if(tSig[iS]==NULL)
				tSig[iS] = (double*) MT->AllocMemory(NtData*sizeof(double));
		
		for(int iS=0; iS<NTDI; iS++)
			memcpy( tSig[iS], a.tSig[iS], NtData*sizeof(double));
	}
	
	NofDat = a.NofDat;
	if(!a.NofDat){
		for(int iS=0; iS<NTDI; iS++)
			if(fSig[iS] == NULL)
				fSig[iS] = (dcomplex*) MT->AllocMemory(NfData*sizeof(dcomplex));
		for(int iS=0; iS<NTDI; iS++)
			memcpy( fSig[iS], a.fSig[iS], NfData*sizeof(dcomplex));
	}
	
	NoSNR = a.NoSNR;
	for(int iS=0; iS<NTDI; iS++)
		SNR[iS] = a.SNR[iS];
	
	NoFIM = a.NoFIM;
	FIM4 = a.FIM4;
	FIM6 = a.FIM6;
	
	NoErr = a.NoErr;
	CovM4 = a.CovM4;
	CovM6 = a.CovM6;
	for(int iP=0; iP<NParStudy; iP++){
		Err4[iP] = a.Err4[iP];
		Err6[iP] = a.Err6[iP];
	}
	
	
}


// ***************************************
// * Linking and initialization methods  *
// ***************************************

void LCModSig::AllocBaseMem()
{
	if(NoBaseMem){
		
		fDSig = (dcomplex***) MT->AllocMemory(NParStudy*sizeof(dcomplex**));
		for(int iP=0; iP<NParStudy; iP++)
			fDSig[iP] = (dcomplex**) MT->AllocMemory(NTDI*sizeof(dcomplex*));
		tSig = (double**) MT->AllocMemory(NTDI*sizeof(double*));
		fSig = (dcomplex**) MT->AllocMemory(NTDI*sizeof(dcomplex*));
		tSig2 = (double**) MT->AllocMemory(NTDI*sizeof(double*));
        tSig3 = (double**) MT->AllocMemory(NTDI*sizeof(double*)); // added by Sofiane
        tSig4 = (double**) MT->AllocMemory(NTDI*sizeof(double*)); // added by Sofiane
		NoiseMinLowf = (double*) MT->AllocMemory(NTDI*sizeof(double));
		SNR = (double*) MT->AllocMemory(NTDI*sizeof(double));
		FIM4.init(MT, NParStudy, NParStudy);
		CovM4.init(MT, NParStudy, NParStudy);
		FIM6.init(MT, NParStudy, NParStudy);
		CovM6.init(MT, NParStudy, NParStudy);
		CheckInv.init(MT, NParStudy, NParStudy);
		for(int i=0; i<NTDI; i++){
			tSig[i] = NULL;
			fSig[i] = NULL;
			tSig2[i] = NULL;
            tSig3[i] = NULL; // added by Sofiane
            tSig4[i] = NULL; // added by Sofiane
			for(int iP=0; iP<NParStudy; iP++)
				fDSig[iP][i] = NULL;
			NoiseMinLowf[i] = LC::DBLMAXALLOW;
			SNR[i] = 0.;
		}
		
		
	}
	
	NoBaseMem = false;
}


void LCModSig::CleanBaseMem()
{
	if(tSig != NULL){
		for(int i=0; i<NTDI; i++)
			if(tSig[i] != NULL)
				MT->Free(tSig[i], NtData*sizeof(double));
		MT->Free(tSig, NTDI*sizeof(double*));
	}
	
	if(tSig2 != NULL){
		for(int i=0; i<NTDI; i++)
			if(tSig2[i] != NULL)
				MT->Free(tSig2[i], NtData*sizeof(double));
		MT->Free(tSig2, NTDI*sizeof(double*));
	}
    
    if(tSig3 != NULL){
		for(int i=0; i<NTDI; i++)
			if(tSig3[i] != NULL)
				MT->Free(tSig3[i], NtData*sizeof(double));
		MT->Free(tSig3, NTDI*sizeof(double*));
	}  // added by Sofiane
    
    if(tSig4 != NULL){
		for(int i=0; i<NTDI; i++)
			if(tSig4[i] != NULL)
				MT->Free(tSig4[i], NtData*sizeof(double));
		MT->Free(tSig4, NTDI*sizeof(double*));
	}  // added by Sofiane
	
	if(fSig != NULL){
		for(int i=0; i<NTDI; i++)
			if(fSig[i] != NULL)
				MT->Free(fSig[i], NfData*sizeof(dcomplex));
		MT->Free(fSig, NTDI*sizeof(dcomplex*));
	}
	
	if(fDSig != NULL){
		for(int iP=0; iP<NParStudy; iP++){
			if(fDSig[iP] != NULL){
				for(int iTDI=0; iTDI<NTDI; iTDI++){
					if(fDSig[iP][iTDI] != NULL)
						MT->Free(fDSig[iP][iTDI], NfDataLoc*sizeof(dcomplex));
					fDSig[iP][iTDI] = NULL;
				}
				MT->Free(fDSig[iP], NTDI*sizeof(dcomplex*));
				fDSig[iP] = NULL;
			}
		}
		MT->Free(fDSig, NParStudy*sizeof(dcomplex**));
		fDSig = NULL;
	}
	
	if(SNR != NULL)
		MT->Free(SNR, NTDI*sizeof(double));
	
	if(NoiseMinLowf != NULL)
		MT->Free(NoiseMinLowf, NTDI*sizeof(double));
	
	tSig = NULL;
	tSig2 = NULL;
    tSig3 = NULL;  // added by Sofiane
    tSig4 = NULL;  // added by Sofiane
	fSig = NULL;
	fDSig = NULL;
	SNR = NULL;
	NoiseMinLowf = NULL;
	
}


void LCModSig::CleanNDatMem()
{
	if(tSig != NULL){
		for(int i=0; i<NTDI; i++){
			if(tSig[i] != NULL)
				MT->Free(tSig[i], NtData*sizeof(double));
			tSig[i] = NULL;
		}
	}
	
	if(tSig2 != NULL){
		for(int i=0; i<NTDI; i++){
			if(tSig2[i] != NULL)
				MT->Free(tSig2[i], NtData*sizeof(double));
			tSig2[i] = NULL;
		}
	}
	
    if(tSig3 != NULL){
		for(int i=0; i<NTDI; i++){
			if(tSig3[i] != NULL)
				MT->Free(tSig3[i], NtData*sizeof(double));
			tSig3[i] = NULL;
		}
	}
    
    if(tSig4 != NULL){
		for(int i=0; i<NTDI; i++){
			if(tSig4[i] != NULL)
				MT->Free(tSig4[i], NtData*sizeof(double));
			tSig4[i] = NULL;
		}
	}
    
	if(fSig != NULL){
		for(int i=0; i<NTDI; i++){
			if(fSig[i] != NULL)
				MT->Free(fSig[i], NfData*sizeof(dcomplex));
			fSig[i] = NULL;
		}
	}
	
	if(fDSig != NULL){
		for(int iP=0; iP<NParStudy; iP++){
			if(fDSig[iP] != NULL){
				for(int iTDI=0; iTDI<NTDI; iTDI++){
					if(fDSig[iP][iTDI] != NULL)
						MT->Free(fDSig[iP][iTDI], NfDataLoc*sizeof(dcomplex));
					fDSig[iP][iTDI] = NULL;
				}
			}
		}
	}
}


void LCModSig::Reset(bool InitToo)
{
	AllocBaseMem();
	
	for(int iP=0; iP<Np; iP++)
		pS[iP] = false;
	if(InitToo)
		NoIni = true;
	NotDat = true;
	NofDat = true;
	NoSNR = true;
	NoFIM = true;
	NoErr = true;
}


int LCModSig::PrepareSim(bool RequireInit)
{	
	if(Sim != NULL)
		delete Sim;
	
	if(NoIni){
		
		AllocBaseMem();
		
		
		double FreqMax(0.);
		bool ChangedtMes(false);
		bool AllowRandom(true);
		
		dtMes = dtMesMin;
		
		if(SrcType==SPINBBHHHARM){
			
			//! ***** Because the time frame id always changing and the computation time of FFT is negligible compare to the time data 
			MT->unsetBestFTFwd();
			
			//! ***** Change spin amplitude to be valide with the model 
			if(p[6]>0.98)
				p[6] = 0.98;
			if(p[7]>0.98)
				p[7] = 0.98;
			
			//! ***** Change symmetric mass ratio to be able to compute FIM : \f$ \eta < 0.25 - 2 \times 10^{-6} \f$ then \f$ < 0.994359 \f$ because \f$ \eta = {q \over (1+q)^2} \f$
			Cout << "Parameters : ";
            for(int iP=0; iP<Np; iP++)
				Cout << " " << p[iP];
			Cout << Endl;
			double M(0.);
            //M = p[3]+p[2];
            M = p[15]/pow(p[17],3./5.);
			if(p[17]>0.25-2.e-6)
				p[17] = 0.25-2.e-6;
			
			//! ***** Set the time frame
			if (M >= 1.e5 && M <= 5.e5)
				dtMes = 4.;    
			if (M > 5.e5)
				dtMes = 16.; 
            //dtMes = 16.; 
            
            Cout << "M = " << M << " ==> dtMes = " << dtMes << Endl;
			
		}
        
        // added by Sofiane **********************
        
        if(SrcType==NRWAVE){
            Cout << "Parameters : ";
            for(int iP=0; iP<Np; iP++)
				Cout << " " << p[iP];
            
            double M(p[3]);
            
            //! ***** Set the time frame
			if (M >= 1.e5 && M <= 5.e5)
				dtMes = 4.;    
			if (M > 5.e5)
				dtMes = 16.;  
            
            Cout << "M = " << M << " ==> dtMes = " << dtMes << Endl;
            
            tDur = tDurMax;
			if(p[4]+tEndSecur < tDurMax)
				tDur = p[4]+tEndSecur ;
            
        }
        
        // end added by Sofiane ******************
        
        
		
		
			
		ChangedtMes = false;
		AllowRandom = true;
		
		//! *** Update the time step regarding the maximal frequency
		do{
			Cout << " Choose time : time step = " << dtMes << " s ,  duration = " << tDur << " s" << Endl; 
			
			//! *** Create the simulator and configure it using the configuration files
			Sim  = new LISACode(MT);
			for(int i=0; i<NfNCfg; i++)
				Sim->config(fNCfg[i]);
			
			
			
			Sim->setTimeInfo(dtMes, t0, tDur);
			Sim->setdtPhy(dtMes);
			
			if(GWmod == NULL){
				if(SrcType==SPINBBHHHARM)
					Sim->AddGW("SpinBBHHighHarm");
				if(SrcType==GALBIN)
					Sim->AddGW("GalacticBinary");
				if(SrcType==NRWAVE)
					Sim->AddGW("NRwave");
				GWmod = Sim->GWget();
			}else{
				Sim->GWLinkExt(GWmod);
			}
			
			//! *** Set parameter values : first randomize everything and after set using the catalogue
			if(AllowRandom){
				AllowRandom = false;
				Sim->GWRandAllParam();
			}else{  
				for(int iP=0; iP<Np; iP++)
					if( (SrcType==NRWAVE)&&((iP<8)||(iP>=16)) )   
						Sim->GWsetParam(0, iP, p[iP]);
			}
			
			for(int iP=0; iP<Np; iP++)
				if(pS[iP])
					Sim->GWsetParam(0, iP, p[iP]); 
			
			//! *** Recover the parameter 
			for(int iP=0; iP<Np; iP++){
				p[iP] = Sim->GWgetParam(0, iP);
				//Cout << "\t\t + P" << iP << " = " << p[iP] << Endl;
			}
			for(int iSP=0; iSP<NGWSpePar; iSP++)
				GWmod->setSpecialParam(GWSpeParInd[iSP], GWSpeParVal[iSP]);
			
			if(SrcType==NRWAVE){
				tDur = tDurMax;
				if(p[4]+7.e6+tEndSecur < tDurMax)
					tDur = p[4]+7.e6+tEndSecur ;
				Cout << "Change duration from " << tDurMax << " to " << tDur << " ( tc = " << p[4] << " , safe addition = " << 7.e6+tEndSecur << " )" << Endl;  
				Sim->setTimeInfo(dtMes, t0, tDur);
			}			
			
			MT->MemDisplay();
			
			//! *** Initialization of simulation
			if(Sim->init())
				return 1;
			
			
			FreqMax = Sim->GWgetFreqMax(0);
			Cout << "FreqMax = " << FreqMax << " , dtMes = " << dtMes << " , 1./(4.*dtMes) = " << 1./(4.*dtMes) <<  Endl; 
			
			//! *** Update time step according to the maximal frequency 
			// ChangedtMes = true;  
			if(!ChangedtMes){
				while(FreqMax < 1./(4.*dtMes)){
					ChangedtMes = true;
					dtMes *= 2.;
					delete Sim;
					Sim = NULL;
				}
				if(dtMes>512.)
					dtMes = 512.;
			}else{
				ChangedtMes = false;
			}
			
			
		}while(ChangedtMes);
		
		
	}else{
		
		Sim  = new LISACode(MT);
		for(int i=0; i<NfNCfg; i++)
			Sim->config(fNCfg[i]);
		Sim->setTimeInfo(dtMes, t0, tDur);
		Sim->setdtPhy(dtMes);
		
		if(GWmod == NULL){
			if(SrcType==SPINBBHHHARM)
				Sim->AddGW("SpinBBHHighHarm");
			if(SrcType==GALBIN)
				Sim->AddGW("GalacticBinary");
			if(SrcType==NRWAVE)
				Sim->AddGW("NRwave");
		}else{
			Sim->GWLinkExt(GWmod);
		}
		
		//! *** Set only the parameter which are necessary to define the waveform
		if(SrcType==SPINBBHHHARM){
			
			for(int iP=0; iP<Np; iP++)
				if( (iP!=2)&&(iP!=3)&&(iP!=16)&&(iP!=18) )
					Sim->GWsetParam(0, iP, p[iP]); 
		}
		if(SrcType==GALBIN){
			for(int iP=0; iP<Np; iP++)
				if( (iP<8) )
					Sim->GWsetParam(0, iP, p[iP]); 
		}
		if(SrcType==NRWAVE){
			for(int iP=0; iP<Np; iP++){
				/*
				if(NGWSpePar==1){
					if( (iP<6) || (iP>=16) )
						Sim->GWsetParam(0, iP, p[iP]);
				}else{
					if( (iP<8) )
                 Sim->GWsetParam(0, iP, p[iP]);
                 }*/
				if( (iP<8) || (iP>15) )  // changed by Sofiane    
					Sim->GWsetParam(0, iP, p[iP]);
			}
		}
		
		if(RequireInit)
			if(Sim->init())
				return 1;
	}
    if(SrcType==NRWAVE)
        Sim->GWsetSpecialParam(0, 1, 0);
    
	NoIni = false;
	
	return 0;
}


// *********************
// *  Running methods  *
// *********************

int LCModSig::RunSim()
{
	if(Sim==NULL)
		if(PrepareSim(true))
			return 1;
	
	//! *** Look if we have to allocate
	bool tAllocMem(true);
	if(tSig==NULL){
		Cout << "ERROR in LCModSig::PrepareSim : The model has not been configured !" << Endl;
		throw std::invalid_argument("ERROR in LCModSig::PrepareSim : The model has not been configured !");
	}
	if(tSig[0]!=NULL){
		if(NtData == Sim->NDatExpected(1))
			tAllocMem = false;
		else
			CleanNDatMem();
	}
	
	//! *** Add output series with allocation of the series. Warning : free should be done outside, it will not be done in the simulator
	for(int i=0; i<NTDI; i++){
		Sim->AddSerieOut(TDIName[i], 1, NtData, tAllocMem);
	}
	//! *** Link output series
	for(int i=0; i<NTDI; i++)
		Sim->LinkSerieOut(i, tSig[i]);
	
	
	//! *** Set parameters for the Fourrier data, if time data has been allocated = it's the first allocation or #NtData has been modified
	if(tAllocMem){
		NfData = MT->getNfFTreal(NtData);
		df = 1./(dtMes*NtData);
		
		double FreqMaxLoc(Sim->GWgetFreqMax(0));
		double FreqMinLoc(Sim->GWgetFreqMin(0));
		
		double tsFreqMaxLoc(0.), tMFreqMaxLoc(0.);
		if((SrcType==SPINBBHHHARM)||(SrcType==NRWAVE)){
			double MLoc(0.), EtaLoc(0.), tcLoc(0.);
			if(SrcType==SPINBBHHHARM){
				MLoc = LC::TSUN * p[15] * pow(p[17],-3./5.);
				EtaLoc = p[17] ;
				tcLoc = p[4] ;
			}
			if(SrcType==NRWAVE){
				MLoc = LC::TSUN * p[3];
				EtaLoc = p[8] / ( (1.+p[8])*(1.+p[8]) ); 
				tcLoc = p[4] ;
			}
			//! *** Uncomment if you want to stop at separation = 10M 

			
			//FreqMaxLoc = 1./(M_PI*MLoc*pow(10.,3./2.));
			//FreqMinLoc =1./(M_PI*MLoc*pow(10.,3./2.));
			
			tMFreqMaxLoc = pow(M_PI*MLoc*FreqMaxLoc,-2./3.);
			tsFreqMaxLoc = tcLoc - pow(tMFreqMaxLoc,4.) * (5.*MLoc)/(EtaLoc*256);
			Cout << "Mtot = " << MLoc << " s, Eta = " << EtaLoc ;
		}
		Cout << " >>> t = " << tsFreqMaxLoc <<  " s = " << tMFreqMaxLoc << " M ==> FreqMax = " << FreqMaxLoc << Endl;
		
		
		iFmin = MAX((FreqMinLoc/df), 0)+2;
		iFmax = MIN((FreqMaxLoc/df), NfData)-2;
		
		
		NfDataLoc = iFmax-iFmin ;
		Cout << "Frequency frame : df = " << df << " Hz ,  Fmin = " << iFmin*df << " Hz ,  Fmax = " << iFmax*df << " Hz" << Endl;
	}
	
	if(DispDetails)
		Sim->DispInfo("");
	
	//! *** Compute the time series : Run the simulator
	Sim->Run();
	
	
	
	//! *** Delete the simulator
	delete Sim;
	Sim = NULL;
	
	//! *** If it's required, the results are recorded for checking 
	if(CheckOut){
		char tmpN[100000];
		sprintf(tmpN, "%s_tSig.txt", CurrentNameOut);
		Cout << "Write time checking output in " << tmpN << " ..." << Endl;
		std::ofstream fChecktSig(tmpN);
		fChecktSig.precision(12);
		for(int i=0; i<NtData; i++){
			fChecktSig << i*dtMes;
			for(int iS=0; iS<NTDI; iS++)
				fChecktSig << " " <<  tSig[iS][i];
			fChecktSig << Endl;
		}
		fChecktSig.close();
	}
	
	NotDat = false;
	
	return 0;
	
}


int LCModSig::ComputeFT()
{
	if(NotDat)
		if(RunSim())
			return 1;
	/*
	 if((tSig==NULL)&&((tSig!=NULL)&&(tSig[0]!=NULL))){
	 Cout << "ERROR in LCModSig::ComputeFT : You need to compute time data before doing a Fourrier transform ! " << Endl;
	 throw std::invalid_argument("ERROR in LCModSig::ComputeFT : You need to compute time data before doing a Fourrier transform ! ");
	 }
	 */
	
	
	
	//! *** Allocate data if needed
	for(int i=0; i<NTDI; i++)
		if(fSig[i] == NULL)
			fSig[i] = (dcomplex*) MT->AllocMemory(NfData*sizeof(dcomplex));
	
	//! *** Compute Fourrier transform
	MT->FTMakeFwdMulti(tSig, fSig, NtData, NTDI);
	
	//! *** Normalize Fourrier transform
	for(int iF=0; iF<NfData; iF++)
		for(int iS=0; iS<NTDI; iS++)
			fSig[iS][iF] *= dtMes;
	
	//! *** If it's required, the results are recorded for checking 
	if(CheckOut){
		char tmpN[10000];
		sprintf(tmpN, "%s_fSig.txt", CurrentNameOut);
		Cout << "Write frequency checking output in " << tmpN << " ..." << Endl;
		std::ofstream fCheckfSig(tmpN);
		fCheckfSig.precision(12);
		for(int i=0; i<NfData; i++){
			fCheckfSig << i*df; 
			for(int iS=0; iS<NTDI; iS++)
				fCheckfSig << " " <<  norm(fSig[iS][i]) << " " <<  arg(fSig[iS][i]) ;
			fCheckfSig << Endl;
		}
		fCheckfSig.close();
	}
	
	NofDat = false;
	
	return 0;
}


void LCModSig::ComputeSNR()
{
	dcomplex tmpRes;
	double tmpSP;
	
	strcpy(CurrentNameOut,BaseNameOut);
	
	ComputeFT();	//!< Compute frequency data 
	
	//! *** For checking : recording the cumulative SNR
	std::ofstream fCheckCumul;
	if(CheckOut){
		char tmpN[100000];
		sprintf(tmpN, "%s_Cumul.txt", CurrentNameOut);
		Cout << "Write SNR checking output in " << tmpN << " ..." << Endl;
		fCheckCumul.open(tmpN);
		fCheckCumul.precision(12);
	}
	
	//! *** Compute individual SNR
	for(int iS=0; iS<NTDI; iS++){
		tmpRes = 0.;
		for(int iF=iFmin; iF<iFmax; iF++){
			tmpRes += fSig[iS][iF] * conj(fSig[iS][iF]) / PSDNoise(iF*df, iS);
			if(CheckOut)
				fCheckCumul << iF*df << " " <<  4.0 * df * tmpRes.real() << " " << PSDNoise(iF*df, iS) << Endl;
		}
		if(CheckOut)
			fCheckCumul << "#" << Endl;
		tmpSP = 4.0 * df * tmpRes.real();
		SNR[iS] = sqrt(tmpSP);
	}
	if(CheckOut)
		fCheckCumul.close();
	
	
	NoSNR = false;
}




void LCModSig::ComputeFIM()
{
	
	//! ***** Set swaping
	
	char ** fNSwapMem(NULL);
	FILE * fSwapMem;
	
	double MemGW(0.);
	if(GWmod!=NULL){
		MemGW = GWmod->getMemUse()/(1024.*1024.*1024.);
	}
	double MemReq( MemGW + ( (NTDI*(3.*NtData + NParStudy*2.*NfDataLoc) + 2*NtData)*sizeof(double)/(1024.*1024.*1024.) ) );
	double MemMinReq( MemGW + ( (NTDI*(2.*NtData + MAX(2*NfDataLoc, NtData)) + 2*NtData)*sizeof(double) /(1024.*1024.*1024.) ) );
	Cout << Endl << "\t Required memory = " << MemReq << " Go (max = " << MaxMemRAM << " Go , minimal required = " << MemMinReq << " Go , waveform = " << MemGW << " Go )"; 
	if(MemReq > MaxMemRAM ){
		//! *** Use the disk
		Cout << " ==> Swap on the disk (if FIM needed)." << Endl;
		SwapMem = true;
		fNSwapMem = (char**) MT->AllocMemory(NParStudy*sizeof(char*));
		int RndNum0(MT->RandUniform(0, 10)), RndNum1(MT->RandUniform(0, 10)), RndNum2(MT->RandUniform(0, 10));
		time_t RawTime;
		time( & RawTime);
		int iTime((int)(RawTime));
		char * NameOut;
		NameOut = (char*) MT->AllocMemory(10000*sizeof(char));
		MT->pathExtractName(CurrentNameOut, NameOut);
		for(int iP=0; iP<NParStudy; iP++){
			fNSwapMem[iP] = (char*) MT->AllocMemory(16000*sizeof(char));
			sprintf(fNSwapMem[iP], "%s/tmpSwapMem-%s-%d-%d%d%d-%d.bin" , DirSwapMem, NameOut, iP, RndNum0, RndNum1, RndNum2, iTime );
		}
		MT->Free(NameOut, 10000*sizeof(char));
	}else{
		//! *** Allocation of memory in RAM
		Cout << " ==> Use RAM only ." << Endl;
		SwapMem = false;
		for(int iP=0; iP<NParStudy; iP++)
			for(int iS=0; iS<NTDI; iS++)
				if(fDSig[iP][iS] == NULL)
					fDSig[iP][iS] = (dcomplex*) MT->AllocMemory(NfDataLoc*sizeof(dcomplex));
	} 
	
	for(int iP=0; iP<Np; iP++)
		Cout << "   " << iP << ":" << p[iP];
	Cout << Endl;
	
	//! ***** Compute analytical derivative just with the waveform as it is : 
    // loop on the parameters to study considering only parameters for which one the derivative is analytic
	
	for(int iP=0; iP<NParStudy; iP++){
		
		sprintf(CurrentNameOut,"%s-P%d",BaseNameOut,iP);
		
		//! ** For spinbbh and NR waveform : Distance \f$ {\partial h \over \partial D } = - { h \over D} \f$ or \f$ {\partial h \over \partial \ln D } = - h  \f$
		if(((SrcType==SPINBBHHHARM)||(SrcType==NRWAVE))&&(iParStudy[iP]==5)){
            
			if(NofDat)
				ComputeFT();
            
			for(int iS=0; iS<NTDI; iS++){
				int iFF(0);
				for(int iF=iFmin; iF<iFmax; iF++, iFF++){
					if(lnParStudy[iP])
						fSig[iS][iF] *= -1.;
					else 
						fSig[iS][iF] /= -p[5];
					if(!SwapMem)
						fDSig[iP][iS][iFF] = fSig[iS][iF];
				}
			}
			
			if(SwapMem){
				fSwapMem = fopen(fNSwapMem[iP], "wb");
				for(int iS=0; iS<NTDI; iS++){
					fwrite(fSig[iS]+iFmin, sizeof(dcomplex), NfDataLoc, fSwapMem);
				}
				fclose(fSwapMem);
			}
		}
        
        
        //! ** For galactic binary : Amplitude \f$ {\partial h \over \partial A } = h \f$ or \f$ {\partial h \over \partial \ln A } = A h  \f$
		if((SrcType==GALBIN)&&(iParStudy[iP]==4)){
			if(NofDat)
				ComputeFT();
			for(int iS=0; iS<NTDI; iS++){
				int iFF(0);
				for(int iF=iFmin; iF<iFmax; iF++, iFF++){
					if(!lnParStudy[iP])
						fSig[iS][iF] /= p[4];
					if(!SwapMem)
						fDSig[iP][iS][iFF] = fSig[iS][iF];
				}
			}
			
			if(SwapMem){
				if(MT->Disp())
					Cout << "Write to swap " << fNSwapMem[iP] << " ..." << Endl;
				fSwapMem = fopen(fNSwapMem[iP], "wb");
				for(int iS=0; iS<NTDI; iS++){
					fwrite(fSig[iS]+iFmin, sizeof(dcomplex), NfDataLoc, fSwapMem);
				}
				fclose(fSwapMem);
			}
		}
    }
    
    
    //! ***** Compute semi-analytical derivative with a modified waveform : total mass and time at coalescence for NR (hybrids)
    for(int iP=0; iP<NParStudy; iP++){
        
        //! ** For NR waveform : \f$ {\partial h \over \partial tc } = {\partial h \over \partial tau }    \f$
		if((SrcType==NRWAVE)&&(iParStudy[iP]==4)){
            
            PrepareSim(false);
            
            Sim->GWsetSpecialParam(0, 1, 1);
            
            Sim->init();
            
            RunSim();
            
            if(CheckOut){
                CheckOut = false;
                ComputeFT();
                CheckOut = true;
            }else{
                ComputeFT();
            }
                
            if(CheckOut){
				char tmpN[10000];
				sprintf(tmpN, "%s_fDSig-%d.txt", CurrentNameOut, iParStudy[iP]);
				std::ofstream fCheckfSig(tmpN);
				fCheckfSig.precision(12);
				for(int i=0; i<NfData; i++){
					fCheckfSig << i*df; 
					for(int iS=0; iS<NTDI; iS++)
						fCheckfSig << " " <<  norm(fSig[iS][i]) << " " <<  arg(fSig[iS][i]) << " " <<  real(fSig[iS][i]) << " " <<  imag(fSig[iS][i]) ;
					fCheckfSig << Endl;
				}
				fCheckfSig.close();
			}
			
            for(int iS=0; iS<NTDI; iS++){
				int iFF(0);
				for(int iF=iFmin; iF<iFmax; iF++, iFF++){
					if(lnParStudy[iP])
						fSig[iS][iF] *=  p[4];
					if(!SwapMem)
						fDSig[iP][iS][iFF] = fSig[iS][iF];
				}
			}
			
			if(SwapMem){
				fSwapMem = fopen(fNSwapMem[iP], "wb");
				for(int iS=0; iS<NTDI; iS++){
					fwrite(fSig[iS]+iFmin, sizeof(dcomplex), NfDataLoc, fSwapMem);
				}
				fclose(fSwapMem);
			}
            
		}
        
        
        //! ** For NR waveform ( \f$ \tau = t/M \f$ ): \f$ {\partial h \over \partial M } = A ( h / M - {t \over M^2} {\partial h \over \partial \tau }     \f$
		if((SrcType==NRWAVE)&&(iParStudy[iP]==3)){

            PrepareSim(false);
            
            Sim->GWsetSpecialParam(0, 1, 2);
            
            Sim->init();
            
            RunSim();
            
            if(CheckOut){
                CheckOut = false;
                ComputeFT();
                CheckOut = true;
            }else{
                ComputeFT();
            }

            
            if(CheckOut){
				char tmpN[10000];
				sprintf(tmpN, "%s_fDSig-%d.txt", CurrentNameOut, iParStudy[iP]);
				std::ofstream fCheckfSig(tmpN);
				fCheckfSig.precision(12);
				for(int i=0; i<NfData; i++){
					fCheckfSig << i*df; 
					for(int iS=0; iS<NTDI; iS++)
						fCheckfSig << " " <<  norm(fSig[iS][i]) << " " <<  arg(fSig[iS][i]) << " " <<  real(fSig[iS][i]) << " " <<  imag(fSig[iS][i]) ;
					fCheckfSig << Endl;
				}
				fCheckfSig.close();
			}
			
            for(int iS=0; iS<NTDI; iS++){
				int iFF(0);
				for(int iF=iFmin; iF<iFmax; iF++, iFF++){
					if(lnParStudy[iP])
						fSig[iS][iF] *=  p[3];  // added Sofiane why ?
					if(!SwapMem)
						fDSig[iP][iS][iFF] = fSig[iS][iF];
				}
			}
			
			if(SwapMem){
				fSwapMem = fopen(fNSwapMem[iP], "wb");
				for(int iS=0; iS<NTDI; iS++){
					fwrite(fSig[iS]+iFmin, sizeof(dcomplex), NfDataLoc, fSwapMem);
				}
				fclose(fSwapMem);
			}
		}
		
	}
	
	
	
	//! ***** Compute numerical derivative : loop on the parameters to study considering only parameters for which one the derivative has to be computed numerically
	
	for(int iP=0; iP<NParStudy; iP++){
		
		
		//! **** Exept the distance which have already been studied
		if( (!(((SrcType==SPINBBHHHARM)||(SrcType==NRWAVE))&&(iParStudy[iP]==5))) && (!((SrcType==NRWAVE)&&(iParStudy[iP]==4))) && (!((SrcType==NRWAVE)&&(iParStudy[iP]==3))) && (!((SrcType==GALBIN)&&(iParStudy[iP]==4))) ){
		//if( (!(((SrcType==SPINBBHHHARM)||(SrcType==NRWAVE))&&(iParStudy[iP]==5))) && (!((SrcType==NRWAVE)&&(iParStudy[iP]==4))) && (!((SrcType==GALBIN)&&(iParStudy[iP]==4))) ){

            
			Cout << Endl << ">>>>>>>>>>>> Compute derivative for parameter " << iP << " --> " << iParStudy[iP] << " ..." << Endl;
			
			//double pSave(p[iParStudy[iP]]);
			double DeltaPar;
			
			//! ***** Compute central derivative : Two steps : one for plus Delta and one for minus Delta
			for(int iD=0; iD<4; iD++){
				// changed by Sofiane iD<2 to iD<4
				sprintf(CurrentNameOut,"%s-P%d-D%d",BaseNameOut,iP,iD);
				
				//! ** Prepare the simulator without initialization because we need to put Delta on the parameters before starting the simulation
				PrepareSim(false);
				
				//! ** Apply the delta
				if(DeltaParSpecInd == iParStudy[iP])
					DeltaPar = DeltaParSpecVal;
				else
					DeltaPar = Sim->GWgetDeltaParam(iParStudy[iP]);
				
                DeltaPar *= FactStepD;
                
				if(iD == 0)
					Sim->GWAddExtDeltaParam(iParStudy[iP], -2.*DeltaPar);
                if(iD == 1)
					Sim->GWAddExtDeltaParam(iParStudy[iP], -DeltaPar);
                if(iD == 2)
					Sim->GWAddExtDeltaParam(iParStudy[iP], +DeltaPar);
                if(iD == 3)
					Sim->GWAddExtDeltaParam(iParStudy[iP], +2.*DeltaPar);
			//	else
			//		Sim->GWAddExtDeltaParam(iParStudy[iP], +DeltaPar);
				ParStuDelta[iP] = DeltaPar;
				
				//! ** Now (the Delta has been putted) we can intialize the simulation and run it
				Sim->init();
				
				RunSim();
				
				if(iD==0){
					if(tSig2 == NULL)
						throw std::invalid_argument("ERROR in ComputeFIM : tSig2 has not been allocated at all !");
					if(tSig2[0] == NULL)
						for(int iS=0; iS<NTDI; iS++)
							tSig2[iS] = (double*) MT->AllocMemory(NtData*sizeof(double));
					for(int iS=0; iS<NTDI; iS++)
						memcpy( tSig2[iS], tSig[iS], NtData*sizeof(double));
				}
                if(iD==1){
					if(tSig3 == NULL)
						throw std::invalid_argument("ERROR in ComputeFIM : tSig3 has not been allocated at all !");
					if(tSig3[0] == NULL)
						for(int iS=0; iS<NTDI; iS++)
							tSig3[iS] = (double*) MT->AllocMemory(NtData*sizeof(double));
					for(int iS=0; iS<NTDI; iS++)
						memcpy( tSig3[iS], tSig[iS], NtData*sizeof(double));
				}
                
                if(iD==2){
					if(tSig4 == NULL)
						throw std::invalid_argument("ERROR in ComputeFIM : tSig4 has not been allocated at all !");
					if(tSig4[0] == NULL)
						for(int iS=0; iS<NTDI; iS++)
							tSig4[iS] = (double*) MT->AllocMemory(NtData*sizeof(double));
					for(int iS=0; iS<NTDI; iS++)
						memcpy( tSig4[iS], tSig[iS], NtData*sizeof(double));
				}
				
			}
			
			//! *** Compute the centered derivative
			for(int iTDI=0; iTDI<NTDI; iTDI++){
				for(int i=0; i<NtData; i++){
					//if(i<20)
					//	Cout << iTDI << " " << i << " " << tSig[iTDI][i] << " " << tSig2[iTDI][i] ;
					tSig[iTDI][i] = ( tSig2[iTDI][i] - 8.*tSig3[iTDI][i]  + 8.*tSig4[iTDI][i] - tSig[iTDI][i]  )/(12.*DeltaPar); 
                //    tSig[iTDI][i] = ( - tSig3[iTDI][i]  + tSig4[iTDI][i]  )/(2.*DeltaPar); 
					//if(i<20)
					//	Cout << " " << tSig[iTDI][i] << Endl;
				}
			}
			
			Cout << "Delta of parameter " << iP << " (" << iParStudy[iP] << ") = " << DeltaPar << Endl;
			
			//! *** If it's required, the results are recorded for checking 
			sprintf(CurrentNameOut,"%s-P%d",BaseNameOut,iP);
			if(CheckOut){
				char tmpN[100000];
				sprintf(tmpN, "%s_tDSig-%d.txt", CurrentNameOut, iParStudy[iP]);
				std::ofstream fChecktSig(tmpN);
				fChecktSig.precision(12);
				for(int i=0; i<NtData; i++){
					fChecktSig << i*dtMes;
					for(int iS=0; iS<NTDI; iS++)
						fChecktSig << " " <<  tSig[iS][i];
					fChecktSig << Endl;
				}
				fChecktSig.close();
			}
			
			//! *** Compute and normalize Fourrier transform
            if(CheckOut){
                CheckOut = false;
                ComputeFT();
                CheckOut = true;
            }else{
                ComputeFT();
            }
			
			
			if(CheckOut){
				char tmpN[10000];
				sprintf(tmpN, "%s_fDSig-%d.txt", CurrentNameOut, iParStudy[iP]);
				std::ofstream fCheckfSig(tmpN);
				fCheckfSig.precision(12);
				for(int i=0; i<NfData; i++){
					fCheckfSig << i*df; 
					for(int iS=0; iS<NTDI; iS++)
						fCheckfSig << " " <<  norm(fSig[iS][i]) << " " <<  arg(fSig[iS][i]) << " " <<  real(fSig[iS][i]) << " " <<  imag(fSig[iS][i]) ;
					fCheckfSig << Endl;
				}
				fCheckfSig.close();
			}
			
			//! *** Multiply by the corresponding value for the Log case
			if(lnParStudy[iP]){
				for(int iS=0; iS<NTDI; iS++){
					int iFF(0);
					for(int iF=iFmin; iF<iFmax; iF++, iFF++)
						fSig[iS][iF] *= p[iParStudy[iP]];
				}
			}
			
			if(!SwapMem){
				for(int iS=0; iS<NTDI; iS++){
					int iFF(0);
					for(int iF=iFmin; iF<iFmax; iF++, iFF++)
						fDSig[iP][iS][iFF] = fSig[iS][iF];
				}
			}else{
				if(MT->Disp())
					Cout << "Write to swap " << fNSwapMem[iP] << " ..." << Endl;
				fSwapMem = fopen(fNSwapMem[iP], "wb");
				for(int iS=0; iS<NTDI; iS++){
					fwrite(fSig[iS]+iFmin, sizeof(dcomplex), NfDataLoc, fSwapMem);
					//Coutm << iP << " : " << iS << " : " << fSig[iS][0] << " " <<  fSig[iS][NfDataLoc/2] << " " <<  fSig[iS][NfDataLoc-1] << Endl;
				}
				fclose(fSwapMem); 
			}
		}
	}
	
	//! *** Deallocattion of memory
	if(SwapMem){
		if(tSig2 != NULL){
			for(int i=0; i<NTDI; i++){
				if(tSig2[i] != NULL)
					MT->Free(tSig2[i], NtData*sizeof(double));
				tSig2[i] = NULL;
			}
		}
        
        if(tSig3 != NULL){
			for(int i=0; i<NTDI; i++){
				if(tSig3[i] != NULL)
					MT->Free(tSig3[i], NtData*sizeof(double));
				tSig3[i] = NULL;
			}
		}// added by Sofiane
        
        if(tSig4 != NULL){
			for(int i=0; i<NTDI; i++){
				if(tSig4[i] != NULL)
					MT->Free(tSig4[i], NtData*sizeof(double));
				tSig4[i] = NULL;
			}
		} // added by Sofiane
	}	
	
	
	
	//! ****** Computation of Fisher Information Matrix 
	dcomplex ** pfDSig1(NULL);
	dcomplex ** pfDSig2(NULL);
	dcomplex tmpRes;
	
	//! ** Allocation of memory for the data temporary stored on disk
	if(SwapMem){
		pfDSig1 = (dcomplex**) MT->AllocMemory( NTDI * sizeof(dcomplex*) );
		pfDSig2 = (dcomplex**) MT->AllocMemory( NTDI * sizeof(dcomplex*) );
		for(int iS=0; iS<NTDI; iS++){
			pfDSig1[iS] = (dcomplex*) MT->AllocMemory( NfDataLoc * sizeof(dcomplex) );
			pfDSig2[iS] = (dcomplex*) MT->AllocMemory( NfDataLoc * sizeof(dcomplex) );
		}
		MT->MemDisplay();
	}
	
	FIMnonan = true;
	for(int iP1=0; iP1<NParStudy; iP1++){
		for(int iP2=0; iP2<NParStudy; iP2++){
			double tmpFIM4(0.), tmpFIM6(0.);
			
			if(SwapMem){
				if(MT->Disp())
					Cout << "Read from swap " << fNSwapMem[iP1] << " ..." << Endl;
				fSwapMem = fopen(fNSwapMem[iP1],"rb");
				for(int iS=0; iS<NTDI; iS++){
					fread(pfDSig1[iS], sizeof(dcomplex), NfDataLoc, fSwapMem);
					//Coutm << iP1 << " : " << iS << " : " << pfDSig1[iS][0] << " " <<  pfDSig1[iS][NfDataLoc/2] << " " <<  pfDSig1[iS][NfDataLoc-1] << Endl;
				}
				fclose(fSwapMem);
				if(MT->Disp())
					Cout << "Read from swap " << fNSwapMem[iP2] << " ..." << Endl;
				fSwapMem = fopen(fNSwapMem[iP2],"rb");
				for(int iS=0; iS<NTDI; iS++){
					fread(pfDSig2[iS], sizeof(dcomplex), NfDataLoc, fSwapMem);
					//Coutm << iP2 << " : " << iS << " : " << pfDSig2[iS][0] << " " <<  pfDSig2[iS][NfDataLoc/2] << " " <<  pfDSig2[iS][NfDataLoc-1] << Endl;
				}
				fclose(fSwapMem);
				
			}else{
				pfDSig1 = fDSig[iP1];
				pfDSig2 = fDSig[iP2];
			}
			
			//! *** For checking : recording the cumulative SNR
			std::ofstream fCheckCumul;
			if(CheckOut){
				char tmpN[100000];
				sprintf(tmpN, "%s_CumulFIM%d%d.txt", CurrentNameOut, iP1, iP2);
				Cout << "Write FIM cumuml checking output in " << tmpN << " ..." << Endl;
				fCheckCumul.open(tmpN);
				fCheckCumul.precision(12);
			}
			
			// ! *** Computing inner products
			for(int iS=0; iS<NTDI; iS++){
				tmpRes = 0.;
				for(int iF=0; iF<NfDataLoc; iF++){
					tmpRes += pfDSig1[iS][iF] * conj(pfDSig2[iS][iF]) / PSDNoise((iFmin+iF)*df, iS);
					if(CheckOut)
						fCheckCumul << iF*df << " " << iS << " " <<  4.0 * df * tmpRes.real() << " " << real(pfDSig1[iS][iF]) << " " << imag(pfDSig1[iS][iF]) << " " << real(pfDSig2[iS][iF]) << " " << imag(pfDSig2[iS][iF]) << " " << PSDNoise((iFmin+iF)*df, iS) <<  Endl;
				}
				if(MT->wcmp(TDIName[iS], "X"))
					tmpFIM4 += 4.0 * df * tmpRes.real();
				if( (MT->wcmp(TDIName[iS], "Am")) || (MT->wcmp(TDIName[iS], "Em")) )
					tmpFIM6 += 4.0 * df * tmpRes.real();
				if(CheckOut)
					fCheckCumul << "# tmpFIM4 = " << tmpFIM4 << "  tmpFIM6 = " << tmpFIM6 << Endl;
				
			}
			if(std::isnan(tmpFIM4)||std::isnan(tmpFIM6))
				FIMnonan = false;
			FIM4(iP1,iP2,tmpFIM4);
			FIM6(iP1,iP2,tmpFIM6);
			if(CheckOut)
				fCheckCumul.close();
		}
	}
	
	//! ** Free of memory for the data temporary stored on disk
	if(SwapMem){
		for(int iS=0; iS<NTDI; iS++){
			MT->Free(pfDSig1[iS], NfDataLoc * sizeof(dcomplex));
			MT->Free(pfDSig2[iS], NfDataLoc * sizeof(dcomplex));
		}
		MT->Free(pfDSig1, NTDI * sizeof(dcomplex*));
		MT->Free(pfDSig2, NTDI * sizeof(dcomplex*));
	}
	pfDSig1 = NULL;
	pfDSig2 = NULL;
	
	
	if(SwapMem){
		char cmd[10010];
		//! ** Remove data file
		for(int iP=0; iP<NParStudy; iP++){
			sprintf(cmd, "rm %s", fNSwapMem[iP]);
			system(cmd);
		}
		
		//! ** Free memory for data file name 
		if(fNSwapMem != NULL){
			for(int iP=0; iP<NParStudy; iP++)
				if(fNSwapMem[iP] != NULL)
					MT->Free(fNSwapMem[iP], 10000*sizeof(char));
			MT->Free(fNSwapMem, NParStudy*sizeof(char*));
			fNSwapMem = NULL;
		}
	}
	
	
	NoFIM = false;
	
	if(DispDetails){
		Cout << FIM4 << Endl;
		Cout << FIM6 << Endl;
	}
}



void LCModSig::ComputeErr()
{
	int InvFIMMethod(1);
	double InvDiffUnit;
	double MinInvDiffUnit;
	
	if(NoFIM)
		ComputeFIM();
	
	
	bool Use4links(false), Use6links(false);
	for(int iS=0; iS<NTDI; iS++){
		if( (MT->wcmp(TDIName[iS], "X"))||(MT->wcmp(TDIName[iS], "Y"))||(MT->wcmp(TDIName[iS], "Z")) )
			Use4links = true;
		if( (MT->wcmp(TDIName[iS], "Am")) || (MT->wcmp(TDIName[iS], "Em")) )
			Use6links = true;
	}
	
	
	if(Use4links){
		
		if(FIMnonan){
			
			//! ********** Inverse 4 links FIM **********
			
			MinInvDiffUnit = LC::DBLMAXALLOW;
			//! ***** Test invertion of Fisher Matrix 4 links with Gauss-Jordan method 
			Cout << "Test invertion of Fisher Matrix 4 links with Gauss-Jordan method : " << Endl;
			CovM4 = FIM4.Inv(0);
			CheckInv = CovM4 * FIM4;
			InvDiffUnit = CheckInv.DiffWithUnit();
			if(InvDiffUnit < MinInvDiffUnit){
				MinInvDiffUnit = InvDiffUnit;
				InvFIMMethod = 0;
			}
			Cout << "--> Results = " << CovM4 << Endl;
			Cout << "--> Checking product (should be unit matrix) : difference with unit = " << InvDiffUnit << CheckInv << Endl;
			
			//! ***** Test invertion of Fisher Matrix 4 links with SVD method 
			Cout << "Test invertion of Fisher Matrix 4 links with SVD method : " << Endl;
			CovM4 = FIM4.Inv(1);
			CheckInv = CovM4 * FIM4;
			InvDiffUnit = CheckInv.DiffWithUnit();
			if(InvDiffUnit < MinInvDiffUnit){
				MinInvDiffUnit = InvDiffUnit;
				InvFIMMethod = 1;
			}
			Cout << "--> Results = " << CovM4 << Endl;
			Cout << "--> Checking product (should be unit matrix) : difference with unit = " << InvDiffUnit << CheckInv << Endl;
			
			//! ***** Test invertion of Fisher Matrix 4 links with LU method 
			Cout << "Test invertion of Fisher Matrix 4 links with LU method : " << Endl;
			CovM4 = FIM4.Inv(2);
			CheckInv = CovM4 * FIM4;
			InvDiffUnit = CheckInv.DiffWithUnit();
			if(InvDiffUnit < MinInvDiffUnit){
				MinInvDiffUnit = InvDiffUnit;
				InvFIMMethod = 2;
			}
			Cout << "--> Results = " << CovM4 << Endl;
			Cout << "--> Checking product (should be unit matrix) : difference with unit = " << InvDiffUnit << CheckInv << Endl;
			
			
			if(MinInvDiffUnit<1.){
				//! ***** Invert Fisher matrix
				Cout << " ==> Used inversion method = " << InvFIMMethod << Endl;
				CovM4 = FIM4.Inv(InvFIMMethod);
				Cout << "============================ Covariance Matrix for 4 links :  ============================" << Endl;
				Cout << CovM4 << Endl;
				
				
				
				//! ***** Normalization of covariance matrix
				for(int iP=0; iP<NParStudy; iP++){
					for(int jP=0; jP<NParStudy; jP++){
						if(iP!=jP)
							CovM4(iP, jP, CovM4(iP,jP)/sqrt(CovM4(iP,iP)*CovM4(jP,jP)) );
					}
				}
			}else{
				Cout << "WARNING: None of the methods seems to give a decent inverse for FIM 4 links so we cannot compute the errors (we set covariance matrix at 0)" << Endl;
				for(int iP=0; iP<NParStudy; iP++)
					for(int jP=0; jP<NParStudy; jP++)
						CovM4(iP, jP, 0.);
			}
		}else{
			Cout << "WARNING: The FIM for 4 links have some problem so we cannot compute the errors (we set covariance matrix at 0)" << Endl;
			for(int iP=0; iP<NParStudy; iP++)
				for(int jP=0; jP<NParStudy; jP++)
					CovM4(iP, jP, 0.);
		}
		
		Cout << "============================ Normalized covariance Matrix for 4 links ============================" << Endl;
		Cout << CovM4 << Endl;
		
		
		//! ***** Errors on parameters
		for(int iP=0; iP<NParStudy; iP++){
			if(lnParStudy[iP])
				Err4[iP] = sqrt(CovM4(iP,iP))*p[iParStudy[iP]];
			else
				Err4[iP] = sqrt(CovM4(iP,iP));
		}
	}
	
	
	if(Use6links){
		
		if(FIMnonan){
			//! ********** Inverse 6 links FIM **********
			
			MinInvDiffUnit = LC::DBLMAXALLOW;
			//! ***** Test invertion of Fisher Matrix 6 links with Gauss-Jordan method 
			Cout << "Test invertion of Fisher Matrix 6 links with Gauss-Jordan method : " << Endl;
			CovM6 = FIM6.Inv(0);
			CheckInv = CovM6 * FIM6;
			InvDiffUnit = CheckInv.DiffWithUnit();
			if(InvDiffUnit < MinInvDiffUnit){
				MinInvDiffUnit = InvDiffUnit;
				InvFIMMethod = 0;
			}
			Cout << "--> Results = " << CovM6 << Endl;
			Cout << "--> Checking product (should be unit matrix) : difference with unit = " << InvDiffUnit << CheckInv << Endl;
			
			//! ***** Test invertion of Fisher Matrix 6 links with SVD method 
			Cout << "Test invertion of Fisher Matrix 6 links with SVD method : " << Endl;
			CovM6 = FIM6.Inv(1);
			CheckInv = CovM6 * FIM6;
			InvDiffUnit = CheckInv.DiffWithUnit();
			if(InvDiffUnit < MinInvDiffUnit){
				MinInvDiffUnit = InvDiffUnit;
				InvFIMMethod = 1;
			}
			Cout << "--> Results = " << CovM6 << Endl;
			Cout << "--> Checking product (should be unit matrix) : difference with unit = " << InvDiffUnit << CheckInv << Endl;
			
			//! ***** Test invertion of Fisher Matrix 6 links with LU method 
			Cout << "Test invertion of Fisher Matrix 6 links with LU method : " << Endl;
			CovM6 = FIM6.Inv(2);
			CheckInv = CovM6 * FIM6;
			InvDiffUnit = CheckInv.DiffWithUnit();
			if(InvDiffUnit < MinInvDiffUnit){
				MinInvDiffUnit = InvDiffUnit;
				InvFIMMethod = 2;
			}
			Cout << "--> Results = " << CovM6 << Endl;
			Cout << "--> Checking product (should be unit matrix) : difference with unit = " << InvDiffUnit << CheckInv << Endl;
			
			
			//! ***** Invert Fisher matrix
			Cout << " ==> Used inversion method = " << InvFIMMethod << Endl;
			CovM6 = FIM6.Inv(InvFIMMethod);
			Cout << "============================ Covariance Matrix for 6 links :  ============================" << Endl;
			Cout << CovM6 << Endl;
			
			
			if(MinInvDiffUnit<1.){
				//! ***** Normalization of covariance matrix
				for(int iP=0; iP<NParStudy; iP++){
					for(int jP=0; jP<NParStudy; jP++){
						if(iP!=jP)
							CovM6(iP, jP, CovM6(iP,jP)/sqrt(CovM6(iP,iP)*CovM6(jP,jP)) );
					}
				}
			}else{
				Cout << "WARNING: None of the methods seems to give a decent inverse for FIM 6 links so we cannot compute the errors (we set covariance matrix at 0)" << Endl;
				for(int iP=0; iP<NParStudy; iP++)
					for(int jP=0; jP<NParStudy; jP++)
						CovM6(iP, jP, 0.);
			}
			
		}else{
			Cout << "WARNING: The FIM for 6 links have some problem so we cannot compute the errors (we set covariance matrix at 0)" << Endl;
			for(int iP=0; iP<NParStudy; iP++)
				for(int jP=0; jP<NParStudy; jP++)
					CovM6(iP, jP, 0.);
		}
		
		Cout << "============================ Normalized covariance Matrix for 6 links ============================" << Endl;
		Cout << CovM6 << Endl;
		
		
		//! ***** Errors on parameters
		for(int iP=0; iP<NParStudy; iP++){
			if(lnParStudy[iP])
				Err6[iP] = sqrt(CovM6(iP,iP))*p[iParStudy[iP]];
			else
				Err6[iP] = sqrt(CovM6(iP,iP));
		}
		
	}
	
	NoErr = false;
}




double LCModSig::LogLikelihood(dcomplex ** fDat, int NSigDat, int NfDatDat, int iFminDat, double dfDat)
{
	if(NofDat)
		if(ComputeFT())
			return(-1.1e30);
	
	if(NSigDat != NTDI){
		Cout << "ERROR in LCModSig::LogLikelihood : The number of channel in data (" << NSigDat << ") and signal ( " << NTDI << "are not the same !" << Endl;
		throw std::invalid_argument("ERROR in LCModSig::LogLikelihood : The number of channel in data and signal are not the same !");
	}
	if(NfDatDat != NfDataLoc){
		Cout << "ERROR in LCModSig::LogLikelihood : The number of data in data (" << NfDatDat << ") and signal ( " << NfDataLoc << "are not the same !" << Endl;
		throw std::invalid_argument("ERROR in LCModSig::LogLikelihood : The number of data in data and signal are not the same !");
	}
	
	double sSPh(0.), Norm(0.);
	dcomplex tmpsSPh(0.), tmpNorm(0.);
	
	for(int iS=0; iS<NTDI; iS++){
		tmpsSPh = 0.;
		tmpNorm = 0.;
		for(int iF=0; iF<NfDataLoc; iF++){
			double PSDNoiseVal( PSDNoise((iFmin+iF)*df, iS) );
			tmpsSPh += fDat[iS][iF]       * conj(fSig[iS][iFmin+iF]) / PSDNoiseVal;
			tmpNorm += fSig[iS][iFmin+iF] * conj(fSig[iS][iFmin+iF]) / PSDNoiseVal;
		}
		sSPh += 4.0 * df * tmpsSPh.real();
		Norm += 4.0 * df * tmpNorm.real();
	}
	
    //! WARNING : TO BE CHECKED : to have a coherent used of likelihood in MCMC we should divide by NTDI
    
	return(sSPh - 0.5*Norm);
}


void LCModSig::MCMCjump(double heat, double scale)
{
	/*for(int iP=0; iP<Np; iP++)
	 Cout << "\t- p["<< iP << "] init = " << p[iP] << Endl;
	 Cout << Endl;*/
	
	LCMatrix EigenVects, EigenVals;
	
	//! ***** Compute eigen vectors and eigen values
	if(NTDI==1)
		FIM4.EigenVsymCompute(EigenVects, EigenVals);
	if(NTDI==2)
		FIM6.EigenVsymCompute(EigenVects, EigenVals);
	
	
	if(MT->Disp()){
		Cout << "MCMCjump : Eigen values  : ";
		for(int iP=0; iP<NParStudy; iP++)
			Cout << " \t" << EigenVals(iP,0);
		Cout << Endl << "MCMCjump : Eigen vectors : " << Endl;
		for(int iP=0; iP<NParStudy; iP++){
			for(int ic=0; ic<NParStudy; ic++)
				Cout << " \t\t" << EigenVects(iP,ic);
			Cout << Endl ;
		}
		Cout << Endl;
	}
	
	
	//! *** Apply temperature 
	EigenVals *= (1./heat);
	
	std::vector<double> Jump(NParStudy);
	
	//! ***** Normal distribution and compute the jump
	std::vector<double> Amps(NParStudy);
	double xdraw;
	double snp ( sqrt((double)NParStudy) );
    double chi ( (double)fabs(MT->RandGaussian( 0., 1.)) );
    double norm ( 0. );
    for(int i = 0 ; i < NParStudy ; i++){
	    xdraw = MT->RandGaussian( 0., 1.); 
	    Amps[i] = xdraw;
        norm += xdraw*xdraw;
	}
    norm = sqrt(norm);
	for(int i = 0 ; i < NParStudy ; i++){
        Amps[i] = snp * Amps[i]/norm;
		//Cout << "\t- Amps["<< i << "] = " << Amps[i] << Endl;
	}
	//Cout << Endl;
	
    norm = 0.0;
    for(int i = 0 ; i < NParStudy ; i++){
        norm += Amps[i]*Amps[i];
	    Amps[i] = chi*Amps[i]/sqrt(EigenVals(i,0));
	}
	/*Cout << "\tNorm of dZ = " << norm  << "  chi =  " << chi << Endl;
	 for(int i = 0 ; i < NParStudy ; i++)
	 Cout << "\t- Amps["<< i << "] = " << Amps[i] << Endl;
	 Cout << Endl;*/
	
	for (int j =0; j<NParStudy; j++){
        for (int i = 0; i<NParStudy; i++){
            Jump[j] += EigenVects(j,i)*Amps[i];
        }
    }
	/*for(int i = 0 ; i < NParStudy ; i++)
	 Cout << "\t+ Jump["<< i << "] = " << Jump[i] << Endl;
	 Cout << Endl;*/
	
	//! **** Apply the jump
	for(int iP=0; iP<NParStudy; iP++){
		if(lnParStudy[iP]){
			p[iParStudy[iP]] = exp(log(p[iParStudy[iP]]) + scale*Jump[iP]);
		}else{
			p[iParStudy[iP]] = p[iParStudy[iP]] + scale*Jump[iP];
		}
	}
	
	
	//! **** Convert spin/L angles back to the proper range
	if((SrcType==SPINBBHHHARM)||(SrcType==NRWAVE)||(SrcType==GALBIN)){
		double betaS(p[0]);
		double lambdaS(p[1]);
		double x = cos(betaS)*cos(lambdaS);
		double y = cos(betaS)*sin(lambdaS);
		double z = sin(betaS);
		betaS = asin(z);              
		lambdaS = atan2(y,x);
		if(lambdaS < 0.0) 
			lambdaS += 2.*M_PI;
		p[0] = betaS;
		p[1] = lambdaS;
	}
	
	//! **** Convert spin angles back to the proper range
	if(SrcType==SPINBBHHHARM){
		BackRangeColatLong(p[8],p[10]);
		BackRangeColatLong(p[9],p[11]);
		BackRangeColatLong(p[13],p[14]);
	}
	
	
	for (int iP=0; iP<NParStudy; iP++) {
		
		//! *** Check angle which should be between 0 and 2 pi  
		if( ((SrcType==SPINBBHHHARM)&&((iParStudy[iP]==10)||(iParStudy[iP]==11)||(iParStudy[iP]==12)||(iParStudy[iP]==14)||(iParStudy[iP]==18))) || ((SrcType==GALBIN)&&((iParStudy[iP]==2)||(iParStudy[iP]==7))) || ((SrcType==NRWAVE)&&((iParStudy[iP]==2)||(iParStudy[iP]==7))) ){
			if(p[iParStudy[iP]] < 0.0)
			{
				do{
					p[iParStudy[iP]] += 2.*M_PI;
				} while(p[iParStudy[iP]] < 0.0);
			}
			else {
				if(p[iParStudy[iP]] > 2.*M_PI){
					do{
						p[iParStudy[iP]] -= 2.*M_PI;
					} while(p[iParStudy[iP]] > 2.*M_PI);
				}
			}
		}
		
		//! *** Check angle which should be between 0 and pi  
		if( ((SrcType==SPINBBHHHARM)&&((iParStudy[iP]==8)||(iParStudy[iP]==9)||(iParStudy[iP]==13)||(iParStudy[iP]==18))) || ((SrcType==GALBIN)&&((iParStudy[iP]==3))) || ((SrcType==NRWAVE)&&((iParStudy[iP]==6))) ){
			if(p[iParStudy[iP]] < 0.0)
			{
				do{
					p[iParStudy[iP]] += 2.*M_PI;
				} while(p[iParStudy[iP]] < 0.0);
			}
			else {
				if(p[iParStudy[iP]] > M_PI){
					p[iParStudy[iP]] = acos(cos(p[iParStudy[iP]]));
				}
			}
		}
		
		//! *** Check spin amplitude
		if( (SrcType==SPINBBHHHARM) && ((iParStudy[iP]==6)||(iParStudy[iP]==7)) ){
			p[iParStudy[iP]] = fabs(p[iParStudy[iP]]);
			if(p[iParStudy[iP]]>0.98)
				p[iParStudy[iP]] = 0.98;
		}
		
		//! *** Check Eta in [0.01,0.25]
		if( (SrcType==SPINBBHHHARM) && ((iParStudy[iP]==17)) ){
			if(p[iParStudy[iP]]<0.01)
				p[iParStudy[iP]] = 0.0100001;
			if(p[iParStudy[iP]]>0.25)
				p[iParStudy[iP]] = 0.249999;
		}
		
		//! *** Check Eta in [0.01,0.25]
		if( (SrcType==SPINBBHHHARM) && ( ((iParStudy[iP]>=2)&&(iParStudy[iP]<=7)) || ((iParStudy[iP]>=15)&&(iParStudy[iP]<=17)) ) ){
			double Pmin,Pmax;
			GWmod->getRange(iParStudy[iP], Pmin, Pmax);
			if(p[iParStudy[iP]]<Pmin)
				p[iParStudy[iP]] = (1.+1.e-10)*Pmin;
			if(p[iParStudy[iP]]>Pmax)
				p[iParStudy[iP]] = (1.-1.e-10)*Pmax;
		}
		
	}
	
	
	Cout << "New jump : ";
	for(int iP=0; iP<Np; iP++)
		Cout << " " << p[iP];	
	Cout << Endl;
	
	
	NotDat = true;
	NofDat = true;
	
}

void LCModSig::BackRangeColatLong(double & thI, double & phI)
{
	double th(thI), ph(phI);
	double x = sin(th)*cos(ph);
	double y = sin(th)*sin(ph);
	double z = cos(th);
	th = acos(z);              
	ph = atan2(y,x);
	if(ph < 0.0) 
		ph += 2.*M_PI;
	thI = th;
	phI = ph;
}



// ******************************
// *  Detector related methods  *
// ******************************



double LCModSig::PSDNoise(double f, int iS)
{
	double Sn;
	double L0s;
	double phiL;
	double Sxy;
	
	//! *** Set armlength in second
	switch (Detector){
		case LISA:
			L0s = 5.0e9/LC::c_SI;
			break;
		case ELISA:
			L0s = 1.0e9/LC::c_SI;
			break;
		case C1:
			L0s = 1.0e9/LC::c_SI;
			break;
		case C2:
			L0s = 1.0e9/LC::c_SI;
			break;
		case C3:
			L0s = 1.0e9/LC::c_SI;
			break;
		case C4:
			L0s = 3.0e9/LC::c_SI;
			break;
		case C5:
			L0s = 2.0e9/LC::c_SI;
			break;
		case C6:
			L0s = 1.0e9/LC::c_SI;
			break;
		case SGOF1:
			L0s = 2.5875e11/LC::c_SI;
			break;
		case SGOM1:
			L0s = 2.1022e10/LC::c_SI;
			break;
		default:
			throw std::invalid_argument("ERROR in gPSDNoise : There is no armlength associated to this detector !");
			break;
	}
	
	phiL = 2.0*M_PI*L0s*f;
	
	//! *** PSD of noise for X, Y, Z
	Sn = 16.0*sin(phiL)*sin(phiL)*(SOpticalNoise(f) + (3.0+cos(2.0*phiL))*SAccelerationNoise(f));
	
	//if(iS==1)
	//	Cout << f << " " << Sn ;
	
	//! *** Add cross-term for A, E (in MLDC convention)
	if((MT->wcmp(TDIName[iS],"Am")) || (MT->wcmp(TDIName[iS],"Em"))){
		Sxy = -4.0*sin(2.0*phiL)*sin(phiL)*(4.0*SAccelerationNoise(f)+SOpticalNoise(f));
		Sn = (2.0/3.0)*(Sn-Sxy);
	}
	//if(iS==1)
	//	Cout << " " << Sxy << " " << Sn << Endl;
	
	//! Add galactic confusion noise
	if(GalacticNoise)
		Sn += SxGal(f, L0s);
	
	
	if((Detector==SGOF1)||(Detector==SGOM1)){
		double NoiseMin(1.e-42);
		if(Detector==SGOF1)
			NoiseMin = 1.e-40;
		if(Detector==SGOM1)
			NoiseMin = 1.e-41;
		if(Sn<NoiseMin)
			Sn = NoiseMin;
	}else{
		
		//! *** Clipping using similar system as Michele
		double fLim(0.25/L0s);
		if(f<fLim){
			//! Look for the mininmum
			if(Sn<NoiseMinLowf[iS]){
				NoiseMinLowf[iS] = Sn;
			}
		}else{
			//! Apply for the mininmum
			if(NoiseMinLowf[iS]>Sn)
				Sn = NoiseMinLowf[iS];
		}
	}
	return(Sn);
}



double LCModSig::SOpticalNoise(double f)
{
	switch (Detector) {
		case LISA:
			return (1.63124e-37*(f*f));
			break;
		case ELISA:
			return (5.07e-38*(f*f));
			break;
		case C1:
			return (7.73e-38*(f*f));
			break;
		case C2:
			return (2.87e-38*(f*f));
			break;
		case C3:
			return (5.12e-38*(f*f));
			break;
		case C4:
			return (23.51e-38*(f*f));
			break;
		case C5:
			return (4.87e-38*(f*f));
			break;
		case C6:
			return (8.45e-38*(f*f));
			break;
		case SGOF1:
			return (1.33e-34*(f*f));
			break;
		case SGOM1:
			return (9.88e-36*(f*f));
			break;
		default:
			throw std::invalid_argument("ERROR in MFLSignal::SOpticalNoise : The detector doesn't have any optical path noise !");
			break;
	}
	
}

double LCModSig::SAccelerationNoise(double f)
{
	if (f < LC::PRECISION)
		return (1.0e-60);
	
	switch (Detector) {
		case LISA:
			return (2.5281e-48*(1.0/(f*f)+1.0e-8/(f*f*f*f)));
			break;
		case ELISA:
			return (6.00e-48*(1+1.e-4/f)/(f*f));
			break;
		case C1:
			return (8.17e-48*pow(1.0/f+1.8e-4/(f*f),2.));
			break;
		case C2:
			return (6.00e-48/(f*f));
			break;
		case C3:
			return (6.00e-48/(f*f));
			break;
		case C4:
			return (6.00e-48/(f*f));
			break;
		case C5:
			return (6.00e-48/(f*f));
			break;
		case C6:
			return (6.00e-48/(f*f));
			break;
		case SGOF1:
			return (2.82e-49 * pow(f,-1.5)/(f*f));
			break;
		case SGOM1:
			return (5.45e-46 * pow((f/0.001),-1.5)/(f*f));
			break;
		default:
			throw std::invalid_argument("ERROR in SAccelerationNoise : The detector doesn't have any acceleration noise !");
			break;
	}
	
}


double LCModSig::SxGal(double f, double L0s)
{
	double phiL, sphiL;
	double SxG(0.0);	
	
	phiL = 2.0*M_PI*f*L0s;							
	sphiL = sin(phiL);
	
	switch (Detector){
		case LISA:
			// ** LISA 
			if((f>1.0e-4)&&(f<=1.0e-3))
				SxG = pow(10.0,-44.62) * pow(f,-2.3);
			if((f>1.0e-3)&&(f<=pow(10.0,-2.7)))
				SxG = pow(10.0,-50.92) * pow(f,-4.4);
			if((f>pow(10.0,-2.7))&&pow(10.0,-2.4))
				SxG = pow(10.0,-62.8) * pow(f,-8.8);
			if((f>pow(10.0,-2.4))&&(f<=1.0e-2))
				SxG = pow(10.0,-89.68) * pow(f,-20.0);
			
			if((f>1.0e-4)&&(f<=1.0e-2))
				SxG *= 16.0 * L0s*L0s * (2.0*M_PI*f)*(2.0*M_PI*f) * sphiL*sphiL ;
			
			break;
		
		case ELISA:
			if((f>4.5e-4)&&(f<=5.3e-4))
				SxG = 1e-13 * pow(f,7.);
			if((f>5.3e-4)&&(f<=2.2e-3))
				SxG = 2.9174e-47 * pow(f,-3.235);
			if((f>2.2e-3)&&(f<=4.e-3))
				SxG = 1.517e-51 * pow(f,-4.85);
			if((f>4.0e-3)&&(f<=5.88e-3))
				SxG = 6.706e-58 * pow(f,-7.5);
			
			if((f>4.5e-4)&&(f<=5.88e-3))
				SxG *= 4.* L0s*L0s * (2.*M_PI*f)*(2.*M_PI*f) * 3./5. * sphiL*sphiL;
			
			break;
			
			
		case C1:
			if((f>4.5e-4)&&(f<=5.3e-4))
				SxG = 1e-13 * pow(f,7.);
			if((f>5.3e-4)&&(f<=2.2e-3))
				SxG = 2.9174e-47 * pow(f,-3.235);
			if((f>2.2e-3)&&(f<=4.e-3))
				SxG = 1.517e-51 * pow(f,-4.85);
			if((f>4.0e-3)&&(f<=5.88e-3))
				SxG = 6.706e-58 * pow(f,-7.5);
			
			if((f>4.5e-4)&&(f<=5.88e-3))
				SxG *= 4.* L0s*L0s * (2.*M_PI*f)*(2.*M_PI*f) *3./5. * sphiL*sphiL;
			
			break;
			
		case C2:
			if((f>4.5e-4)&&(f<=5.3e-4))
				SxG = 1e-13 * pow(f,7.);
			if((f>5.3e-4)&&(f<=2.2e-3))
				SxG = 2.9174e-47 * pow(f,-3.235);
			if((f>2.2e-3)&&(f<=4.e-3))
				SxG = 1.517e-51 * pow(f,-4.85);
			if((f>4.0e-3)&&(f<=5.88e-3))
				SxG = 6.706e-58 * pow(f,-7.5);
			
			if((f>4.5e-4)&&(f<=5.88e-3))
				SxG *= 4.* L0s*L0s * (2.*M_PI*f)*(2.*M_PI*f) *3./5. * sphiL*sphiL;
			
			break;
			
		case C3:
			if((f>4.5e-4)&&(f<=5.3e-4))
				SxG = 1e-13 * pow(f,7.);
			if((f>5.3e-4)&&(f<=2.2e-3))
				SxG = 2.9174e-47 * pow(f,-3.235);
			if((f>2.2e-3)&&(f<=4.e-3))
				SxG = 1.517e-51 * pow(f,-4.85);
			if((f>4.0e-3)&&(f<=5.88e-3))
				SxG = 6.706e-58 * pow(f,-7.5);
			
			if((f>4.5e-4)&&(f<=5.88e-3))
				SxG *= 4.* L0s*L0s * (2.*M_PI*f)*(2.*M_PI*f) *3./5. * sphiL*sphiL;
			
			break;
			
		case C4:
			if((f>1.e-4)&&(f<=5.01e-4))
				SxG = 1.3516e-43 * pow(f,-2.1);
			if((f>5.01e-4)&&(f<=2.07e-3))
				SxG = 1.4813e-47  * pow(f,-3.3);
			if((f>2.07e-3)&&(f<=3.4e-3))
				SxG = 1.17757e-52 * pow(f,-5.2);
			if((f>3.4e-3)&&(f<=5.2e-3))
				SxG = 2.7781e-62 * pow(f,-9.1);
			
			if((f>1.e-4)&&(f<=5.2e-3))
				SxG *= 4.* L0s*L0s * (2.*M_PI*f)*(2.*M_PI*f) *3./5. * sphiL*sphiL;
			
			break;
			
		case C5:
			if((f>1.8e-4)&&(f<=2.4e-4))
				SxG = 5.4e-36;
			if((f>2.4e-4)&&(f<=5.01e-4))
				SxG = 1.3516e-43 * pow(f,-2.1);
			if((f>5.01e-4)&&(f<=2.07e-3))
				SxG = 1.4813e-47  * pow(f,-3.3);
			if((f>2.07e-3)&&(f<=3.4e-3))
				SxG = 1.17757e-52 * pow(f,-5.2);
			if((f>3.4e-3)&&(f<=5.2e-3))
				SxG = 2.7781e-62 * pow(f,-9.1);
			
			if((f>1.8e-4)&&(f<=5.2e-3))
				SxG *= 4.* L0s*L0s * (2.*M_PI*f)*(2.*M_PI*f) *3./5. * sphiL*sphiL;
			
			break;
			
		case C6:
			if((f>4.5e-4)&&(f<=5.3e-4))
				SxG = 1e-13 * pow(f,7.);
			if((f>5.3e-4)&&(f<=2.2e-3))
				SxG = 2.9174e-47 * pow(f,-3.235);
			if((f>2.2e-3)&&(f<=4.e-3))
				SxG = 1.517e-51 * pow(f,-4.85);
			if((f>4.0e-3)&&(f<=5.88e-3))
				SxG = 6.706e-58 * pow(f,-7.5);
			
			if((f>4.5e-4)&&(f<=5.88e-3))
				SxG *= 4.* L0s*L0s * (2.*M_PI*f)*(2.*M_PI*f) *3./5. * sphiL*sphiL;
			
			break;
			
		default:
			throw std::invalid_argument("ERROR in SxGal : There is no confusion noise for this detector !");
			break;
	}
	return(SxG);
}






// end of LISACODE-ModSig.cpp

