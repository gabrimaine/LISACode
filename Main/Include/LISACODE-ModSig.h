/*
 *  LISACODE-ModSig.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 28/05/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \ingroup main 
 * \{
 */

#ifndef __LCMODSIG_H
#define __LCMODSIG_H

#include "LISACODE-Matrix.h"
#include "LISACODE-LISACode.h"

enum DetectorConfig{
	LISA,
	ELISA,
	C1,
	C2,
	C3,
	C4,
	C5,
	C6,
	SGOF1,
	SGOM1,
	UNKNOWNDC
};


enum GWSrcType{
	SPINBBHHHARM,
	GALBIN,
	NRWAVE,
	UNKNOWNGW
};


/** \brief  Class for managing evaluation of signal based on LISACode.
 * \author A. Petiteau
 * \version 2.0
 * \date 28/06/2011
 *
 * This class is the base class for managing evaluation of signal based on LISACode like SNR, FIM, ... .
 * 
 *
 */
class LCModSig
{
	protected:
	
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	
	/********** Results internal variables **********/
	//! NOTE: The main internal variables are protected because
	
	/** \brief List of SNR of each TDI channel  */
	double * SNR;
	
	/** \brief Fisher Information Matrix for 4 links */
	LCMatrix FIM4;
	
	/** \brief Fisher Information Matrix for 6 links */
	LCMatrix FIM6;
	
	/** \brief Covariance matrix for 4 links */
	LCMatrix CovM4;
	
	/** \brief Covariance matrix for 6 links */
	LCMatrix CovM6;
	
	/** \brief List of parameters errors for 4 links. Size = #Np */
	double * Err4;
	
	/** \brief List of parameters errors for 6 links. Size = #Np */
	double * Err6;
	
	/********** Main internal variables **********/
	
	/** \brief List of parameters values. Size = #Np */
	double * p;
	
	/** \brief True if the corresponding parameters has been set and should not be randomize. Size = #Np */
	bool * pS;
	
	/** \brief Number of parameters */
	int Np;
	
	/** \brief Pointer on simulator */
	LISACode * Sim;
	
	/** \brief Type of source */
	GWSrcType SrcType;
	
	/** \brief Name of TDI channels. Size : NTDI x 16 */
	char ** TDIName;
	
	/** \brief Number of TDI channels */
	int NTDI;
	
	
	/********** Variables for time data **********/
	
	/** \brief Time evolution of signals. Size: NTDI x NtData */
	double ** tSig;
	
	/** \brief Time evolution of signals. Size: NTDI x NtData */
	double ** tSig2;
    
    /** \brief Time evolution of signals. Size: NTDI x NtData  added by Sofiane*/
	double ** tSig3;
    
    /** \brief Time evolution of signals. Size: NTDI x NtData added by Sofiane */
	double ** tSig4;
	
	/** \brief Number of time data in one TDI channel */
	int NtData;
	
	/** \brief Time step  */
	double dtMes;
	
	/** \brief Duration */
	double tDur;
	
	/** \brief Additional time added to the waveform for been sure that we are not missing a part because the real end (for spinbbh) is defined at MECO  */
	double tEndSecur;
	
	
	/********** Variables for frequency data **********/ 
	
	/** \brief Frequency evolution of signals. Size: NTDI x NtData/2+1 */
	dcomplex ** fSig;
	
	/** \brief Derivative of frequency evolution of signals. Size: NTDI x NtData/2+1 */
	dcomplex *** fDSig;
	
	/** \brief Number of frequency data in one TDI channel */
	int NfData;
	
	/** \brief Number of frequency data in the integration range of one TDI channel */
	int NfDataLoc;
	
	/** \brief Frequency step */
	double df;
	
	/** \brief Index of minimum frequency of integration range */
	int iFmin;
	
	/** \brief Index of maximum frequency of integration range */
	int iFmax;
	
	
	/********** Internal variables for parameter studies **********/ 
	
	/** \brief Index of parameters to be studied with the FIM */
	int * iParStudy;
	
	/** \brief Number of parameters to be studied with the FIM */
	int NParStudy;
	
	/** \brief True if we need to consider the log of the parameters to be studied with the FIM */
	bool * lnParStudy;
	
	/** \brief Delta of parameters to be studied with the FIM */
	double * ParStuDelta;
	
	/********** Internal variables storing the change needs **********/ 
	
	/** \brief If true, we need to change the base memories (clean then reallocation) */
	bool NoBaseMem;
	
	/** \brief If true, we need to initialize the simulation with #PrepareSim() */
	bool NoIni;
	
	/** \brief If true, we need to compute time data with #RunSim() */
	bool NotDat;
	
	/** \brief If true, we need to compute frequency data with #ComputeFT() */
	bool NofDat;
	
	/** \brief If true, we need to recompute the SNR */
	bool NoSNR;
	
	/** \brief If true, we need to recompute the FIM */
	bool NoFIM;
	
	/** \brief If true, we need to recompute the errors */
	bool NoErr;
	
	
	/********** For internal swaping **********/ 
	
	/** \brief Name of directory for swaping */
	char DirSwapMem[10000];
	
	/** \brief Maximum memory in RAM (in Go) */
	double MaxMemRAM;
	
	/** \brief True if we need to swap for storing fDSig data  */
	bool SwapMem;
	
	
	/********** Others internal variables **********/ 
	
	/** \brief Minimal value of the noise in each TDI channel. Size = NTDI */
	double * NoiseMinLowf;
	
	/** \brief Matrix use for checking the inversion */
	LCMatrix CheckInv;
	
	/** \brief Number of character allocated for configuration file name */
	int NCharfNCfg;
	
	/** \brief Number of character allocated for TDI name */
	int NCharTDIName;
	
	/** \brief True if there is no nan in FIM */
	bool FIMnonan;
	
	/** \brief Pointer on gravitational wave model */
	LCGW * GWmod;
	
	/** \brief List of index of specific parameters */
	int * GWSpeParInd;
	
	/** \brief List of value of specific parameters */
	double * GWSpeParVal;
	
	/** \brief Number of specific parameters */
	int NGWSpePar;
	
	
public:
	
	/********** Variables defining the configurations **********/
	
	/** \brief Type of detector */
	DetectorConfig Detector;
	
	/** \brief 1 if we consider the galactic noise */
	int GalacticNoise;
	
	/** \brief Minimal time step  */
	double dtMesMin;
	
	/** \brief Maximal duration */
	double tDurMax;
	
	/** \brief Initial time */
	double t0;
    
    /** \brief Factor on step used for numerical derivative */
    double FactStepD;
	
	/** \brief List of names of configurations files. Size : NfCfg x 2048  */
	char ** fNCfg;
	
	/** \brief Number of configuration files */
	int NfNCfg;
	
	/** \brief True if we have to record checking output */
	bool CheckOut;
	
	/** \brief Base name for swaping and checking output */
	char BaseNameOut[10000];
	
	/** \brief Base name for swaping and checking output */
	char CurrentNameOut[10000];
		
	/** \brief Index of parameter for which one we have a special value of Delta */
	double DeltaParSpecInd;
	
	/** \brief Specified value of Delta on parameter. */
	double DeltaParSpecVal;

	/** \brief True if the waveform is not valid and return 0 as results */
	bool NotValidPar;
	
	/** \brief True if we want to display details in log */
	bool DispDetails;
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCModSig();
	
	/*! \brief Standard constructor */
	LCModSig(LCTools * MT_n);
	
	/*! \brief Destructor */
	~LCModSig();
	
	/*! \brief Initialization of all at NULL */
	void initNULL(bool CleanMem);
	
	
	/**********  Access methods  **********/
	
	/*! \brief Set the parameter value 
	 *	@param[in]	iP		Index of parameter
	 *	@param[in]	Val		New value 
	 */
	void sP(int iP, double Val);
	
	/*! \brief Add specific parameter for graviational wave model 
	 *	@param[in]	iSP		Index of specific parameter
	 *	@param[in]	Val		New value 
	 */
	void AddGWSpePar(int iSP, double Value);
	
	/*! \brief Set the source type 	*/
	void sSrcType(GWSrcType SrcType_n);
	
	
	/*! \brief Set TDI generators from 
	 *	@param[in]	TDIAllName	Char string containing all the name separated by coma 
	 */
	void sTDI(char * TDIAllName);
	
	/*! \brief Set the swap system :
	 *	@param[in]	DirSwapMem_n	New path of directory where we will swap
	 *	@param[in]	MaxMemRAM_n		Maximal available size in RAM 
	 */
	void sSwap(const char * DirSwapMem_n, double MaxMemRAM_n);
	
	/*! \brief Set list of parameters to study 
	 *	@param[in]	CharIndPar	Index of parameters to study separated by coma. L before the index means that we have to consider the log. 
	 */
	void sParStudy(char * CharIndPar);
	
	
	GWSrcType gSrcType() const { return(SrcType); };
	double gP(int iP) const {return(p[iP]);};
	double gNp() const {return(Np);};
	int gNTDI() const {return(NTDI);};
	char * gTDIName(int iTDI) const {return(TDIName[iTDI]);};
	int giParStudy(int iP) const {return(iParStudy[iP]);};
	int gNParStudy() const {return(NParStudy);};
	double gMaxMemRAM() const {return(MaxMemRAM); };
	char * gDirSwapMem() {return(DirSwapMem); };
	bool gFIMnonan() {return(FIMnonan);};
	
	void gPStudy(double * & pStudy);
	void gfSig(dcomplex ** &fSigCpy, int & NfDataLocCpy, int & NTDICpy, int & iFminCpy, double & dfCpy);
	double gSNR(int iTDI);
	double gSNR4();
	double gSNR6();
	LCMatrix gFIM4();
	LCMatrix gFIM6();
	LCMatrix gCovM4();
	LCMatrix gCovM6();
	double gErr4(int iP);
	double gErr6(int iP);
	double gErr4Sky();
	double gErr6Sky();
	
	
	/*! \brief Get the error on sky :
	 @param[in]		Nlinks		Number of links 4 or 6
	 @param[out]	ErrMajAxis	Major axis of the ellipse
	 @param[out]	ErrMajAxis	Minor axis of the ellipse
	 @return	Error in steradians
	 */
	double gErrSky(int Nlinks, double & ErrMajAxis, double & ErrMinAxis);
	

	/*! \brief Copy the base informations */
	void CopyBase(const LCModSig &a);
	
	/*! \brief Copy the parameters only*/
	void CopyParam(const LCModSig &a);
	
	
	/********** Operators **********/
	void operator=(const LCModSig &a);
	
	
	
	/**********  Configuration methods  **********/
	
	/*! \brief Add a configuration file name */
	void addfCfg(char * fCfgNew);
	
	
	/**********  Linking and initalization methods  **********/
	
	/*! \brief Reallocate base memory */
	void AllocBaseMem();
	
	/*! \brief Clean base memory */
	void CleanBaseMem();
	
	/*! \brief Clean all memory depending on the number of time step #NtData */
	void CleanNDatMem();
	
	/*! \brief Reset the parameters ... */
	void Reset(bool InitToo);
	
	/*! \brief Initialization of simulation : check the parameter, set the time frame ... */
	int PrepareSim(bool RequireInit);
	
	
	/**********  Running methods  **********/
	
	/*! \brief Running simulation : filling time data. If fNCheck is not NULL, the data are recorded in file fNCheck_tSig.txt . */
	int RunSim();
	
	/*! \brief Compute Fourrier transform of time data. If fNCheck is not NULL, the data are recorded in file fNCheck_fSig.txt . */
	int ComputeFT();
	
	/*! \brief Compute SNR . If fNCheck is not NULL, the cumulative SNR data are recorded in file fNCheck_fSig.txt . */
	void ComputeSNR();
	
	/*! \brief Compute Fisher Information Matrix*/
	void ComputeFIM();
		
	/*! \brief Compute covariance matrix and errors */
	void ComputeErr();
	
	/*! \brief Compute Likelihodd */
	double LogLikelihood(dcomplex ** fDat, int NSigDat, int NfDatDat, int iFminDat, double dfDat);
	
	
	/*! \brief Make an jump for Monte-Carlo Markov Chain from inital point X */
	void MCMCjump(double heat, double scale);
	
	
	/*! \brief Convert colatitude / longitude back to the proper range */
	void BackRangeColatLong(double & thI, double & phI);
	
	
	/**********  Detector related methods  **********/	
	/*! \brief Return the value of noise PSD */ 
	double PSDNoise(double f, int iS);
	
	/*! \brief Return the value of optical path noise */ 
	double SOpticalNoise(double f);
	
	/*! \brief Return the value of acceleration */ 
	double SAccelerationNoise(double f);
	
	/*! \brief Return the value of galactic confusion noise */ 
	double SxGal(double f, double L0s);
	
};

#endif //__LCMODSIG_H

/**\}*/

// end of LISACODE-ModSig.h

