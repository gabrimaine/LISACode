/*
 *  LISACODE-LISACode.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 12/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \ingroup main 
 * \{
 */


#ifndef __LISACODE_H
#define __LISACODE_H

#include <unistd.h>
#include "LISACODE-Tools.h"
#include "ezxml.h"
#include "LISACODE-GW.h"
#include "LISACODE-GWFile.h"
#include "LISACODE-GWStochastic.h"
#include "LISACODE-GWGalBin.h"
#include "LISACODE-GWCosmicString.h"
#include "LISACODE-GWSpinBBHHHarm1.h"
#include "LISACODE-GWSpinBBHNR1.h"
#include "LISACODE-Detector.h"
#include "LISACODE-Orbits.h"
#include "LISACODE-OrbitsAnaLISA.h"
#include "LISACODE-OrbitsAnaMLDC.h"
#include "LISACODE-OrbitsAnaOctahedron.h"
#include "LISACODE-OrbitsData.h"
#include "LISACODE-TDIInt.h"
#include "LISACODE-TDIIntStd2002.h"
#include "LISACODE-TDIGen.h"
#include "LISACODE-DataFileWrite.h"


/** \brief Base class for managing the detector.
 * \author A. Petiteau
 * \version 2.0
 * \date 10/04/2011
 *
 * This class is the base class for managing the detector.
 * 
 *
 */
class LISACode
{
protected:
	
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	/** \brief List of pointers on gravitational waves (extern allocation) */
	LCGW ** GWs;
	
	/** \brief Number of pointers on gravitational waves */
	int NGWs;
	
	/** \brief Pointer on the detector */
	LCDetector * LISA;
	
	/** \brief Pointer on orbits */
	LCOrbits * Orb;
	
	/** \brief List of pointer on delays used in TDI
	 *	ind	: iD	:	link
	 *	0	:	1	:	3'->2
	 *	1	:	2	:	1'->3
	 *	2	:	3	:	2'->1
	 *	3	:	1'	:	2 ->3'
	 *	4	:	2'	:	3 ->1'
	 *	5	:	3'	:	1 ->2'
	 */
	LCSerie2 ** DelaysTDI;
	
	/** \brief Last values of delays  */
	double ** LastDelays;
	
	/** \brief Name of delays : NDelaysTDI x 4  */
	char ** NameDelays;
	
	/** \brief Number of pointer on used delays */
	int NDelaysTDI;
	
	/** \brief Pointer on intermediate TDI */
	LCTDIInt ** TDII;
	
	/** \brief Number of pointer on intermediate TDI */
	int NTDII;
	
	/** \brief Pointer on TDI generator */
	LCTDIGen ** TDIG;
	
	/** \brief Number of pointer on TDI generator */
	int NTDIG;
	
	
	/** \brief Time offset */
	double t0;
	
	/** \brief Measurements time step */
	double dtMes;
	
	/** \brief Duration of the simulation */
	double tDur;
	
	/** \brief Duration of the storage for one delay. Fix at 30 by hard at the moment because it's already almost 2 times 9e9 m for a arm */
	double tStoreDelay;
	
	/** \brief Time shift between the computation of detector measurements and the computation of TDI intermediate data */
	double tShiftTDII;
	
	/** \brief Time shift between the computation of TDI intermediate data and the computation of TDI final data (generators) */
	double tShiftTDIG;
	
	/** \brief Number of data in TDI output : \f$ N_{Data}^{TDIG} = t_{Dur} / \Delta_{mes} \f$ : number of step in the phase 3 */
	int NDataTDIG;
		
	/** \brief Number of data to be computed in advanced in TDI intermediate in order to be able to compute TDI output : \f$ N_{AdvData}^{TDII} = N_{MaxDelay}^{TDIG} * t_{StoreDelay} / \Delta_{mes} \f$ : number of step in the phase 2 */
	int NDataReqTDIG;

	/** \brief Number of data to be computed in advanced in the detector in order to be able to compute TDI intermediate : \f$ N_{AdvData}^{Det} = N_{MaxDelay}^{TDII} * t_{StoreDelay} / \Delta_{mes} \f$ : number of step in the phase 1 */
	int NDataReqTDII;

	
	/************ For file outputs ********** */
	
	/** \brief List of pointer on file output */
	LCDataFileWrite ** fOut;
	
	/** \brief Number of output files */
	int NfOut;
	
	/** \brief Type of output :
	 *	0 : Phasemeter
	 *	1 : TDI
	 *	2 : Delays (armlength)
	 */
	std::vector<int> fOutType;
	
	/** \brief List of names of obsevables for each output file : char[16] */
	char ***  fOutObsNames;
	
	/** \brief List of names of obsevables for each output file : char[16] */
	int * NfOutObsNames;
	
    /** \brief Global header : essentially used for xml output*/
    std::ofstream * GlobalHeader;
	
	/************ For serie outputs ********** */
	
	/** \brief List of pointer on time vector output (local allocation but extern deallocation) */
	double ** sOut;
	
	/** \brief List of size of time vector output (local allocation) */
	int * NDatsOut; 
	
	/** \brief Number of time vector output */
	int NsOut;
	
	
	/** \brief Type of serie output (same coding as for #fOutType ) */
	std::vector<int> sOutType;
	
	/** \brief Pointers on data to be store in series during record phase (local allocation) */
	double ** sOutDat;
	
	/** \brief List of index giving the current position where we are currently writing in the serie (local allocation) */
	int * sOutInd;
	
	/** \brief If true display progression of simulation */
	bool DispProg;
	
	/** \brief True if GW are defined only locally and therefore will be deleted with the simulator */
	bool GWLocal;
	
	/** \brief True if we have to optimize the time for gravitational wave : have a good estimation of the first time and last time  
	 *	at which one we will have to compute \f$ h_{+} \f$ and \f$ h_{\times} \f$ 
	 */
	bool GWOptimizeTime;
	
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LISACode();
	
	/*! \brief Standard constructor */
	LISACode(LCTools * MT);
	
	
	/*! \brief Destructor */
	~LISACode();
	
	/*! \brief Initialization of all at NULL */
	void initNULL(bool CleanMem);
	
	
	
	/**********  Configuration methods  **********/
	
	/*! \brief Configuration from XML file */
	void config(char * fNInXML);
	
	/*! \brief Configuration of intermediate TDI */
	void configTDIIntermediate(ezxml_t xmlbloc);
	
	/*! \brief Configuration of orbits from a bloc : <XSIL Name="..." Type="Output"> ... </XSIL> */
	void configOuput(ezxml_t xmlbloc);
	
	
	/* \brief Add output file   (allocating 16 char for each name)
	 *
	 *	@param[in]	TmpfOutType			Type of output : 0->Phasemeter, 1->TDI, 2->Delays
	 *	@param[in]	TmpFileName			File name
	 *	@param[in]	TmpFileType			File type : ASCII, Binary, XML
	 *	@param[in]	TmpfOutObsNames		List of observables
	 *	@param[in]	TmpNfOutObsNames	Number of observables
	 */ 
	void AddOutputFile(int TmpfOutType, char * TmpFileName, TypeFile TmpFileType, char ** TmpfOutObsNames, int TmpNfOutObsNames);
	
	
	/* \brief Read the name of the observables in InAll and store them in OutName (allocating 16 char for each name)
	 *
	 *	@param[in]	InAll		String of characters containing all names in bloc with ',' between each name
	 *	@param[in]	ipIA		Current position of reading in input
	 *	@param[out]	OutName		List of read names
	 *	@param[out]	NOutName	Number of read names in the list
	 */ 
	//void ReadNameObs(const char * InAll, int & posInAll, char ** & OutName, int & NOutName);
	
	/**********  Linking and initalization methods  **********/
	
	/*! \brief Initialization */
	int init();
	
	/*! \brief Return the number of data expected
	 *	@param[in]	ObsName		Name of observable
	 */
	int NDatExpected(int TypeObs);
	
	/*! \brief Add a serie output 
	 *	@param[in]	sOut_n		Pointer on the series
	 *	@param[in]	TypeOut		Type of output
	 *	@param[in]	ObsName		Name of observable
	 */
	void AddSerieOut(char * ObsName, int TypeObs, int & NDatExpect, bool Allocate);
	
	/*! \brief Add a serie output 
	 *	@param[in]	iS			Index of the series
	 *	@param[out]	psOut		Pointer on the serie
	 */
	void LinkSerieOut(int iS, double * &psOut);
	
	
	/**********  Access methods  **********/
	
	double gett0() { return(t0); };
	double getdtMes() { return(dtMes); };
	double gettDur() { return(tDur); };
	void sett0(double t0_n) {t0 = t0_n; };
	void settDur(double tDur_n) {tDur = tDur_n; };
	void setdtMes(double dtMes_n) {dtMes = dtMes_n; };
	
	/*! \brief Set default time info (useful if ther are not define in any xml file) */
	void setTimeInfo(double dtMes_n, double t0_n, double tDur_n);
	/*! \brief Set physical time step in detector */
	void setdtPhy(double dtPhy) { LISA->setdtPhy(dtPhy); };
	
	/*! \brief Set true or false the progression display  */ 
	void setDispProg(bool DispProg_n) {DispProg = DispProg_n; };
	/*! \brief Set the global xml header */
	void setGlobalXMLHeader(char * GlobXMLHeadName);
	
	/**********  Running methods  **********/
	
	/*! \brief Main running */
	void Run();
	
	/*! \brief Record output regarding the running phase (1,2 or 3) */
	void RecordOutputs(int RunPhase, double t, double tTDII, double tTDIG);
	
	/*! \brief Store the delay for a given time in the corresponding list */
	void StoreDelaysTDI(double t);
	
	
	
	/**********  Methods for GWs  **********/
	
	/*! \brief Add a gravitational wave of type :
	 *	@param[in] GWTypeName   Type of GW : SampledPlaneWave , Stochastic , GalacticBinary , SpinBBHHighHarm
	 */
	void AddGW(char * GWTypeName);
	
	/*! \brief Set parameter iP of GW iGW */
	void GWsetParam(int iGW, int iParam, double ParValue);
	
	/*! \brief get parameter iP of GW iGW */
	double GWgetParam(int iGW, int iParam);
	
	/*! \brief Choose randomly the specified GW's parameter */
	void GWRandParam(int iParamRand);
	
	/*! \brief Choose randomly all GW's parameter */
	void GWRandAllParam();
	
	
	/*! \brief Return Delta value of the parameter iP of the first GW :
	 *	@param	iP			Index of parameter to be modified
	 *	@return				Value of delta
	 */
	double GWgetDeltaParam(int iP);
	/*! \brief Add Delta x FactDelta to value of the parameter iP of the first GW :
	 *	@param	iP			Index of parameter to be modified
	 *	@param	FactDelta	Factor apply on delta
	 *	@return				Value of delta
	 */
	double GWAddDeltaParam(int iP, double FactDelta);
	
	/*! \brief Add exterm Delta to value of the parameter iP of the first GW :
	 *	@param	iP			Index of parameter to be modified
	 *	@param	FExtDelta	external value of delta
	 */
	void GWAddExtDeltaParam(int iP, double ExtDelta);
	
	/*! \brief Display all the GW's parameters */
	void GWDispParam(std::ostream * out);
	
	/*! \brief Display all the name of GW's parameters */
	void GWDispParamName(std::ostream * out);
	
	/*! \brief Return minimal frequency of GW iGW */
	double GWgetFreqMin(int iGW);
	
	/*! \brief Return maximal frequency of GW iGW */
	double GWgetFreqMax(int iGW);
	
    /*! \brief Set special parameter iSP of GW iGW */
	void GWsetSpecialParam(int iGW, int iSParam, double SpeParValue);
    
	/*! \brief Return pointer on first GW and set GW as defined externally (it won't be deleted with the simulator) */
	LCGW * GWget();
	
	/*! \brief Link to an extern GW */
	void GWLinkExt(LCGW * ExtGW);
	
	
	/**********  Others methods  **********/
	
	/*! \brief Display information */ 
	void DispInfo(char * BTab);
	
};

#endif //__LISACODE_H

/**\}*/

// end of LISACODE-LISACode.h
