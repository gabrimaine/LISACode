/*
 *  LISACODE-TDI.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 20/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \defgroup TDI TDI
 * This group contains all the things related to Time Delay Interferometry
 * \{
 */

#ifndef __LCTDIGEN_H
#define __LCTDIGEN_H

#include "LISACODE-Tools.h"
#include "ezxml.h"
#include "LISACODE-Serie2.h"
#include "LISACODE-TDIInt.h"

/*! \brief Description of a pack : one bloc of TDI */
struct TDIPack {
	//! \brief Pointer on signal on which one the delay  
	LCSerie2 * Sig;
	
	//! \brief List of pointers on delays  
	LCSerie2 ** Del;
	
	//! \brief Number of pointers on delays  
	int NDel;
	
	//! \brief Multiplicative factor  
	double Fact ;
	
	
	/* ***** Variable used for constructing the pack ***** */
	
	//! \brief Coded pack
	int Code;
	
};



/** \brief Class for applying the Time Delay Interferometry methods .
 * \author A. Petiteau
 * \version 2.0
 * \date 20/04/2011
 *
 * This class apply the Time Delay Interferometry methods
 * 
 *
 */
class LCTDIGen
{
protected:
	
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	/** \brief List of packs */
	TDIPack * p;
	
	/** \brief Number of packs */
	int Np;
	
	/** \brief Pointer to the result */
	double * Result;
	
	/** \brief True if local allocation for the value of the result */
	bool ResultLocalAlloc;
	
	/** \brief Time shift between the computation of the delay (synchronize with detector measurements) and the computation of TDI generators */
	double tShiftDelay;
	
	/** \brief Time shift between the computation of the signals (intermediate data or detector measurements) and the computation of TDI generator */
	double tShiftSig;
	
	/** \brief Difference between the time shift for delay and time shift for signal. */
	double DtShiftSigDelay;

	/** \brief Type of interpolation */
	INTERP InterpType;
	
	/** \brief Value used in interpolation */
	double InterpUtilValue;

	/** \brief Type of interpolation */
	INTERP InterpTypeDelay;
	
	/** \brief Value used in interpolation */
	double InterpUtilValueDelay;
	
	/** \bief Maximal number of delays */
	int NbMaxDelays;
	
	/** \brief Name of the TDI generator */
	char Name[8];
	
	//! \brief 10 if the code use one digit, 100 if they use 2 digits 
	int powDigit;
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Standard constructor */
	LCTDIGen(LCTools * MT_n);
	
	/*! \brief Destructor */
	~LCTDIGen();
	
	/*! \brief Initialization at NULL */
	void initNULL(bool CleanMem);
	
	
	/**********  Configuration methods  **********/
	
	/*! \brief Configuration of TDI from XML bloc <XSIL Name="..." Type="TDIGenerator"> ... </XSIL> */
	void config(ezxml_t xmlbloc);
	
	/*! \brief Add a pack with the coding and the factor specified in arguments */
	void addPack(int CodePack, double Fact);
	
	/*! \brief Load pack form the name for the pre-regiestred TDI generator */
	void PreRegistred(const char * generatorname);
	
	
	/**********  Linking and initalization methods  **********/
	
	/*! \brief Initialization */
	void init();
	
	/*! \brief Return the code corresponding to the signal of pack iP on which one the delays will be applied
	 */
	int iSigCodePack(int iP);
	
	/*! \brief Link the signal of packs to the detector outputs */
	void LinkSigDetector(LCDetector * LISA);
	
	/*! \brief Link the signal of packs to the intermediate TDI */
	void LinkSigTDIInt(LCTDIInt ** TDII, int NTDII);
	
	/*! \brief Link delays of packs to the stored delays */ 
	void LinkDelays(LCSerie2** AllDelays, int NAllDelays, char ** NameDelays);
	
	/*! \brief Return true if the argument matches the name of the generator */
	bool MatchName(char * ObsName);
	
	/*! \brief Link result value to pointer on output data (usually the pointer on the current data in the output or in series for link LISACode ) */
	void LinkResult2Output(double * pOutData);
	
	
	
	/**********  Access methods  **********/
	
	/*! \brief Return the number of data needed in advance for the interpolation */
	int getNDatInterp();
	
	/*! \brief Return the maximum number of combined delay */
	int getNMaxDelay();
	
	/*! \brief Set time informations */
	void setTimeInfo(double tShiftSig_n, double tShiftDelay_n);
	
	
	/**********  Running methods  **********/
	
	/*! \brief Compute one step of TDI generator */
	double RunStep(double t);
	
	
	/**********  Other methods  **********/
	
	/*! \brief Display informations */
	void DispInfo(char * BTab);
	
	
};

#endif //__LCTDIGEN_H

/**\}*/

// end of LISACODE-OBPM.h