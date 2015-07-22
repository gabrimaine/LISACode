/*
 *  LISACODE-TDIInt.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 20/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \ingroup TDI 
 * \{
 */

#ifndef __LCTDIINT_H
#define __LCTDIINT_H

#include "LISACODE-Tools.h"
#include "ezxml.h"
#include "LISACODE-Serie2.h"
#include "LISACODE-Detector.h"

/** \brief Base class for computing intermediate Time Delay Interferometry.
 * \author A. Petiteau
 * \version 2.0
 * \date 20/04/2011
 *
 * This base class for computing intermediate Time Delay Interferometry.
 * 
 *
 */
class LCTDIInt
{
protected:
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	/** \brief List of pointer on used input data (phasemeter measurements) [extern allocation] */
	LCSerie2 ** Sig;
	
	/** \brief Number pointer on used input data */
	int NSig;
	
	/** \brief Pointer on output data [local allocation] */
	LCSerie2 * Dat;
	
	/** \brief List of pointer on used delays [extern allocation] */
	LCSerie2 ** Delays;
	
	/** \brief Number of pointer on used delays [extern allocation] */
	int NDelays;
	
	/** \brief Type of interpolation */
	INTERP InterpType;
	
	/** \brief Value used in interpolation */
	double InterpUtilValue;
	
	/** \brief Type of interpolation */
	INTERP InterpTypeDelay;
	
	/** \brief Value used in interpolation */
	double InterpUtilValueDelay;
	
	/** \brief Name of the intermediate TDI */
	char Name[8];
	
	
	/******** Equivalent position of the optical bench and phasemeter  ********/
	
	/*! \brief  Equivalent index of spacecraft corresponding to the photodetector-phasemeter */
	int iSC; 
	
	/*!\brief  Equivalent direction flag of the optical bench :  0 if the optical bench is in the direct direction (1->2->3->1), else 1 (1->3->2->1) */
	int IndirectDir;
	
	/** ********** Variables use in initialization but not after ********** */
	
	/** \brief Measurement time step [use in initialization but not after] */
	double dtMes;
	
	/** \brief Number of data in data serie [given by LISACode] */
	int NData;
	
	/** \brief Time shift between the computation of detector measurements and the computation of TDI intermediate data */
	double tShift;
	
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Standard constructor */
	LCTDIInt(LCTools * MT_n);
	
	/*! \brief Destructor */
	virtual ~LCTDIInt();
	
	/*! \brief Initialization at NULL */
	void initNULLBase(bool CleanMem);
	
	
	/***************  Required methods (to be present in derived classes)  ***************/
	
	/*! \brief Initialization at NULL */
	virtual void initNULL(bool CleanMem);
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="TDIIntmediate"> ... </XSIL> */
	virtual void config(ezxml_t xmlbloc);
	
	/*! \brief Link to detector : link the signal and maybe other stuff */
	virtual void LinkDetector(LCDetector * LISA);

	/*! \brief Link to the delays */
	virtual void LinkDelays(LCSerie2** AllDelays, int AllNDelays);

	/*! \brief Initialization */
	virtual void init();
	
	/*! \brief Return the maximum number of combined delay */
	virtual int getNMaxDelay();
	
	/*! \brief Running one step */
	virtual void RunStep(double t);
	
	/*! \brief Display informations */
	virtual void DispInfo(char * BTab);
	
	/***************  Local methods  ***************/
	
	/**********  Configuration methods  **********/
	

	/**********  Access methods  **********/
	
	/*! \brief Return the number of data needed in advance for the interpolation */
	int getNDatInterp();
	
	/*! \brief Set name and equivalent position of the measurments system */
	void setNameEqPos(const char * NameEqPos);
	
	/*! \brief Set time informations */
	void setTimeInfo(double dtMes_n, int NDataTDII_n, double tShiftTDII_n);
	
	double getLastRes(){return(Dat->getBinValue(0));};

	int getiSC() {return(iSC);};
	int getIndirectDir() {return(IndirectDir);};
	LCSerie2 * getDat() {return(Dat);};
	
	
	/**********  Linking and initalization methods  **********/
	
	/*! \brief Initialization */
	void initBase();	
	
	
	/**********  Running methods  **********/
	
	
	/**********  Other methods  **********/
	
	/** \brief Return the noise value corresponding to iN and the delay .
	 *	It return 0 if the noise is NULL
	 *	@param[in] iN		index of the noise in the list of noise
	 *	@param[in] tDelay	time delay refering to the global time 
	 */
	double gN(int iN, double tDelay);
	
	/** \brief Display the basic information */
	void DispInfoBase(char * BTab);
	
};

#endif //__LCTDIInt_H

/**\}*/

// end of LISACODE-TDIInt.h