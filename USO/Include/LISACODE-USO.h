/*
 *  LISACODE-USO.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 12/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \defgroup USO USO clock
 * This group contains all the things related to clock (Utra Stable Oscillator).
 * \{
 */


#ifndef __LCUSO_H
#define __LCUSO_H

#include "LISACODE-Tools.h"
#include "LISACODE-Serie2.h"
#include "ezxml.h"


/** \brief Base class for managing the clocks (USO) of spacecrafts.
 * \author A. Petiteau
 * \version 2.0
 * \date 12/04/2011
 *
 * This class is the base class for managing the clocks (USO) of spacecrafts..
 * 
 *
 */
class LCUSO
{
protected:
	
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	/** \brief Index of the spacecraft */
	int iSC;
	
	/** \brief Time serie containing the difference between the time of the USO and the global time : \f$ t_{Shift} = t_{USO} - t_{ref} \f$*/
	LCSerie2 * tShift ;
	
	/** \brief Reference (global) time corresponding to 0 delay in #tShift serie */
	double t0tShift;
	
	/** \brief Fact for converting time shift data (serie) into phase data */
	double tShift2Phase;
	
	/** \brief Type of interpolation */
	INTERP InterpType;
	
	/** \brief Value used in interpolation */
	double InterpUtilValue;
	
	/** \brief Number of data to add at each run step (measurement time step) */
	int N2Add;
	
	/** \brief Type of output data : if true frequency data else phase data */ 
	bool IsFreqOutput;
	
	/** ********** Variables use in initialization but not after ********** */
	
	/** \brief Measurement time step [use in initialization but not after] */
	double dtMes;
	
	/** \brief Physical time step [use in initialization but not after] */
	double dtPhy;
	
	/** \brief Time of the first data [use in initialization but not after] */
	double tFirst;
	
	/** \brief Time of the last data [use in initialization but not after] */
	double tLast;
	
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCUSO();
	
	/*! \brief Standard constructor */
	LCUSO(LCTools * MT_n);
	
	/*! \brief Destructor */
	virtual ~LCUSO();
	
	/*! \brief Initialization at NULL of all base of noise */
	void initNULLBase(bool CleanMem);
	
	/***************  Required methods (to be present in derived classes)  ***************/
	
	/*! \brief Initialization at NULL */
	virtual void initNULL(bool CleanMem);
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="USO"> ... </XSIL> */
	virtual void config(ezxml_t noisexmlbloc);
	
	/*! \brief Configuration of individual parameter */
	virtual void config(int iParam, double ParamVal);
	
	/*! \brief Initialization */
	virtual void init();
	
	/*! \brief Specific computation of the travel time  */
	virtual void genertShift(int iStartBin);
	
	/*! \brief Display information */ 
	virtual void DispInfo(char * BTab);
	
	
	/***************  Local methods  ***************/
	
	/**********  Configuration methods  **********/
	
	/**********  Linking and initalization methods  **********/
	
	/*! \brief Initialization of base of noise */
	void initBase();
	
	
	/**********  Access methods  **********/
	
	void setdtMes(double dtMes_n)	{dtMes=dtMes_n;};
	void setdtPhy(double dtPhy_n)	{dtPhy=dtPhy_n;};
	void settFirst(double tFirst_n)	{tFirst=tFirst_n;};
	void settLast(double tLast_n)	{tLast=tLast_n;};
	void setTimeInfo(double dtMes_n, double dtPhy_n, double tFirst_n, double tLast_n);
	void setFreqOutput() {IsFreqOutput = true;};
	void setPhaseOutput() {IsFreqOutput = false;};
	
	int getiSC(){return(iSC);};
	
	/*! \brief Set the frequency used for defining #tShift2Phase the conversion between time shift and phase or frequency noise */
	void RefFrequency(double RefFreq);
	
	/**********  Running methods  **********/
	
	/*! \brief Run one measurement step : generation of noise */
	void RunStep(double t);
	
	/*! \brief Get time of the USO computed as : \f$ t_{USO} = t + ... \f$  */
	double gT(double t);
	
	/*! \brief Get noise (phase or frequency) : Return the noise value corresponding to the delay relative to the global time of the simulation */
	double gN(double tDelay);
	
	/**********  Others methods  **********/
	
};

#endif //__LCUSO_H

/**\}*/

// end of LISACODE-USO.h