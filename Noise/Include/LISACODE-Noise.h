/*
 *  LISACODE-Noise.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 05/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


/** \defgroup Noise Noise
 * This group contains all the things related to noises.
 * \{
 */


#ifndef __LCNOISE_H
#define __LCNOISE_H

#include "LISACODE-Tools.h"
#include "LISACODE-Serie2.h"
#include "ezxml.h"


/** \brief Base class for managing noises.
 * \author A. Petiteau
 * \version 2.0
 * \date 06/04/2011
 *
 * This class is the base class for managing noises.
 * 
 *
 */
class LCNoise
{
protected:
	
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	/** \brief Serie which contains noise data */
	LCSerie2 * NoiseData;
	
	/** \brief Type of interpolation */
	INTERP InterpType;
	
	/** \brief Value used in interpolation */
	double InterpUtilValue;
	
	
	/** \brief Name of the noise refering the localization of the noise. ex : pm1s, c2, ... */ 
	char Name[4];
	
	/*! \brief  Index of spacecraft corresponding to the noise */
	int iSC; 
	
	/*!\brief  Direction flag of the optical bench containing the noise */
	int IndirectDir;
	
	/** \brief Number of data to add at each run step (measurement time step) \f$ N_{2Add} = \Delta t_{Mes} / \Delta t_{Phy} \f$  */
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
	LCNoise();
	
	/*! \brief Standard constructor */
	LCNoise(LCTools * MT);
	
	
	/*! \brief Destructor */
	virtual ~LCNoise();
	
	
	/*! \brief Initialization at NULL of all base of noise */
	void initNULLBase(bool CleanMem);
	
	/***************  Required methods (to be present in derived classes)  ***************/
	
	/*! \brief Initialization at NULL */
	virtual void initNULL(bool CleanMem);
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="Noise"> ... </XSIL> */
	virtual void config(ezxml_t noisexmlbloc);
	
	/*! \brief Configuration of individual parameter */
	virtual void config(int iParam, double ParamVal);
	
	/*! \brief Initialization */
	virtual void init();
	
	/*! \brief Generation of noise corresponding from bin iStartBin until bin 0 */
	virtual void generNoise(int iStartBin);
	
	/*! \brief Display information about the noise */ 
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
	
	/*! \brief Set the name of the noise */
	void setName(const char * NamePos);
	
	void setFreqOutput() {IsFreqOutput = true;};
	void setPhaseOutput() {IsFreqOutput = false;};
	
	/*! \brief Compare the nchar first characters of Name_cmp with the Name and return true if there are the same. */
	bool CmpName(const char * Name_cmp, int nchar);
	
	int getiSC() {return(iSC);};
	int getIndirectDir() {return(IndirectDir);};
	int getInterpUtilValue() {return(InterpUtilValue);};
	
	void setInterpolation(INTERP InterpType_n, double InterpUtilValue_n);
	
	
	/**********  Running methods  **********/
	
	/*! \brief Run one measurement step : generation of noise */
	void RunStep();
	
	/*! \brief Return the noise value corresponding to the delay relative to the global time of the simulation */
	double gN(double tDelay); 

	/*! \brief Return the noise value corresponding to the bin */
	double gB(int iBin);
	
	
	/**********  Others methods  **********/
	
	/*! \brief Display the basis */
	void DispInfoBase(char * BTab);
	
	void DispData() {NoiseData->DispData();};
	
};

#endif //__LCNOISE_H

/**\}*/

// end of LISACODE-Noise.h