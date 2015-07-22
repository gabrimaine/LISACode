/*
 *  LISACODE-NoiseWhite.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 05/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \ingroup Noise Noise
 * \{
 */


#ifndef __LCNOISEWHITE_H
#define __LCNOISEWHITE_H

#include "LISACODE-Noise.h"

/** \brief Class for generating white noise.
 * \author A. Petiteau
 * \version 2.0
 * \date 08/04/2011
 *
 * This class generates white noise.
 * 
 *
 */
class LCNoiseWhite : public LCNoise
{
protected:
	
	/*! \brief White noise standard deviation : Root square Power Spectral Density (in unit.Hz-1) over root square 2 */
	double Sigma;
	 
	
	/** ********** Variables use in initialization but not after ********** */
	
	/*! \brief Square root of Power Spectral Density */
	double SqPSD;
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCNoiseWhite();
	
	/*! \brief Standard constructor */
	LCNoiseWhite(LCTools * MT);
	
	/*! \brief Standard constructor */
	LCNoiseWhite(LCTools * MT, double SqPSD);
	
	
	/*! \brief Destructor */
	~LCNoiseWhite();
	
	
	/***************  Required methods  ***************/
	
	/*! \brief Initialization at NULL */
	void initNULL(bool CleanMem);
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="Noise"> ... </XSIL> */
	void config(ezxml_t noisexmlbloc);
	
	/*! \brief Configuration of parameters individually :
	 *	\arg 0 : square root of PSD
	 */
	void config(int iParam, double ParamVal);
	
	/*! \brief Initialization */
	void init();
	
	/*! \brief Generation of noise corresponding from bin iStartBin until bin 0 */
	void generNoise(int iStartBin);
	
	/*! \brief Display information about the noise */ 
	void DispInfo(char * BTab);
	
	
	/***************  Local methods  ***************/
	
	/**********  Configuration methods  **********/
	
	/**********  Linking and initalization methods  **********/
	
	/**********  Access methods  **********/
	
	/*! \brief Return PSD */
	double getPSD();
	
	/*! \brief Return square root of PSD */
	double getSqPSD();
	
	/*! \brief Sets #Sigma(\f$\sigma\f$) by giving the root square of PSD */
	void setSqPSD(double SqPSD);
	
	/**********  Running methods  **********/
	
	/**********  Others methods  **********/
	
	
	
};

#endif //__LCNOISEWHITE_H

/**\}*/

// end of LISACODE-NoiseWhite.h