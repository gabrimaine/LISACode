/*
 *  LISACODE-NoiseFile.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 14/03/13.
 *  Copyright 2013 Laboratory AstroParicles and Cosmology - University Paris-Diderot (Paris, France). All rights reserved.
 *
 */



/** \ingroup Noise Noise
 * \{
 */


#ifndef __LCNOISEFILE_H
#define __LCNOISEFILE_H

#include "LISACODE-Noise.h"
#include "LISACODE-Filter.h"
#include "LISACODE-DataFileRead.h"

/** \brief Class for generating noise from a file
 * \author A. Petiteau
 * \version 2.0
 * \date 08/04/2011
 *
 * This class generates noise from data read in a file. 
 * This data can be filtered for convertion, for example to convert displacement noise in relative frequency noise.
 * 
 *
 */
class LCNoiseFile : public LCNoise
{
protected:
    
    /************* For managing input data ************/
    
    /*! \brief Pointer on the data file class */
	LCDataFileRead * fDat;
	
	/** \brief Type of interpolation */
	INTERP InterpType;
	
	/** \brief Value used in interpolation */
	double InterpUtilValue;
    
    /** \brief Time in data reading */
    double tReadData;

    /** \brief Starting time in data reading (consider points needed for interpolation) */
    double tStartReadData;
    
    /** \brief Maximal time in data reading (consider points needed for interpolation) */
    double tMaxReadData;
    
    /** \brief Index of record to read : 1 to maximum numbre of records */
    int IndexRecToRead;
	
    /*************** For filtering *******************/
    
	/*! \brief Raw data : white noise */
	LCSerie2 ** RawData;
	
	/*! \brief Filter */
	LCFilter ** Filter;
	
	/*! \brief Filtered data (used when there are more than one filter) */
	LCSerie2 ** TmpFilData;
	
	/*! \brief Number of filter and white noise associated */
	int NFil;
	
	/*! \brief Factor multiplying the input data  */
	double Factor;
	
	
	/************ Variables use in initialization but not after ********** */
	
	/** \brief Number of data for stabilization. */
	int NbDataStab;

	
	/******* The following informations are used for computing the coefficient alpha and beta */
	
	/** \brief Type of Filtering : frequency shape. Coding: 
	 *	\arg  0 : None
	 *	\arg  1 : Filter f (\f$ y(f) = x(f) f^2 \f$)
	 *  \arg  2 : Filter 1/f (\f$ y(f) = x(f) f^{-2} \f$)
	 */
	int TypeFilter;
	
	
	/*! \brief List of specific frequency values used for constructing the filter (used by PreStabLaserNoiseFreq )*/
	std::vector<double> fSpec;
	
	
	/*! \brief Value used for debug/checking */
	double tct;
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCNoiseFile();
	
	/*! \brief Standard constructor */
	LCNoiseFile(LCTools * MT);
	
	
	/*! \brief Destructor */
	~LCNoiseFile();
	
	
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
	
	/**********  Running methods  **********/
	
	/**********  Others methods  **********/
	
	
	
	
};

#endif //__LCNOISEFILE_H

/**\}*/

// end of LISACODE-NoiseFile.h