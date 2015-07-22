/*
 *  LISACODE-NoiseFilter.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 05/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \ingroup Noise Noise
 * \{
 */


#ifndef __LCNOISEFILTER_H
#define __LCNOISEFILTER_H

#include "LISACODE-Noise.h"
#include "LISACODE-Filter.h"

/** \brief Class for generating filtered noise.
 * \author A. Petiteau
 * \version 2.0
 * \date 08/04/2011
 *
 * This class generates filtered noise. 
 * It can also generate white noise.
 * 
 *
 */
class LCNoiseFilter : public LCNoise
{
protected:
	
	/*! \brief Raw data : white noise */
	LCSerie2 ** RawData;
	
	/*! \brief Filter */
	LCFilter ** Filter;
	
	/*! \brief Filtered data (used when there are more than one filter) */
	LCSerie2 ** TmpFilData;
	
	/*! \brief Number of filter and white noise associated */
	int NFil;
	
	/*! \brief Factor multiplying the output : \f$ \sigma = \sqrt{PSD} / \sqrt{2 \Delta t } \f$  */
	double Sigma;
	
	
	/************ Variables use in initialization but not after ********** */
	
	/** \brief Number of data for stabilization. */
	int NbDataStab;

	
	/******* The following informations are used for computing the coefficient alpha and beta */
	
	/** \brief Type of noise : frequency shape. Coding : 
	 *	\arg -1 : Undefined 
	 *	\arg  0 : White
	 *	\arg  1 : Filter f (\f$ PSD(f) = PSD_0 f^2 \f$)
	 *  \arg  2 : Filter 1/f (\f$ PSD(f) = PSD_0 f^{-2} \f$)
	 *  \arg  3 : Filter with \f$ PSD(f) = PSD_0 {\left({1 \over f} + {f_{knee} \over f^{2}} \right)}^2 \f$
	 *  \arg  4 : Filter f^a (\f$ PSD(f) = PSD_0 f^{a} \f$)
	 *	\arg 10 : White phase = Blue frequency
	 *	\arg 11 : White frequency
	 *	\arg 12 : Red frequency
	 *	\arg 13 : Pink frequency = Pink acceleration	
	 *	\arg 14 : Filter for pre-stabilized laser frequency noise  
	 */
	int TypeNoise;
	
	/** \brief Power spectral density */
	double SqPSD;
	
	/*! \brief Knee frequency (used for pink acceleration, f^a) */
	double fknee;
	
	/*! \brief Low frequency for starting the slope (used for pink acceleration, f^a) */
	double fLow;
	
	/*! \brief Slope : \f$ a \f$ for a filter in , \f$ f^a \f$  (used for pink acceleration, f^a) */
	double slope;
	
	/*! \brief List of specific frequency values used for constructing the filter (used by PreStabLaserNoiseFreq )*/
	std::vector<double> fSpec;
	
	
	/*! \brief Value used for debug/checking */
	double tct;
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCNoiseFilter();
	
	/*! \brief Standard constructor */
	LCNoiseFilter(LCTools * MT);
	
	
	/*! \brief Destructor */
	~LCNoiseFilter();
	
	
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
	
	/*! \brief Computation of coefficient of filter for a filter with shape \f$ f^{-\alpha} \f$ using
	 *	@param[out]	alpha_n		\f$ \alpha \f$ coefficients of output filter
	 *	@param[out]	beta_n		\f$ \beta \f$ coefficients of output filter
	 *	@param[in]	fmin		Minimal frequency for the range where the filter have the good shape
	 *	@param[in]	fmax		Maximal frequency for the range where the filter have the good shape
	 *	@param[in]	fpow		Power of f   --> Shape : \f$ f^{f_{pow}} \f$
	 *	@param[in]	SqrPSDf0	Square root
	 */
	void CoefOofFilter(std::vector< std::vector<double> > & alpha_n, std::vector< std::vector<double> > & beta_n, double fmin, double fmax, double fpow, double SqrPSDf0 );
	
	/**********  Access methods  **********/
	
	/**********  Running methods  **********/
	
	/**********  Others methods  **********/
	
	
	
	
};

#endif //__LCNOISEFILTER_H

/**\}*/

// end of LISACODE-NoiseFilter.h