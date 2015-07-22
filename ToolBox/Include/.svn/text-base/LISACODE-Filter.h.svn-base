// $Id:  Exp $
/*
 *  LISACode-Filter.h
 *  LISACode V 2.0
 *
 *  Created on 08/04/11 by Antoine PETITEAU (AEI)
 *  Last modification on 08/04/11 by Antoine PETITEAU (AEI)
 *
 */

/** \ingroup ToolBox
 * \defgroup Filter Class Filter
 * (See class #LCFilter for a detailed description)
 * \{
 */


#ifndef __LCFILTER_H
#define __LCFILTER_H

#include <stdexcept>
#include <fstream>
#include <math.h>
#include <complex>
#include <stdlib.h>
#include <vector>
#include "LISACODE-Tools.h"
#include "LISACODE-Serie2.h"
#include "LISACODE-EllipticFilter.h"

/** \brief Class for managing filter.
 * \author A. Petiteau
 * \version 2.0
 * \date 13/05/2011
 *
 * This class manage filter. It can apply mutliple filters with multiple cells on #LCSerie2 data.
 * 
 *
 */
class LCFilter 
{
protected:
	
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	/** \brief Pointer on the raw data (allocated outside) */
	LCSerie2 * RawData;
	
	/** \brief Pointer on the filtered data (allocated outside) */
	LCSerie2 * FilData;
	
	/** \brief List of temporay data (allocated inside using #RawData format) */
	LCSerie2 ** TmpData;
	
	/** \brief Number of filter = number of intermediate series + 1 */
	int NFil;
	
	/** \brief List of list of alpha coefficient. Size (NIntDat-1) x NAlpha[i] */
	double ** Alpha;
	
	/** \brief List of size of list of alpha coefficient. Size (NIntDat-1) */
	int * NAlpha;
	
	/** \brief List of list of beta coefficient. Size (NIntDat-1) x NBeta[i] */
	double ** Beta;
	
	/** \brief List of size of list of beta coefficient. Size (NIntDat-1) */
	int * NBeta;
	
	/** \brief Number of data for stabilization. */
	int NbDataStab;
	
public:
	
	/********** Constructor **********/
	/** \brief Default constructor */
	LCFilter();
	
	/** \brief Standard constructor */
	LCFilter(LCTools * MT_n);
	
	/** \brief Destructor */
	~LCFilter();
	
	
	
	/********** Initalization methods **********/
	/*! \brief Initialization of the variables at NULL */
	void initNULL(bool CleanMem);
	
	/*! \brief Initialization of filter. 
	 *	@param[in] Alpha : 2D list of alpha coefficients
	 *	@param[in] Beta  : 2D list of beta coefficients 
	 *	@param[in] NbDataStabilization_n : IN : Number of data to compute for the stabilization 
	 */
	void init(std::vector< std::vector<double> > Alpha_n,
			  std::vector< std::vector<double> > Beta_n,
			  int NbDataStabilization_n);
	
	
	/*! \brief Initialization of filter. 
	 *	@param[in] at : attenuation [dB]
	 *	@param[in] fe : oscillations in bandwidth [dB]
	 *	@param[in] fb : low transition frequency [Hz]
	 *	@param[in] fa : high transition frequency [Hz]
	 */
	void init(double at, 
			  double bp, 
			  double fb,
			  double fa);    
	
	/********** Access methods **********/
	void setRawdata(LCSerie2 * PtrRawData) {RawData = PtrRawData; };
	void setFilteredData(LCSerie2 * PtrFilData) {FilData = PtrFilData; };
	void setInOutData(LCSerie2 * PtrRawData, LCSerie2 * PtrFilData) {RawData = PtrRawData; FilData = PtrFilData;};
	int getNbDataStab() { return(NbDataStab); }; 
	int getDepth();
	
	/********** Running methods **********/
	/*! \brief Application of the filter from the bin StartBin until 0 
	 *	@param[in] StartBin : Starting bin for filtering
	 */
	void App(int StartBin);

	/*! \brief Display information */
	void DispInfo(char * BTab );
	
};

#endif // __LCFILTER_H

/**\}*/

// end of LISACODE-LCFilter.h


