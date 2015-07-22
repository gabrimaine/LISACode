/*
 *  LISACODE-GWFile.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 23/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


/** \ingroup GW
 * \{
 */


#ifndef __LCGWFILE_H
#define __LCGWFILE_H

#include "LISACODE-Tools.h"
#include "ezxml.h"
#include "LISACODE-GW.h"
#include "LISACODE-DataFileRead.h"

/** \brief Class for managing gravitational waves read in a file.
 * \author A. Petiteau
 * \version 2.0
 * \date 23/04/2011
 *
 * This class manages graviational waves read in a file.
 * 
 *
 */
class LCGWFile : public LCGW
{
protected:
	
	/*! \brief Pointer on the data file class */
	LCDataFileRead * fDat;
	
	/** \brief Type of interpolation */
	INTERP InterpType;
	
	/** \brief Value used in interpolation */
	double InterpUtilValue;
	
	
	/** ********** Variables use in initialization but not after ********** */
	
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCGWFile();
	
	/*! \brief Standard constructor */
	LCGWFile(LCTools * MT_n);
	
	
	/*! \brief Destructor */
	~LCGWFile();
	
	
	/***************  Required methods (to be present in derived classes)  ***************/
	
	/*! \brief Initialization at NULL */
	void initNULL(bool CleanMem);
	
	/**********  Configuration methods  **********/
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="GWFile"> ... </XSIL> */
	void config(ezxml_t xmlbloc);
	
	/*! \brief Configuration of individual parameter */
	void config(int iParam, double ParamVal);
	
	
	/**********  Access methods  **********/
	
	/*! \brief Return value of parameter @iP :
	 * - 0 : Beta
	 * - 1 : Lambda
	 */
	double getParam(int iP);
	
	/*! \brief Set value of parameter @iP at @Param_n */
	void setParam(int iP, double Param_n);
	
	/*! \brief Set \Pmin and \Pmax at range of parameter @iP */
	void getRange(int iP, double &Pmin, double &Pmax);
	
	/*! \brief Return small variation delta used in computation of derivative for parameter @iP */
	double getDelta(int iP);
	
	/*! \brief Set value of specific parameter @iPS at @SpecParam_n */
	void setSpecialParam(int iPS, double SpecParam_n);
	
	/*! \brief Get value of specific parameter @iPS */
	double getSpecialParam(int iPS);
	
	/*! \brief Set time informations : see GW::setTimeInfo(double,double,double,double,double) for details */
	void setTimeInfo(double t0_n, double dt_n, double TObs_n, double tAskMin_n, double tDOrbMax_n,  double tMaxDiff);
	
	/*! \brief Choose randomly parameter iP */
	void RandParam(int iP);
	
	
	/**********  Linking and initalization methods  **********/
	
	/*! \brief Initialization 
	 @	return error code (0 if everything ok, 1 if there is a problem) 
	 */
	int init();
	
	
	/**********  Running methods  **********/
	
	/*! \brief Compute gravitational strain in barycentric reference frame.  
	 * If polarization is null, computation in barycentric reference frame.
	 * Else computation in source reference frame
	 */
	void Computehpc(double t); 
	
	
	/**********  Others methods  **********/
	
	/*! \brief Display information about the GWFile */ 
	void DispInfo(char * BTab);
	
	/*! \brief Display all parameters */
	void DispAllParam(std::ostream * out);
	
	/*! \brief Display name of all parameters */
	void DispAllParamName(std::ostream * out);
	
	
	
	
	
	/***************  Local methods  ***************/
	
	
	/**********  Configuration methods  **********/
	
	
	/**********  Linking and initalization methods  **********/
	
	
	/**********  Access methods  **********/
	
	
	/**********  Running methods  **********/
	
	
	/**********  Others methods  **********/
	
	
	
};

#endif //__LCGWFILE_H

/**\}*/

// end of LISACODE-GWFile.h