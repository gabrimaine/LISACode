/*
 *  LISACODE-USONoiseDeriv.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 12/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \ingroup USO 
 * \{
 */


#ifndef __LCUSONOISEDERIV_H
#define __LCUSONOISEDERIV_H

#include "LISACODE-Noise.h"
#include "LISACODE-NoiseWhite.h"
#include "LISACODE-NoiseFilter.h"
#include "ezxml.h"
#include "LISACODE-USO.h"


/** \brief Derived class modeling a clock (USO) with a noise + a constant deriv.
 * \author A. Petiteau
 * \version 2.0
 * \date 12/04/2011
 *
 * This class is the derived class for modeling a clock (USO) with a noise + a constant deriv.
 * 
 *
 */
class LCUSONoiseDeriv : public LCUSO
{
protected:
	
	/*! \brief Time Offset */
	double tOffset;
	
	/*! \brief Time slope : Linear coefficient of the derivative */
	double tSlope;
	
	/*! \brief Noise */
	LCNoise * noise;
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCUSONoiseDeriv();
	
	/*! \brief Standard constructor */
	LCUSONoiseDeriv(LCTools * MT_n);
	
	/*! \brief Destructor */
	virtual ~LCUSONoiseDeriv();
	
	
	/***************  Required methods (to be present in derived classes)  ***************/
	
	/*! \brief Initialization at NULL */
	void initNULL(bool CleanMem);
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="USO"> ... </XSIL> */
	void config(ezxml_t usodata);
	
	/*! \brief Configuration of individual parameter */
	void config(int iParam, double ParamVal);
	
	/*! \brief Initialization */
	void init();
	
	/*! \brief Specific computation of the travel time  */
	void genertShift(int iStartBin); 
	
	/*! \brief Display information */ 
	void DispInfo(char * BTab);
	
	
	
	/***************  Local methods  ***************/
	
	/**********  Linking and initalization methods  **********/
	
	/**********  Access methods  **********/
	
	/**********  Running methods  **********/
	
	/**********  Others methods  **********/
		
	
};

#endif //__LCUSONOISEDERIV_H

/**\}*/

// end of LISACODE-USONoiseDeriv.h