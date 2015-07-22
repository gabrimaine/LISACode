/*
 *  LISACODE-TDIIntStd2002.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 20/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \ingroup TDI 
 * \{
 */

#ifndef __LCTDIINTSTD2002_H
#define __LCTDIINTSTD2002_H

#include "LISACODE-Tools.h"
#include "ezxml.h"
#include "LISACODE-Serie2.h"
#include "LISACODE-Detector.h"
#include "LISACODE-TDIInt.h"

/** \brief Class for computing intermediate Time Delay Interferometry in standard 2002.
 * \author A. Petiteau
 * \version 2.0
 * \date 21/04/2011
 *
 * This class computes intermediate Time Delay Interferometry in standard 2002.
 * 
 *
 */
class LCTDIIntStd2002 : public LCTDIInt 
{
protected:
	
	
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Standard constructor */
	LCTDIIntStd2002(LCTools * MT_n);
	
	/*! \brief Destructor */
	~LCTDIIntStd2002();
	
	
	/***************  Required methods (to be present in derived classes)  ***************/
	
	/*! \brief Initialization at NULL */
	void initNULL(bool CleanMem);
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="TDIIntmediate"> ... </XSIL> */
	void config(ezxml_t xmlbloc);
	
	/*! \brief Link to detector : link the signal and maybe other stuff */
	void LinkDetector(LCDetector * LISA);

	/*! \brief Link to the delays */
	void LinkDelays(LCSerie2** AllDelays, int NAllDelays);

	/*! \brief Initialization */
	void init();
	
	/*! \brief Return the maximum number of combined delay */
	int getNMaxDelay();
	
	/*! \brief Running one step */
	void RunStep(double t);
	
	/*! \brief Display informations */
	void DispInfo(char * BTab);
	
	
	/***************  Local methods  ***************/
		
	
	/**********  Other methods  **********/
	
	
};

#endif //__LCTDIINTSTD2002_H

/**\}*/

// end of LISACODE-TDIIntStd2002.h