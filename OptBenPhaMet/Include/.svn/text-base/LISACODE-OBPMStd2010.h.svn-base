/*
 *  LISACODE-OBPMStd2010.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 15/11/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \ingroup OBPM
 * \{
 */

#ifndef __LCOBPMSTD2010_H
#define __LCOBPMSTD2010_H

#include "LISACODE-Tools.h"
#include "ezxml.h"
#include "LISACODE-OBPM.h"


/** \brief Class for managing the measurement system (optical bench + photodiode + phasemeter) based on standard LISA design 2010.
 * \author A. Petiteau
 * \version 2.0
 * \date 13/04/2011
 *
 * This class is the derived class for managing the measurement system (optical bench + photodiode + phasemeter) based on standard LISA design 2010.
 * The design and related measurements correspond to the ones described in M. Otto & al. paper (2011)
 * 
 *
 */
class LCOBPMStd2010 : public LCOBPM
{
protected:
	
	/** \brief Type of optical bench path + phasemeter 
	 *	0 : 'scl'  : science measurement including clock noise : \f$ s_i^c  \f$     : incoming beam - local beam
	 *	1 : 'ssb' : science sideband measurement               : \f$ s_i^sb \f$     : incoming sideband beam - local sideband beam
	 *	2 : 'tau'   : test mass measurement                    : \f$ \tau_i \f$     : join optical bench bounced on local test mass beam - local beam 
	 *	3 : 'eps'   : reference measurement                    : \f$ \epsilon_i \f$ : join optical bench beam - local beam 
	 */
	int TypeOBPM;
	
	
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCOBPMStd2010();
	
	/*! \brief Standard constructor */
	LCOBPMStd2010(LCTools * MT);
	
	
	/*! \brief Destructor */
	~LCOBPMStd2010();
	
	
	/***************  Required methods (to be present in derived classes)  ***************/
	
	/*! \brief Initialization at NULL */
	void initNULL(bool CleanMem);
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="Orbits"> ... </XSIL> */
	void config(ezxml_t xmlbloc);
	
	/*! \brief Configuration of individual parameter */
	void config(int iParam, double ParamVal);
	
	/*! \brief Link the noises */ 
	void LinkNoise(LCNoise ** AllNoises, int NAllNoises);
	
	/*! \brief Link the factor to apply on shot noise times delay (in telescope) */ 
	void LinkFactShotNoise(double * FactShotNoise_n);
	
	/*! \brief Initialization */
	void init();
	
	/*! \brief Make a measurements using :
	 *	@param[in] TimeLocal Local time (time when the clock "decides" to take a measurements)
	 *	@param[in] TShiftEff Shift between global time and USO time
	 */ 
	double MeasurePho(double TimeLocal, double TShiftEff);
	
	/*! \brief Display information */ 
	void DispInfo(char * BTab, bool LinkedModule);
	
	/***************  Local methods  ***************/
	
	
	
	
};

#endif //__LCOBPMSTD2010_H

/**\}*/

// end of LISACODE-OBPMStd2010.h