/*
 *  LISACODE-OBPMOcta1.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 13/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \ingroup OBPM
 * \{
 */


#ifndef __LCOBPMOCTA1_H
#define __LCOBPMOCTA1_H

#include "LISACODE-Tools.h"
#include "ezxml.h"
#include "LISACODE-OBPM.h"


/** \brief Class for managing the measurement system (optical bench + photodiode + phasemeter) for the first test on Octahedron.
 * \author A. Petiteau
 * \version 2.0
 * \date 21/11/2011
 *
 * This class is the derived class for managing the measurement system (optical bench + photodiode + phasemeter) for the first test on Octahedron.
 * 
 *
 */
class LCOBPMOcta1 : public LCOBPM
{
protected:
	
	/** \brief Type of optical bench path + phasemeter 
	 *	0 : sci : scientific measurement : incoming beam - local
	 */
	int TypeOBPM;
	
	
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCOBPMOcta1();
	
	/*! \brief Standard constructor */
	LCOBPMOcta1(LCTools * MT);
	
	
	/*! \brief Destructor */
	~LCOBPMOcta1();
	
	
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

#endif //__LCOBPMOCTA1_H

/**\}*/

// end of LISACODE-OBPMOcta1.h