/*
 *  LISACODE-OrbitsOctaAnalytic.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 10/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \ingroup Orbits Orbits
 * \{
 */


#ifndef __LCORBITSANALYTICOCTAHEDRON_H
#define __LCORBITSANALYTICOCTAHEDRON_H

#include "LISACODE-Orbits.h"
#include "ezxml.h"


/** \brief Class computing the analytic orbits of spacecrafts.
 * \author A. Petiteau
 * \version 2.0
 * \date 10/04/2011
 *
 * This class computes the analytic orbits of spacecrafts.
 * 
 *
 */
class LCOrbitsAnaOctahedron : public LCOrbits
{
protected:
	
	/** \brief Distance between the Octahedron barycenter and the Sun. */
	double Rgc;
	
	/** \brief Angular velocity of the Octahedron barycenter around the Sun. */
	double omega;
	
	/** \brief Offset on spacecraft relative to the nominal position, i.e. deformation factor : 0 = equal arm */
	double pSCOff[6];
	
	
	/************  Internal variables  ************/
	
	/** \brief Pre-computed L*sqrt(2)/2 */
	double Lo2, Lsq2o2;
	
	/** \brief Position of spacecraft relative to L1 */ 
	LCVector posRL1[6];
	
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCOrbitsAnaOctahedron();
	
	/*! \brief Standard constructor */
	LCOrbitsAnaOctahedron(LCTools * MT);
	
	
	/*! \brief Destructor */
	~LCOrbitsAnaOctahedron();
	
	
	/***************  Required methods (to be present in derived classes)  ***************/
	
	/*! \brief Initialization at NULL */
	void initNULL(bool CleanMem);
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="Orbits"> ... </XSIL> */
	void config(ezxml_t noisexmlbloc);
	
	/*! \brief Configuration of individual parameter */
	void config(int iParam, double ParamVal);
	
	/*! \brief Initialization */
	void init();
	
	/*! \brief Specific computation of the travel time  */
	double ArmSpecific(int em, int rec, double trec);
	
	/*! \brief Returns the position of the spacecraft in the barycentric frame for the time t (s) as argument and spacecraft number (1, 2 or 3) */
	LCVector position(int iSC, double t);
	
	/*! \brief  Returns the position of the spacecraft in the barycentric frame for the time t (s) as argument and spacecraft number (1, 2 or 3) */
	LCVector velocity(int iSC, double t);
	
	/*! \brief Display information */ 
	void DispInfo(char * BTab);
	
	/***************  Local methods  ***************/
	
	/**********  Configuration methods  **********/
	
	/**********  Linking and initalization methods  **********/
	
	/**********  Access methods  **********/
	
	/**********  Running methods  **********/
	
	/**********  Others methods  **********/
	
	
};

#endif //__LCORBITSANALYTICOCTAHEDRON_H

/**\}*/

// end of LISACODE-OrbitsAnaOctahedron.h