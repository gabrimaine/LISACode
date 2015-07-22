/*
 *  LISACODE-LCOrbitsData.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 10/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \ingroup Orbits Orbits
 * \{
 */


#ifndef __LCORBITSDATA_H
#define __LCORBITSDATA_H

#include "LISACODE-Orbits.h"
#include "LISACODE-DataFileRead.h"
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
class LCOrbitsData : public LCOrbits
{
protected:
	
	/*! \brief Pointer on the data file class */
	LCDataFileRead * fDat;
	
	/** \brief Type of interpolation */
	INTERP InterpType;
	
	/** \brief Value used in interpolation */
	double InterpUtilValue;
	
	
	/** \brief Table containing index of column for the position components with [iSC-1][x->0,y->1,z->2] */
	int iColPos[3][3];
	
	/** \brief Table containing index of column for the velocity components with [iSC-1][x->0,y->1,z->2] */
	int iColVel[3][3];
	
	/** \brief True if the velocity is also read in file */
	bool VelocityInFile;
	
	/** \brief Factor to apply on the read time to matched the unit of the simulation (second) */
	double ConvTime;
	
	/** \brief Factor to apply on the read position to matched the unit of the simulation (meter) */
	double ConvPos;
	
	/** \brief Factor to apply on the read velocity to matched the unit of the simulation (meter/second) */
	double ConvVel;
	
	/** \brief Time shift between the beginning of the data and the beginning of the simulation */
	double TimeShift;
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCOrbitsData();
	
	/*! \brief Standard constructor */
	LCOrbitsData(LCTools * MT);
	
	
	/*! \brief Destructor */
	~LCOrbitsData();
	
	
	/***************  Required methods (to be present in derived classes)  ***************/
	
	/*! \brief Initialization at NULL */
	void initNULL(bool CleanMem);
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="Orbits"> ... </XSIL> */
	void config(ezxml_t noisexmlbloc);
	
	/*! \brief Configuration of individual parameter */
	void config(int iParam, double ParamVal);
	
	/*! \brief Initialization */
	void init();
	
	/*! \brief Display information */ 
	void DispInfo(char * BTab);
	
	/*! \brief Specific computation of the travel time  */
	double ArmSpecific(int em, int rec, double trec);
	
	/*! \brief Returns the position of the spacecraft in the barycentric frame for the time t (s) as argument and spacecraft number (1, 2 or 3) */
	LCVector position(int iSC, double t);
	
	/*! \brief  Returns the position of the spacecraft in the barycentric frame for the time t (s) as argument and spacecraft number (1, 2 or 3) */
	LCVector velocity(int iSC, double t);
	
	
	/***************  Local methods  ***************/
	
	/**********  Configuration methods  **********/
	
	/**********  Linking and initalization methods  **********/
	
	/**********  Access methods  **********/
	
	/** \brief Set the name of the input file */
	void setFile(char * NewFileName);
	
	/**********  Running methods  **********/
	
	/**********  Others methods  **********/
	
	
	
	
	
};

#endif //__LCORBITSDATA_H

/**\}*/

// end of LISACODE-LCOrbitsData.h