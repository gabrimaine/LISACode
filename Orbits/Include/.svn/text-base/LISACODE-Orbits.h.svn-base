/*
 *  LISACODE-Orbits.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 10/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \defgroup Orbits Orbits
 * This group contains all the things related to orbits.
 * \{
 */


#ifndef __LCORBITS_H
#define __LCORBITS_H

#include "LISACODE-Tools.h"
#include "LISACODE-Vector.h"
#include "ezxml.h"

#define NSCMAX 6
#define NARMMAX 30

/** \brief Base class for managing the orbits of spacecrafts.
 * \author A. Petiteau
 * \version 2.0
 * \date 10/04/2011
 *
 * This class is the base class for managing the orbits of spacecrafts.
 * 
 *
 */
class LCOrbits
{
protected:
	
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	
	/** \brief Number of spacecraft (Max NSC = 6 --> size of #SCposStore)*/
	int NSC;
	
	/** \brief Moving (1) or static(0) orbits   */
	int move;
	
	/** \brief Initial time (seconds) */
	double t0;
	
	/** \brief  Nominal distance between two satellites (meter) */
	double L0m;
	
	
	/** Order for computation of time travel along arm
	 * > 0 : Compute the time travel based on a spacecraft position
	 * -  0 : order 0 
	 * -  1 : order 1/2 
	 * -  2 : order 1 
	 * < 0 : use specific function for computing the travel time (defined in derived class) : tArmSpecific
	 */
	int OrderArm;
	
	/****** For keeping last position in memory */
	
	/** \brief Stored positions : size = max of NSC */	
	LCVector SCposStore[NSCMAX];
	
	/** \brief last position computation time. */	
	double tStorePos; 
	
	/** \brief last position storage time. */	
	double tRangeStorePos;
	
	
	/****** For keeping last armlength in memory */
	
	/** \brief Stored time travel for each one of the arms : size = maxNSC*(maxNSC-1) .
	 *
	 * Each arm corresponds to a pair of emitter (em) and receiver (rec):
	 * For the case with 3 spacecrafts :
	 * - arm 2 : arm 2  : (em,rec)=(1,3)
	 * - arm 3 : arm 3  : (em,rec)=(2,1)
	 * - arm 1 : arm 1  : (em,rec)=(3,2)
	 * - arm 6 : arm 3' : (em,rec)=(1,2)
	 * - arm 5 : arm 2' : (em,rec)=(2,3)
	 * - arm 4 : arm 1' : (em,rec)=(3,1)
	 * For the general case (number of spacecrafts different 3):
	 * - arm i : em=1+floor(i/(NSC-1)) , rec = i%(NSC-1) + (i%(NSC-1)>=em)
	 * for example with 6 spacecrafts :
	 * - arm  1 : (em,rec)=(1,2)
	 * - arm  2 : (em,rec)=(1,3)
	 * - arm  3 : (em,rec)=(1,4)
	 * - arm  4 : (em,rec)=(1,5)
	 * - arm  5 : (em,rec)=(1,6)
	 * - arm  6 : (em,rec)=(2,1)
	 * - arm  7 : (em,rec)=(2,3)
	 * - arm  8 : (em,rec)=(2,4)
	 * - arm  9 : (em,rec)=(2,5)
	 * - arm 10 : (em,rec)=(2,6)
	 * - arm 11 : (em,rec)=(3,1)
	 * - arm 12 : (em,rec)=(3,2)
	 
	 */
	double ArmStore[NARMMAX];
	
	/** \brief Last time travel computation time. */	
	double tStoreArm; 
	
	/** \brief Precision range on the delay computation. */	
	double tRangeStoreArm;
	 
	
	
	
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCOrbits();
	
	/*! \brief Standard constructor */
	LCOrbits(LCTools * MT);
	
	
	/*! \brief Destructor */
	virtual ~LCOrbits();
	
	/*! \brief Initialization at NULL of all base of noise */
	void initNULLBase(bool CleanMem);
	
	/***************  Required methods (to be present in derived classes)  ***************/
	
	/*! \brief Initialization at NULL */
	virtual void initNULL(bool CleanMem);
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="Orbits"> ... </XSIL> */
	virtual void config(ezxml_t noisexmlbloc);
	
	/*! \brief Configuration of individual parameter */
	virtual void config(int iParam, double ParamVal);
	
	/*! \brief Initialization */
	virtual void init();
	
	/*! \brief Display information */ 
	virtual void DispInfo(char * BTab);
	
	/*! \brief Specific computation of the travel time  */
	virtual double ArmSpecific(int em, int rec, double trec);
	
	/*! \brief Returns the position of the spacecraft in the barycentric frame for the time t (s) as argument and spacecraft number (1, 2 or 3) */
	virtual LCVector position(int iSC, double t);
	
	/*! \brief  Returns the position of the spacecraft in the barycentric frame for the time t (s) as argument and spacecraft number (1, 2 or 3) */
	virtual LCVector velocity(int iSC, double t);
	
	
	/***************  Local methods  ***************/
	
	/**********  Configuration methods  **********/
	
	void configBase(ezxml_t xmlbloc);
	
	/**********  Linking and initalization methods  **********/
	
	/*! \brief Initialization of base of noise */
	void initBase();
	
	
	/**********  Access methods  **********/
	
	/*! \brief Return the nominal travel time (in second) */
	double getNominalArm() { return(L0m/LC::c_SI); } ; 
	
	/*! \brief Return the number of spacecraft */
	int getNSC() {return(NSC);};
	
	/*! \brief Return the number of arm */
	int getNArm() {return(NSC*(NSC-1));};
	
	/**********  Methods for time travel along arm  **********/
	
	/*! \brief Compute the travel time :
	 *	@param[in] emi	Index of the emitter (1, 2 or 3)
	 *	@param[in] rec	Index of the receiver (1, 2 or 3)
	 *	@param[in] trec	Time at the receiver
	 */
	double ArmCompute(int em, int rec, double trec);
	
	/*! \brief  Compute the contribution of the specified order to the travel time :
	 *	@param[in] emi	 Index of the emitter (1, 2 or 3)
	 *	@param[in] rec	 Index of the receiver (1, 2 or 3)
	 *	@param[in] trec	 Time at the receiver
	 *	@param[in] order Specified order
	 */
	double ArmOrderContrib(int order, double tij, LCVector ri, LCVector rj, LCVector vi, LCVector vj, LCVector rij, LCVector nij, LCVector n);
	
	/*! \brief Return the time travel of the light along the specified arm (recompute the value if needed)
	 *	@param[in] emi	Index of the emitter (1, 2 or 3)
	 *	@param[in] rec	Index of the receiver (1, 2 or 3)
	 *	@param[in] trec	Time at the receiver
	 */
	double Arm(int em, int rec, double trec);
	
	/*! \brief Return the time travel of the light along the specified arm (recompute the value if needed)
	 *	@param[in] emi	Index of the emitter (1, 2 or 3)
	 *	@param[in] rec	Index of the receiver (1, 2 or 3)
	 *	@param[in] trec	Time at the receiver
	 */
	LCVector unitArm(int em, int rec, double trec);

	/**********  Running methods  **********/
	
	/*! \brief Return the position of the spacecraft iSC (recompute the value if needed) */
	LCVector Pos(int iSC, double t);
	
	/*! \brief  Return spacecrafts relative velocity along an arm : can be used for computing the Doppler*/
	double ArmVelocity(int em, int rec, double trec);
	
	/*! \brief  Return  vector normal to the constellation */
	LCVector VectNormal(double t);
	
	/**********  Others methods  **********/
	
	/*! \brief  Estimate the minimal and maximal at which one we will have to compute GW components (for GW time optimization) */
	void tGWMinMax(double GWbet, double GWlam, double t0_l, double tDur_l, double dt_l,  double & tGWmin, double & tGWmax);
	
	/*! \brief Display informations */
	void DispInfoBase(char * BTab);
};

#endif //__LCORBITS_H

/**\}*/

// end of LISACODE-Orbits.h