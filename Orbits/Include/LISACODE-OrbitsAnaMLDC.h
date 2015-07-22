/*
 *  LISACODE-OrbitsAnalytic.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 10/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \ingroup Orbits Orbits
 * \{
 */


#ifndef __LCORBITSANAMLDC_H
#define __LCORBITSANAMLDC_H

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
class LCOrbitsAnaMLDC : public LCOrbits
{
protected:
	
    /*!\brief 0 for staticLISA, 1 for rigid orbits and 2 for eccentric orbits */
    double Approx;
    
    /*!\brief Eccentricity */
	double e_mldc;
	
    /*!\brief Armlength */
	double L0s;
	
	
	
	/* \brief Satellites phase. */	
	std::vector<double> rot;
	
	/*! \brief Satellites phase cosinus. */	
	std::vector<double> crot;
	
	/*! \brief Satellites phase sinus. */	
	std::vector<double> srot;	
	
	
	/************ Variables use in initialization but not after ********** */
	
    /*! \brief Phase between satellites : constant,  \f$ arot=\frac{2 \cdot \pi}{3} \f$. */	
	double arot;
    
	/** \brief Initial rotation (radians) */
	double rot0;
	
	/** \brief Distance between the LISA barycenter and the Sun. */
	double Rgc;
	
	/** \brief Angular velocity of the LISA barycenter around the Sun. */
	double omega;
    
	/** \brief Constant */
    double sqrt_3;

    /** \brief Constant */
    double i32;

    /** \brief Constant */
    double r1532;
    
    /** \brief Constant */
    double pi3o2;
    
    /** \brief Constant */
    double pi2o3;
    
    /** \brief Constant */
    double pi4o3;
	
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCOrbitsAnaMLDC();
	
	/*! \brief Standard constructor */
	LCOrbitsAnaMLDC(LCTools * MT);
	
	
	/*! \brief Destructor */
	~LCOrbitsAnaMLDC();
	
	
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

#endif //__LCORBITSANAMLDC_H

/**\}*/

// end of LISACODE-OrbitsAnaMLDC.h