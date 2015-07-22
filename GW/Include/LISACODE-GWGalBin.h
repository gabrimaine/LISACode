/*
 *  LISACODE-GWGalBin.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 06/05/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


/** \ingroup GW
 * \{
 */


#ifndef __LCGWGALBIN_H
#define __LCGWGALBIN_H

#include "LISACODE-Tools.h"
#include "ezxml.h"
#include "LISACODE-GW.h"

/** \brief Class for managing gravitational waves emitted by galactic binaries.
 * \author A. Petiteau
 * \version 2.0
 * \date 23/04/2011
 *
 * This class manages gravitational waves emitted by galactic binaries.
 * The waveform is :
 *	\f[ h^S_{+}(t) = A (1+ \cos^2 \iota ) \cos [ 2 \pi (f t + \dot{f} t^2 / 2) + \phi_0 ] = h{+,0} \cos [ \omega t + \dot{\omega} t^2 ] \f]
 *	\f[ h^S_{\times}(t) = - 2 A (\cos \iota) \sin [ 2 \pi (f t + \dot{f} t^2 / 2) + \phi_0 ] = h{\times,0} \sin [ \omega t + \dot{\omega} t^2 ] \f]
 *
 */
class LCGWGalBin : public LCGW
{
protected:
	
	/*! \brief Inclination \f$ iota \f$ (in radian) */
	double inc;
	
	/*! \brief Amplitude \f$ A \f$ (in radian)*/
	double Amp;
	
	/*! \brief Frequency \f$ f \f$ (in Hz)*/
	double Freq;
	
	/*! \brief First time derivative of frequency \f$ \dot{f} \f$ (in 10^-16 Hz /s)*/
	double DFreq;
	
	/*! \brief Initial phase \f$ \phi_0 \f$ (in radian) */
	double Phi0;
	
	
	/** ********** Variables use in initialization but not after ********** */
	
	/*! \brief Mass 1 \f$ m_1 \f$ (in solar mass) : used for computing amplitude #Amp */
	double m1;
	
	/*! \brief Mass 2 \f$ m_2 \f$ (in solar mass) : used for computing amplitude #Amp */
	double m2;
	
	/*! \brief Luminosity distance \f$ D_L \f$ (in kpc) : used for computing amplitude #Amp */
	double DL;
	
	/*! \brief Eccentricity \f$ D_L \f$ (in kpc) : not use at the moment */
	double ecc;
	
	/*! \brief First period time derivative \f$ D\dot{P} \f$ (in kpc) : used for computing first time derivative of frequency  #DFreq */
	double DPeriod;
	
	/*! \brief True if amplitude is directly set : Avoid any recomputation of amplitude */
	bool DontChangeAmp;
	
	/*! \brief True if first time derivative of frequency is directly set : Avoid any recomputation */
	bool DontChangeDFreq;
	
	/** ********** Variables computed during the initialization ********** */
	
	/*! \brief Amplitude of \f$ h_+ \f$ , \f$ h{+,0} = A (1+ \cos^2 \iota) \f$ */
	double hp0;
	
	/*! \brief Amplitude of \f$ h_{times} \f$ , \f$ h{\times,0} = - 2 A (\cos \iota) \f$ */
	double hc0;
	
	/*! \brief Pulsation : \f$ \omega = 2 \pi f \f$ */
	double om;
	
	/*! \brief Derivative of pulsation : \f$ \dot{\omega} = \pi \dot{f} \f$ */
	double Dom;
	
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCGWGalBin();
	
	/*! \brief Standard constructor */
	LCGWGalBin(LCTools * MT_n);
	
	
	/*! \brief Destructor */
	~LCGWGalBin();
	
	
	/***************  Required methods (to be present in derived classes)  ***************/
	
	/*! \brief Initialization at NULL */
	void initNULL(bool CleanMem);
	
	/**********  Configuration methods  **********/
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="GWSpinBBH2"> ... </XSIL> */
	void config(ezxml_t xmlbloc);
	
	/*! \brief Configuration of individual parameter */
	void config(int iParam, double ParamVal);
	
	
	/**********  Access methods  **********/
	
	/*! \brief Return value of parameter 
	 *	@param[in]	iP : Index of parameter
	 *	@return		Value of parameter
	 *	Coding for index iP :
	 *	- 0  : \f$ \beta \f$	: #Beta			: Ecliptic latitude (unit: rad)
	 *	- 1  : \f$ \lambda \f$  : #Lambda		: Ecliptic longitude (unit: rad)
	 *	- 2  : \f$ \psi \f$		: #Polarization	: Polarization (unit: rad)
	 *	- 3  : \f$ \iota \f$	: #inc			: Inclination (unit: rad)
	 *	- 4  : \f$ A \f$		: #Amp			: Amplitude
	 *	- 5  : \f$ f \f$		: #Freq			: Frequency (unit: Hz/s)
	 *	- 6  : \f$ \dot{f} \f$	: #DFreq		: First time derivative of frequency (unit: 10^-16 Hz/s)
	 *	- 7  : \f$ \phi0 \f$	: #Phi0			: Initial phase (unit: rad)
	 *	- 8  : \f$ m_1 \f$		: #m1			: Mass1 (used for computing amplitude ;  unit: solar mass)
	 *	- 9  : \f$ m_2 \f$		: #m2			: Mass2 (used for computing amplitude ;  unit: solar mass)
	 *	- 10 : \f$ D_L \f$		: #DL			: Luminosity distance (used for computing amplitude ;  unit: kpc)
	 *	- 11 : \f$ ecc \f$		: #ecc			: Eccentricity (not use at the moment)
	 *	- 12 : \f$ \dot{P} \f$	: #DPeriod		: First time derivative of period (unit: s/s)
	 */
	double getParam(int iP);
	
	/*! \brief Set value of parameter  
	 *	@param[in]	iP		: Index of parameter (see #getParam for coding)
	 *	@param[in]	param_n :Value of parameter
	 */
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
	
	/*! \brief Display informations */
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

#endif //__LCGWGALBIN_H

/**\}*/

// end of LISACODE-GWSpinBBH2.h