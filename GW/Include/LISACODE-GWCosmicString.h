/*
 *  LISACODE-GWGalBin.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 06/05/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

//tool maispas toolbox, constantes ->  LC::nomconstante

/** \ingroup GW
 * \{
 */


#ifndef __LCGWCosmicString_H 
#define __LCGWGCosmicString_H

#include "LISACODE-Tools.h" /* outils tf etc */
#include "ezxml.h"
#include "LISACODE-GW.h"

/** \brief Class for managing gravitational waves emitted by Cosmic String's cusps.
 * \author A. Petiteau
 * \version 2.0
 * \date 23/04/2011
 *
 * This class manages gravitational waves emitted by Cosmic String's Cusps.
 * The waveform is :
 *	\f[ h^S_{+}(t) = A (1+ \cos^2 \iota ) \cos [ 2 \pi (f t + \dot{f} t^2 / 2) + \phi_0 ] = h{+,0} \cos [ \omega t + \dot{\omega} t^2 ] \f]
 *	\f[ h^S_{\times}(t) = - 2 A (\cos \iota) \sin [ 2 \pi (f t + \dot{f} t^2 / 2) + \phi_0 ] = h{\times,0} \sin [ \omega t + \dot{\omega} t^2 ] \f]
 *
 */


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! up : a changer !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
class LCGWCosmicString : public LCGW 
{
protected:
       
	
	/** ********** Variables use in initialization but not after ********** */


        // !!!
        //PARAMETRES USED IN THE COMPUTATION OF A (amplitude)
	// !!!

	// mu : masse lineique kg*m-1 valeur moyenne : 6e21 kg*m-1 (masse terre pour L=1km)
	double Gmu;
	
	// Dl: luminosity distance to the source : longueur (Dl= 10**(((m-M)/5)+1) où m & M : magnitude relative et absolue.
	double DL;

	// alpha : l=alpha*t where l is the length of the loop.
	double alpha;

	//flow 
	double flow;

	//temps cosmologique au moment de la formation de la boucle ?? ou au moment de la creation de la cusp plutot ?
	double tcosmo;
	
	// !!! END !!!
    
    /*! \brief Number of bursts */
    int Nb;

	/*! \brief List of bursts' central time of arrival. Size : #Nb */
	double * centralTime;
	
	/*! \brief List of initial amplitude of h+. Size : #Nb */
	double * amplitudeInit;
    
    /*! \brief List of cos(2 x polarization). Size : #Nb */
	double * cos2Pol;
    
    /*! \brief List of -sin(2 x polarization). Size : #Nb */
	double * msin2Pol;
    
    /*! \brief q=1/3 for Cusps and q=2/3 for kinks / */
	double q;
    
    /* \brief Burst type : Cusp or Kink */
    char BurstType[10];
    
    /* \brief Polarization */
    double PolarizationDelta;
    
    /*! \brief Name of input file containing the list of bursts */
    char FileNameBurstList[1024];
	
	/*! \brief minimal and maximal frequencies of the burst / */
	double fhigh;

	/*! \brief time variables describing the observation : t0 and duration of the observation. / */
	double t0Real;

	double dtReal;

	double TObsReal;

	double fminSim;
    double fmaxSim;
	double TObsSim;
	double dtSim;


	
	/*! \brief True if amplitude is directly set : Avoid any recomputation of amplitude / */
	bool DontChangeAmp;

	bool RecomputeAmp;
       
	/*! \brief True if first time derivative of frequency is directly set : Avoid any recomputation /
	bool DontChangeDFreq;
	
	/** ********** Variables computed during the initialization ********** */


	/*! \brief frequencies used in the computation, may differ from the "real" frequency of the burst. / */
	double fmin;

	double fmax;
	
	/*! \brief number of computed points of h_+(t) & h_+(f) / */
	int N;

	int Nf;

	double df;


	// L : Size of the feature that produces the cusp : longueur
	double ParentSize;


	double* hp_time;

	dcomplex* hp_freq;
	
	
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! up : qu'est ce ? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCGWCosmicString();
	
	/*! \brief Standard constructor */
	LCGWCosmicString(LCTools * MT_n);
	
	
	/*! \brief Destructor */
	~LCGWCosmicString();
	
	
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
	
	
	
	
	
};

#endif //__LCGWGALBIN_H

/**\}*/

// end of LISACODE-GWCosmicString.h
