/*
 *  LISACODE-GWSpinBBHHHarm1.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 16/08/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


/** \ingroup GW
 * \{
 */


#ifndef __LCGWSPINBBHNR1_H
#define __LCGWSPINBBHNR1_H

#include "LISACODE-Tools.h"
#include "ezxml.h"
#include "LISACODE-GW.h"
#include "LISACODE-Vector.h"
#include "LISACODE-Matrix.h"
#include "LISACODE-DataFileRead.h"
#include "LISACODE-GWSpinBBHHHarm1.h"



const int LCGWSpinBBHNR1_ORDERLAGINTER = 10;


/** \brief Class for managing gravitational waves emitted by black hole binaries with spin reading each harmonic in file generated from numerical relativity.
 * \author A. Petiteau
 * \version 2.0
 * \date 16/08/2011
 *
 * This class manages gravitational waves emitted by black hole binaries with spin 
 * reading each harmonic in file generated by numerical relativity or hybrid waveform.
 * The waveform is based on [<a href="Bibliography.html#Arun-2009">Arun &nbsp; 2009</a>]
 *
 */
class LCGWSpinBBHNR1 : public LCGW
{
protected:
	
    /*! \brief Type of output :
     \arg 0 standard : \f$ h_{+,\times} \f$ 
     \arg 1 for tc derivative : \f$ \partial h_{+,\times} / \partial \tau \f$
     \arg 2 for total mass derivative : \f$ \partial h_{+,\times} / \partial \tau \f$
     */
	int OutputType; 
    
    
	// *************** Parameters of system *************** 
	
	// ********** Free parameters ********** 
	
	/*! \brief Total mass (in solar mass) */
	double Mtot;             
	
	/*! \brief Time of coalescence (in seconds) */
	double tc;
	
	/*! \brief Distance between source and detector  (in kiloparsec) */
	double DL;
	
	/*! \brief Polar angle of the detector seen from the source (in radians) */
	double Thd;
	
	/*! \brief Azimuthal angle of the detector seen from the source (in radians) */
	double Phd;
    
    /*! \brief Polar angle of the total angular momentum in SSB (in radians) */
	//  double thetaJ;
	
	/*! \brief Azimuthal angle of the total angular momentum in SSB (in radians) */
	//  double phiJ;
    
    /*! \brief Random polarization angle (in radians) */
	double RandPolarization;
    
    /*! \brief NR polarization angle (in radians) */
	double NRPolarization;
	
	
	// ********** Fixed parameters by simulation ********** 
	
	/*! \brief Mass ratio (>1) */
	double mrat;
	
	/*! \brief Spin of object 1  */
	double chi1;
	
	/*! \brief Spin of object 2 */
	double chi2;
	
	/*! \brief Initial polar angle of spin of object 1 in source reference frame */
	double thS1;

	/*! \brief Initial polar angle of spin of object 2 in source reference frame */
	double thS2;

	/*! \brief Initial azimuthal angle of spin of object 1 in source reference frame */
	double phS1;
	
	/*! \brief Initial azimuthal angle of spin of object 2 in source reference frame */
	double phS2;
	
	/*! \brief Initial phase */
	double Phi0;
	
	/*! \brief Initial polar angle of orbital angular momentum in source reference frame  */
	double thL;
	
	/*! \brief Initial azimuthal angle of orbital angular momentum in source reference frame  */
	double phL;
	
	/*! \brief Cartesian coordinates of diretion of Ln, S1 and S2 in source (NR) reference frame */
	LCVector LnN, S1N, S2N;
	//double S1x, S1y, S1z, S2x, S2y, S2z, LNx, LNy, LNz;
	
	// *************** Rescaled/Derived parameters ***************
	
	/*! \brief Total mass (in second) */
	double M;
	
	/*! \brief symetric mass ratio (in second) */
	double eta;
	
	/*! \brief Distance between source and detector (in seconds) */
	double dist;
	
	/*! \brief Initial polar angle of orbital angular momentum in solar system barycentric (SSB) reference frame */
	double thBL;
	
	/*! \brief Initial azimuthal angle of orbital angular momentum in solar system barycentric (SSB) reference frame */
	double phBL;
	
	/*! \brief Initial polar angle of spin of object 1 in solar system barycentric (SSB) reference frame */
	double thBS1;
	
	/*! \brief Initial polar angle of spin of object 2 in solar system barycentric (SSB) reference frame */
	double thBS2;
	
	/*! \brief Initial azimuthal angle of spin of object 1 in solar system barycentric (SSB) reference frame */
	double phBS1;
	
	/*! \brief Initial azimuthal angle of spin of object 2 in solar system barycentric (SSB) reference frame */
	double phBS2;
	
	
	// *************** For computation of extra parameters ***************
	
	/*! \brief Cartesian coordinates of Ln, S1 and S2 in source reference frame */
	LCVector LnB, S1B, S2B;
	//double S1Bx, S1By, S1Bz, S2Bx, S2By, S2Bz, LnBx, LnBy, LnBz;
	
	/*! \brief True if extra parameters (as thLB, phLB, etc) have to be compute */
	bool NeedExtraParamCompute;
	
	/*! \brief Type of extra parameter that we have to compute : 1->(#thLB, #phLB, #thS1, #phS1, #thS2 and #phS2) or 2->(#Polarization, #Thd and #Phd)*/
	int TypeExtraParamCompute;
	
	/*! \brief Mass of BH 1 (biggest black hole) */
	double m1;
	
	/*! \brief Mass of BH 2 (smallest black hole) */
	double m2;
	
	/*! \brief Amplitude of orbital the angular momentum */
	double AmpL;
	
	/*! \brief Amplitude of the spin of the object 1 (biggest object) */
	double AmpS1;
	
	/*! \brief Amplitude of the spin of the object 2 (smallest object) */
	double AmpS2;
	
	/*! \brief Cartesian coordinates of total angular momentum J in source (NR) reference frame */
	LCVector JN;
	
    
	// *************** Storing of numerical relativity data *************** 
	
	/*! \brief Name of file containing numerical relativity data to load */
	char * NRFileName;
	
	/*! \brief Number of harmonics */
	int NHarm;
	
	/*! \brief List of l index of harmonics. Size = [NHarm] */
	int * lHarm;
	
	/*! \brief List of m index of harmonics. Size = [NHarm] */
	int * mHarm;
	
	/*! \brief List of table ofdata containing harmonics from numerical relativity data. Size = [NHarm][2=(h+,hx)][NRdataN]  */
	double *** NRdata;
	
	/*! \brief Number of data for each harmonics*/
	int NRdataN;
	
	/*! \brief Initial time of data */
	double NRdatat0;
	
	/*! \brief End time of data */
	double NRdatatend;
	
	/*! \brief Time step of data */
	double NRdatadt;
	
	/*! \brief Observation time in seconds corresponding where initial condition are defined (computed) */
	double tInitObs;
	
	/*! \brief Hybrid time in t/M corresponding where initial condition are defined (read in NRdata file) */
	double toMInitHyb;
	
	/*! \brief Initial orbital frequency in \f$ omega \times M \f$ */
	double omMInit;
	
	/*! \brief Time shift between observation time (in sec) and hybrid time (in sec) : \f$ \Delta t^{obs \rightarrow hyb} = t^{obs} - t^{hyb} \f$ (computed)*/
	double tShiftHybObs;
	
	/*! \brief List of value of spherical harmonics : [NHarm] */
	dcomplex * SpherHarmVal;
	
	/*! \brief Maximal frequency times toatl mass (in seconds) : FreqMax = FreqMaxM/M  */
	double FreqMaxM;
	
	/**   Set the \f$l_{max} \f$ */
	int lmax;
    
    /** \brief Space for storing \f$ \partial h / \partial \tau \f$   */
    double * dhidtau;
    
    
	// ********** For taper and inspiral **********
	
	/*! \brief Pointer on GW SpinBBHHHarm1 */
	LCGWSpinBBHHHarm1 * GWSBHHH;
	
	/*! \brief If true, the taper will be applied; set via special parameter 0  */
	bool ApplyTaper;

	/*! \brief Reference time for applying the taper */
	double tTaper;
	
	/*! \brief Inspiral end time of inspiral used for applying the taper */
	double tendTaper;
	
	/*! \brief Maximal frequency which is used to parametrized the taper and which have to be used as maximal frequency when taper is applied */
	double FreqMaxTaper;
	
	/*! \brief Keeping in memory some time informations for constructing for example interne inspiral */
	double tInfot0, tInfotMaxDiff;
	
	
	// ********** For integration **********
	
	/*! \brief Firt time for which one the waveform is defined (in seconds) */
	double tStartWave;
		
	/*! \brief Last time for which one the waveform is defined (in seconds) */
	double tEndWave;
	
	/*! \brief Real firt time  : start of the observation */
	double t0Start;
	
	/*!  \brief Duration of observation */
	double maxDur;
	
	/*!  \brief Time step of the integration */
	double dt;
	
	
	
	
	// *************** Other internal variables ***************

	/*! \brief Amplitude */
	double Amp;
	
	/*! \brief Colattitude and longitude (in degrees) */
	double thetaS, phiS;
	
	/*! \brief Cosinus and sinus of colattitude */
	double cThS, sThS;
	
	/*! \brief Cosinus and sinus of 2 times the polarization */
	double c2psi, s2psi;

    double PSIN;

	
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCGWSpinBBHNR1();
	
	/*! \brief Standard constructor */
	LCGWSpinBBHNR1(LCTools * MT_n);
	
	
	/*! \brief Destructor */
	~LCGWSpinBBHNR1();
	
	
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
	 *	Coding for index iP of free parameters :
	 *	- 0  : \f$ \beta \f$		: #Beta						: Ecliptic latitude (in radians)
	 *	- 1  : \f$ \lambda \f$		: #Lambda					: Ecliptic longitude (in radians)
	 *	- 2  : \f$ \Psi \f$			: #Polarization				: Polarization (in radians)
	 *	- 3  : \f$ M_{tot} \f$		: #Mtot						: Total mass (in solar mass)
	 *	- 4  : \f$ tc \f$			: #tc						: Time at coalescence (in seconds)
	 *	- 5  : \f$ D_L \f$			: #DL						: Distance (in kpc)
	 *	- 6  : \f$ \theta_d  \f$	: #Thd						: Polar angle of the detector seen from the source (in radians)
	 *	- 7  : \f$ \phi_d \f$		: #Phd						: Azimuthal angle of the detector seen from the source (in radians)
	 *	Coding for index iP of fixed parameters :
	 *	- 8  : \f$ q \f$			: #mrat						: Mass ratio (>1)
	 *	- 9  : \f$ chi1 \f$			: #chi1						: Amplitude of spin1
	 *	- 10 : \f$ chi2 \f$			: #chi2						: Amplitude of spin2
	 *	- 11 : \f$ \theta_{S,S1} \f$: #thS1						: Initial polar angle of spin of object 1 in source reference frame
	 *	- 12 : \f$ \theta_{S,S2} \f$: #thS2						: Initial polar angle of spin of object 2 in source reference frame
	 *	- 13 : \f$ \phi_{S,S1} \f$	: #phS1						: Initial azimuthal angle of spin of object 1 in source reference frame
	 *	- 14 : \f$ \phi_{S,S2} \f$	: #phS2						: Initial azimuthal angle of spin of object 2 in source reference frame
	 *	- 15 : \f$ \phi_{0} \f$		: #Phi0						: Initial phase
	 *	Coding for index iP of parameters which can be use instead of defining Thd, Phd and Polarization : 
	 *	- 16 : \f$ \theta_{B,S1} \f$: #thBS1					: Initial polar angle of spin of object 1 in solar system barycentric (SSB) reference frame
	 *	- 17 : \f$ \theta_{B,S2} \f$: #thBS2					: Initial polar angle of spin of object 2 in solar system barycentric (SSB) reference frame
	 *	- 18 : \f$ \phi_{B,S1} \f$	: #phBS1					: Initial azimuthal angle of spin of object 1 in solar system barycentric (SSB) reference frame
	 *	- 19 : \f$ \phi_{B,S2} \f$	: #phBS2					: Initial azimuthal angle of spin of object 2 in in solar system barycentric (SSB) reference frame
	 *	- 20 : \f$ \theta_{B,L} \f$	: #thBL						: Initial polar angle of orbital angular momentum in solar system barycentric (SSB) reference frame
	 *	- 21 : \f$ \phi_{B,L} \f$	: #phBL						: Initial azimuthal angleof orbital angular momentum in solar system barycentric (SSB) reference frame
	 
	 */
	double getParam(int iP);
	
	/*! \brief Set value of parameter (same for ) 
	 *	@param[in]	iP		: Index of parameter (see #getParam for coding)
	 *	@param[in]	param_n : Value of parameter
	 */
	void setParam(int iP, double Param_n);
	
	/*! \brief Set \Pmin and \Pmax at range of parameter @iP */
	void getRange(int iP, double &Pmin, double &Pmax);
	
	/*! \brief Return small variation delta used in computation of derivative for parameter @iP */
	double getDelta(int iP);
	
	/*! \brief Get value of specific parameter @iPS 
	 *	Coding for index iP :
	 *	- 0  : ApplyTaper	: If > 0, the taper will be aplied
	 */
	double getSpecialParam(int iPS);
	
	/*! \brief Set value of specific parameter @iPS at @SpecParam_n */
	void setSpecialParam(int iPS, double SpecParam_n);
	
	/*! \brief Set time informations */
	void setTimeInfo(double dt);
	
	/*! \brief Choose randomly parameter iP */
	void RandParam(int iP);
	
	/*! \brief Set time informations :
	 *	@param[in]	t0_n		: Initial time
	 *	@param[in]	dt_n		: Time step
	 *	@param[in]	tLast_n		: Maximum time shift between detector and GW (reference for GW is the barycenter of the solar system)
	 *	@param[in]	TObs_n		: Duration of observation
	 *	@param[in]	tMaxDiff	: Maximal time difference between estimation of GW : typically 2 times armlength is OK 
	 */
	void setTimeInfo(double t0_n, double dt_n, double TObs_n, double tAskMin_n, double tDOrbMax_n,  double tMaxDiff);
	
	
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
	
	/** \brief Load numerical relativity data */
	void LoadNRdata();

	/** \brief Function to read data in information line of NRdata file */
	double ReadInfoNRdata(std::vector<char*> Words);
	
	/** \brief Return \f$ h_+ \f$ of harmonic l,m at t/M  */
	double hpharm(int l, int m, double toM);
	
	/** \brief Return \f$ h_{times} \f$ of harmonic l,m at t/M  */
	double hcharm(int l, int m, double toM);
    
    /** \brief Return \f$ \partial h_+ / \partial t \f$ of harmonic l,m at t/M  */
	double dhpharmdtau(int l, int m, double toM);
	
	/** \brief Return \f$ \partial h_{\times} / \partial t \f$ of harmonic l,m at t/M  */
	double dhcharmdtau(int l, int m, double toM);
    
    /** \brief Return \f$ \partial h_{+ or \times} / \partial t \f$ of harmonic l,m at t/M  */
	double dhharmdtau(int ihphc, int l, int m, double toM);
	
	/** \brief Find the index of the harmonic corresponding to the l and m in argument. If don't find return -1 */
	int IndHarm(int l, int m);
	
	/** \brief Return the value corresponding to time toM for h+ (ihphc=0) or hx (ihphc=1) of harmonic iH using linear interpolation */
	double gHarmDataLIN(int iH, int ihphc, double toM);
	
    /** \brief Return Lagrange interpolation on hybrids :
        @param[in]  iH    : Index of harmonic
        @param[in]  ihphc : 0 for h+ and 1 for hx
        @param[in]  toM   : \f$ \tau = t / M \f$ 
        @param[in]  order : Order of the interpolation 
        @return Result of interpolation
     */
    double InterLagrangeNR(int iH, int ihphc, double toM, int order);
    
	/** \brief Return the spherical harmonics \f$Y_{lm}(\theta, \phi)$    */
	dcomplex ComputeSpherHarm( int l, int m, double theta, double phi);
	
	/** \brief Return the spherical harmonics form stored data */
	dcomplex SpherHarm(int l, int m);
	
	/**  \brief Return  n! : compute the factorial function  */
	int factorial (int num);
	
	/** \brief Compute extra parameters */
	void ComputeExtraParam();
	
	/** \brief Compute the total angular momentum */
	void ComputeTotalAngMom();
	
	/**  \brief Compute angle describing the direction of spin of index iSpin in Solar System Barycentric frame  */
	void ConvertSpinDirToSSB(int iSpin);
	
	
	
	/*! \brief Compute approximation of orbital frequency using PN model */
	/*double ComputeApprOrbFreq(double t, 
							  double Tc,
							  double x1,
							  double x2,
							  double LS1, 
							  double LS2,
							  double S1S2,
							  double m1M2,
							  double m2M2,
							  double eta);
	 */
};

#endif //__LCGWSPINBBHNR1_H

/**\}*/

// end of LISACODE-LCGWSPINBBHNR1.h