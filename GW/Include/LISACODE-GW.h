/*
 *  LISACODE-GW.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 23/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


/** \defgroup GW GW
 * This group contains all the things related to gravitational wave modeling.
 * \{
 */


#ifndef __LCGW_H
#define __LCGW_H

#include "LISACODE-Tools.h"
#include "LISACODE-Vector.h"
#include "LISACODE-Matrix.h"
#include "ezxml.h"


/** \brief Base class for managing gravitational waves.
 * \author A. Petiteau
 * \version 2.0
 * \date 23/04/2011
 *
 * This class is the base class for managing graviational waves.
 * 
 *
 */
class LCGW
{
protected:
	
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	/** Last value of \f$ h_{B+}(t) \f$ in barycentric reference frame computed at #tLast */
	double hBpLast;
	
	/** Last value of \f$ h_{B \times}(t) \f$ in barycentric reference frame computed at #tLast */
	double hBcLast;
	
	/** Last time of computation of \f$ h_{B+}(t) \f$  and  \f$ h_{B \times}(t) \f$ */
	double tLast;
	
	/** Minimal (smallest or last) time of computation of \f$ h_{B+}(t) \f$  and  \f$ h_{B \times}(t) \f$ */
	double tAskMin;
	
	/** Maximal (largest or first) time of computation of \f$ h_{B+}(t) \f$  and  \f$ h_{B \times}(t) \f$ */
	double tDOrbMax;
    
    /** Real t0 for correcting time */
    double t0Real;
	
	/** Maximal duration before recomputing \f$ h_{B+}(t) \f$  and  \f$ h_{B \times}(t) \f$ */
	double tStore;
	
	/** ********** Parameters ********** */
	
	/** \brief Number of parameters */
	int NParams;
	
	/** \brief Ecliptic latitude (declination) of sky position */
	double Beta;
	
	/** \brief Ecliptic longitude  of sky position */
	double Lambda;
    
    /** \brief Ecliptic latitude (declination) f the total Angular momentum */
	double thetaJ;
	
	/** \brief Ecliptic longitude  of the total Angular momentum */
	double phiJ;
	
	/** \brief Polarization angle : if 0 no computation at all */
	double Polarization;
	
	/** \brief If true, make computation of polarization */
	bool ComputePolarization;
	
	/** \brief Cosinus of polarization */
	double c2Pol;
	
	/** \brief Sinus of polarization */
	double s2Pol;
	
	/*! \brief Minimal frequency of GW (-1.0 --> undefined) */
	double FreqMin;
	
	/*! \brief Maximal frequency of GW (-1.0 --> undefined) */
	double FreqMax;
	
	/*! \brief Name of the source : fix size at 128 */
	char * Name;
	
	/*! Size of memory used for GW modelling (in o) */
	double MemUse;
	
    /* \brief Flag to detect if a parameter has been found or not */
    bool ParameterReaded;
    
    
#ifdef _DEBUG_GW_     
	/*! \brief Pointer on a file used for checking */
    std::ofstream * DEBUGfCheck;
    bool DEBUGDispAll;
#endif	
	
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCGW();
	
	/*! \brief Standard constructor */
	LCGW(LCTools * MT_n);
	
	
	/*! \brief Destructor */
	virtual ~LCGW();
	
	/*! \brief Initialization at NULL of all base of GW */
	void initNULLBase(bool CleanMem);
	
	
	/***************  Required methods (to be present in derived classes)  ***************/
	
	/*! \brief Initialization at NULL */
	virtual void initNULL(bool CleanMem);
	
	/**********  Configuration methods  **********/
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="GW"> ... </XSIL> */
	virtual void config(ezxml_t xmlbloc);
	
	/*! \brief Configuration of individual parameter */
	virtual void config(int iParam, double ParamVal);
	
	
	/**********  Access methods  **********/
	
	/*! \brief Return value of parameter @iP */
	virtual double getParam(int iP);
	
	/*! \brief Set value of parameter @iP at @Param_n */
	virtual void setParam(int iP, double Param_n);
	
	/*! \brief Set \Pmin and \Pmax at range of parameter @iP */
	virtual void getRange(int iP, double &Pmin, double &Pmax);
	
	/*! \brief Return small variation delta used in computation of derivative for parameter @iP */
	virtual double getDelta(int iP);
	
	/*! \brief Set value of specific parameter @iPS at @SpecParam_n */
	virtual void setSpecialParam(int iPS, double SpecParam_n);
	
	/*! \brief Get value of specific parameter @iPS */
	virtual double getSpecialParam(int iPS);
	
	/*! \brief Choose randomly parameter iP */
	virtual void RandParam(int iP);
	
	/*! \brief Set time informations :
	 *	@param[in]	t0_n		: Initial time
	 *	@param[in]	dt_n		: Time step
	 *	@param[in]	TObs_n		: Duration of observation
	 *	@param[in]	tAskMin_n	: Smallest time at which one we will ask to compute \f$ h_{+}\f$ and/or \f$ h_{\times}\f$
	 *	@param[in]	tDOrbMax_n	: Largest time at which one we will ask to compute \f$ h_{+}\f$ and/or \f$ h_{\times}\f$ 
	 *	@param[in]	tMaxDiff	: Maximal time difference between estimation of GW : typically 2 times armlength is OK 
	 */
	virtual void setTimeInfo(double t0_n, double dt_n, double TObs_n, double tAskMin_n, double tDOrbMax_n,  double tMaxDiff);
	
	
	/**********  Linking and initalization methods  **********/
	
	/*! \brief Initialization */
	virtual int init();
	
	
	/**********  Running methods  **********/
	
	/*! \brief Compute gravitational strain in barycentric reference frame.  
	 * If polarization is null, computation in barycentric reference frame.
	 * Else computation in source reference frame
	 */
	virtual void Computehpc(double t); 
		
	
	/**********  Others methods  **********/
	
	/*! \brief Display information about the GW */ 
	virtual void DispInfo(char * BTab);
	
	/*! \brief Display all parameters */
	virtual void DispAllParam(std::ostream * out);
	
	/*! \brief Display name of all parameters */
	virtual void DispAllParamName(std::ostream * out);
	
	
	/***************  Local methods  ***************/
	
	
	/**********  Configuration methods  **********/
	
	/*! \brief Basic configuration : read sky postion, etc */
	void configBase(ezxml_t xmlbloc);
	
	
	/**********  Linking and initalization methods  **********/
	
	/*! \brief Initialization of base of GW */
	void initBase();
	
	
	/**********  Access methods  **********/
	
	/*! \brief Return the number of parameters */
	int getNParam() {return(NParams);};
	
	/*! \brief Return the ecliptic latitude (declination) of sky position */
	double getSkyBet() {return(Beta);};
	
	/*! \brief Return the ecliptic longitude of sky position */
	double getSkyLam() {return(Lambda);};
	
	
	/*! \brief Return minimal frequency */
	double getFreqMin() {return(FreqMin);} ;
	
	/*! \brief Return maximal frequency */
	double getFreqMax() {return(FreqMax);} ;
    
	
	/*! \brief Return the size of memory used for GW modelling (in o)*/
	double getMemUse() {return(MemUse);};
	
	
	/*! \brief Choose randomly all parameters */
	void RandAllParams();
	
	
	/*! \brief Return pointer on source Name */
	char * getName() {return(Name);};
	
	/*! \brief Set the source Name */
	void setName(const char * Name_n);
	
	
	/**********  Running methods  **********/
	
	/*! \brief Compute \f$ h_{+}(t) \f$ and \f$ h_{\times}(t) \f$ in barycentric reference frame :
	 * If polarization is null, direct computation in barycentric reference frame.
	 * Else computation in source reference frame then conversion in barycentric reference frame
	 */
	void ComputehBpc(double t);
	
	/*! \brief Return gravitational strain \f$ h_{B+}(t) \f$ in barycentric reference frame  */
	double hBp(double t); 
	
	/*! \brief Return gravitational strain \f$ h_{B\times}(t) \f$ in barycentric reference frame  */
	double hBc(double t); 
	
	/*! \brief Return gravitational strain \f$ h_{+}(t) \f$ and \f$ h_{\times}(t) \f$ in barycentric reference frame */
	void hBpc(double t, double & hBpV, double & hBcV); 
		
	
	/********** Tools methods **********/
	
	/*! \brief Compute half of the Hanning window */
	double halfhann(double t, double t0, double t1);
	
	
	/**  \brief Compute angle describing the direction of spin of index iSpin in Solar System Barycentric frame  
	 * @param[in]	Thd_i			Polar angle of the detector seen from the source (in radians)
	 * @param[in]	Phd_i			Azimuthal angle of the detector seen from the source (in radians)
	 * @param[in]	Pol_i			Polarization 
	 * @param[in]	JN_i			Vector in Numerical relativity frame of total angular
	 * @param[in]	LnN_i			Unit vector in Numerical relativity frame of orbital angular momentum 
	 * @param[in]	S1N_i			Unit vector in Numerical relativity frame of spin 1
	 * @param[in]	S2N_i			Unit vector in Numerical relativity frame of spin 2 
	 * @param[in]	AmpL_i			Amplitude of orbital angular momentum 
	 * @param[in]	AmpS1_i			Amplitude of spin 1
	 * @param[in]	AmpS2_i			Amplitude of spin 2
	 * @param[out]	LnB_o			Unit vector in SSB frame of orbital angular momentum
	 * @param[out]	S1B_o			Unit vector in SSB frame of spin 1
	 * @param[out]	S2B_o			Unit vector in SSB frame of spin 2
	 * @param[out]	PolRecompute	Recomputed polarization for checking
	 */
	void ConvertLS1S2dirNR2SSB(double Thd_i, 
							   double Phd_i,
							   double Pol_i,
							   LCVector JN_i,
							   LCVector LnN_i,
							   LCVector S1N_i,
							   LCVector S2N_i,
							   double AmpL_i,
							   double AmpS1_i,
							   double AmpS2_i,
							   LCVector & LnB_o,
							   LCVector & S1B_o,
							   LCVector & S2B_o,
							   double & PolRecompute);
    
    
    void ConvertLS1S2dirNR2SSBSofVer(double & Phd_o, 
                                           double ThS_i,
                                           double PhS_i,
                                           double thetaJ_i,
                                           double phiJ_i,
                                           LCVector LnN_i,
                                           LCVector S1N_i,
                                           LCVector S2N_i,
                                           double AmpL_i,
                                           double AmpS1_i,
                                           double AmpS2_i,
                                           LCVector & LnB_o,
                                           LCVector & S1B_o,
                                           LCVector & S2B_o);
    
	
	/*! \brief Compute angle describing the angles \f$ \theta_d \f$, \f$ \phi_d \f$ and polarization from other parameters
	 * @param[in]	LnB_i		Unit vector in SSB frame of orbital angular momentum
	 * @param[in]	S1B_i		Unit vector in SSB frame of spin 1
	 * @param[in]	S2B_i		Unit vector in SSB frame of spin 2
	 * @param[in]	LnN_i		Unit vector in Numerical relativity frame of orbital angular momentum 
	 * @param[in]	S1N_i		Unit vector in Numerical relativity frame of spin 1
	 * @param[in]	S2N_i		Unit vector in Numerical relativity frame of spin 2
	 * @param[in]	AmpL_i		Amplitude of orbital angular momentum 
	 * @param[in]	AmpS1_i		Amplitude of spin 1
	 * @param[in]	AmpS2_i		Amplitude of spin 2
	 * @param[out]	Thd_o		Polar angle of the detector seen from the source (in radians)
	 * @param[out]	Phd_o		Azimuthal angle of the detector seen from the source (in radians)
	 * @param[out]	Pol_o		Polarization
	 *
	 */
	void ComputeThdPhdPol(LCVector LnB_i,
						  LCVector S1B_i,
						  LCVector S2B_i,
						  LCVector LnN_i,
						  LCVector S1N_i,
						  LCVector S2N_i,
						  double AmpL_i,
						  double AmpS1_i,
						  double AmpS2_i,
						  double & Thd_o, 
						  double & Phd_o, 
						  double & Pol_o);
    
    
    
    void ComputeThdPhdPolSofVer(LCVector LnB_i,
                                      LCVector S1B_i,
                                      LCVector S2B_i,
                                      LCVector LnN_i,
                                      LCVector S1N_i,
                                      LCVector S2N_i,
                                      double AmpL_i,
                                      double AmpS1_i,
                                      double AmpS2_i,
                                      double & Thd_o, 
                                      double & Phd_o, 
                                      double & Pol_o,
                                      double  thetaJ_i,
                                      double  phiJ_i);
	
	/** \brief Function to read NR parameters in  NRdata file */
	int ReadNRParam(char * FileName_i,
						  double & mrat_o,
						  double & eta_o,
						  double & chi1_o,
						  double & chi2_o,
						  double & Phi0_o,
						  double & FreqMaxM_o,
						  double & omMInit_o,
						  double & toMInitHyb_o,
						  LCVector & S1N_o,
						  LCVector & S2N_o,
						  LCVector & LnN_o);
	
	/** \brief Function to read data in information line of NRdata file */
	double ReadInfoNRdata(std::vector<char*> Words);
	
	
	/**********  Others methods  **********/
	
	/*! \brief Add Delta x FactDelta to value of the parameter iP o:
	 *	@param[in]	iP			Index of parameter to be modified
	 *	@param[in]	FactDelta	Factor apply on delta
	 *	@return				Value of delta
	 */
	double AddDeltaPar(int iP, double FactDelta);
	
	
	/*! \brief Display the basis */
	void DispInfoBase(char * BTab);
	
	
		
	
};

#endif //__LCGW_H

/**\}*/

// end of LISACODE-GW.h