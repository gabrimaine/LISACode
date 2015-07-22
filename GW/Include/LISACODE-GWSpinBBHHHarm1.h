/*
 *  LISACODE-GWSpinBBHHHarm1.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 31/05/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


/** \ingroup GW
 * \{
 */


#ifndef __LCGWSpinBBHHHarm1_H
#define __LCGWSpinBBHHHarm1_H

#include "LISACODE-Tools.h"
#include "ezxml.h"
#include "LISACODE-GW.h"

/** \brief Class for managing gravitational waves emitted by black hole binaries with spin considering inspiral and higher harmonics.
 * \author A. Petiteau
 * \version 2.0
 * \date 23/04/2011
 *
 * This class manages gravitational waves emitted by galactic binaries.
 * The waveform is based on [<a href="Bibliography.html#Arun-2009">Arun &nbsp; 2009</a>]
 *
 */
class LCGWSpinBBHHHarm1 : public LCGW
{
protected:
	
	
	// *************** Parameters of system *************** 
	
	/*! \brief Mass of object 1 (in solar masses) */
	double mass1;                
	
	/*! \brief Mass of object 2 (in solar masses) */
	double mass2;
	
	/*! \brief Time of coalescence (in seconds) */
	double tc;
	
	/*! \brief Distance between source and detector  (in kiloparsec) */
	double DL;
	
	/*! \brief Spin of object 1  */
	double chi1;
	
	/*! \brief Spin of object 2 */
	double chi2;
	
	/*! \brief Polar angle of spin S1 */
	double PolarAngleOfSpin1;
	
	/*! \brief Polar angle of spin S2 */
	double PolarAngleOfSpin2; 
	
	/*! \brief Azimuthal angle of spin S1 */   
	double AzimuthalAngleOfSpin1;
	
	/*! \brief Azimuthal angle of spin S1 */
	double AzimuthalAngleOfSpin2;
	
	/*! \brief Phase at coalescence */
	double phic, Phi0;
	
	/*! \brief Initial polar angle of orbital momentum \brief  */
	double InitialPolarAngleL;
	
	/*! \brief Initial azimuthal angle of orbital momentum \brief  */
	double InitialAzimuthalAngleL;
	
	
	// *************** Mass variables ***************
	
	/*! \brief Mass ratio and it's inverse */
	double eta;
	
	/*! \brief Chirp mass (in solar mass) */
	double Mchirp;
	
	/*! \brief Reduced mass (in solar mass) */
	double mum;
	
	/*! \brief Total mass (in solar mass) */
	double Mtot;
	
	/*! \brief  \f$ m_1^2/M^2 \f$  */
	double m1M2;
	
	/*! \brief  \f$ m_2^2/M^2 \f$ */
	double m2M2;
	
	/*! \brief Mass of object 1 (in second) */
	double m1; 
	
	/*! \brief Mass of object 2 (in second) */
	double m2;
	
	/*! \brief Total mass (in second) */
	double M;
	
	/*! \brief Reduced mass (in second) */
	double mu;
	
	/*! \brief Difference of mass (in second) */
	double dm;
	
	/*! \brief Distance between source and detector (in seconds) */
	double dist;
	
	
	
	// *************** Precession variables ***************
	
	/* \brief Spin amplitude of BH1 and BH2 (equivalent to chi1 and chi2) */
	double x1, x2;
	
	/* \brief  */
	double iota, alpha;
	
	/* \brief Components of unit vector of orbital angular momentum */
	double Lnx, Lny, Lnz;
	
	/* \brief Components of unit vector of spin of BH1  */
	double S1x, S1y, S1z;
	
	/* \brief Components of unit vector of spin of BH2  */
	double S2x, S2y, S2z;
	
	/* \brief Scalar product */
	double LS1, LS2, S1S2;
	
	/* \brief Polarization (NR parameter) */
	double PolNR;
	
	/* \brief Polar angle of the detector seen from the source in NR frame (in radians) (NR parameter) */
	double ThdNR;
	
	/* \brief Azimuthal angle of the detector seen from the source in NR frame (in radians) (NR parameter) */
	double PhdNR;
	
	
	// ********** For taper **********
	
	/*! \brief Q-factor for the new taper : if > 0.0 apply the MECO, else (f[i]<f[i-1])&&(x)  */
	double taperQ;

	/*! \brief Reference time for applying the taper */
	double tTaper;
	
	/*! \brief True if we should apply taper : if merger during the observation */
	bool ApplyTaper;
	
	
	// ********** For integration **********
	
	/*!  \brief Initial time for starting the integration */
	double t0Intg;
	
	/*!  \brief Initial time for starting the integration in t/M unit if > -1.0e30, it will be use */
	double t0IntgtoM;
	
	/*!  \brief Initial omegaM, orbital ferquency times total mass in seconds.  If > -1.0e30, it will be use */
	double omM0Intg;
	
	/*! \brief Real firt time  : start of the observation */
	double t0Start;
	
	/*!  \brief Duration of observation */
	double maxDur;
	
	/*!  \brief True if backward integration */
	bool back;
	
	/*!  \brief Time step of the integration */
	double dt;
	
	/*! \brief Current time (last point computed) of integration */
	double tIc;
	
	/*! \brief Number of last steps kept in memory */
	int NMemItg;
	
	/*! \brief Last values (last point computed) of integration storing the #NMemItg last step. Size : NMemItg x 13 */
	double ** coordI;
	
	/*! \brief End time of the inspiral waveform (in seconds)*/
	double tend;
	
	
	
	// *************** Internal variables about precession etc ***************
	
	/*! \brief True if extra parameters (as thLB, phLB, etc) have to be compute */
	bool NeedExtraParamCompute;
	
	/*! \brief Type of extra parameter that we have to compute : 1->(#thLB, #phLB, #thS1, #phS1, #thS2 and #phS2) or 2->(#Polarization, #Thd and #Phd).
     So 1 for NR parameter and 2 for standard parameters 
     */
	int TypeExtraParamCompute;
	
	/*! \brief Amplitude of orbital the angular momentum */
	double AmpL;
	
	/*! \brief Amplitude of the spin of the object 1 (biggest object) */
	double AmpS1;
	
	/*! \brief Amplitude of the spin of the object 2 (smallest object) */
	double AmpS2;
	
	/*! \brief Cartesian coordinates of Ln, S1 and S2 in source reference frame */
	LCVector LnB, S1B, S2B;
	
	/*! \brief Cartesian coordinates of diretion of Ln, S1 and S2 in source (NR) reference frame */
	LCVector LnN, S1N, S2N;
	
	/*! \brief Cartesian coordinates of total angular momentum J in source (NR) reference frame */
	LCVector JN;
	
	double NRmrat, NReta, NRchi1, NRchi2, NRPhi0, NRFreqMaxM, NRomMInit, NRtoMInitHyb;
	
	// *************** Other internal variables ***************
	/*! \brief Fixed amplitude */
	double Amp;
	
	/*! \brief Fixed amplitude with taper applied */
	double wk, v;
	
	/*! \brief  */
	double theta;
	
	/*! \brief Colattitude and longitude (in degrees) */
	double thetaS, phiS;
	
	/*! \brief Cosinus and sinus of colattitude */
	double cThS, sThS;
	
	/*! \brief  */
	double psi;
	
	/*! \brief Last value of energy */
	double En_prev;
	
	/*! \brief Last value of angular velocity */
	double om_prev;
	
	/*! \brief Initial value of angular velocity */
	double om0;
	
	/*! \brief Current (last out of the interpolation) of angular velocity */
	double om;
	
	/*! \brief Current value (last out of the interpolation) of phase */
	double phi;
	
	/*! \brief Current value (last out of the interpolation) of corrected phase \f$ \Psi = \phi - 2 M \omega \log {\omega \over \omega_0} \f$ */
	double Psi;
	
	/*! \brief True if we are considering the non-spinning case */
	bool nonspin;
	
	/*! \brief \f$ h_+ \f$ and \f$ h_{\times} \f$ in source frame */
	double hp, hc;
	
	/*!  \brief Unit vectors for switching between Source frame and Barycentric frame  */
	double ex[3];
	double ey[3];
	double ez[3];
	
	/*! \brief Direction of the source */
	double n[3];
	
	double LnBx, LnBy, LnBz; 
	double PhLnB, ThLnB;
	double PolB, c2PB, s2PB;
	
	
	double h0PNp, h0PNc, h05PNp, h05PNc, h1PNp, h1PNc;
	double beta, sigma, tau, tau38, tau14, tau58, x;
	
	
	/*! \brief Cosninus and sinus of some angles */ 
	double ci, si, cth, sth, ciby2, siby2;
	double s2th, c2th, s3th, c3th, c4th, s4th, cth2, sth2, cth3;
	double cthsth2, sthcth3;
	
	
	double ChpC1th, ChpS1th, ChpC1thE2, ChpS1thE2, ChpC1thE3, ChpC1thE4, ChpC2th, ChpC3th, ChpC4th;
	double ChpS2th, ChpS3th, ChpS4th, Chp5S1thPS3th, Chp3PC2th, Chp1P3C2th, ChpS1thM3S3th, Chp5p4C2thP7C4th, ChpC1thXS1th; 
	
	std::ofstream * fOutCheckPrec;
	double tCheckPrecOld;
	
	
	bool CheckData;
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCGWSpinBBHHHarm1();
	
	/*! \brief Standard constructor */
	LCGWSpinBBHHHarm1(LCTools * MT_n);
	
	
	/*! \brief Destructor */
	~LCGWSpinBBHHHarm1();
	
	
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
	 *	- 0  : \f$ \beta \f$		: #Beta						: Ecliptic latitude (in radians)
	 *	- 1  : \f$ \lambda \f$		: #Lambda					: Ecliptic longitude (in radians)
	 *	- 2  : \f$ m_1 \f$			: #mass1					: Mass1 (in solar mass) redshifted (what we see)
	 *	- 3  : \f$ m_2 \f$			: #mass2					: Mass2 (in solar mass) redshifted (what we see)
	 *	- 4  : \f$ tc \f$			: #tc						: Time at coalescence (in seconds)
	 *	- 5  : \f$ D_L \f$			: #DL						: Distance (in kpc)
	 *	- 6  : \f$ \chi_1 \f$		: #chi1						: Amplitude of spin of BH1
	 *	- 7  : \f$ \chi_2 \f$		: #chi2						: Amplitude of spin of BH2
	 *	- 8  : \f$ \theta_{S_1} \f$	: #PolarAngleOfSpin1		: Initial polar angle of spin of BH1 (in radians)
	 *	- 9  : \f$ \theta_{S_2} \f$	: #PolarAngleOfSpin2		: Initial polar angle of spin of BH2 (in radians)
	 *	- 10 : \f$ \varphi_{S_1} \f$: #AzimuthalAngleOfSpin1	: Azimuthal angle of spin of BH1 (in radians)
	 *	- 11 : \f$ \varphi_{S_2} \f$: #AzimuthalAngleOfSpin2	: Azimuthal angle of spin of BH2 (in radians)
	 *	- 12 : \f$ \Phi_0  \f$		: #Phi0						: Initial phase (in radians)
	 *	- 13 : \f$ \theta_L \f$		: #InitialPolarAngleL		: Initial polar angle of orbital momentum (in radians)
	 *	- 14 : \f$ \varphi_L  \f$	: #InitialAzimuthalAngleL	: Initial azimuthal angle of orbital momentum (in radians)
	 *	- 15 : \f$ M_{c}  \f$		: #Mchirp					: Chirp mass (in solar mass) redshifted (what we see)
	 *	- 16 : \f$ \mu  \f$			: #mum						: Reduced mass (in solar mass) redshifted (what we see)
	 *	- 17 : \f$ \eta  \f$		: #eta						: Symetric mass ratio
	 *	- 18 : \f$ \phi_c  \f$		: #phic						: Phase at coalescence (in radians)
	 *	- 19 : \f$ M_{tot}  \f$		: #Mtot					    : Total mass (in solar mass) redshifted (what we see)
	 *	- 20 : \f$ \Psi \f$		    : #PolNR					: Polarization (NR parameter)
	 *	- 21 : \f$ \theta_d \f$		: #ThdNR					: Polar angle of the detector seen from the source in NR frame (in radians) (NR parameter)
	 *	- 22 : \f$ \phi_d \f$		: #PhdNR					: Azimuthal angle of the detector seen from the source in NR frame (in radians) (NR parameter)
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
	
	/*! \brief Get value of specific parameter @iPS 
	 *	Coding for index iP :
	 *	- 0  : #t0Intg		: Initial integration time (in seconds)
	 *	- 1  : #t0Intg		: Initial separation (in M distance. ex : 9 for 9M)
	 *	- 2  : #t0IntgtoM	: Initial integration time (in unit of t/M )
	 *	- 3  : #omM0Intg	: Initial omegaM, orbital ferquency times total mass in seconds.
	 */
	double getSpecialParam(int iPS);
	
	/*! \brief Set value of specific parameter @iPS at @SpecParam_n */
	void setSpecialParam(int iPS, double SpecParam_n);
	
	/*! \brief Set time informations */
	void setTimeInfo(double dt);
	
	/*! \brief Choose randomly parameter iP */
	void RandParam(int iP);
	
	/*! \brief Set time informations : see LCGW::setTimeInfo */
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
	void getTaperInfo(double & tTaper_s, double & tendTaper_s, double & FreqMaxTaper_s);
	
	/**********  Running methods  **********/
	
	
	/**********  Others methods  **********/
	
	
	/*! \brief Compute inspiral */
	void ComputeInspiral();
	
	/*! \brief Compute approximation of orbital frequency */
	double ComputeApprOrbFreq(double t, double Tc);
	
	/*! \brief Compute approximation of orbital Phase */
	double ComputeApprOrbPhase(double t, double Tc);
	
	/*! \brief Integrator */
	double Integrator(double* ystart, double* dydx, double* yend, double tstart, double dt, int n);
	
	/*! \brief Compute derivative */
	void Derivs(double x, double* y, double* dydx, int n);
	
	
	/*! \brief Compute amplitude of \f$ h_+ \f$ and \f$ h_{\times} \f$ at order 0 PN */
	void Computeh0PN(double& h0PNp, double& h0PNc);

	/*! \brief Compute amplitude of \f$ h_+ \f$ and \f$ h_{\times} \f$ at order 0 PN and for harmonic 22*/
	void Computeh0PNl2(double& h0PNp, double& h0PNc);
	
	/*! \brief Compute amplitude of \f$ h_+ \f$ and \f$ h_{\times} \f$ at order 0.5 PN */
	void Computeh05PN(double& h05PNp, double& h05PNc);
	
	/*! \brief Compute amplitude of \f$ h_+ \f$ and \f$ h_{\times} \f$ at order 1 PN */
	void Computeh1PN(double& h1PNp, double& h1PNc);
	
	/*! \brief Compute amplitude of \f$ h_+ \f$ and \f$ h_{\times} \f$ at all PN orders (0, 0.5 and 1 PN) */
	void ComputehAllPN(double& h0PNp, double& h0PNc, double& h05PNp, double& h05PNc, double& h1PNp, double& h1PNc);
	
	/*! \brief Compute amplitude of \f$ h_+ \f$ and \f$ h_{\times} \f$ at all PN orders (0, 0.5 and 1 PN) : second version */
	void ComputehAllPN2(double& h0PNp, double& h0PNc, double& h05PNp, double& h05PNc, double& h1PNp, double& h1PNc);
	
	
	/*! \brief Compute instantaneous energy of the system */
	double ComputeEnergy();
	
	
};

#endif //__LCGWSpinBBHHHarm1_H

/**\}*/

// end of LISACODE-LCGWSpinBBHHHarm1.h