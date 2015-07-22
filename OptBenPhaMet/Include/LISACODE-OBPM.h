/*
 *  LISACODE-OBPM.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 12/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \defgroup OBPM Optical Bench Phasemeter
 * This group contains all the things related to the measurement system : optical bench path + photodiode + phasemeter 
 * \{
 */


#ifndef __LCOBPM_H
#define __LCOBPM_H

#include "LISACODE-Tools.h"
#include "ezxml.h"
#include "LISACODE-ArmResp.h"
#include "LISACODE-Noise.h"
#include "LISACODE-Orbits.h"
#include "LISACODE-USO.h"
#include "LISACODE-Filter.h"


/** \brief Base class for managing the measurement system : optical bench path + photodiode + phasemeter .
 * \author A. Petiteau
 * \version 2.0
 * \date 10/04/2011
 *
 * This class is the base class for managing the measurement system : optical bench path + photodiode + phasemeter 
 * 
 *
 */
class LCOBPM
{
protected:
	
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	/** \brief Name of the optical bench path + phasemeter */
	char Name[8];
	
	/** \brief Pointer on module computing arm response to GWs (extern allocation : in Detector) */
	LCArmResp * ArmGW;
	
	/** \brief List of pointers on used noises */
	LCNoise ** pn;
	
	/** \brief Number of pointers on used noises */
	int Npn;
	
	/** \brief Pointer on local clock */
	LCUSO * clockl;
	
	/** \brief Pointer on orbits */
	LCOrbits * orb;
	
	/** \brief Pointer on measurement of photodiode : raw results */
	LCSerie2 * MesPho;
	
	/** \brief Pointer on measurement of the  phasemeter : filtered results at measurement time step */
	LCSerie2 * MesPha;
	
	/** \brief Phasemeter filter  */
	LCFilter * PMFilter;
	
	/** \brief Pointer on filtered data (allocation in the filter) */
	LCSerie2 * MesFilter;
	
	/** \brief Pointer on the value of the last measurement */
	double * LastMes;
	
	/** \brief True if local allocation for the value of the last measurement */
	bool LastMesLocalAlloc;
	
	/** \brief Current value of the graviational wave signal */
	double GWSignal;
	
	/** \brief Pointer on the factor to be applied on shot noise times delays */
	double * FactShotNoise;
	
	
	/******** Position of the optical bench and phasemeter  ********/
	
	/*! \brief  Index of spacecraft corresponding to the photodetector-phasemeter */
	int iSC; 
	
	/*!\brief  Direction flag of the optical bench :  0 if the optical bench is in the direct direction (1->2->3->1), else 1 (1->3->2->1) */
	int IndirectDir;
	
	/*! \brief  Index of the linked spacecraft corresponding to the photodetector-phasemeter */
	int iSCLink;
	
	/*!\brief  Direction flag of the linked optical bench :  0 if the optical bench is in the direct direction (1->2->3->1), else 1 (1->3->2->1) */
	int IndirectDirLink;
	
	
	/************ Variables use in initialization but not after ********** */
	
	/******** Parameters of filter  ********/
	
	//! ** Attenuation (in dB)
	double PMFAttenuation; 
	
	//! ** Oscillations in bandwidth (in dB)
	double PMFOscillation;
	
	//! ** Low transition frequency in factor of frequency of measurment (factor of measurement frequency)
	double PMFLowFreqFact;
	
	//! ** High transition frequency in factor of frequency of measurment (factor of measurement frequency)
	double PMFHighFreqFact;
	
	//! ** Alpha coefficients of the filter
	std::vector< std::vector<double> > PMFalpha;
	
	//! ** Alpha coefficients of the filter
	std::vector< std::vector<double> > PMFbeta;
	
	/** ********** Variables use in initialization but not after ********** */
	
	/** \brief Time offset */
	double t0;
	
	/** \brief Measurement time step [use in initialization but not after] */
	double dtMes;
	
	/** \brief Physical time step [use in initialization but not after] */
	double dtPhy;
	
	/** \brief Time of the first data [use in initialization but not after] */
	double tFirst;
	
	/** \brief Time of the last data [use in initialization but not after] */
	double tLast;
	
	/** \brief Number of data to add at each run step (measurement time step) \f$ N_{2Add} = \Delta t_{Mes} / \Delta t_{Phy} \f$  */
	int N2Add;
	
	/** \brief Number of data in phasemeter measurements serie [given by LISACode] */
	int NDataPha;
	
	//std::ofstream fCheck;
	
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Default constructor */
	LCOBPM();
	
	/*! \brief Standard constructor */
	LCOBPM(LCTools * MT);
	
	
	/*! \brief Destructor */
	virtual ~LCOBPM();
	
	/*! \brief Initialization at NULL of all base of noise */
	void initNULLBase(bool CleanMem);
	
	/***************  Required methods (to be present in derived classes)  ***************/
	
	/*! \brief Initialization at NULL */
	virtual void initNULL(bool CleanMem);
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="Orbits"> ... </XSIL> */
	virtual void config(ezxml_t xmlbloc);
	
	/*! \brief Configuration of individual parameter */
	virtual void config(int iParam, double ParamVal);
	
	/*! \brief Link the noises */ 
	virtual void LinkNoise(LCNoise ** AllNoises, int NAllNoises);
	
	/*! \brief Link the factor to apply on shot noise (in telescope) */ 
	virtual void LinkFactShotNoise(double * FactShotNoise_n);
	
	/*! \brief Initialization */
	virtual void init();
	
	/*! \brief Make a photodiode measurements using :
	 *	@param[in] t Time
	 */ 
	virtual double MeasurePho(double TimeLocal, double TShiftEff);
	
	/*! \brief Make a phasemeter measurements using :
	 *	@param[in] t Time
	 */ 
	virtual double MeasurePha(double tGlob);
	
	/*! \brief Display information */ 
	virtual void DispInfo(char * BTab, bool LinkedModule);
	
	
	/***************  Local methods  ***************/
	
	/**********  Configuration methods  **********/
	
	/*! \brief Configuration of filter */
	void configFilter(ezxml_t xmlbloc);
	
	/**********  Linking and initalization methods  **********/
	
	/*! \brief Link the orbits to the detector */
	void LinkOrbits(LCOrbits * Orb_n) { orb = Orb_n; } ;
	
	/*! \brief Link arm response to GWs */
	void LinkArmResp(LCArmResp * ArmGW_n) { ArmGW = ArmGW_n; } ;
	
	/*! \brief Link the clock */ 
	void LinkClock(LCUSO ** Clocks, int NClocks);
	
	/*! \brief Initialization of base of noise */
	void initBase();
	
	/*! \brief Link last measurement value to RecLoc (usually the pointer on the current data in the output) */
	void LinkLastMesOutput(double * RecLoc);
	
	
	/**********  Access methods  **********/
	
	/*! \brief Set name and position of the measurments system */
	void setNamePos(const char * NamePos);
	
	/*! \brief Set time informations */
	void setTimeInfo( double t0_n, double dtMes_n, double dtPhy_n, double tFirst_n, double tLast_n, int NDataPha_n);
	
	/*! \brief Return true if the name in argument matches the name of the measurement system */
	bool MatchName(const char * Cmpname);
	
	/*! \brief Compare the nchar first characters of Name_cmp with the Name and return true if there are the same. */
	bool CmpPosName(int iSC_cmp, int IndirectDir_cmp, const char * Name_cmp, int nchar);
	
	/*! \brief Return index of spacecraft */
	int getiSC() {return(iSC);};
	/*! \brief Return index of the linked spacecraft */
	int getiSCLink() {return(iSCLink);};
	
	/*! \brief Return the pointer on phasemeter measurements (example : for TDI) */
	LCSerie2 * getPhaMes() {return(MesPha); };
	
	
	/**********  Running methods  **********/
	
	/**********  Other methods  **********/
	
	/** \brief Return the signal corresponding to GWs .
	 *	It return 0 if the noise is NULL
	 *	@param[in] t	time  
	 */
	double sGW(double t);
	
	/** \brief Return the noise value corresponding to iN and the delay .
	 *	It return 0 if the noise is NULL
	 *	@param[in] iN		index of the noise in the list of noise
	 *	@param[in] tDelay	time delay refering to the global time 
	 */
	double gN(int iN, double tDelay);
	
	/** \brief Return the factor to be applied on the shot noise.
	 *	@param[in] trec		Time when the beam is received on the spacecraft
	 */
	double gFSN(double trec);
	
	/** \brief Display the basic information */
	void DispInfoBase(char * BTab, bool LinkedModule);
	
};

#endif //__LCOBPM_H

/**\}*/

// end of LISACODE-OBPM.h