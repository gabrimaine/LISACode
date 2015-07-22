/*
 *  LISACODE-ArmResp.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 24/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */



/** \defgroup ArmResp ArmResp
 * This group contains all the things related to the computation of arm response
 * \{
 */

#ifndef __LCARMRESP_H
#define __LCARMRESP_H

#include "LISACODE-Tools.h"
#include "LISACODE-Vector.h"
#include "LISACODE-Orbits.h"
#include "LISACODE-GW.h"



/** \brief Class for computing the arm response.
 * \author A. Petiteau
 * \version 2.0
 * \date 24/04/2011
 *
 * This class compute the arm response to the gravitational waves.\n
 * It is based on the computation of the delay time \f$ t - t_0(t) \f$.
 * It computes the signal either 
 * in relative frequency : \f[ { \delta \nu (t) \over \nu_0 } = { d (t - t_0)  \over dt } \f]
 * or in phase : \f[ \phi(t) = 2 \pi \nu_0 (t - t_0) \f] . \n
 * 
 * For the first order approximation in relative frequency the formulation is simple :
 * \f[
 *    { \delta \nu \over \nu_0 } (t) = { -1 \over 2(1+\hat{k}.\hat{n}(t))}
 *			\left[ H \left( t+{\hat{k}.\vec{r}_{rec}(t) \over c} \right) - H \left( t+{\hat{k}.\vec{r}_{em}(t) \over c} - {L \over c} \right) \right] 
 * \f]
 *	with \f$ \hat{k} \f$ the propagation direction of the gravitational wave, 
 *		\f$ \hat{n}(t) \f$ the direction of the arm,
 *		\f$ \vec{r}_{em}(t) \f$ the position of the emitter,
 *		\f$ \vec{r}_{rec}(t) \f$ the position of the receiver,
 *		\f$ L \f$ the arm length and 
 *	\f[ H(t) = h_{B+} (t) \xi_{+} (\hat{\theta},\hat{\phi},\hat{n}(t)) + h_{B\times} (t) \xi_{\times} (\hat{\theta},\hat{\phi},\hat{n}(t)) \f]
 *  with \f[ \left\{ \begin{array}{lll} 
 *		\xi_{+}(\vec{u}, \vec{v}, \vec{n})  & = & {\left(\vec{u}. \vec{n} \right)}^{2} - {\left(\vec{v}. \vec{n} \right)}^{2} \\
 *		\xi_{\times}(\vec{u}, \vec{v}, \vec{n}) & = & 2 \left(\vec{u}. \vec{n} \right) \left(\vec{v}. \vec{n} \right) \end{array} \right.
 *	\f]
 *
 * For the other cases (phase signal or relative frequency signal for order higher than 1) we have to solve the following equation for \f$ t - t_0 (t) \f$:
 * \f[ t - t_0 = { L \over c} - {1\over2} \int^t_{t_0} dt' 
 *		\left[ H(t_k) + {3\over4} H^2(t_k) + {5\over8} H^3(t_k) + ... {{1\over2} ({1\over2}-1) ({1\over2}-2) ... ({1\over2}-n) \over n!} H^n(t) + o(H^n) \right] 
 *	\f]
 * with :
 * \f{eqnarray*}{ t_k(t) & = & t + {\hat{k}.\vec{r}(t) \over c} \\
 *						 & = & t (1 + \hat{k}.\hat{n}(t)) + { \hat{k}.\vec{r}_{em}(t) \over c} - t_0(t) \hat{k}.\hat{n}(t)
 * \f}
 *
 *	NB : The secon part has not yet been implemented.
 */
class LCArmResp
{
protected:
	
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	/** \brief List of pointers on gravitational waves (extern allocation) */
	 LCGW ** GWs;
	
	/** \brief Number of pointers on gravitational waves */
	int NGWs;
	
	/** \brief Pointer on orbit (extern allocation) */
	LCOrbits * Orb;
	
	
	/** \brief Order of precision in response */
	int order;
	
	/** \brief Type of output data : if true frequency data else phase data */ 
	bool IsFreqOutput;
	
	/** \brief True if we can use the simplied formulation until order 1 for relative frequency signal */
	bool O1Freq;
	
	
	/*********** Variables computed during the initialization and used during the computation  *********/
	
	/** \brief List of direction of the gravitational waves. Size NGWs.
	 *	\f[ \hat{k} = - \hat{w} =	\left( \begin{array}{c} 
	 *			\sin \theta_k \cos \phi_k \\ 
	 *			\sin \theta_k  \sin \phi_k \\ 
	 *			\cos \theta_k
	 *	\end{array} \right) = \left( \begin{array}{c} 
	 *			- \cos \beta \cos \lambda \\ 
	 *			- \cos \beta  \sin \lambda \\ 
	 *			- \sin \beta
	 *	\end{array} \right) \f]
	 */
	LCVector ** k;
	
	/** \brief List of unit vector \f$ \hat{u} \f$ for describing the barycentric frame \f$ (\hat{u}, \hat{v}, \hat{k}) \f$ associated to the gravitational waves. Size NGWs. 
	 * \f[ \hat{u} =  \hat{\theta} = {\partial \hat{k}\over \partial \theta_k} = -{\partial \hat{w}\over \partial \beta} = \left( \begin{array}{c} 
	 *		\cos \theta_k \cos \phi_k \\ 
	 *		\cos \theta_k  \sin \phi_k \\ 
	 *		- \sin \theta
	 *	\end{array} \right) = \left( \begin{array}{c} 
	 *		\sin \beta \cos \lambda \\ 
	 *		\sin \beta  \sin \lambda \\ 
	 *		- \cos \beta
	 *	\end{array} \right) \f]
	 */
	LCVector ** u;
	
	/** \brief List of unit vector \f$ \hat{v} \f$ for describing the barycentric frame \f$ (\hat{u}, \hat{v}, \hat{k}) \f$ associated to the gravitational waves. Size NGWs 
	 * \f[ \hat{v} = \hat{\phi} = {1 \over \sin \theta_k} {\partial \hat{k}\over \partial \phi_k} = {-1  \over \cos \beta} {\partial \hat{w}\over \partial \lambda} = \left( \begin{array}{c} 
	 *		- \sin \phi_k\\ 
	 *		\cos \phi_k\\ 
	 *		0
	 *	\end{array} \right) = \left( \begin{array}{c} 
	 *		\sin \lambda\\ 
	 *		- \cos \lambda\\ 
	 *		0
	 *	\end{array} \right) \f]
	 */
	LCVector ** v;
	
#ifdef _DEBUG_ARMRESP_     
	std::ofstream * DEBUGfCheck;
#endif	
    
public:
	
	/*************** Constructor ***************/
	
	/*! \brief Standard constructor */
	LCArmResp(LCTools * MT_n);
	
	/*! \brief Destructor */
	~LCArmResp();
	
	/*! \brief Initialization at NULL */
	void initNULL(bool CleanMem);
	
	
	/**********  Configuration methods  **********/
	
	
	/**********  Linking and initalization methods  **********/
	
	/*! \brief Link the list of gravitationnal waves */
	void LinkGWs(LCGW** GWs_n, int NGWs_n);
	
	/*! \brief Link the orbits */
	void LinkOrb(LCOrbits* Orb_n);
	
	/*! Initialization */
	void init();
	
	/**********  Running methods  **********/
	
	/*! \brief Compute and return the gravitaional waves signal on the arm   
	 *	@param[in] em	index of emitter [1,2,3]
	 *	@param[in] em	index of receiver [1,2,3]
	 *	@param[in] trec	time at reception
	 */
	double gS(int em, int rec, double trec);
	
	/*! \brief For order 1 and relative frequency signal, compute and return the gravitaional waves signal on the arm   
	 *	@param[in] em	index of emitter [1,2,3]
	 *	@param[in] em	index of receiver [1,2,3]
	 *	@param[in] trec	time at reception
	 */
	double SigO1Freq(int em, int rec, double trec);
	
	/*! \brief For computation at any order for phase or relative frequency signal, compute and return the gravitaional waves signal on the arm   
	 *	@param[in] em	index of emitter [1,2,3]
	 *	@param[in] em	index of receiver [1,2,3]
	 *	@param[in] trec	time at reception
	 */
	double SigOAll(int em, int rec, double trec);
	
	
	/**********  Other methods  **********/
	
	/*! \brief Display informations */
	void DispInfo(char * BTab);
	
	
};

#endif //__LCARMRESP_H

/**\}*/

// end of LISACODE-ArmResp.h