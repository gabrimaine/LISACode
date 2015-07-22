/*
 *  LISACODE-ArmResp.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 24/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-ArmResp.h"


// *********************
// ***  Constructor  ***
// *********************

LCArmResp::LCArmResp(LCTools * MT_n)
{
	MT = MT_n;
#ifdef _DEBUG_ARMRESP_    
    DEBUGfCheck = NULL;
#endif 
	initNULL(false);
}


LCArmResp::~LCArmResp()
{
	initNULL(true);
}


void LCArmResp::initNULL(bool CleanMem)
{
	if(CleanMem){
		if(k!=NULL)
			MT->Free(k, NGWs*sizeof(LCVector*));
		if(u!=NULL)
			MT->Free(u, NGWs*sizeof(LCVector*));
		if(v!=NULL)
			MT->Free(v, NGWs*sizeof(LCVector*));
	}
	
#ifdef _DEBUG_ARMRESP_
    if (DEBUGfCheck != NULL) {
        DEBUGfCheck->close();
        delete DEBUGfCheck;
    }
	DEBUGfCheck = NULL;
#endif    
    
	k = NULL;
	u = NULL;
	v = NULL;
	
	GWs = NULL;
	NGWs = 0;
	Orb = NULL;
	
	order = 1;
	IsFreqOutput = true;
	O1Freq = true;
    
	
}


// ***************************
// *  Configuration methods  *
// ***************************



// ***************************************
// * Linking and initialization methods  *
// ***************************************


void LCArmResp::LinkGWs(LCGW** GWs_n, int NGWs_n)
{
	GWs = GWs_n;
	NGWs = NGWs_n;
}


void LCArmResp::LinkOrb(LCOrbits* Orb_n)
{
	Orb = Orb_n;
}


void LCArmResp::init()
{
	double bet,lam;
	double cb,sb,cl,sl;
	
	
	//! **** Allocation and initialization of barycentric frame
	k = (LCVector**) MT->AllocMemory(NGWs*sizeof(LCVector*));
	u = (LCVector**) MT->AllocMemory(NGWs*sizeof(LCVector*));
	v = (LCVector**) MT->AllocMemory(NGWs*sizeof(LCVector*));
	for(int i=0; i<NGWs; i++){
		k[i] = new LCVector(MT);
		u[i] = new LCVector(MT);
		v[i] = new LCVector(MT);
	}
	
	for(int iGW=0; iGW<NGWs; iGW++){
		bet = GWs[iGW]->getSkyBet();
		lam = GWs[iGW]->getSkyLam();
		cb = cos(bet);
		sb = sin(bet);
		cl = cos(lam);
		sl = sin(lam);
		
		//! *** Computation of the barycentric reference frame associated from the source direction
		
		//! ** Computation of direction of propagation (see #k definition for details)
		k[iGW]->p[0] = -cb*cl;
		k[iGW]->p[1] = -cb*sl;
		k[iGW]->p[2] = -sb;
		
		//! ** Computation of first unit vector (see #u definition for details)
		u[iGW]->p[0] = sb*cl;
		u[iGW]->p[1] = sb*sl;
		u[iGW]->p[2] = -cb;
		
		//! ** Computation of first unit vector (see #v definition for details)
		v[iGW]->p[0] = sl;
		v[iGW]->p[1] = -cl;
		v[iGW]->p[2] = 0;
		
		/*
		Cout << "In GW constructor : lambda = " <<  lam << "  ,  beta = " <<  bet << Endl;
		Cout << "  - k = " << k[iGW]->p[0] << " " << k[iGW]->p[1] << " " << k[iGW]->p[2] << " " << Endl;
		Cout << "  - u = " << u[iGW]->p[0] << " " << u[iGW]->p[1] << " " << u[iGW]->p[2] << " " << Endl;
		Cout << "  - v = " << v[iGW]->p[0] << " " << v[iGW]->p[1] << " " << v[iGW]->p[2] << " " << Endl;	
		*/
	}
    
#ifdef _DEBUG_ARMRESP_    
	char fNCheck[512];
    int iRCheck(MT->ifloor(MT->RandUniform(0, 10000.)));
    sprintf(fNCheck,"CheckArm_%d.txt",iRCheck);
    std::cerr << "DEBUG:ArmResp : File = " << fNCheck << Endl;
    DEBUGfCheck = new std::ofstream(fNCheck);
#endif
    
}


// ********************
// *  Access methods  *
// ********************



// *********************
// *  Running methods  *
// *********************

double LCArmResp::gS(int em, int rec, double trec)
{
	if(O1Freq)
		return(SigO1Freq(em,rec,trec));
	else 
		return(SigOAll(em,rec,trec));
}






double LCArmResp::SigO1Freq(int em, int rec, double trec)
{
	/*! First order approximation in relative frequency : computation of the following formula :
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
	*		\xi_{+}(\vec{u}, \vec{v}, \vec{n})  & = & {\left(\vec{u}. \vec{n} \right)}^2 - {\left(\vec{v}. \vec{n} \right)}^2 \\
	*		\xi_{\times}(\vec{u}, \vec{v}, \vec{n}) & = & 2 \left(\vec{u}. \vec{n} \right) \left(\vec{v}. \vec{n} \right) \end{array} \right.
	*	\f]
	*/
	
	//! **** Compute quantities associated to the detector
	
	LCVector rrec(MT), rem(MT), vArm(MT), n(MT);
	double tem, tr, te;
	double hpr, hcr, hpe, hce;
	
	//! ** Time at emission : \f$ t_{em} = t_{rec} - L/c \f$
	tem = trec + Orb->Arm(em, rec, trec);
	
	//! ** Position of receiver at reception time
	rrec = Orb->Pos(rec, trec)/LC::c_SI;
	
	//! ** Position of emitter at reception time (is it OK ?)
	rem = Orb->Pos(em, trec)/LC::c_SI;  
	
	//! ** Unit vector of the arm
	vArm = rrec - rem;
	n = vArm.unit();
	
	//! **** Loop on gravitaional waves 
	double un, vn, kn;
	double xipo2, xico2, krrec, krem; 
	double sig(0.);
	for(int iGW=0; iGW<NGWs; iGW++){
		
		//! *** The tree scalar products between the frame and the arm unit vector : \f$ \hat{u}.\hat{n} \f$ , \f$ \hat{v}.\hat{n} \f$ and \f$ \hat{k}.\hat{n} \f$
		un = (*u[iGW])*n;
		vn = (*v[iGW])*n;
		kn = (*k[iGW])*n; 
		
		//! *** Compute \f[ { \xi_{+} \over 2} (t_r)  = { {\left(\vec{u}. \vec{n}(t_r) \right)}^2 - {\left(\vec{v}. \vec{n}(t_r) \right)}^2 \over 2 } \f]
		xipo2 = 0.5*(un*un-vn*vn);
		
		//! *** Compute \f[ {\xi_{\times} \over 2} (t_r) = \left(\vec{u}. \vec{n}(t_r) \right) \left(\vec{v}. \vec{n}(t_r) \right) \f]
		xico2 = un*vn;
		
		//! *** Scalar product between the direction of propagation and the spacecraft positions \f$ \hat{k}.\vec{r}_{rec} \f$ and \f$ \hat{k}.\vec{r}_{rec} \f$ 
		//! for computing the time \f$ t_r, t_e \f$ when the gravitational wave act on each spacecraft 
		krrec = (*k[iGW])*rrec;
		krem = (*k[iGW])*rem;
		tr = trec + krrec;
		te = tem + krem;
		
		//! *** Compute the 2 components of the gravitational wave strain at time when it act on each spacecraft
		GWs[iGW]->hBpc(tr, hpr, hcr);
		GWs[iGW]->hBpc(te, hpe, hce);
		
		/*
		if(((tr>1.8e6)&&(tr<1.9e6))||((tr>0.8e6)&&(tr<0.9e6))||((tr>0.2e6)&&(tr<0.3e6))){
			MT->o->precision(12);
			Cout << tr << " " << hpr << " " << hcr << " " << te << " " << hpe << " " << hce << Endl; 
		}
		 */
		
		
		/*! *** Compute the signal which is :
		 * \f[ s^{GW} = { ( h_{B+}(t_e) - h_{B+}(t_r) ) {\xi_{+} \over 2} (t_r) + ( h_{B\times}(t_e) - h_{B\times}(t_r) ) {\xi_{\times} \over 2} (t_r)  \over 1 + \hat{k}.\hat{n}  }
		 * \f]
		 *	NB : We consider \f$ \xi_{+,\times} (t_e) \sim \xi_{+,\times} (t_r) \f$ which corresponds to \f$ \hat{n}(t_e) \sim \hat{n} (t_r) \f$
		 */
		sig += ( (hpe-hpr) * xipo2 + (hce-hcr) * xico2 ) / (1. + kn) ;
		
		//if((rec==1)&&(em==2))
		//	DEBUGfCheck << trec << " " << tem << " " << un << " " << vn << " " << kn << " " << xipo2 << " " << xico2 << " " << krrec << " " << krem << " " << tr << " " << te << " " << hpr << " " << hcr << " " << hpe << " " << hce << " " << sig << Endl; 
		
	}
    
#ifdef _DEBUG_ARMRESP_     
    (*DEBUGfCheck) << trec << " " << sig << " " << em << " " << rec << Endl;
#endif
	
	return(sig);
	
}


double LCArmResp::SigOAll(int em, int rec, double trec)
{
	/*! Cases of phase signal or relative frequency signal for order higher than 1 : solve the following equation for \f$ (t-t_0) \f$:
	 * \f[ t-t_0 = { L \over c} - {1\over2} \int^t_{t_0} dt' 
	 *		\left[ H(t_k) + {3\over4} H^2(t_k) + {5\over8} H^3(t_k) + ... {{1\over2} ({1\over2}-1) ({1\over2}-2) ... ({1\over2}-n) \over n!} H^n(t) + o(H^n) \right] 
	 *	\f]
	 * with :
	 * \f{eqnarray*}{ t_k(t) & = & t + {\hat{k}.\vec{r}(t) \over c} \\
	 *						 & = & t (1 + \hat{k}.\hat{n}(t)) + { \hat{k}.\vec{r}_{em}(t) \over c} - t_0 \hat{k}.\hat{n}(t)
	 * \f}
	 */
	
	//! \todo TODO ... and it won't be easy !
	
	return(0.);
}


// *******************
// *  Other methods  *
// *******************


// end of LISACODE-ArmResp.cpp