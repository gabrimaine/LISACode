/*
 *  LISACODE-OrbitsAnaLISA.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 10/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-OrbitsAnaLISA.h"

// *********************
// ***  Constructor  ***
// *********************

/**
 *
 */
LCOrbitsAnaLISA::LCOrbitsAnaLISA()
: LCOrbits()
{
	MT = new LCTools;
	MT->LocTools = true;
	initNULL(false);
}


LCOrbitsAnaLISA::LCOrbitsAnaLISA(LCTools * MT_n)
: LCOrbits(MT_n)
{
	MT = MT_n;
	initNULL(false);
}


LCOrbitsAnaLISA::~LCOrbitsAnaLISA()
{
	initNULL(true);
}




// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCOrbitsAnaLISA::initNULL(bool CleanMem)
{	
	NSC = 3;
	alpha = 0.;
	nu = 0.;
	tmu = 0.;
	cmu = 0.;
	smu = 0.;
	e = 1.;
	sqrtee = 0.;
	arot = 0.;
	
	Rgc = LC::au_m;
	omega = 2.*M_PI/LC::Yr_SI;
}


void LCOrbitsAnaLISA::config(ezxml_t xmlbloc)
{
	configBase(xmlbloc);
	ezxml_t param;
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		if(MT->wcmp(ezxml_attr(param,"Name"),"InitialRotation"))
			rot0 = MT->gXMLAngle(param);
		if(MT->wcmp(ezxml_attr(param,"Name"),"InitialPosition")){
			double tmpInitPos(0.0);
			tmpInitPos = MT->gXMLAngle(param);
			t0 = tmpInitPos*LC::Yr_SI/(2.0*M_PI);
		}
		if(strcmp(ezxml_attr(param,"Name"),"Armlength")==0){
			L0m = MT->gXMLLength(param);
		}
		if(strcmp(ezxml_attr(param,"Name"),"OrderDelay")==0){
			OrderArm = atoi(ezxml_txt(param));
		}
		if(strcmp(ezxml_attr(param,"Name"),"OrbitApproximation")==0){
			char * tmpStr(NULL);
			MT->stripcopy((*param).txt, tmpStr);
			
			if(MT->wcmp(tmpStr,"Static"))
				move = 0;
			if(MT->wcmp(tmpStr,"Moving"))
				move = 1;
			if(MT->wcmp(tmpStr,"Eccentric"))
				move = 1;
			if(tmpStr != NULL)
				MT->Free(tmpStr, (strlen(tmpStr)+1) * sizeof(char));
		}
	}
}


void LCOrbitsAnaLISA::config(int iParam, double ParamVal)
{
	
}


void LCOrbitsAnaLISA::init()
{
	initBase();
	
	/*! *** Initialization of orbital parameters :
	 *  \f{eqnarray*}{
	 *	\alpha & = & { L_0 \over 2 R } \\
	 *	\nu & = & {\pi \over 3} + {5 \over 8} \alpha \\
	 *	t_{\mu} & = & \alpha { \sin{\nu} \over \sin{\pi/3} } + \alpha \cos{\nu} \\
	 *	c_{\mu} & = & { 1 \over 1 + t_{\mu}^2 } \\
	 *	s_{\mu} & = & t_{\mu} \ c_{\mu} \\
	 *	e & = & \sqrt{ 1 + {4 \over \sqrt{3}} \cos{\nu} } \\
	 *  sqrtee & = & \sqrt{1-e^2} \\
	 *	a_{rot} & = & { 2 \pi \over 3 }
	 *	\f}
	 */
	
	alpha = L0m/(2*Rgc);
	nu  = (M_PI/3.)+(5./8.)*alpha;
	tmu = alpha*sin(nu)/(sin(M_PI/3.)+alpha*cos(nu));
	cmu = 1./sqrt(1.+tmu*tmu);
	smu = tmu*cmu;
	e  = sqrt(1.+(4./sqrt(3.)*cos(nu)*alpha)+(4./3.*alpha*alpha))-1.;
	sqrtee = sqrt(1.-e*e);
	arot = 2./3.*M_PI;
	
	
	/*! *** Computation of spacecraft phase \f$ rot_i = i \ a_{rot} + rot_0 \f$ 
	 * and associated cosinus and sinus */
	
	rot.resize(3);
	crot.resize(3);
	srot.resize(3);
	for(int iS=0; iS<3; iS++){
		rot[iS]  = iS * arot + rot0;
		crot[iS] = cos(rot[iS]);
		srot[iS] = sin(rot[iS]);
	}
	
	
	/*! *** Store the value corresponding at t = 0 */
	
	tStorePos = 0.;
	for(int i=1; i<4; i++)
		SCposStore[i-1] = position(i,0.0);
	
	tStoreArm = 0.;
	for(int i=1; i<4; i++){
		ArmStore[i-1] = ArmCompute( i, (i+1)%3+1, 0.0);
		ArmStore[i+2] = ArmCompute( i, i%3+1, 0.0);
	}	
}


double LCOrbitsAnaLISA::ArmSpecific(int em, int rec, double trec)
{
	//*! Not used
	return(L0m/LC::c_SI);
}


LCVector LCOrbitsAnaLISA::position(int iSC, double t)
{
	LCVector r(MT);
	double cpsi,spsi;
	
	exanom(iSC, t, cpsi, spsi);
	
	/*! *** Computation of position :
	 *	\f$ \vec{r}	= \left( \begin{array}{l}
	 *					Rgc \cdot (cpsi-e) \cdot cmu \cdot crot_{nb-1} - sqrtee  \cdot spsi \cdot srot_{nb-1} \\
	 *					Rgc \cdot (cpsi-e) \cdot cmu \cdot srot_{nb-1} - sqrtee  \cdot spsi \cdot crot_{nb-1} \\ 
	 *					- Rgc \cdot smu \cdot (cpsi-e) 
	 *				\end{array} \right) 
	 *				= \left( \begin{array}{l}
	 *					R (\cos{\psi}-e) \ c_{\mu} \ \cos{rot_i} - sqrtee  \sin{\psi} \sin{rot_i} \\
	 *					R (\cos{\psi}-e) \ c_{\mu} \ \sin{rot_i} - sqrtee  \sin{\psi} \cos{rot_i} \\ 
	 *					- R s_{\mu} \cdot (\cos{\psi} - e) 
	 *			\end{array} \right) 
	 *  \f$
	 */
	
	r.p[0] = Rgc*((cpsi - e)*cmu*crot[iSC-1] - sqrtee*spsi*srot[iSC-1]);
	r.p[1] = Rgc*((cpsi - e)*cmu*srot[iSC-1] + sqrtee*spsi*crot[iSC-1]);
	r.p[2] = -Rgc*smu*(cpsi - e);
	
	return r;
}


LCVector LCOrbitsAnaLISA::velocity(int iSC, double t)
{
	LCVector v(MT);
	double cpsi,spsi;
	double psidot;
	
	exanom(iSC, t, cpsi, spsi);
	
	
	/*! *** Computation of velocity :
	 *	\f$  \dot{\psi} = { \Omega R \over 1 - e \cos{\psi} }\f$
	 *	and
	 *	\f$ \vec{v}	= \left( \begin{array}{l}
	 *				\dot{\psi} \cdot (-spsi \cdot cmu \cdot crot_{nb-1} - sqrtee  \cdot cpsi \cdot srot_{nb-1}) \\
	 *				\dot{\psi} \cdot (-spsi \cdot cmu \cdot srot_{nb-1} + sqrtee  \cdot cpsi \cdot crot_{nb-1}) \\ 
	 *				\dot{\psi} \cdot smu \cdot spsi 
	 *			\end{array} \right) 
	 *			= \left( \begin{array}{l}
	 *				\dot{\psi} (-spsi  c_{\mu}  \cos{rot_{iSC-1}} - sqrtee  \cos{\psi}  \sin{rot_{iSC-1}} ) \\
	 *				\dot{\psi} (-spsi  c_{\mu}  \sin{rot_{iSC-1}} + sqrtee  \cos{\psi}  \cos{rot_{iSC-1}} ) \\ 
	 *				\dot{\psi} s_{mu}  \sin{\psi}
	 *			\end{array} \right) 
	 *	\f$
	 */
	
	
	psidot = omega * Rgc / (1.-e*cpsi);
	
	v.p[0] = psidot*(-spsi*cmu*crot[iSC-1] - sqrtee*cpsi*srot[iSC-1]);
	v.p[1] = psidot*(-spsi*cmu*srot[iSC-1] + sqrtee*cpsi*crot[iSC-1]);
	v.p[2] = psidot*smu*spsi;
	
	
	return v;
}

void LCOrbitsAnaLISA::DispInfo(char * BTab)
{
	if(MT->Disp()){
		Cout << BTab << "Analytic LISA orbits :" << Endl;
		DispInfoBase(BTab);
		Cout << "    - Initial rotation      = " << rot0  << Endl;
		Cout << "    - Armlength             = " << L0m << " m" << Endl;
	}
}	


// ***********************
// ***  Local mehtods  ***
// ***********************




// ********************
// *  Access methods  *
// ********************




// ********************
// *  Others methods  *
// ********************

void LCOrbitsAnaLISA::exanom(int iSC, double ts, double &cpsi, double &spsi)
{
	double t(t0 + move*ts); 
	double psi;
	double exold, exnew, dex;
	
	/*! *** Computation of initial quantities :
	 *  \f{eqnarray*}{
	 *	\psi & = & \Omega t - rot_i = {2 \pi \over T} t - rot_i \\
	 *	E_{old} & = & \psi + e \sin{\psi} \left( 1 + e \left( \cos{\psi} - e  { 1 - 3 \cos^2{\psi} \over 2 } \right) \right)
	 *	\f}
	 */
	
	psi = omega*t - rot[iSC-1];
	cpsi=cos(psi);
	spsi=sin(psi);
	exold = psi + e * spsi * (1.+e*(cpsi-e*(0.5-1.5*cpsi*cpsi)));
	
	/*! *** Iterative resolution in 3 steps :
	 *  \f{eqnarray*}{
	 *	\Delta E & = & { E_{old} - e sin{\psi} - \psi /over 1 - e \cos{\psi} } \\
	 *	E_{new} & = & E_{old} - \Delta E \\
	 *	E_{old} & = & E_{new}
	 *	\f}
	 */
	
	for(int k=0; k<3; k++){
		cpsi = cos(exold);
		spsi = sin(exold);
		dex = (exold-e*spsi-psi)/(1.-e*cpsi);
		exnew = exold-dex;
		exold = exnew;
	}
}


// end of LISACODE-OrbitsAnaLISA.cpp