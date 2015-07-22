/*
 *  LISACODE-Orbits.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 05/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


#include "LISACODE-Orbits.h"


// *********************
// ***  Constructor  ***
// *********************

/**
 *
 */
LCOrbits::LCOrbits()
{
	MT = new LCTools;
	MT->LocTools = true;
	initNULLBase(false);
}


LCOrbits::LCOrbits(LCTools * MT_n)
{
	MT = MT_n;
	for(int i=0; i<NSCMAX; i++)
		SCposStore[i].setTools(MT);
	initNULLBase(false);
}


LCOrbits::~LCOrbits()
{
	initNULLBase(true);
}



// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCOrbits::initNULL(bool CleanMem)
{	
	
}


void LCOrbits::config(ezxml_t xmlbloc)
{
	
}


void LCOrbits::config(int iParam, double ParamVal)
{
	
}


void LCOrbits::init()
{
	
}


double LCOrbits::ArmSpecific(int em, int rec, double trec)
{
	return(0.0);
}


void LCOrbits::DispInfo(char * BTab)
{

}	


LCVector LCOrbits::position(int iSC, double t)
{
	LCVector r(MT);
	return(r);
}


LCVector LCOrbits::velocity(int iSC, double t)
{
	LCVector r(MT);
	return(r);
}



// ***********************
// ***  Local mehtods  ***
// ***********************


// **************************
// * Configuration methods  *
// **************************

void LCOrbits::configBase(ezxml_t xmlbloc)
{
	ezxml_t param;
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		if(MT->wcmp(ezxml_attr(param,"Name"),"TimeStoringPosition"))
			tRangeStorePos = MT->gXMLTime(param);
		if(MT->wcmp(ezxml_attr(param,"Name"),"TimeStoringArmlength"))
			tRangeStoreArm = MT->gXMLTime(param);
	}
}


// ***************************************
// * Linking and initialization methods  *
// ***************************************


void LCOrbits::initNULLBase(bool CleanMem)
{
	// Nothing to clean in memory at the moment
	//if(CleanMem){
	//}
	
	move = 1;
	t0 = 0.;
	OrderArm = 2;
	
	
	tStorePos = LC::DBLMINALLOW;
	tStoreArm = LC::DBLMINALLOW;
	
	tRangeStorePos = 0.1;
	tRangeStoreArm = 0.1;
	
}



void LCOrbits::initBase()
{
	
}




// ********************
// *  Access methods  *
// ********************



// ***************************************
// *  Methods for time travel along arm  *
// ***************************************

double LCOrbits::ArmCompute(int em, int rec, double trec)
{
	double tA(0.), tij;
	
	if ((em<1)||(em>NARMMAX)||((rec<1)&&(rec>NARMMAX))){
		std::cerr << "ERROR in LCOrbits::ArmCompute : Spacecraft number must be in [0," << NARMMAX << "]  ! " << Endl;
		throw std::invalid_argument ("ERROR in LCOrbits::ArmCompute : Spacecraft number must be in [0,NARMMAX]  ! ");
	}
	/*
	if(OrderArm != -10)
		OrderArm = order_default;
	if ((OrderArm<-2)||(OrderArm>2))
		throw std::invalid_argument ("ERROR in LCOrbits::ArmCompute : Order must be 0 (for 0), 1 (for 1/2), 2 (for 1) or -1 (for analytical MLDC eccentric formulation)");
	 */
	
	if(OrderArm >= 0){
		LCVector ri(MT), rj(MT), vi(MT), vj(MT), rij(MT), nij(MT), n(MT);
		
		/*! **** Read positions and velocities and compute related quantities
		 *	Index i refers to receiver and j to emitter
		 *  \f{eqnarray*}{
		 *	\vec{r}_{ij} & = & \vec{r}_j - \vec{r}_i = \vec{r}_{em} - \vec{r}_{rec}  \\
		 *	\vec{n}_{ij} & = & \vec{r}_{ij} / | \vec{r}_{ij} | \\
		 *	\hat{n}  & = & - \vec{n}_{ij} \\
		 *	t_{ij}   & = & - | \vec{r}_{ij} | / c = - | \vec{r}_{em} - \vec{r}_{rec} | / c
		 *  \f}
		 */ 
		
		ri = position(rec, trec);
		vi = velocity(rec, trec);
		vi = vi/LC::c_SI;
		
		rj = position(em, trec);
		vj = velocity(em, trec);
		vj = vj/LC::c_SI;
		
		/*
		ri.p[0] = 5.e9;
		ri.p[1] = 0.;
		ri.p[2] = 0.;
		rj.p[0] = -1.e-5;
		rj.p[1] = 0.;
		rj.p[2] = 0.;
		rij = rj-ri;
		*/
		
		 
		rij = rj-ri;
		nij = rij.unit();
		n = nij*(-1.);
		
		tij = rij.norm()/LC::c_SI;
		tij = -tij;  
		
		/*
		MT->o->precision(16);
		Cout << "em = " << em << "  -> rec = " << rec << " : " << rij.norm() << " ==> " << tij << Endl;
		for(int i=0; i<3; i++)
			Cout << "\t\t" << ri.p[i] << "\t\t" << rj.p[i] << "\t\t" << rij.p[i] << Endl;   
		*/
		 
		if (OrderArm>=0)
			tA += ArmOrderContrib(0, tij, ri, rj, vi, vj, rij, nij, n);
		
		if (OrderArm>=1)
			tA += ArmOrderContrib(1, tij, ri, rj, vi, vj, rij, nij, n);
		
		if (OrderArm>=2)
			tA += ArmOrderContrib(2, tij, ri, rj, vi, vj, rij, nij, n);
		
		return (tA);
			
	}else{
		return( ArmSpecific(em, rec, trec) );
	}
}



double LCOrbits::ArmOrderContrib(int order, double tij, LCVector ri, LCVector rj, LCVector vi, LCVector vj, LCVector rij, LCVector nij, LCVector n)
{
	double cr(0.);
	double c21, c22, c23, num, denom;
	switch (order) {
		case 0:
			//! *** Order 0 : \f{eqnarray*}{ t^{(0)} & = & t_{ij} \f}
			cr = tij;
			break;
		case 1:
			//! *** Order 1/2 : \f{eqnarray*}{ t^{(1/2)} & = & t_{ij} \; \vec{n} . \vec{v}_j \f}
			cr = tij*(n*vj);
			break;
		case 2:
			/*! *** Order 1 :
			 *  \f{eqnarray*}{
			 *	c_{21} & = & { t_{ij} \over 2 } \left( \vec{v}_j . \vec{v}_j + { \left( \hat{n}.\vec{v}_j \right) }^2  \right) \\
			 *  c_{22} & = & - { t_{ij}^2 \over 2 } \ c \ R_{Schw} \ { \hat{n} . \vec{r}_j \over (\vec{r}_j . \vec{r}_j ) ^{3/2} } \\
			 *  c_{23} & = &  { c \over R_{Schw} } (1 + \Gamma) \log{ \left( \sqrt{ {( c t_{ij} + \hat{n}.\vec{r}_i )}^2 } + \vec{r}_i . \vec{r}_i - {(\vec{r}_i . \hat{n})}^2 \over \hat{n}.\vec{r}_i + \sqrt{ \vec{r}_i . \vec{r}_i } \right) }\\
			 *	t^{(1)} & = & c_{21} + c_{22} + c_{23}
			 *	\f}
			 */
			c21 = 0.5 * tij * ( (vj*vj) + pow((n*vj),2.) );
			c22 = -0.5 * tij*tij * LC::c_SI * LC::RSchw * (n*rj) / pow((rj*rj),1.5);
			num = sqrt( pow(LC::c_SI*tij+(n*ri),2.) + (ri*ri)-pow((ri*n),2.) ) + LC::c_SI*tij+(n*ri);
			denom=(n*ri)+sqrt((ri*ri));
			c23 = (LC::RSchw/LC::c_SI)*(1.+LC::gamma_u)*log(num/denom);
			cr=c21+c22+c23;
			break;
		default:
			throw std::invalid_argument("ERROR in LCOrbits::ArmOrderContrib : Unknow order !");
			break;
	}
	return (cr);
}


double LCOrbits::Arm(int em, int rec, double trec)
{
	//! *** If the difference between the required received time #trec and the last computation time of the time travel is superior than the precision range ... 
	if( fabs( trec*move - tStoreArm ) > tRangeStoreArm){
		if(NSC==3){
			//! ... recompute the time travel
			tStoreArm = trec;
			//! - for arm 2 : arm 2  : (em,rec)=(1,3)
			ArmStore[0] = ArmCompute(1,3,trec);
			//! - for arm 3 : arm 3  : (em,rec)=(2,1)
			ArmStore[1] = ArmCompute(2,1,trec);
			//! - for arm 1 : arm 1  : (em,rec)=(3,2)
			ArmStore[2] = ArmCompute(3,2,trec);
			//! - for arm 6 : arm 3' : (em,rec)=(1,2)
			ArmStore[3] = ArmCompute(1,2,trec);
			//! - for arm 4 : arm 1' : (em,rec)=(2,3)
			ArmStore[4] = ArmCompute(2,3,trec);
			//! - for arm 5 : arm 2' : (em,rec)=(3,1)
			ArmStore[5] = ArmCompute(3,1,trec);
		}else{
			for(int i=0; i<(NSC*(NSC-1)); i++){
				int emC(1+floor(i/(NSC-1)));
				int reC((i+emC)%NSC+((i+emC)%NSC>=emC?1:0));
				//Coutm << i << " : " << emC << " --> " << reC << Endl;
				ArmStore[i]= ArmCompute(emC, reC,trec);
			}
		}
	}
	
		
	//! *** Return the required time travel : index = \f$  b + i_{em} - 1 \quad \textrm{with} \quad b = 0 , \textrm{if} \; mod(i_{em},3) + 1 = i_{rec}, 0 else  \f$ 
	//Cout << "===>>>" << (em%3+1==rec) * 3 + em - 1 << Endl;
	if(NSC==3)
		return(ArmStore[ (em%3+1==rec) * 3 + em - 1]);
	return(ArmStore[(em-1)*(NSC-1)+(rec<em?rec-1:rec-2)]);
}



// *********************
// *  Running methods  *
// *********************

LCVector LCOrbits::Pos(int iSC, double t)
{
	//! *** If the difference between the required received time #trec and the last computation time of the time travel is superior than the precision range ... 
	if( fabs( t*move - tStorePos ) > tRangeStorePos){
		//! ... recompute the
		tStorePos = t;
		for(int i=1; i<=NSC; i++){
			SCposStore[i-1] = position(i,t);
		}
	}
	return(SCposStore[iSC-1]);
}


double LCOrbits::ArmVelocity(int em, int rec, double trec)
{
	LCVector ri(MT), rj(MT), vi(MT), vj(MT), rij(MT), nij(MT), n(MT);	
	
	ri = position(rec, trec);
	vi = velocity(rec, trec);
	
	rj = position(em, trec);
	vj = velocity(em, trec);
	
	rij = rj-ri;
	nij = rij.unit();
	n   = nij * (-1.);
	
	return(vj*n-vi*n);	
}


LCVector LCOrbits::VectNormal(double t)
{
	LCVector n1(MT), n2(MT), r(MT);
	
	n1 = position(2,t)-position(3,t);
	n2 = position(3,t)-position(1,t);
	
	r.p[0] = n1.p[1]*n2.p[2]-n1.p[2]*n2.p[1];
	r.p[1] = n1.p[2]*n2.p[0]-n1.p[0]*n2.p[2];
	r.p[2] = n1.p[0]*n2.p[1]-n1.p[1]*n2.p[0];
	
	
	return (r.unit());
}


// ********************
// *  Others methods  *
// ********************

void LCOrbits::tGWMinMax(double GWbet, double GWlam, double t0_l, double tDur_l, double dt_l,  double & tGWmin, double & tGWmax)
{
	LCVector kGW(MT, -cos(GWbet)*cos(GWlam), -cos(GWbet)*sin(GWlam), -sin(GWbet));
	double dtf(100.*dt_l), t(t0_l);
	double dmin(1.e30), dmax(-1.e30);
	double dtmp;
	do{
		for(int iSC=1; iSC<=NSC; iSC++){
			dtmp = kGW * position(iSC, t);
			if ((dtmp)<dmin)
				dmin = dtmp;
			if (dtmp>dmax)
				dmax = dtmp;
		}
		t += dtf;
	}while(t<tDur_l);
	
	tGWmin = (dmin-2.*L0m)/LC::c_SI;
	tGWmax = dmax/LC::c_SI;
}



void LCOrbits::DispInfoBase(char * BTab)
{
	if(MT->Disp()){
		Cout.precision(13);
		Cout << "    - Number of spacecrafts = " << NSC << Endl; 
		Cout << "    - Initial time          = " << t0 << " s" << Endl;
		Cout << "    - Motion                = " << move << Endl;
		Cout << "    - Arm time travel order = " << OrderArm << Endl;
		Cout << "    - Time range for computing position  = " << tRangeStorePos << " s" << Endl;
		Cout << "    - Time range for computing armlength = " << tRangeStoreArm << " s" << Endl;
	}
}

// end of LISACODE-Orbits.cpp



