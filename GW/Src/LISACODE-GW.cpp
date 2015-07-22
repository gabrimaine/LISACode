/*
 *  LISACODE-GW.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 23/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-GW.h"


// *********************
// ***  Constructor  ***
// *********************

LCGW::LCGW()
{
	MT = new LCTools;
	MT->LocTools = true;
#ifdef _DEBUG_GW_    
    DEBUGfCheck = NULL;
    DEBUGDispAll = false;
#endif 
	initNULLBase(false);
}


LCGW::LCGW(LCTools * MT_n)
{
	MT = MT_n;
#ifdef _DEBUG_GW_    
    DEBUGfCheck = NULL;
    DEBUGDispAll = false;
#endif 
	initNULLBase(false);
}


LCGW::~LCGW()
{
	initNULLBase(true);
}


void LCGW::initNULLBase(bool CleanMem)
{
	if(CleanMem){
		
		if(Name != NULL)
			MT->Free(Name, 128*sizeof(char));
	}
	
	Name = NULL;
	
	hBpLast = 0.;
	hBcLast = 0.;
	tLast = LC::DBLMINALLOW;
	tStore = LC::PRECISION;
	
	
	tAskMin = 0.;
	tDOrbMax = 0.;
	
	NParams = 2;
	Beta = 0.;
	Lambda = 0.;
    thetaJ = 0.;
    phiJ = 0.;
	
	Polarization = 0.;
	ComputePolarization = false;
	
	FreqMin = 1.e-5;
	FreqMax = 1.;
	
	MemUse = 0.;
    
    ParameterReaded = false;
    
	
#ifdef _DEBUG_GW_
	if(DEBUGfCheck!=NULL){
		DEBUGfCheck->close();
		delete DEBUGfCheck;
	}
	DEBUGfCheck = NULL;
#endif     

	
    
}





// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCGW::initNULL(bool CleanMem)
{
	
}


// ***************************
// *  Configuration methods  *
// ***************************


void LCGW::config(ezxml_t xmlbloc)
{
	
}


void LCGW::config(int iParam, double ParamVal)
{
	
}


// ********************
// *  Access methods  *
// ********************


double LCGW::getParam(int iP)
{
	return(0.);
}


void LCGW::setParam(int iP, double Param_n)
{
	
}


void LCGW::getRange(int iP, double &Pmin, double &Pmax)
{
	
}


double LCGW::getDelta(int iP)
{
	//! \todo TODO properly
	return(1.);
}


void LCGW::setSpecialParam(int iPS, double SpecParam_n)
{
	
}


double LCGW::getSpecialParam(int iPS)
{
	return(0.);
}


void LCGW::setTimeInfo(double t0_n, double dt_n, double TObs_n, double tAskMin_n, double tDOrbMax_n,  double tMaxDiff)
{
    t0Real = t0_n; 
	tAskMin = tAskMin_n;
	tDOrbMax = tDOrbMax_n;
}


void LCGW::RandParam(int iP)
{
	
}


// ***************************************
// * Linking and initialization methods  *
// ***************************************

int LCGW::init()
{
	initBase();
	
	return 0;
}


// *********************
// *  Running methods  *
// *********************


void LCGW::Computehpc(double t)
{
	
}




// *******************
// *  Other methods  *
// *******************

void LCGW::DispInfo(char * BTab)
{
	if(MT->Disp()){
		DispInfoBase(BTab);
	}
}

void LCGW::DispAllParam(std::ostream * out)
{
	for(int iP=0; iP<NParams; iP++)
		(*out) << " " << getParam(iP);
}


void LCGW::DispAllParamName(std::ostream * out)
{
	(*out) << "Beta Lambda Polarization";
}



// ***********************
// ***  Local mehtods  ***
// ***********************

// ***************************
// *  Configuration methods  *
// ***************************

void LCGW::configBase(ezxml_t xmlbloc)
{
	//! *** Set the name
	setName(ezxml_attr(xmlbloc,"Name"));
	
	ezxml_t param;
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
		
        ParameterReaded = false;
        
		//! *** Read ecliptic lattitude 
		if(MT->wcmp(ezxml_attr(param,"Name"),"EclipticLatitude")){
            ParameterReaded = true;
			Beta = MT->gXMLAngle(param);
        }
		
		//! *** Read ecliptic colattitude 
		if(MT->wcmp(ezxml_attr(param,"Name"),"EclipticColatitude")){
            ParameterReaded = true;
			Beta = M_PI/2. - MT->gXMLAngle(param);
        }
		
		//! *** Read ecliptic longitude
        if(MT->wcmp(ezxml_attr(param,"Name"),"EclipticLongitude")){
            ParameterReaded = true;
			Lambda = MT->gXMLAngle(param);
        }
		
		//! *** Read polarization
        if((MT->wcmp(ezxml_attr(param,"Name"),"Polarization"))||(MT->wcmp(ezxml_attr(param,"Name"),"Polarisation"))){
            ParameterReaded = true;
			Polarization = MT->gXMLAngle(param);
        }
        
        if(ParameterReaded)
            Cout << " " << ezxml_attr(param,"Name") ; 
	}
}


// ********************
// *  Access methods  *
// ********************


void LCGW::RandAllParams()
{
	for(int iP=0; iP<NParams; iP++)
		RandParam(iP);
}


void LCGW::setName(const char * Name_n)
{	
	if(Name == NULL)
		Name = (char *) MT->AllocMemory(128*sizeof(char));
	
	int NChar(MIN(strlen(Name_n), 127));
	for(int i=0; i<NChar; i++)
		Name[i] = Name_n[i];
	Name[NChar] = '\0';	
}


// ***************************************
// * Linking and initialization methods  *
// ***************************************

void LCGW::initBase()
{
	if(fabs(Polarization)>LC::PRECISION){
		ComputePolarization = true;
		c2Pol = cos(2.*Polarization);
		s2Pol = sin(2.*Polarization);
	}
}


// *********************
// *  Running methods  *
// *********************


void LCGW::ComputehBpc(double t)
{
	Computehpc(t);
	
	tLast = t;
	
	if(ComputePolarization){
		double hSp, hSc;
		hSp = hBpLast;
		hSc = hBcLast;
        //! WARNING : TO BE CALRIFIED IF THE 2 NEXT LINES ARE NOT MAKING A BUG FOR SPINBBH SOURCE : THERE CANCELED BY SOFIANE BUT DECANCELLED BY ANTOINE FOR GALBIN
		hBpLast =   c2Pol * hSp + s2Pol * hSc;     
		hBcLast = - s2Pol * hSp + c2Pol * hSc;    
	}
	
}

double LCGW::hBp(double t)
{
	if(fabs(t-tLast)>tStore)
		ComputehBpc(t);
	
	return(hBpLast);
}


double LCGW::hBc(double t)
{
	if(fabs(t-tLast)>tStore)
		ComputehBpc(t);
	
	return(hBcLast);
}


void LCGW::hBpc(double t, double & hBpV, double & hBcV)
{
	if(fabs(t-tLast)>tStore){
#ifdef _DEBUG_GW_     
        (*DEBUGfCheck) << t;
#endif

		ComputehBpc(t);

#ifdef _DEBUG_GW_     
        (*DEBUGfCheck) << " " << hBpLast << " " << hBcLast << Endl;
#endif
    }
	
	hBpV = hBpLast;
	hBcV = hBcLast;
}



// *******************
// *  Tools methods  *
// *******************


double LCGW::halfhann(double t, double t0, double t1){
	
 	double w;
	
 	if(t > t0 && t > t1){
 		if(t0 > t1) return 0;
 		else return 1;
 	} else if(t < t0 && t < t1){
 		if(t0 < t1) return 0;
 		else return 1;
 	} else {
 		w = cos(M_PI*(t-t1)/(2.0*(t0-t1)));
 		return w*w;
 	}
	
}

// *******************************
//     Added by Sofiane

void LCGW::ConvertLS1S2dirNR2SSBSofVer(double & Phd_o, 
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
								 LCVector & S2B_o
								 )
{
	
	LCVector kB(MT, -sin(ThS_i)*cos(PhS_i), -sin(ThS_i)*sin(PhS_i), - cos(ThS_i) );
    LCVector ZB(MT, 0., 0., 1. );
    LCVector JB(MT, sin(thetaJ_i)*cos(phiJ_i), sin(thetaJ_i)*sin(phiJ_i), cos(thetaJ_i) );
	
	//! *** Compute the Source frame (xN, yN, JN) expressed in in SSB  
	LCVector ZBcJB(MT);
	LCVector xNB(MT);
	LCVector yNB(MT);
	ZBcJB =  JB.cross(ZB);
	yNB =  ZBcJB / ZBcJB.norm();
	xNB = yNB.cross(JB);
    
    double kx, ky;
    
    kx = kB.scalar(xNB);
    ky = kB.scalar(yNB);
    
    Phd_o = acos( kx/sqrt(kx*kx + ky*ky)  );
    if(Phd_o < 0.)
        Phd_o += 2.*M_PI;
	
		
	//! *** Apply the rotation on L, S1, S2
    

    
    
    
	LnB_o = ( xNB * LnN_i.p[0] ) + ( yNB * LnN_i.p[1] ) + ( JB * LnN_i.p[2] );
	S1B_o = ( xNB * S1N_i.p[0] ) + ( yNB * S1N_i.p[1] ) + ( JB * S1N_i.p[2] );
    S2B_o = ( xNB * S2N_i.p[0] ) + ( yNB * S2N_i.p[1] ) + ( JB * S2N_i.p[2] );
    
    

    
	
	}




//     End added by Sofiane
// *******************************




void LCGW::ConvertLS1S2dirNR2SSB(double Thd_i, 
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
								 double & PolRecompute)
{
	
	LCVector kN(MT, sin(Thd_i)*cos(Phd_i), sin(Thd_i)*sin(Phd_i), cos(Thd_i) );
	
	//! *** Compute {k,p,q} expressed in NR frame 
	LCVector kNcJN(MT);
	LCVector pN(MT);
	LCVector qN(MT);
	kNcJN = kN.cross(JN_i);
	pN = kNcJN / kNcJN.norm();
	qN = kN.cross(pN);
	
	
	//! *** Compute rotation matrix of {k,p,q} expressed in NR frame
	LCMatrix Rpn(MT,3,3);
	Rpn.s(0,0, kN(0) );
	Rpn.s(0,1, kN(1) );
	Rpn.s(0,2, kN(2) );
	Rpn.s(1,0, pN(0) );
	Rpn.s(1,1, pN(1) );
	Rpn.s(1,2, pN(2) );
	Rpn.s(2,0, qN(0) );
	Rpn.s(2,1, qN(1) );
	Rpn.s(2,2, qN(2) );
	
	//! *** Compute rotation matrix of {k,p,q} expressed in NR frame
	LCMatrix Rnp(MT,3,3);
	Rnp = Rpn.Inv(1);
	
	
	//! *** Compute rotation matrix of SSB unit vectors expressed in {k,p,q}
	double th, ph ;
	th = M_PI/2. - Beta;
	ph = Lambda;
	LCVector kB(MT, -(cos(ph)*sin(th)), -(sin(th)*sin(ph)), -cos(th) );
	LCMatrix Rpb(MT,3,3);
	Rpb.s(0,0, -sin(th)*cos(ph) );
	Rpb.s(0,1, -sin(th)*sin(ph) );
	Rpb.s(0,2, -cos(th) );
	Rpb.s(1,0, cos(th)*cos(ph)*cos(Polarization) - sin(ph)*sin(Pol_i) );
	Rpb.s(1,1, cos(th)*cos(Pol_i)*sin(ph) + cos(ph)*sin(Pol_i) );
	Rpb.s(1,2, -(cos(Pol_i)*sin(th)) );
	Rpb.s(2,0, cos(Pol_i)*sin(ph) + cos(th)*cos(ph)*sin(Pol_i) );
	Rpb.s(2,1, -(cos(ph)*cos(Pol_i)) + cos(th)*sin(ph)*sin(Pol_i) );
	Rpb.s(2,2, -(sin(th)*sin(Pol_i)) );
	
	LCMatrix Rbp(MT,3,3);
	Rbp.s(0,0, -(cos(ph)*sin(th)) );
	Rbp.s(0,1, cos(th)*cos(ph)*cos(Pol_i) - sin(ph)*sin(Pol_i) );
	Rbp.s(0,2, cos(Pol_i)*sin(ph) + cos(th)*cos(ph)*sin(Pol_i) );
	Rbp.s(1,0, -(sin(th)*sin(ph)) );
	Rbp.s(1,1, cos(th)*cos(Pol_i)*sin(ph) + cos(ph)*sin(Pol_i) );
	Rbp.s(1,2, -(cos(ph)*cos(Pol_i)) + cos(th)*sin(ph)*sin(Pol_i) );
	Rbp.s(2,0, -cos(th) );
	Rbp.s(2,1, -(cos(Pol_i)*sin(th)) );
	Rbp.s(2,2, -(sin(th)*sin(Pol_i)) );
	
	
	//! *** Compute rotation matrix of SSB unit vectors expressed in NR frame
	LCMatrix Rnb(MT,3,3);
	Rnb = Rnp * Rpb ;
	
	LCVector pB1(MT,Rpb(1,0),Rpb(1,1),Rpb(1,2));
	LCVector qB1(MT,Rpb(2,0),Rpb(2,1),Rpb(2,2));
	LCVector pB2(MT), qB2(MT);
	pB2 = Rnb.T() * pN;
	qB2 = Rnb.T() * qN;
	
	if(MT->DispDet()){
		//! ** For checking
		Cout << "kN = " << kN << Endl;
		Cout << "JN = " << JN_i << Endl;
		Cout << "kN x JN = " << kNcJN << Endl;
		Cout << "pN = " << pN <<  " --> norm=" << pN.norm() << Endl;
		Cout << "qN = " << qN <<  " --> norm=" << qN.norm() << Endl;
		
		LCVector kp(MT, 1.,0.,0.);
		LCVector pp(MT, 0.,1.,0.);
		LCVector qp(MT, 0.,0.,1.);
		Cout << "Rpn = " << Rpn << Endl;
		Cout << "kN = " << kN << Endl;
		Cout << "kN test1 = Rpn.T * kp = Rpn.T * " << kp << " = " << Rpn.T() * kp << Endl;
		Cout << "pN = " << pN << Endl;
		Cout << "pN test1 = Rpn.T * pp = Rpn.T * " << pp << " = " << Rpn.T() * pp << Endl;
		Cout << "qN = " << qN << Endl;
		Cout << "qN test1 = Rpn.T * qp = Rpn.T * " << qp << " = " << Rpn.T() * qp << Endl;
		
		Cout << "Rnp = " << Rnp << Endl;
		Cout << "Rpn = " << Rpn << Endl;
		Cout << "Rnp*Rpn = " << Rnp*Rpn << Endl;
		
		Cout << "Rpb = " << Rpb << Endl;
		Cout << "RpbInv = " << Rbp << Endl;
		Cout << "Rpb*RpbInv = " << Rpb*Rbp << Endl;
		Cout << "Rnb = " << Rnb << Endl;
		Cout << "kB ref = " << kB << Endl; 
		Cout << "kB test2 = Rbn * kN = " << Rnb.T() * kN << Endl;
		Cout << "pB1 = " << pB1 << Endl;
		Cout << "pB2 = " << pB2 << Endl;
		Cout << "qB1 = " << qB1 << Endl;
		Cout << "qB2 = " << qB2 << Endl;
	}	
	
	//! *** Apply the rotation on L, S1, S2
	LnB_o = Rnb.T() * LnN_i;
	S1B_o = Rnb.T() * S1N_i;
	S2B_o = Rnb.T() * S2N_i;
	
	
	if(MT->DispDet()){
		
		//! *** Checking by recomputing polarization
		LCVector JB(MT);
		double thBJ, phBJ;
		double up, down;
		JB = LnB_o*AmpL_i + S1B_o*AmpS1_i + S2B_o*AmpS2_i ;
		Cout << "JB = " << JB << " (unit=" << JB.unit() << ")" << Endl;
		Cout << "L  = " << AmpL_i << "  " << LnB_o << " ] " << Endl;
		Cout << "S1 = " << AmpS1_i << "  " << S1B_o << " ] " << Endl;
		Cout << "S2 = " << AmpS2_i << "  " << S2B_o << " ] " << Endl;
		thBJ = acos(JB.z()/JB.norm()); 
		if(thBJ<0.)
			thBJ += M_PI; 
		if((fabs(thBJ) <= 1.e-6) || (fabs(thBJ - M_PI) <= 1.e-6) )
			phBJ = 0.0;
		else
			phBJ = atan2(JB.y(),JB.x());
		up = cos(th)*sin(thBJ)*cos(ph-phBJ) - cos(thBJ)*sin(th);
		down = sin(thBJ)*sin(ph - phBJ);
		PolRecompute = atan2(up,down);
		PolRecompute += M_PI;
		
		LCVector uB(MT, cos(th)*cos(ph), cos(th)*sin(ph), -sin(th) );
		LCVector vB(MT, sin(ph), -cos(ph), 0. );
		LCVector kBcJB(MT);
		kBcJB = kB.cross(JB);
		
		
		Cout << "tan(PolRef)   = tan(" << Pol_i << ")   = " << tan(Pol_i) << Endl;
		Cout << "tan(PolTest1) = tan(" << PolRecompute << ")   = " << tan(PolRecompute) << Endl;
		Cout << "tan(PolTest2) = v.(kxJ)/u.(kxJ)   = " << -(vB*kBcJB)/(uB*kBcJB) << "  = " << -(vB*kBcJB) << " / " << (uB*kBcJB) << Endl;
		
		//ComputeThdPhdPol();
		
		Cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << Endl;
	}
}







//    added by Sofiane
//  $$$$$$$$$$$$$$$$$$$


void LCGW::ComputeThdPhdPolSofVer(LCVector LnB_i,
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
                            double  phiJ_i
                            )
{

	
	double th, ph ;
	th = M_PI/2. - Beta;
	ph = Lambda;
		
	double up, down;

    Pol_o = Polarization;
    
    Thd_o = acos(-cos(thetaJ_i)*cos(th) - sin(thetaJ_i)*sin(th)*cos(ph - phiJ_i));
    
    LCVector kB(MT, -sin(th)*cos(ph), -sin(th)*sin(ph), - cos(th) );
    LCVector ZB(MT, 0., 0., 1. );
    LCVector JBnew(MT, sin(thetaJ_i)*cos(phiJ_i), sin(thetaJ_i)*sin(phiJ_i), cos(thetaJ_i) );
	
	//! *** Compute the Source frame (xN, yN, JN) expressed in in SSB  
	LCVector ZBcJB(MT);
	LCVector xNB(MT);
	LCVector yNB(MT);
	ZBcJB =  JBnew.cross(ZB);
	yNB =  ZBcJB / ZBcJB.norm();
	xNB = yNB.cross(JBnew);
    
    double kx, ky;
    
    kx = kB.scalar(xNB);
    ky = kB.scalar(yNB);
    
    Phd_o = acos( kx/sqrt(kx*kx + ky*ky)  );
    if(Phd_o < 0.)
        Phd_o += 2.*M_PI;
    
    
    up =  cos(Phd_o - ph);    // introduced by Sofiane
    down = cos(th)*sin(Phd_o - ph);  // introduced by Sofiane
    
    double PolB = atan2(up, down);
    if(PolB<0.)
        PolB += M_PI;  // changed by Sofiane*2
    Polarization = PolB; 
    
}



//  $$$$$$$$$$$$$$$$$$$$$
//   end added by Sofiane





void LCGW::ComputeThdPhdPol(LCVector LnB_i,
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
							double & Pol_o)
{
	double Thd_old(Thd_o), Phd_old(Phd_o), Pol_old(Pol_o);
	
	double th, ph ;
	th = M_PI/2. - Beta;
	ph = Lambda;
	
	
	//! *** Compute rotation matrix for SSB frame to NR using the fact that we have expression of L, S1 and S2 in both
	LCMatrix Rlb(MT,3,3) , Rln(MT,3,3) , Rbl(MT,3,3) , Rbn(MT,3,3) ;
	LCVector kB(MT, -(cos(ph)*sin(th)), -(sin(th)*sin(ph)), -cos(th) );
	LCVector kN(MT);
	
	
	
	Rlb.s(0,0, LnB_i(0) );
	Rlb.s(0,1, LnB_i(1) );
	Rlb.s(0,2, LnB_i(2) );
	Rlb.s(1,0, S1B_i(0) );
	Rlb.s(1,1, S1B_i(1) );
	Rlb.s(1,2, S1B_i(2) );
	Rlb.s(2,0, S2B_i(0) );
	Rlb.s(2,1, S2B_i(1) );
	Rlb.s(2,2, S2B_i(2) );
	
	Rbl = Rlb.Inv(1);
	
	Rln.s(0,0, LnN_i(0) );
	Rln.s(0,1, LnN_i(1) );
	Rln.s(0,2, LnN_i(2) );
	Rln.s(1,0, S1N_i(0) );
	Rln.s(1,1, S1N_i(1) );
	Rln.s(1,2, S1N_i(2) );
	Rln.s(2,0, S2N_i(0) );
	Rln.s(2,1, S2N_i(1) );
	Rln.s(2,2, S2N_i(2) );
	
	Rbn = Rbl * Rln ; 
	
	kN = Rbn.T() * kB;
	
	Thd_o = kN.th();
	Phd_o = kN.ph();
	
	if(MT->DispDet()){	
		//! *** For checking
		Cout << "Rlb = " << Rlb << Endl;
		Cout << "Rbl = " << Rbl << Endl;
		Cout << "Rln = " << Rln << Endl;
		Cout << "Rbn = " << Rbn << Endl;
		Cout << "kB = " << kB << Endl;
		Cout << "kN = " << kN << Endl;
		Cout << "new Thd = " << Thd_o << "  ( old = " << Thd_old << ")" << Endl;
		Cout << "new Phd = " << Phd_o << "  ( old = " << Phd_old << ")" << Endl;
	}
	
	//! *** Compute polarization
	LCVector JB(MT);
	double thBJ, phBJ;
	double up, down;
	JB = LnB_i*AmpL_i + S1B_i*AmpS1_i + S2B_i*AmpS2_i ;
	thBJ = acos(JB(2)/JB.norm()); 
	if(thBJ<0.)
        thBJ += M_PI; 
    if((fabs(thBJ) <= 1.e-6) || (fabs(thBJ - M_PI) <= 1.e-6) )
		phBJ = 0.0;
	else
		phBJ = atan2(JB(1),JB(0));
	up = cos(th)*sin(thBJ)*cos(ph-phBJ) - cos(thBJ)*sin(th);
	down = sin(thBJ)*sin(ph - phBJ);
	Polarization = atan2(up,down);
	Polarization += M_PI;
	
	LCVector uB(MT, cos(th)*cos(ph), cos(th)*sin(ph), -sin(th) );
	LCVector vB(MT, sin(ph), -cos(ph), 0. );
	LCVector kBcJB(MT);
	kBcJB = kB.cross(JB);
	
	Cout << "tan(PolOld) = " << "tan(" << Pol_old << ")   = " << tan(Pol_old) << Endl;
	Cout << "tan(PolNew) = " << "tan(" << Pol_o << ")   = "<< tan(Pol_o) << Endl;
	
	
}


int LCGW::ReadNRParam(char * FileName_i,
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
					  LCVector & LnN_o)
{
	char tmpLine[16384];
	double thS1tmp, phS1tmp, thS2tmp, phS2tmp, thLtmp, phLtmp;
	
	bool UseInitAng(false), UseInitXYZ(false);
	
	std::ifstream fInMain(FileName_i);
	
	if(!fInMain)
		return 1;
	
	std::vector<char*> Words(0);
	
	//! **** Read fixed parameters 
	int iLine(0);
	do{
		iLine++;
		fInMain.getline(tmpLine,16384,'\n');
		MT->wextract(tmpLine, Words);
		
		if(Words.size()<1){
			Cout << "ERROR on LCGW::LoadNRdata : Probelm of reading : we need at least one word on the line " << iLine <<" !" << Endl;
			throw std::invalid_argument("ERROR on LCGW::LoadNRdata : Probelm of reading : we need at least one word on the line !");
		}
		
		if((!MT->wcmp(Words[0],"Harmonics"))&&(Words.size()<2)){
			Cout << "ERROR on LCGW::LoadNRdata : Probelm of reading : we need at least 2 words (ex: 'x = 1.' ) on the line " << iLine << " (except if the first one is 'Harmonics') !" << Endl;
			throw std::invalid_argument("ERROR on LCGW::LoadNRdata : Probelm of reading : we need at least 2 words (ex: 'x = 1.' ) on the line (except if the first one is 'Harmonics') !");
		}
		
		if(MT->wcmp(Words[0],"q")){
			mrat_o = ReadInfoNRdata(Words);
			eta_o = mrat_o/pow(mrat_o+1.,2.);
		}
		
		if((MT->wcmp(Words[0],"a1"))||(MT->wcmp(Words[0],"chi1")))
			chi1_o = ReadInfoNRdata(Words);
		
		if((MT->wcmp(Words[0],"a2"))||(MT->wcmp(Words[0],"chi2")))
			chi2_o = ReadInfoNRdata(Words);
		
		if(MT->wcmp(Words[0],"Phi0"))
			Phi0_o = ReadInfoNRdata(Words);
		
		if(MT->wcmp(Words[0],"phistart"))
			Phi0_o = ReadInfoNRdata(Words);
		
		if(MT->wcmp(Words[0],"FreqMax"))
			FreqMaxM_o = ReadInfoNRdata(Words);
		
		if(MT->wcmp(Words[0],"omega_orb"))
			omMInit_o = ReadInfoNRdata(Words);
		
		if(MT->wcmp(Words[0],"tstart"))
			toMInitHyb_o = ReadInfoNRdata(Words);
		
		
		if((MT->wcmp(Words[0],"thS1"))||(MT->wcmp(Words[0],"thS2"))||(MT->wcmp(Words[0],"thLN"))||(MT->wcmp(Words[0],"phS1"))||(MT->wcmp(Words[0],"phS2"))||(MT->wcmp(Words[0],"phLN"))){
			
			UseInitAng = true;
			if(UseInitXYZ){
				throw std::invalid_argument("ERROR on LCGWSpinBBHNR1::LoadNRdata : You already start to define initial condiftion using x,y,z so you cannot used angles now !" );
				Cout << "ERROR on LCGWSpinBBHNR1::LoadNRdata : You already start to define initial condiftion using x,y,z so you cannot used angles now !" << Endl;
			}
			
			if(MT->wcmp(Words[0],"thS1"))
				thS1tmp = ReadInfoNRdata(Words);
			
			if(MT->wcmp(Words[0],"thS2"))
				thS2tmp = ReadInfoNRdata(Words);
			
			if(MT->wcmp(Words[0],"thLN"))
				thLtmp = ReadInfoNRdata(Words);
			
			if(MT->wcmp(Words[0],"phS1"))
				phS1tmp = ReadInfoNRdata(Words);
			
			if(MT->wcmp(Words[0],"phS2"))
				phS2tmp = ReadInfoNRdata(Words);
			
			if(MT->wcmp(Words[0],"phLN"))
				phLtmp = ReadInfoNRdata(Words);
		}
		
		if((MT->wcmp(Words[0],"S1x"))||(MT->wcmp(Words[0],"S1y"))||(MT->wcmp(Words[0],"S1z"))||(MT->wcmp(Words[0],"S2x"))||(MT->wcmp(Words[0],"S2y"))||(MT->wcmp(Words[0],"S2z"))||(MT->wcmp(Words[0],"LNx"))||(MT->wcmp(Words[0],"LNy"))||(MT->wcmp(Words[0],"LNz"))){
			
			UseInitXYZ = true;
			if(UseInitAng){
				throw std::invalid_argument("ERROR on LCGWSpinBBHNR1::LoadNRdata : You already start to define initial condiftion using angles so you cannot used x,y,z now !");
				Cout << "ERROR on LCGWSpinBBHNR1::LoadNRdata : You already start to define initial condiftion using angles so you cannot used x,y,z now !" << Endl;
			}
			
			if(MT->wcmp(Words[0],"S1x"))
				S1N_o.x(ReadInfoNRdata(Words));
			
			if(MT->wcmp(Words[0],"S1y"))
				S1N_o.y(ReadInfoNRdata(Words));
			
			if(MT->wcmp(Words[0],"S1z"))
				S1N_o.z(ReadInfoNRdata(Words));
			
			
			if(MT->wcmp(Words[0],"S2x"))
				S2N_o.x(ReadInfoNRdata(Words));
			
			if(MT->wcmp(Words[0],"S2y"))
				S2N_o.y(ReadInfoNRdata(Words));
			
			if(MT->wcmp(Words[0],"S2z"))
				S2N_o.z(ReadInfoNRdata(Words));
			
			
			if(MT->wcmp(Words[0],"LNx"))
				LnN_o.x(ReadInfoNRdata(Words));
			
			if(MT->wcmp(Words[0],"LNy"))
				LnN_o.y(ReadInfoNRdata(Words));
			
			if(MT->wcmp(Words[0],"LNz"))
				LnN_o.z(ReadInfoNRdata(Words));
		}	
		
		
	}while((!MT->wcmp(Words[0],"Harmonics"))&&(!fInMain.eof()));
	
	
	fInMain.close();
	fInMain.clear();
	
	//Cout << "S1N = " << S1N << Endl;
	//Cout << "S2N = " << S2N << Endl;
	//Cout << "LnN = " << LnN << Endl;
	
	if(UseInitAng){
		S1N_o.x( sin(thS1tmp)*cos(phS1tmp) );
		S1N_o.y( sin(thS1tmp)*sin(phS1tmp) );
		S1N_o.z( cos(thS1tmp) );
		S2N_o.x( sin(thS2tmp)*cos(phS2tmp) );
		S2N_o.y( sin(thS2tmp)*sin(phS2tmp) );
		S2N_o.z( cos(thS2tmp) );
		LnN_o.x( sin(thLtmp)*cos(phLtmp) );
		LnN_o.y( sin(thLtmp)*sin(phLtmp) );
		LnN_o.z( cos(thLtmp) );
	}
	
	/*
	 if(UseInitXYZ){
	 thS1 = S1N_o.th();
	 phS1 = S1N_o.ph();
	 thS2 = S2N_o.th();
	 phS2 = S2N_o.ph();
	 thL = LnN_o.th();
	 phL = LnN_o.ph();
	 }
	 */
	
	for(int iW=0; iW<Words.size(); iW++)
		MT->Free(Words[iW],256*sizeof(char));
	
	return 0;
}


double LCGW::ReadInfoNRdata(std::vector<char*> Words)
{
	double Res(0.);
	if(Words.size()<2){
		throw std::invalid_argument("ERROR in LCGWSpinBBHNR1::ReadInfoNRdata : We need more than one word in information line of NRdata file.");
		Cout << "ERROR in LCGWSpinBBHNR1::ReadInfoNRdata : We need more than one word in information line of NRdata file." << Endl;
	}
	if(MT->wcmp(Words[1],"=")){
		if(Words.size()<2){
			throw std::invalid_argument("ERROR in LCGWSpinBBHNR1::ReadInfoNRdata : We need more than two words in information line of NRdata file when '=' is used.");
			Cout << "ERROR in LCGWSpinBBHNR1::ReadInfoNRdata : We need more than two words in information line of NRdata file when '=' is used." << Endl;
		}
		Res = atof(Words[2]);
	}else{
		Res = atof(Words[1]);
	}
	return Res;
}




// *******************
// *  Other methods  *
// *******************

double LCGW::AddDeltaPar(int iP, double FactDelta)
{
	double Delta(getDelta(iP));
	double Param(getParam(iP));
	
	setParam(iP, Param + FactDelta*Delta);
	
	return(Delta);
}



void LCGW::DispInfoBase(char * BTab)
{
	if(MT->Disp()){
		Cout << BTab << "\t- Sky position : beta = " << Beta << " rad ,  lambda = " << Lambda << " rad" << Endl;
		Cout << BTab << "\t- Polarization  = " << Polarization << " rad" << Endl;
	}
	//std::cerr << Beta << " " << Lambda << Endl;
}





// end of LISACODE-GW.cpp