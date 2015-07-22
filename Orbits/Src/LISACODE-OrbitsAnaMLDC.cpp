/*
 *  LISACODE-OrbitsAnalytic.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 10/04/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-OrbitsAnaMLDC.h"

// *********************
// ***  Constructor  ***
// *********************

/**
 *
 */
LCOrbitsAnaMLDC::LCOrbitsAnaMLDC()
: LCOrbits()
{
	MT = new LCTools;
	MT->LocTools = true;
	initNULL(false);
}


LCOrbitsAnaMLDC::LCOrbitsAnaMLDC(LCTools * MT_n)
: LCOrbits(MT_n)
{
	MT = MT_n;
	initNULL(false);
}


LCOrbitsAnaMLDC::~LCOrbitsAnaMLDC()
{
	initNULL(true);
}




// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCOrbitsAnaMLDC::initNULL(bool CleanMem)
{	
	NSC = 3;
    
    arot = 2./3.*M_PI;
    
    e_mldc = 0.00964838;
    
	sqrt_3 = sqrt(3.0);
	i32 = 1.0/32.0;
	r1532 = 15.0/32.0;
	pi3o2 = 3.0*M_PI/2.0;
	pi2o3 = 2.0*M_PI/3.0;
	pi4o3 = 4.0*M_PI/3.0;
    
    Rgc = LC::au_m;
	omega = 2.*M_PI/LC::Yr_SI;
    
}


void LCOrbitsAnaMLDC::config(ezxml_t xmlbloc)
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
			
			if(MT->wcmp(tmpStr,"Static")){
				move = 0;
                Approx = 0;
                OrderArm = -1;
            }
			if(MT->wcmp(tmpStr,"Moving")){
				move = 1;
                Approx = 1;
                OrderArm = -1;
            }
			if(MT->wcmp(tmpStr,"Eccentric")){
				move = 1;
                Approx = 2;
                OrderArm = -2;
            }
				
			if(tmpStr != NULL)
				MT->Free(tmpStr, (strlen(tmpStr)+1) * sizeof(char));
		}
	}
}


void LCOrbitsAnaMLDC::config(int iParam, double ParamVal)
{
	
}


void LCOrbitsAnaMLDC::init()
{
	initBase();
    
    L0s = L0m/LC::c_SI;
	
	
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


double LCOrbitsAnaMLDC::ArmSpecific(int em, int rec, double trec)
{
	double xi0(rot0+pi3o2+t0/omega);
	double omegat(omega*trec);
	int signArm;
	double deltaArm(xi0);
	
	// THERE IS STILL A DIFFERENCE WITH SYNTHETICLISA : NEED TO UNDERSTAND IT
	
	if(em == 1){
		if(rec == 2){
			deltaArm += pi4o3;
			signArm = 1.0;
		}
		if(rec == 3){
			deltaArm += pi2o3;
			signArm = -1.0;
		}
	}
	if(em == 2){
		if(rec == 3){
			signArm = 1.0;
		}
		if(rec == 1){
			deltaArm += pi4o3;
			signArm = -1.0;
		}
	}
	if(em == 3){
		if(rec == 1){
			deltaArm += pi2o3;
			signArm = 1.0;
		}
		if(rec == 2){
			signArm = -1.0;
		}
	}		
	
	if(OrderArm == -1){
		// ** Rigid
		//cout << signArm << " " << signArm*omega*Rgc/c_SI << " " << omegat - deltaArm  << " " << -L0s * (1.0 + (signArm*omega*Rgc/c_SI) * sin(omegat - deltaArm)) << endl;
		return ( -L0s * (1.0 + (signArm*omega*Rgc/LC::c_SI) * sin(omegat - deltaArm) )  );
	}else{
		// ** Eccentric
		//if((trec>-335)&&(trec<-295))
		//	cout << trec << " " << em << "->" << rec << " : " << signArm << " " << signArm*omega*Rgc/c_SI << " " << omegat - deltaArm  << " " << L0s << " " << i32*e_mldc*sin(3.0*(omegat-xi0)) << " " << (signArm*omega*Rgc/c_SI - r1532*e_mldc) << " " << -L0s * (1.0 + i32*e_mldc*sin(3.0*(omegat-xi0)) + (signArm*omega*Rgc/c_SI - r1532*e_mldc) * sin(omegat - deltaArm) )  << endl;
		return ( -L0s * (1.0 + i32*e_mldc*sin(3.0*(omegat-xi0)) + (signArm*omega*Rgc/LC::c_SI - r1532*e_mldc) * sin(omegat - deltaArm) )  );
	}
}


LCVector LCOrbitsAnaMLDC::position(int iSC, double t)
{
	LCVector r(MT);
	double beta;
	double alpha;
	double c_alpha,s_alpha,c_beta,s_beta;
	double sq_e_mldc;
	
	alpha=omega*(t0+move*t);
	c_alpha=cos(alpha); 
	s_alpha=sin(alpha);
	beta=rot[iSC-1];
	c_beta=crot[iSC-1] ;
	s_beta=srot[iSC-1];
	sq_e_mldc=pow(e_mldc,2);
	
	// TO DO : IMPLEMENT Static(0) and Rigid(1) (see SyntheticLISA doc)
	
	if(move<=1){
		//cout << "1st order in excentricity" << endl;
		// ** 1st order in excentricity  : Eccentric
		
		r.p[0] = LC::au_m*(c_alpha+e_mldc*(s_alpha*c_alpha*s_beta-(1+s_alpha*s_alpha)*c_beta));
		r.p[1] = LC::au_m*(s_alpha+e_mldc*(s_alpha*c_alpha*c_beta-(1+c_alpha*c_alpha)*s_beta));
		r.p[2] = -LC::au_m*e_mldc*sqrt_3*cos(alpha-beta);
		
	}else{
		// ** 2nd order in excentricity  : Eccentric2
		//cout << "2nd order in excentricity " << endl;		
		
		/*
		 r.p[0] =   0.5 * au_m * e_mldc * ( cos(2.0*alpha-beta) - 3.0*cos(beta) )
		 + 0.125 * au_m * sq_e_mldc * ( 3.0*cos(3.0*alpha-2.0*beta) - 5.0*( 2.0*cos(alpha)+cos(alpha-2.0*beta) ) )
		 + au_m * cos(alpha);
		 
		 r.p[1] =   0.5 * au_m * e_mldc * ( sin(2.0*alpha-beta) - 3.0*sin(beta) )
		 + 0.125 * au_m * sq_e_mldc * ( 3.0*sin(3.0*alpha-2.0*beta) - 5.0*( 2.0*sin(alpha)-sin(alpha-2.0*beta) ) )
		 + au_m * sin(alpha);
		 
		 r.p[2] = - sqrt_3 * au_m * e_mldc * cos(alpha-beta) 
		 + sqrt_3 * au_m * sq_e_mldc * ( cos(alpha-beta)*cos(alpha-beta) + 2.0*sin(alpha-beta)*sin(alpha-beta) );
         */
		
		r.p[0] =  LC::au_m * ( 0.5 * e_mldc * ( cos(2.0*alpha-beta) - 3.0*c_beta )
                          + 0.125 * sq_e_mldc * ( 3.0*cos(3.0*alpha-2.0*beta) - 5.0*( 2.0*c_alpha+cos(alpha-2.0*beta) ) )
                          +  c_alpha ) ;
		
		r.p[1] =   LC::au_m * ( 0.5 * e_mldc * ( sin(2.0*alpha-beta) - 3.0*s_beta)
                           + 0.125 * sq_e_mldc * ( 3.0*sin(3.0*alpha-2.0*beta) - 5.0*( 2.0*s_alpha-sin(alpha-2.0*beta) ) )
                           +  s_alpha ) ;
		
		r.p[2] = - LC::au_m * ( sqrt_3 * e_mldc * cos(alpha-beta) 
                           + sqrt_3 * sq_e_mldc * ( cos(alpha-beta)*cos(alpha-beta) + 2.0*sin(alpha-beta)*sin(alpha-beta) ) );
	}
	
	return r;
}


LCVector LCOrbitsAnaMLDC::velocity(int iSC, double t)
{
	LCVector v(MT);
    double beta;
	double alpha;
	double c_alpha,s_alpha,c_beta,s_beta,d_alpha;
	
	alpha=omega*(t0+move*t);
	d_alpha=omega;
	c_alpha=cos(alpha);  
	s_alpha=sin(alpha);
	beta=rot[iSC-1];  
	c_beta=crot[iSC-1] ; 
	s_beta=srot[iSC-1];
	
	
	// TO DO : IMPLEMENT Static(Approx=0, OriginalLISA()) and Rigid(Approx=1, CircularRotating(0.0,1.5*pi,-1)) (see SyntheticLISA doc)
	
	
	if(Approx<=1){
		//cout << "1st" << endl;
		// ** 1st order in excentricity : Eccentric
	    v.p[0] = LC::au_m*(-s_alpha+e_mldc*((c_alpha*c_alpha-s_alpha*s_alpha)*s_beta-2*s_alpha*c_alpha*c_beta))*d_alpha;
	    v.p[1] = LC::au_m*( c_alpha+e_mldc*((c_alpha*c_alpha-s_alpha*s_alpha)*c_beta+2*c_alpha*s_alpha*s_beta))*d_alpha;
	    v.p[2] = LC::au_m*e_mldc*sqrt_3*sin(alpha-beta)*d_alpha;
	    //cout << "Vel old : t,nb,x,y,z =" << t << "  "<<nb<<"  "<<v.p[0]<<"  "<<v.p[1]<<"  "<<v.p[2]<< endl;
	    //vx_old=v.p[0] ; vy_old=v.p[1] ; vz_old=v.p[2] ;
	}else{
		// ** 2nd order in excentricity  : Eccentric2
		//cout << "2nd" << endl;
		v.p[0] = 0.5*LC::au_m*e_mldc*(-2*sin(2.0*alpha-beta))
	    +0.125*LC::au_m*pow(e_mldc,2)*(-9.0*sin(3.0*alpha-2.0*beta)-5.0*(-2.0*sin(alpha)-sin(alpha-2.0*beta)))
	    -LC::au_m*sin(alpha);
		v.p[0] = v.p[0]*d_alpha ;
		
		v.p[1] = 0.5*LC::au_m*e_mldc*(2*cos(2.0*alpha-beta))
	    +0.125*LC::au_m*pow(e_mldc,2)*(9.0*cos(3.0*alpha-2.0*beta)-5.0*(2.0*cos(alpha)-cos(alpha-2.0*beta)))
	    +LC::au_m*cos(alpha);
		v.p[1] = v.p[1]*d_alpha ;
		
		v.p[2] = +sqrt_3*LC::au_m*e_mldc*sin(alpha-beta)
	    +sqrt_3*LC::au_m*pow(e_mldc,2)*(-2*cos(alpha-beta)*sin(alpha-beta)+4.0*sin(alpha-beta)*cos(alpha-beta));
		v.p[2] = v.p[2]*d_alpha ;
	}
	
	/*	cout << "Vel : t,nb,x,y,z =" << t << "  "<<nb<<"  "<<v.p[0]<<"  "<<v.p[1]<<"  "<<v.p[2]<< endl;
     if (t > 1000.){
     throw invalid_argument("GeometryAnalytic::velocity : On arrete !");
     }*/
    
	return v;
}

void LCOrbitsAnaMLDC::DispInfo(char * BTab)
{
	if(MT->Disp()){
		Cout << BTab << "Analytic MLDC orbits :" << Endl;
		DispInfoBase(BTab);
		Cout << "    - Initial rotation      = " << rot0  << Endl;
		Cout << "    - Armlength             = " << L0m << " m" << Endl;
        Cout << "    - Approximation         = " << Approx << Endl;
        
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


// end of LISACODE-OrbitsAnalytic.cpp