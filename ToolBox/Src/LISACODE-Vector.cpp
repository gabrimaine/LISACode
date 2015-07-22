// $Id:  Exp $
/*
 *  LISACode-Vector.cpp
 *  LISACode V 2.0
 *
 *  Created on 08/04/11 by Antoine PETITEAU (AEI)
 *  Last modification on 08/04/11 by Antoine PETITEAU (AEI)
 *
 */

#include "LISACODE-Vector.h"

// *****************
// *  Constructor  *
// *****************

LCVector::LCVector()
{
	MT = new LCTools;
	MT->LocTools = true;
	for(int k=0; k<3; k++)
		p[k] = 0.;
}


LCVector::LCVector(LCTools * MT_n)
{
	MT = MT_n;
	for(int k=0; k<3; k++)
		p[k] = 0.;
}


LCVector::LCVector(LCTools * MT_n, double x, double y, double z)
{
	MT = MT_n;
	p[0] = x;
	p[1] = y;
	p[2] = z;
}


LCVector::LCVector(LCTools * MT_n, double pRef[3])
{
	MT = MT_n;
	for(int k=0; k<3; k++)
		p[k] = pRef[k];
}


LCVector::~LCVector()
{
	
}

// ********************
// *  Initialization  *
// ********************

void LCVector::setTools(LCTools * MT_n)
{
	if(MT->LocTools){
		delete MT;
		MT = MT_n;
	}
}

// ***************
// *  Operators  *
// ***************


void LCVector::operator=(const LCVector &a)
{
	for(int k=0; k<3; k++)
		p[k] = a.p[k];
}


LCVector LCVector::operator+(const LCVector &a)
{
	LCVector r(MT);
	for(int k=0; k<3; k++)
		r.p[k] = p[k] + a.p[k];
	return r;
}


LCVector LCVector::operator-(const LCVector &a)
{
	LCVector r(MT);
	for(int k=0; k<3; k++)
		r.p[k] = p[k] - a.p[k];
	return r;
}


LCVector LCVector::operator*(const double &v)
{
	LCVector r(MT);
	for(int k=0; k<3; k++)
		r.p[k] = p[k] * v;
	return r;
}


LCVector LCVector::operator/(const double &v)
{
	LCVector r(MT);
	for(int k=0; k<3; k++)
		r.p[k] = p[k] / v;
	return r;
}


void LCVector::operator+=(const LCVector &a)
{
	for(int k=0; k<3; k++)
		p[k] +=  a.p[k];
}


void LCVector::operator-=(const LCVector &a)
{
	for(int k=0; k<3; k++)
		p[k] -=  a.p[k];
}


void LCVector::operator*=(const double &v)
{
	for(int k=0; k<3; k++)
		p[k] *=  v;
}


void LCVector::operator/=(const double &v)
{
	for(int k=0; k<3; k++)
		p[k] /=  v;
}


double LCVector::operator*(const LCVector &a)
{
	return(scalar(a));
}


double LCVector::operator()(const int k)
{
	return p[k];
}


void LCVector::operator()(const int k, const double v)
{
	p[k] = v;
}


std::ostream &operator<<(std::ostream &output, const LCVector &a)
{
	output << "[ " << a.p[0] << " , " << a.p[1] << " , " << a.p[2] << " ]";
	return output;
}


double LCVector::ph()
{
	double Res;
	Res = atan2( y(), x() );
	while(Res < 0.)
		Res += 2.*M_PI;
	return(Res);
}


// *****************
// * Other methods *
// *****************

double LCVector::norm()
{
	double r(0.);
	for(int k=0;k<3;k++)
		r += p[k]*p[k];
	//Cout << r << Endl;
	return(sqrt(r));
}


LCVector LCVector::unit()
{
	double vn(norm());
	LCVector r(MT);
	for(int k=0;k<3;k++)
		r.p[k] = p[k]/vn;
	return(r);
}


double LCVector::scalar(const LCVector &a)
{
	double r(0.);
	for(int k=0; k<3; k++)
		r += p[k] * a.p[k];
	return r;
}


LCVector LCVector::cross(const LCVector &a)
{
	LCVector r(MT);
	
	r.p[0] = p[1]*a.p[2] - p[2]*a.p[1];
	r.p[1] = p[2]*a.p[0] - p[0]*a.p[2];
	r.p[2] = p[0]*a.p[1] - p[1]*a.p[0];
	
	return(r);
}


void LCVector::DispMathematica()
{
	Cout << "{";
	for(int i=0; i<3; i++){
		if(i==0)
			Cout <<  p[i];
		else
			Cout << "," <<  p[i];
	}
	Cout << "}" << Endl;
}










// end of LISACODE-LCVector.cpp



