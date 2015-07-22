// $Id:  Exp $
/*
 *  LISACode-Vector.h
 *  LISACode V 2.0
 *
 *  Created on 08/04/11 by Antoine PETITEAU (AEI)
 *  Last modification on 08/04/11 by Antoine PETITEAU (AEI)
 *
 */

/** \ingroup ToolBox
 * \defgroup Vector Class Vector
 * (See class #LCVector for a detailed description)
 * \{
 */


#ifndef __LCVECTOR_H
#define __LCVECTOR_H

#include <stdexcept>
#include <fstream>
#include <math.h>
#include <complex>
#include <stdlib.h>
#include <vector>
#include "LISACODE-Tools.h"



/** \brief Class for managing vector.
 * \author A. Petiteau
 * \version 2.0
 * \date 08/04/2011
 *
 * This class manage vector
 * 
 *
 */
class LCVector 
{
public:
	
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	/** \brief Components x, y and z */
	double p[3];
	
	
	/********** Constructor **********/
	/** \brief Default constructor */
	LCVector();
	
	/** \brief Standard constructor */
	LCVector(LCTools * MT);
	
	/** \brief Standard constructor with the values of coordinates individually */
	LCVector(LCTools * MT, double x, double y, double z);
	
	/** \brief Standard constructor with the values of coordinates in table */
	LCVector(LCTools * MT_n, double pRef[3]);
	
	/** \brief Destructor */
	~LCVector();
	
	
	/********** Initialization *******/
	void setTools(LCTools * MT_n);
	
	/********** Operators **********/
	
	/*! \brief Operator equality */
	void operator=(const LCVector &a);
	
	/*! \brief Operator add */
	LCVector operator+(const LCVector &a);
	
	/*! \brief Operator substract */
	LCVector operator-(const LCVector &a);
	
	/*! \brief Operator multiply by scalar */
	LCVector operator*(const double &v);
	
	/*! \brief Operator divide by scalar */
	LCVector operator/(const double &v);
	
	/*! \brief Operator quick add */
	void operator+=(const LCVector &a);
	
	/*! \brief Operator quick substract */
	void operator-=(const LCVector &a);

	/*! \brief Operator quick multiply by a scalar */
	void operator*=(const double &v);
	
	/*! \brief Operator quick divied by a scalar */
	void operator/=(const double &v);
	
	/*! \brief Operator scalar product */
	double operator*(const LCVector &a);
	
	/*! \brief Return the value of the k component */
	double operator()(const int k);
	
	/*! \brief Set the value of the k component */
	void operator()(const int k, const double v);
	
	/*! \brief Display in specified output */
	friend std::ostream & operator<<(std::ostream &output, const LCVector &a);
	
	/********** Access methods  **********/
	
	/** \brief Return x component */
	double x() {return(p[0]);};
	
	/** \brief Return y component */
	double y() {return(p[1]);};
	
	/** \brief Return z component */
	double z() {return(p[2]);};
	
	/** \brief Return polar angle, $\theta$, of spherical coordinates with value in [0,pi] */
	double th() {return( acos(z()) );};
	
	/** \brief Return azimuthal angle, $\phi$, of spherical coordinates with value in [0,2pi] */
	double ph();
	
	/** \brief Set x component to argument */
	void x(double v) {p[0]=v;};
	
	/** \brief Set y component to argument */
	void y(double v) {p[1]=v;};
	
	/** \brief Set z component to argument */
	void z(double v) {p[2]=v;};
	
	/** \brief Set x,y,z components to arguments */
	void xyz(double vx, double vy, double vz) {p[0]=vx; p[1]=vy; p[2]=vz;};
	
	/********** Other methods  **********/
	
	/*! \brief Return norm */
	double norm();
	
	/*! \brief Return norm */
	LCVector unit();
	
	/*! \brief Return result of scalar product between vector local and vector a */
	double scalar(const LCVector &a);
	
	/*! \brief Return result of cross product between vector local and vector a */
	LCVector cross(const LCVector &a);
	
	/*! \brief Display in Mathematica form */
	void DispMathematica();
	
	
	
};

#endif // __LCVECTOR_H

/**\}*/

// end of LISACODE-LCVector.h


