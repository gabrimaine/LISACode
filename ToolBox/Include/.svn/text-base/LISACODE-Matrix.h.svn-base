// $Id:  Exp $
/*
 *  LISACode-Matrix.h
 *  LISACode V 2.0
 *
 *  Created on 08/04/11 by Antoine PETITEAU (AEI)
 *  Last modification on 08/04/11 by Antoine PETITEAU (AEI)
 *
 */

/** \ingroup ToolBox
 * \defgroup Matrix Class Matrix
 * (See class #LCMatrix for a detailed description)
 * \{
 */

#ifndef __LCMATRIX_H
#define __LCMATRIX_H

#include <stdexcept>
#include <fstream>
#include <math.h>
#include <complex>
#include <stdlib.h>
#include <vector>
#include "LISACODE-Tools.h"
#include "LISACODE-Vector.h"

//#ifdef _GSLLINKED_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
//#endif



/** \brief Class for managing matrix.
 * \author A. Petiteau
 * \version 2.0
 * \date 13/05/2011
 *
 * This class manage matrix
 * 
 *
 */
class LCMatrix
{
protected:
	
	LCTools * MT;
	
	/*! \brief Matrix : [irow][icolumn]*/
	double ** m; 
	
	/*! \brief Number of columns */
	int nc;
	
	/*! \brief Number of rows */
	int nr;
	
public:
	
	/********** Constructor **********/
	
	/*! \brief Default */
	LCMatrix();
	
	/** \brief Standard constructor for an empty matrix */
	LCMatrix(LCTools * MT_n, int nc_n, int nr_n);
	
	/** \brief Standard constructor for a matrix with all value at v */
	LCMatrix(LCTools * MT_n, int nc_n, int nr_n, double v);
	
	/** \brief Constructor copying an array */
	LCMatrix(LCTools * MT_n, int nc_n, int nr_n, double** n);
	
	/** \brief Constructor copying */
	LCMatrix(const LCMatrix &a);
	
	/** \brief Destructor */
	~LCMatrix();
	
	
	/********** Initalization methods **********/
	
	/*! \brief  */
	void init(LCTools * MT_n, int nc_n, int nr_n);
	
	
	/********** Access methods **********/
	
	/*! \brief  */
	double g(int ir, int il);
	
	/*! \brief  */
	void s(int ir, int il, double v);
	
	/*! \brief  */
	void add(int ir, int il, double v);
	
	int getnrows(){return(nr);};
	
	/*! \brief  */
	int getncolumns(){return(nc);};
	
	/*! \brief  */
	int getMinDim();
	
	/*! \brief  */
	void setNull();
	
	/*! \brief  */
	void setIdentity();
	
	
	/********** Operators **********/
	
	/*! \brief Operator equality */
	void operator=(const LCMatrix &a);
	
	/*! \brief Operator add */
	LCMatrix operator+(const LCMatrix &a);
	
	/*! \brief Operator substract */
	LCMatrix operator-(const LCMatrix &a);
	
	/*! \brief Operator for matrix multiplication */
	LCMatrix operator*(const LCMatrix &a);
	
	/*! \brief Operator multiply by scalar */
	LCMatrix operator*(const double &v);
	
	/*! \brief Operator multiply by vector */
	LCVector operator*(const LCVector &v);
	
	/*! \brief Operator quick add */
	void operator+=(const LCMatrix &a);
	
	/*! \brief Operator quick substract */
	void operator-=(const LCMatrix &a);
	
	/*! \brief Operator quick multiply by another matrix */
	void operator*=(const LCMatrix &a);
	
	/*! \brief Operator quick multiply by a scalar */
	void operator*=(const double &v);
	
	/*! \brief Return the value of the component of row ir and column ic */
	double operator()(const int &ir, const int &ic);
	
	/*! \brief Set the value of the component of row ir and column ic at v */
	void operator()(const int &ir, const int &ic, const double &v);
	
	/*! \brief Display in specified output */
	friend std::ostream &operator<<(std::ostream &output, const LCMatrix &a);
	
	
	/********** Other methods standard **********/
	
	/*! \brief Return the transpose */
	LCMatrix T();
	
	/*! \brief Return the inverse 
	 *	@param KeyMethod Key code defining the method used : 0:Gauss, 1:SVD(GSL link required), 2:LU(GSL link required)
	 */
	LCMatrix Inv(int KeyMethod);
	
	/*! \brief Return the inverse using Gauss-Jourdan pivot */
	LCMatrix InvGJ();
	
	/*! \brief Extract part of matrix : from column ic_e to ic_e + nc_e and from row ir_e to ir_e + nr_e*/
	LCMatrix Extract(int ic_e, int ir_e, int nc_e, int nr_e);
	
	/*! \brief Projection on N components */
	LCMatrix Proj(int * iCom, int NCom);
	
	/*! \brief Compute difference with unit matrix */
	double DiffWithUnit();
	
	/*! \brief Return eigen vector and eigen value  */
	void EigenVsym(LCMatrix & EigenVect, LCMatrix & EigenVal);
	
	/*! \brief Display in Mathematica form */
	void DispMathematica();
	
	
	/********** Other methods requiring GSL **********/
	
//#ifdef _GSLLINKED_
	
	/*! \brief Return the inverse using SVD decomposition */
	LCMatrix InvSVD();
	
	/*! \brief Singular value decomposition */
	void SVD(LCMatrix &U, LCMatrix &S, LCMatrix &V);
	
	/*! \brief Return the inverse using LU refine methods (...) */
	LCMatrix InvLURefine();
	
	/*! \brief Compute eigen vector and eigen value  */
	void EigenVsymCompute(LCMatrix & EigenVect, LCMatrix & EigenVal);
	
//#endif
	
};






#endif // __LCMATRIX_H

/**\}*/

// end of LISACODE-LCMatrix.h



