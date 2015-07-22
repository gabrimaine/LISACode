// $Id:  Exp $
/*
 *  LISACode-Serie2.h
 *  LISACode V 2.0
 *
 *  Created on 15/02/11 by Antoine PETITEAU (AEI)
 *  Last modification on 08/04/11 by Antoine PETITEAU (AEI)
 *
 */

/** \ingroup ToolBox
 * \defgroup LCSerie2_ Class LCSerie2
 * (See class #LCSerie2 for a detailed description)
 * \{
 */

#ifndef __LCSERIE2_H
#define __LCSERIE2_H

#include <stdexcept>
#include <fstream>
#include <math.h>
#include <complex>
#include <stdlib.h>
#include <vector>
#include "LISACODE-Tools.h"


/**\brief Interpolation type*/
enum INTERP{TRU, /**<Truncated interpolation*/
	LIN,  /**<Linear interpolation*/
	CUB,  /**<Cubic interpolation*/
	LAG,  /**<Lagrange interpolation*/
	SIN,   /**<Truncated sinc interpolation*/
	UND /**<Undefined*/
};


/**
 * \brief Class for managing each data output.
 * \author A. Petiteau
 * \version 2.0
 * \date 15/02/2011
 *
 * This class stored data directly using pointer.
 * It gives the value \f$y^{I}\f$ (\f$I\f$ refers to the signal) to the corresponding \f$x\f$ 
 * (usually time or frequency) using interpolation methods
 * which can be specified.
 *
 */
class LCSerie2 
{
protected:
	
	/**\brief Pointer on toolbox class */
	LCTools * MT;
	
	/**\brief Reference serie data. */
	double ** ys;
	
	/**\brief Number of data. */
	int N;
	
	/**\brief Maximal number of data. */
	int Nmax;
	
	/**\brief Reference (starting) value : \f$ x_0 \f$ */
	double x0;
	
	/**\brief Reference serie step \f$ {\delta x} \f$ */
	double dx;
	
	/**\brief Table containing Lagrange coefficient (filled at the first used). */
	double * LagPolDen[LC::ORDERMAXLAG];
	
	/**\brief Relative reference of the current \f$x\f$ : \f$ ((x-x_0)/{\delta x}) \f$ */
	double xr;
	
	/**\brief Bin reference of the current \f$x\f$ : \f$ Int_{min}((x-x_0)/{\delta x}) \f$ */
	int bin;
	
	/**\brief True if #xr and #bin for the current \f$x\f$ is already computed  
	 * The goal here is to avoid the multiple computation of \f$ Int_{min}((x-x_0)/{\delta x}) \f$ in 
	 * gData and TruncVal, InterLinear, InterCubic, InterHermite, or InterLagrange for computational time saving.
	 */
	bool xrbinE;
	
	/** Sum, minimal and maximal value */
	double Sum, Min, Max;
	
public:
	// **** Constructor ****
	
	/**\brief Default constructor */
	LCSerie2();
	
	/**\brief Standard constructor 
	 * @param[in]  MT_n   Pointer on toolbox
	 * @param[in]  start  Reference (starting point) \f$ x_0 \f$. 
	 * @param[in]  delta  Step \f$ \delta x \f$.
	 */
	LCSerie2(LCTools * MT_n, double start, double delta);
	
	/**\brief Standard constructor with specification of the length of the serie
	 * @param[in]  MT_n    Pointer on toolbox
	 * @param[in]  start   Reference (starting point) \f$ x_0 \f$. 
	 * @param[in]  delta   Step \f$ \delta x \f$.
	 * @param[in]  length  Length of the data (number of data) \f$ N \f$.
	 */
	LCSerie2(LCTools * MT_n, double start, double delta, int length);
	
	/**\brief Constructor including data to be stored :
	 * @param[in]  MT_n    Pointer on toolbox
	 * @param[in]  start   Reference (starting point) \f$ x_0 \f$. 
	 * @param[in]  delta   Step \f$ \delta x \f$.
	 * @param[in]  ys_n    Initial data,
	 * @param[in]  N_n     Length of the data (number of data) \f$ N \f$,
	 * @param[in]  Nmax_n  Maximal length of the data (number of data) \f$ N_{max} \f$.
	 */
	LCSerie2(LCTools * MT_n, double start, double delta, double ** ys_n, int N_n, int Nmax_n);
	
	/**\brief Constructor including file for downloading data
	 * @param[in] MT_n       Pointer on toolbox
	 * @param[in] FileName   Name of file. 
	 * @param[in] iCol       Index of the read column (the column 1 is used as reference)
	 * @param[in] xMinExtra  Minimal value for reference \f$ x \f$ (reached by linear extrapolation)
	 * @param[in] xMaxExtra  Maximal value for reference \f$ x \f$ (reached by linear extrapolation)
	 */
	LCSerie2(LCTools * MT_n, char* FileName, int iCol, double xMinExtra, double xMaxExtra);
	
	/**\brief Destructor 
	 */
	~LCSerie2();
	
	
	/**********  Access methods  **********/
	
	/**\brief Returns number of values (size of #ys attribute). 
	 */
	int getNbVal() const {return(N);};
	
	/**\brief Returns number of values (size of #ys attribute). 
	 */
	int getNmax() const {return(Nmax);};
	
	/**\brief Returns reference step in \f$x\f$ 
	 */
	double getRefStep() const {return(dx);};
	
	/**\brief Set start value of \f$x\f$ for reference 
	 */
	void setRefStart(double start) {x0 = start;};
	
	/**\brief Set reference step in \f$x\f$ 
	 */
	void setRefStep(double delta) {dx = delta;};
	
	/**\brief Set maximal number of value 
	 */
	void setNmax(int Nmax_n);
	
	/**\brief Return reference value corresponding to the bin 
	 */
	double getRef(int bin) const;
	
	/**\brief Return data value corresponding to the bin 
	 */
	double getBinValue(int bin) const;
	
	/**\brief Set data value corresponding to the bin 
	 */
	void setBinValue(int bin, double x);
	
	/**\brief Return sum of value (computed in function #ComputeSumMinMax()) */
	double getSum() {return(Sum);};
	
	/**\brief Return minmal value (computed in function #ComputeSumMinMax())*/
	double getMin() {return(Min);};
	
	/**\brief Return maximal value (computed in function #ComputeSumMinMax())*/
	double getMax() {return(Max);};
	
	/** \brief Display data (for checking) */
	void DispData();
	
	/**********  Initialization methods  **********/
	void initNULL(bool CleanMem);
	
	
	/**********  Load methods  **********/
	
	/** \brief Add data at the begining and remove the last one if the maximal number of data is reached 
	 * @param[in] y Value which is added
	 */ 
	void addData(double y);
	
	/** \brief Allocation of all the memory with default value 0 for new data */ 
	void allocAll();
	
	/** \brief Load data in file and extrapole 
	 * @param[in]  FileName   Name of file. 
	 * @param[in]  iCol       Index of the read column (the column 1 is used as reference)
	 * @param[in]  xMinExtra  Minimal value for reference \f$ x \f$ (reached by linear extrapolation)
	 * @param[in]  xMaxExtra  Maximal value for reference \f$ x \f$ (reached by linear extrapolation)
	 */ 
	void LoadFile(char* FileName, int iCol, double xMinExtra, double xMaxExtra);
	
	
	
	/**********  Get methods  **********/
	
	/** \brief Return the data for the reference x with the specified interpolation method.
	 * @param[in] x                Reference value \f$x\f$ (function found the corresponding data)   
	 * @param[in] InterpType     Type of interpolation
	 * @param[in] InterpUtilValue  Value used during the interpolation : it can be necessary for making the interpolation (For example for lagrange interpolation, InterpUtilValue is the order of the interpolation)
	 *
	 */
	double gData(double x, INTERP InterpType, double InterpUtilValue);
	
	/**\brief Return truncated value 
	 * @param[in] x  Reference value \f$ x \f$ (function found the corresponding data)
	 */
	inline double TruncVal(double x);
	
	/**\brief Return value obtained by linear interpolation
	 * @param[in] x  Reference value \f$x\f$ (function found the corresponding data)
	 */
	inline double InterLinear(double x);
	
	/**\brief Return value obtained by cubic interpolation
	 * @param[in] x  Reference value \f$x\f$ (function found the corresponding data)
	 */
	double InterCubic(double x);
	
	/**\brief Return value obtained by hermite interpolation 
	 * @param[in] x       Reference value \f$x\f$ (function found the corresponding data)
	 * @param[in] tension tension parameter used by the hermte interpolator
	 * @param[in] bias    bias parameter used by the hermte interpolator
	 */
	double InterHermite(double x, double tension, double bias);
	
	/**\brief Initialize the lagrange coefficient 
	 * @param[in] order Order of Lagrange interpolation
	 */
	void InitLagPolDen(int order);
	
	/**\brief Return data value corresponding to reference \f$x\f$ using Lagrange interpolation at specified order 
	 * @param[in] x     Reference value \f$x\f$ (function found the corresponding data)
	 * @param[in] order Order of Lagrange interpolation
	 */
	inline double InterLagrange(double x, int order);
	
	/**\brief Compute x corresponding the value with the specified interpolation method.
	 * @param[in] Val             Required value
	 * @param[in] InterpType    Type of interpolation
	 * @param[in] InterpUtilValue Value used during the interpolation
	 * @param[in] ReqPrecision    Required precision
	 * \return  Reference x corresponding the value with the specified interpolation method
	 */
	double gxInv(double Val, INTERP InterpType, double InterpUtilValue, double ReqPrecision);
	
	
	
	/**********  Operation methods  **********/
	
	/** \brief Return the sum */
	void ComputeSumMinMax();
	
	/** \brief Normalize the data : divide each one by the sum */
	void Normalize();
	
	/** \brief Transform in cumulative */
	void Cumul();
	
	/** \brief Transform : Absolute value  */
	void Abs();
	
	/** \brief Transform : Put 0 instead of negative values  */
	void Positive();

	
	/**********  Derivative methods  **********/
	
	/*! \brief Compute raw derivative using current bin and bin at -NDbin */
	double DerivRawCur(int NDbin);
	
	/*! \brief Compute backward derivative at second order of precision (use the 3 last bins)  */
	double DerivBackOrder2Cur();
	
	/*! \brief Compute raw derivative using x and x - Dx points,   */
	double DerivRawSpe(double x, double Dx, INTERP InterpType, double InterpUtilValue);
	
	/*! \brief Compute backward derivative et second order of precision (use x, x-Dx and x-2Dx points)  */
	double DerivBackOrder2Spe(double x, double Dx, INTERP InterpType, double InterpUtilValue);


};


#endif // __LCSERIE2_H

/**\}*/

// end of LISACODE-Serie2.h