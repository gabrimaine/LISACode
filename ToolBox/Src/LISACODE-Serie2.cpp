// $Id:  $
/*
 *  LISACODE-Serie2.cpp
 *  LISACode V 2.0
 *
 *
 *  Created on 27/05/09 by Antoine Petiteau (AEI)
 *  Last modification on 26/01/10 by Antoine Petiteau (AEI)
 *
 */

#include "LISACODE-Serie2.h"


// *****************
// *  Constructor  *
// *****************

/** 
 * Constructs an instance and initializes it with default values.
 *
 * \arg  #x0 = 0
 * \arg  #dx = 1
 * \arg  Initialize #LagPolDen at NULL
 */
LCSerie2::LCSerie2()
{
	MT = new(LCTools);
	initNULL(false);
}

/** 
 * Constructs an instance and initializes it with inputs and default values.
 *
 * Some details on running :
 * \arg  Initialize data vector #ys and number of data #N
 * \arg  Set reference by #x0 and #dx using arguments
 * \arg  Initialize pointer on denominator of Lagrange polynom #LagPolDen
 * \arg  No exitstig bin for current \f$ x \f$ (of course, because there is not current \f$ x \f$)
 */
LCSerie2::LCSerie2(LCTools * MT_n, double start, double delta)
{
	MT = MT_n;
	initNULL(false);
	
	x0 = start;
	dx = delta;
}


/** 
 * Constructs an instance and initializes it with inputs and default values.
 *
 * Some details on running :
 * \arg  Initialize data vector #ys and number of data #N
 * \arg  Set reference by #x0 and #dx using arguments
 * \arg  Initialize pointer on denominator of Lagrange polynom #LagPolDen
 * \arg  No exitstig bin for current \f$ x \f$ (of course, because there is not current \f$ x \f$)
 *
 */
LCSerie2::LCSerie2(LCTools * MT_n, double start, double delta, int length)
{
	MT = MT_n;
	initNULL(false);
	ys = NULL;
	Nmax = length;
	x0 = start;
	dx = delta;
	ys = (double**)MT->AllocMemory(Nmax*sizeof(double*));
}


/**  
 * Constructs an instance and initializes it with inputs. 
 *
 * Some details on running :
 * \arg  Initialize data vector #ys and number of data #N
 * \arg  Set reference by #x0 and #dx using arguments
 * \arg  Initialize pointer on denominator of Lagrange polynom #LagPolDen
 * \arg  No exitstig bin for current \f$ x \f$ (of course, because there is not current \f$ x \f$)
 */
LCSerie2::LCSerie2(LCTools * MT_n, double start, double delta, double ** ys_n, int N_n, int Nmax_n)
{
	MT = MT_n;
	initNULL(false);
 	x0 = start;
	dx = delta;
	Nmax = Nmax_n;
	N = N_n;
	ys = (double**)MT->AllocMemory(Nmax*sizeof(double*));
	for(int i=0; i<N; i++){
		ys[i] = (double*)MT->AllocMemory(sizeof(double));
		(*ys[i]) = (*ys_n[i]);
	}
}


/**  
 * Constructs an instance and initializes it with inputs :
 *  
 * Some details on running :
 * \arg  Load file
 * \arg  Initialize pointer on denominator of Lagrange polynom #LagPolDen
 * \arg  No exitstig bin for current \f$ x \f$ (of course, because there is not current \f$ x \f$)
 */
LCSerie2::LCSerie2(LCTools * MT_n, char* FileName, int iCol, double xMinExtra, double xMaxExtra)
{
	MT = MT_n;
	initNULL(false);
	LoadFile(FileName, iCol, xMinExtra, xMaxExtra);
	
}


/** 
 * Destructor. 
 */
LCSerie2::~LCSerie2()
{
	initNULL(true);
}


// **********************
// **  Access methods  **
// **********************

/** 
 * @param[in] Nmax_n  New maximal number of data
 *
 * Details on running :
 *
 */
void LCSerie2::setNmax(int Nmax_n)
{
	// \todo Describe the running
	if(Nmax == 0){
		Nmax = Nmax_n;
		ys = (double**)MT->AllocMemory(Nmax*sizeof(double*));
	}else{
		double ** tmpys(ys);
		ys = NULL;
		ys = (double**)MT->AllocMemory(Nmax_n*sizeof(double*));
		for(int i=0; i<Nmax_n; i++)
			ys[i] = tmpys[i];
		for(int i=Nmax; i<N; i++)
			delete tmpys[i];
		MT->Free(tmpys, Nmax*sizeof(double*));
		Nmax = Nmax_n;
	}
}


/** \brief Gets reference value corresponding to bin input.
 *
 * Input is checked ; it must be positive or null, and lower than #ys attribute size. \n
 * \return
 * \f$ x0 + dx \cdot bin \f$
 */
double LCSerie2::getRef(int bin) const
{
	//if((bin<0)||(bin > N-1))
	//	throw std::invalid_argument("LCSerie2: The required bin does not exist !");  
	return(x0 + dx*(bin));
}

/** \brief Gets reference y value corresponding to bin input.
 *
 * Input is checked ; it must be positive or null, and lower than #ys attribute size. \n
 * \return
 * ys[bin]
 */
double LCSerie2::getBinValue(int bin) const
{
	if((bin<0)||(bin > N-1)){
		Cout << "ERROR in LCSerie2::getBinValue : The required bin " << bin << " is not in correct range [0," << N-1 << "]" << Endl; 
		throw std::invalid_argument("LCSerie2::getBinValue : The required bin does not exist !"); 
	}
	return((*ys[bin]));
}

/** \brief Sets reference y value corresponding to bin and x inputs.
 *
 * Input is checked ; it must be positive or null, and lower than #ys attribute size. \n
 * Then ys is filled :
 * ys[bin]=x .
 */
void LCSerie2::setBinValue(int bin, double x)
{
	if((bin<0)||(bin > N-1)){
		Cout << "ERROR in LCSerie2::setBinValue : The required bin " << bin << " is not in correct range [0," << N-1 << "]" << Endl; 
		throw std::invalid_argument("LCSerie2::setBinValue : The required bin does not exist !");
	}
	(*ys[bin]) = x;
}

void LCSerie2::DispData()
{
	if(MT->Disp()){
		for(int k=0; k<Nmax; k++)
			if(ys[k]!=NULL)
				Cout << "\t>\t" << k << " " << (*ys[k]) << Endl;
	}
}


// ******************************
// **  Initialization methods  **
// ******************************

void LCSerie2::initNULL(bool CleanMem)
{
	if(CleanMem){
		
		for(int i=0; i<N; i++)
			MT->Free(ys[i],sizeof(double));
		MT->Free(ys, Nmax*sizeof(double*));
		
		for(int order=0; order<LC::ORDERMAXLAG; order++)
			if(LagPolDen[order] != NULL)
				MT->Free(LagPolDen[order], (order+1)*sizeof(double) );
		
	}
	
	ys = NULL;
	N = 0;
	x0 = 0.0;
	dx = 1.0;
	for(int i=0; i<LC::ORDERMAXLAG; i++)
		LagPolDen[i] = NULL;
	xrbinE = false;
	Sum = 0.0;
	Min = 1.0e100;
	Max = -1.0e100;
	
}



// ********************
// **  Load methods  **
// ********************

/**  
 * Adds data at the begining of the serie.
 */
void LCSerie2::addData(double y)
{
	double * ysswap;
	
	if(N < Nmax){
		N++;
		for(int k=N-1; k>0; k--)
			ys[k] = ys[k-1];
		ys[0] = (double*) MT->AllocMemory(sizeof(double));
	}else{
		ysswap = ys[N-1];
		for(int k=N-1; k>0; k--)
			ys[k] = ys[k-1];
		ys[0] = ysswap;
	}
	(*ys[0]) = y; 
} 


/**  
 * 
 */
void LCSerie2::allocAll()
{
	for(int i=N; i<Nmax; i++){
		ys[i] = (double*) MT->AllocMemory(sizeof(double));
		(*ys[i]) = 0.0;
	}
	N = Nmax;
}


/**  
 *
 * Details on running : none
 *
 */
void LCSerie2::LoadFile(char* FileName, int iCol, double xMinExtra, double xMaxExtra)
{
	std::vector<double> x(0), Dat(0);
	std::ifstream fin(FileName);
	std::string junk;
	double tmpx, tmpDat;
	int xbin;
	double xMax, xVal;
	
	//Cout << Endl << "In LoadxDataFromFile ... " << FileName << Endl;
	
	
	/** \arg Count the number of columns, */ 
	if (!fin)
		throw std::invalid_argument("ERROR: LCSerie2::LoadFile : Problem in openning data file.");
	
	/** \arg Read the file and store it in temporary memory, */ 
	while(fin.peek() == '#'){
		fin.ignore(16384,'\n');
	};
	while(!fin.eof()){
		fin >>  tmpx; 
		for(int iC=1; iC<iCol; iC++)
			fin >> junk;
		fin >> tmpDat;
		fin.ignore(16384,'\n');
		if((x.size()<1)||(!MT->deq( tmpx, x[x.size()-1] ))){
			x.push_back(tmpx);
			Dat.push_back(tmpDat);
		}
		//Cout << x.size()-1 << " " << tmpx << " " << tmpDat << Endl; 
	};
	fin.close();
	
	
	/** \arg Minimal step found is used for setting \f$ \delta_x \f$,  */
	dx = 1.e30;
	//if(xMinExtra<x[0])
	//	dx = x[0]-xMinExtra;
	for(int i=0; i<x.size()-1; i++)
		dx = MIN(dx, x[i+1]-x[i]);
	
	/** \arg Define boundaries using data and the parameters xMinExtra, xMaxExtra */
	x0 = x[0];
	while (x0 > xMinExtra)
		x0 -= dx;
	xMax = x[x.size()-1];
	while (xMax < xMaxExtra)
		xMax += dx;
	
	
	N = MT->iceil((xMax - x0) / ((double)(dx))) + 1;
	Nmax = N;
	
	/** \arg Allocation of memory */
	ys = (double**)MT->AllocMemory(Nmax*sizeof(double*));
	for(int i=0; i<Nmax; i++){
		ys[i] = (double*) MT->AllocMemory(sizeof(double));
		(*ys[i]) = 0.0;
	}
	
	//Cout << Endl << Endl;
	
	/** \arg Fill data base using linear interpolation */
	xbin = 0;
	for(int i=0; i<N; i++){
		// Cout << *ys[i] << Endl;
		xVal = x0 + i * dx;
		while((xVal<x[xbin])&&(xbin>0)){
			xbin--;
		};
		while((xVal>=x[xbin+1])&&(xbin+1<x.size()-1)){
			xbin++;
		};
		
		*ys[i] = ( (x[xbin+1]-xVal) * (Dat[xbin]) + (xVal-x[xbin]) * (Dat[xbin+1]) ) / (x[xbin+1]-x[xbin]);
		
		//Cout << i << " " << xVal << " " << *ys[i] << Endl;
	}
	
	//Cout << Endl << Endl;
	ComputeSumMinMax();
}


// *******************
// **  Get methods  **
// *******************


/**
 * Returns the exact value if x is a multiple of dx, else interpolated value.
 *
 * Interpolation method depends on InterpType input. \n
 * It is checked. Its expected values are TRU, LIN, CUB and LAG.
 *
 * \return
 * \f[ \left\{ \begin{array}{ll} 
 TruncVal(x) & \textrm{if } modulo(\frac{x-x_0}{dx},1) \approx 0 \\
 TruncVal(x) & \textrm{if } InterpType = TRU \\
 InterLinear(x) & \textrm{if } InterpType = LIN \\
 InterCubic(x) & \textrm{if } InterpType = CUB \\
 InterLagrange(x, int(InterpUtilValue)) & \textrm{if } InterpType =  LAG
 \end{array} \right. \f]
 *
 */
double LCSerie2::gData(double x, INTERP InterpType, double InterpUtilValue)
{
	//Cout << "LCSerie2::gData :" << Endl;
	//Cout << "  fmod(xr,1.0) = " << fmod(xr,1.0) << Endl;
	
	xr = (x-x0)/dx ;
	bin = MT->ifloor(xr);
	xrbinE = true;
	//if(fmod((x-x0)/dx,1.0) < PRECISION){
	//Cout << "  --> Tronque" << Endl;
	if((xr-bin) < LC::PRECISION){
		return(TruncVal(x));
	}else{
		switch(InterpType){
			case TRU :
				return(TruncVal(x));
				break;
			case LIN :
				return(InterLinear(x));
				break;
			case CUB :
				return(InterCubic(x));
				break;
			case LAG :
				return(InterLagrange(x, int(InterpUtilValue)));
				break;
				//case SIN :
				//	return(InterSincTrunc(x, int(InterpUtilValue)));
				//	break;	
			default :
				throw std::invalid_argument("LCSerie2: No interpolation method !");
		}
	}
}



/** \brief Returns truncated value.
 *
 * First, bin index is computed :
 * \f[ bin= floor(\frac{x-x_0}{dx}) \f]
 * Then bin is checked ; it must be positive or null, and lower than #ys attribute size. \n
 *
 * \return
 * \f[ ys[bin] \f]
 */
inline double LCSerie2::TruncVal(double x)
{
	if(!xrbinE)
		bin = MT->ifloor((x-x0)/dx);
	if((bin<0)||(bin > N-1)){
		std::cerr << "ERROR in LCSerie2::TruncVal : bin = " << bin << " (x = " << x << ") is not included in [0," << N-1 << "]." << Endl;  
		throw std::invalid_argument("LCSerie2: TruncVal(): The required bin does not exist !");
	}
	//Cout << "  TruncVal : x0 = " << x0 << " , x = " << x << " , bin = " << bin;
	//Cout << "  , res = " << *ys[bin] << Endl;
	xrbinE = false;
	return((*ys[bin]));
}


// Return value obtained by linear interpolation
/** \brief Returns linear interpolation result.
 *
 * First, bin index is computed :
 * \f[ bin= floor(\frac{x-x_0}{dx}) \f]
 * Then bin is checked ; it must be positive or null, and lower than #ys attribute size. \n
 *
 * \return
 * \f[ (bin+1-\frac{x-x_0}{dx}) \cdot ys[bin] + \frac{x-x_0}{dx-bin} \cdot ys[bin+1] \f]
 */
inline double LCSerie2::InterLinear(double x)
{
	if(!xrbinE){
		xr = (x-x0)/dx ;
		bin = MT->ifloor(xr);
	}
	if((bin<0)||(bin+1 > N-1)){
		std::cerr << "ERROR in LCSerie2::InterLinear : bin = " << bin << " (x = " << x << ") is not included in [0," << N-1 << "]." << Endl;  
		throw std::invalid_argument("LCSerie2: InterLinear(): The required bin does not exist !");
	}
	xrbinE = false;
	return( (bin + 1.0 - xr) * (*ys[bin]) + (xr-bin) * (*ys[bin+1]) ); //Linear interpolation	
}


// Return value obtained by cubic interpolation
/** \brief Returns cubic interpolation result.
 *
 * Indices are computed : \n
 * \f$ bin_1= floor(\frac{x-x_0}{dx}) \f$ \n
 * \f$ bin_2= bin_1+1 \f$ \n
 * \f$ bin_0= bin_1-1 \f$ \n
 * \f$ bin_3= bin_2+1 \f$ \n
 * Indices are checked ; \f$bin_0 \f$ must be positive or null, and \f$bin_3 \f$ must be lower than #ys attribute size.\n
 *
 * \return
 * \f$ {\mu}^3 \cdot (ys[bin_3]-ys[bin_2]-ys[bin_0]+ys[bin_1]) \f$ \n
 * \f$ + {\mu}^2 \cdot (2 \cdot ys[bin_0]-2 \cdot ys[bin_1]+ys[bin_2]-ys[bin_3]) \f$ \n
 * \f$ + {\mu} \cdot (ys[bin_2]-ys[bin_0]) \f$ \n
 * \f$ + ys[bin_3] \f$, \n
 * where :\n
 * \f$ \mu= \frac{x-x_0}{dx}-floor(\frac{x-x_0}{dx}) \f$
 */
double LCSerie2::InterCubic(double x)
{
	if(!xrbinE){
		xr = (x-x0)/dx ;
		bin = MT->ifloor(xr);
	}
	int bin0, bin1, bin2, bin3;
	double y0,y1,y2,y3,a0,a1,a2,a3;
	double mu;
	
	bin1=bin;
	bin2=bin1+1;  
	bin0=bin1-1;
	bin3=bin2+1;
	if((bin0<0)||(bin3 > N-1))
		throw std::invalid_argument("LCSerie2: InterCubic(): The required bin does not exist !");  
	
	y0=(*ys[bin0]);
	y1=(*ys[bin1]);
	y2=(*ys[bin2]);
	y3=(*ys[bin3]);
	
	a0=y3-y2-y0+y1;
	a1=y0-y1-a0;
	a2=y2-y0;
	a3=y1;
	
	mu=(xr)-floor(xr);
	
	xrbinE = false;
	return ((mu*mu*mu)*a0+(mu*mu)*a1+mu*a2+a3);
}

//Return value obtained by hermite interpoilation
/** \brief Returns hermite interpolation result, depending on x, tension and bias inputs.
 *
 * Indices are computed : \n
 * \f$ bin_1= floor(\frac{x-x_0}{dx}) \f$ \n
 * \f$ bin_2= bin_1+1 \f$ \n
 * \f$ bin_0= bin_1-1 \f$ \n
 * \f$ bin_3= bin_2+1 \f$ \n
 * Indices are checked ; \f$bin_0 \f$ must be positive or null, and \f$bin_3 \f$ must be lower than #ys attribute size.\n
 *
 * \return
 * \f$ (\mu^3 - 2 \cdot \mu^2 + \mu) \cdot ys[bin_1] \f$
 * \f$	+( \mu^3 - 2 \cdot \mu^2 + \mu)) \cdot  \f$
 * \f$((ys[bin_1]-ys[bin_0]) \cdot (1+bias) \cdot \f$ \n
 * \f$ \frac{1-tension}{2}+(ys[bin_2]-ys[bin_1]) \cdot (1-bias) \cdot \frac{1-tension}{2}) \f$, \n
 * where :\n
 * \f$  \mu= \frac{x-x_0}{dx}-floor(\frac{x-x_0}{dx}) \f$
 *
 */
double LCSerie2::InterHermite(double x, double tension, double bias)
{
	if(!xrbinE){
		xr = (x-x0)/dx ;
		bin = MT->ifloor(xr);
	}
	int bin0, bin1, bin2, bin3;
	double y0,y1,y2,y3,a0,a1,a2,a3;
	double mu, mu2, mu3;
	double m0, m1;
	
	bin1=bin+1;
	bin2=bin1+1;
	bin0=bin1-1;
	bin3=bin2+1;
	if((bin0<0)||(bin3 > N-1))
		throw std::invalid_argument("LCSerie2: InterHermite(): The required bin does not exist !");  
	
	y0=(*ys[bin0]);
	y1=(*ys[bin1]);
	y2=(*ys[bin2]);
	y3=(*ys[bin3]);
	
	
	mu=(xr)-floor(xr);
	mu2 = mu * mu;
	mu3 = mu2 * mu;
	m0  = (y1-y0)*(1+bias)*(1-tension)/2;
	m0 += (y2-y1)*(1-bias)*(1-tension)/2;
	m1  = (y2-y1)*(1+bias)*(1-tension)/2;
	m1 += (y3-y2)*(1-bias)*(1-tension)/2;
	a0 = 2*mu3 - 3*mu2 + 1;
	a1 = mu3 - 2*mu2 + mu;
	a2 =  mu3 - mu2;
	a3 = -2*mu3 + 3*mu2;
	
	xrbinE = false;
	return(a0*y1+a1*m0+a2*m1+a3*y2); 
}


void LCSerie2::InitLagPolDen(int order)
{
	LagPolDen[order] = (double*) MT->AllocMemory((order+1)*sizeof(double));
	for(int k=0; k<=order; k++){
		LagPolDen[order][k] = 1.0;
		for(int j=0; j<=order; j++){
			if(j!=k)
				LagPolDen[order][k] *= 1.0/(k-j);
		}
	}
}


/** \brief Returns Lagrange interpolation result.
 *
 * Indices are computed : \n
 * \f$ bin= floor(\frac{x-x_0}{dx}) \f$ \n
 * \f$ kmin= bin-ordermin+1  \textrm{ , where } ordermin=floor(\frac{order+1}{2}) \f$ \n
 * \f$ kmax= bin+order+1-ordermin \f$ \n
 * Indices are checked : 
 * \arg bin must be positive or null, and lower than (#ys attribute size -1)
 * \arg kmin must be positive or null\n
 * \arg kmax must be and lower than (#ys attribute size -1)
 *
 * \return
 * \f$ \sum_{k=kmin}^{kmax} ys[k] \cdot P_k
 \textrm{ , where }
 P_k=\prod_{j=kmin, j \ne k}^{kmax} \frac{x-x_0-j \cdot dx}{(k-j) \cdot dx} \f$
 */
inline double LCSerie2::InterLagrange(double x, int order)
{
	//Cout << "////\\\\LCSerie2::InterLagrange : ordre = " << order << "  et x = " << x << Endl;
	if(!xrbinE){
		xr = (x-x0)/dx ;
		bin = MT->ifloor(xr);
	}
	double res(0.0), Pk(0.0);
	int ordermin(MT->ifloor(double(order+1)/2.0));
	//Cout << int(ys.size())-1 << Endl;
	if((bin<0)||(bin+1 > N-1)){
		std::cerr << "ERROR in LCSerie2::InterLagrange : bin = " << bin << " (x = " << x << ") is not included in [0," << N-1 << "]." << Endl;  
		throw std::invalid_argument("LCSerie2::InterLagrange : The required bin does not exist !");  
	}
	int kmin(bin-ordermin+1), kmax(bin+(order+1-ordermin));
	if(kmin < 0){
		std::cerr << "ERROR in LCSerie2::InterLagrange : For bin = " << bin << " (x = " << x << "), kmin = bin-ordermin+1 < 0 " << Endl;
		throw std::invalid_argument("LCSerie2::InterLagrange : The required bin does not exist !");
	}
	if(kmax > N-1){
		std::cerr << "ERROR in LCSerie2::InterLagrange : For bin = " << bin << " (x = " << x << "), kmax = order+1-ordermin > " << N-1 << Endl;
		throw std::invalid_argument("LCSerie2::InterLagrange : The required bin does not exist !");
	}
	//Cout << "  InterLagrange : x0 = " << x0 << " , x = " << x;
	//Cout << " , bin = " << bin << " , kmin = " << kmin << " , kmax = " << kmax;
	
	/*
	 // Old Lagrange
	 for(int k=kmin; k<=kmax; k++){
	 Pk = 1.0;
	 for(int j=kmin; j<=kmax; j++){
	 if(j!=k)
	 Pk *= (x-x0-j*dx)/((k-j)*dx);
	 }
	 res += (*ys[k])*Pk;
	 }
	 */
	
	// ** New Lagrange with tabulation of polynome denominator
	int RealOrder (kmax - kmin);
	if(LagPolDen[RealOrder] == NULL)
		InitLagPolDen(RealOrder);
	
	for(int k=kmin; k<=kmax; k++){
		Pk = LagPolDen[RealOrder][k-kmin];
		for(int j=kmin; j<=kmax; j++){
			if(j!=k)
				Pk *= xr-j ;
		}
		res += (*ys[k])*Pk;
	}
	
	
	//Cout << "  , res = " << res << Endl;
	xrbinE = false;
	return(res);
}


/*double LCSerie2::InterSincTrunc(double x, int Nlenght) const
 {
 //Cout << "////\\\\LCSerie2::InterLagrange : ordre = " << order << "  et x = " << x << Endl;
 double res(0.0);
 int Nhalf((int)(floor((double)(Nlenght-1)/2.0)));
 int bin((int)(round(xr)));
 double D(fmod(xr,1.0));
 if(D>=0.5)
 D -= 1.0;
 //Cout << "LCSerie2::InterSincTrunc : N = " << Nlenght << Endl;
 //Cout << "  x   = " << x << Endl;
 //Cout << "  D   = " << D << Endl;
 if(bin-Nhalf < 0)
 throw std::invalid_argument("LCSerie2::InterSincTrunc : The required bin does not exist !");
 if(bin+Nhalf > int(ys.size())-1)
 throw std::invalid_argument("LCSerie2::InterSincTrunc : The required bin does not exist !");
 //Cout << "  InterLagrange : x0 = " << x0 << " , x = " << x;
 //Cout << " , bin = " << bin << " , kmin = " << kmin << " , kmax = " << kmax;
 
 for(int k=-1*Nhalf; k<=Nhalf; k++){
 res += ys[bin+k]*sin(D-k)/(D-k);
 }
 //Cout << "res = " << res << Endl;
 //Cout << Endl;
 return(res);
 }
 */


/**
 * Use bisection method. 
 */
double LCSerie2::gxInv(double Val, INTERP InterpType, double InterpUtilValue, double ReqPrecision)
{
	double xlow, xhigh, xmean;
	double ylow, ymean;
	
	if(Val<Min){
		std::cerr << "Error : Reuqired value " << Val << " < minimal value " << Min << Endl;
		throw std::invalid_argument("LCSerie2::gxInv: Value too low.");
	}
	if(Val>Max){
		std::cerr << "Error : Reuqired value " << Val << " < minimal value " << Max << Endl;
		throw std::invalid_argument("LCSerie2::gxInv: Value too high.");
	}	
	
	xlow = x0;
	xhigh = getRef(N-1);
	ylow = gData(xlow, InterpType, InterpUtilValue) - Val;
	do{
		xmean = (xhigh+xlow)/2.0;
		ymean = gData(xmean, InterpType, InterpUtilValue) - Val;
		if(ylow*ymean>0){
			xlow = xmean;
			ylow = ymean;
		}else{
			xhigh = xmean;
		}
			
	}while(fabs(ymean)>ReqPrecision);

	return (xmean);
}


// *************************
// **  Operation methods  **
// *************************


/** 
 *
 */
void LCSerie2::ComputeSumMinMax()
{
	Sum = 0.;
	Min = 1.e100;
	Max = -1.0e100;
	for(int i=0; i<N; i++){
		Sum += *ys[i];
		if(*ys[i]<Min)
			Min = *ys[i];
		if(*ys[i]>Max)
			Max = *ys[i];
	}
}



/** 
 *
 */
void LCSerie2::Normalize()
{
	for(int i=0; i<N; i++){
		*ys[i] = *ys[i]/Sum;
	}
	ComputeSumMinMax();
}


/** 
 *
 */
void LCSerie2::Cumul()
{
	double TmpSum(0.0);
	for(int i=0; i<N; i++){
		TmpSum += *ys[i];
		*ys[i] = TmpSum;
	}
	ComputeSumMinMax();
}


/** 
 *
 */
void LCSerie2::Abs()
{
	for(int i=0; i<N; i++){
		*ys[i] = fabs(*ys[i]);
	}
	ComputeSumMinMax();
}



/** 
 *
 */
void LCSerie2::Positive()
{
	for(int i=0; i<N; i++){
		if(*ys[i] < 0.0)
			*ys[i] = 0.0;
	}
	ComputeSumMinMax();
}



// **************************
// **  Derivative methods  **
// **************************


double LCSerie2::DerivRawCur(int NDbin)
{
	return ( (getBinValue(0) - getBinValue(NDbin) ) / (NDbin*dx) );
}


double LCSerie2::DerivBackOrder2Cur()
{
	return ( ( 3.*getBinValue(0) - 4.*getBinValue(1) + getBinValue(2) ) / (2.*dx) );
}


double LCSerie2::DerivRawSpe(double x, double Dx, INTERP InterpType, double InterpUtilValue)
{
	return ( (gData(x, InterpType, InterpUtilValue) - gData(x+Dx, InterpType, InterpUtilValue) ) / Dx);
}


double LCSerie2::DerivBackOrder2Spe(double x, double Dx, INTERP InterpType, double InterpUtilValue)
{
	return ( ( 3.*gData(x, InterpType, InterpUtilValue) - 4.*gData(x+Dx, InterpType, InterpUtilValue) + gData(x+2.*Dx, InterpType, InterpUtilValue) ) / (2.*Dx) );
}




//end of LCSerie2.cpp
