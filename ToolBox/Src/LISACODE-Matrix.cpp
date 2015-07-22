// $Id:  Exp $
/*
 *  LISACode-Matrix.cpp
 *  LISACode V 2.0
 *
 *  Created on 13/05/11 by Antoine PETITEAU (AEI)
 *  Last modification on 13/05/11 by Antoine PETITEAU (AEI)
 *
 */

#include "LISACODE-Matrix.h"

// *****************
// *  Constructor  *
// *****************


LCMatrix::LCMatrix()
{
	MT = new LCTools;
	m = NULL;
	nc = 0;
	nr = 0;
}


LCMatrix::LCMatrix(LCTools * MT_n, int nr_n, int nc_n)
{
	m = NULL;
	init(MT_n, nr_n, nc_n);
}


LCMatrix::LCMatrix(LCTools * MT_n, int nr_n, int nc_n, double v)
{
	m = NULL;
	init(MT_n, nr_n, nc_n);
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nc; ic++)
			m[ir][ic] = v;
}


LCMatrix::LCMatrix(LCTools * MT_n, int nr_n, int nc_n, double ** n)
{
	m = NULL;
	init(MT_n, nr_n, nc_n);
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nc; ic++)
			m[ir][ic] = n[ir][ic];
}


LCMatrix::LCMatrix(const LCMatrix &a)
{
	m = NULL;
	init(a.MT, a.nr, a.nc);
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nc; ic++)
			m[ir][ic] = a.m[ir][ic];
	
	
}


LCMatrix::~LCMatrix()
{
	if(m != NULL){
		for(int i=0; i<nr; i++)
			MT->Free(m[i], nc*sizeof(double));
		MT->Free(m, nr*sizeof(double*));
	}
}


// ********************
// *  Initialization  *
// ********************

void LCMatrix::init(LCTools * MT_n, int nr_n, int nc_n)
{
	if(m != NULL){
		for(int ir=0; ir<nr; ir++)
			MT->Free(m[ir], nc*sizeof(double));
		MT->Free(m, nr*sizeof(double*));
	}
	MT = MT_n;
	nr = nr_n;
	nc = nc_n;
	m = (double**) MT->AllocMemory(nr*sizeof(double*));
	for(int ir=0; ir<nr; ir++){
		m[ir] = (double*) MT->AllocMemory(nc*sizeof(double));
		for(int ic=0; ic<nc; ic++)
			m[ir][ic] = 0.0;
	}
}


// ********************
// *  Access methods  *
// ********************


double LCMatrix::g(int ir, int ic)
{
	return m[ir][ic];
}


void LCMatrix::s(int ir, int ic, double v)
{
	m[ir][ic] = v;
}

void LCMatrix::add(int ir, int ic, double v)
{
	m[ir][ic] += v;
}


int LCMatrix::getMinDim()
{
	if(nc<nr)
		return nc;
	else
		return nr;
}


void LCMatrix::setNull()
{
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nc; ic++)
			m[ir][ic] = 0.0;
}


void LCMatrix::setIdentity()
{
	int n(getMinDim());
	setNull();
	for(int i=0; i<n; i++)
		m[i][i] = 1.0;
}


// ***************
// *  Operators  *
// ***************


void LCMatrix::operator=(const LCMatrix &a)
{
	if((nr!=a.nr)||(nc!=a.nc))
		throw std::invalid_argument("LCMatrix::operator= : The two LCMatrix must have the same size.");
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nc; ic++)
			m[ir][ic] = a.m[ir][ic];
}


LCMatrix LCMatrix::operator+(const LCMatrix &a)
{
	if((nr!=a.nr)||(nc!=a.nc))
		throw std::invalid_argument("LCMatrix::operator+ : The two LCMatrix must have the same size.");
	LCMatrix r(MT, nr, nc);
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nc; ic++)
			r.m[ir][ic] = m[ir][ic] + a.m[ir][ic];
	return r;
}


LCMatrix LCMatrix::operator-(const LCMatrix &a)
{
	if((nr!=a.nr)||(nc!=a.nc))
		throw std::invalid_argument("LCMatrix::operator- : The two LCMatrix must have the same size.");
	LCMatrix r(MT, nr, nc);
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nc; ic++)
			r.m[ir][ic] = m[ir][ic] - a.m[ir][ic];
	return r;
}


LCMatrix LCMatrix::operator*(const LCMatrix &a)
{
	if((nr!=a.nc)||(nc!=a.nr))
		throw std::invalid_argument("LCMatrix::operator* : The two LCMatrix must have the a.columns = b.lrows and a.rows = b.columns .");
	
	LCMatrix r(MT, nr, nr);
	for(int ir=0; ir<nr; ir++){
		for(int ic=0; ic<nr; ic++){
			r.m[ir][ic] = 0.0;
			for(int i=0; i<nc; i++)
				r.m[ir][ic] += m[ir][i] * a.m[i][ic];
		}
	}
	return r;
}

LCMatrix LCMatrix::operator*(const double &v)
{
	LCMatrix r(MT, nr, nc);
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nc; ic++)
			r.m[ir][ic] = m[ir][ic] * v;
	return r;
}

LCVector LCMatrix::operator*(const LCVector &v)
{
	if((nc!=3)||(nr!=3))
		throw std::invalid_argument("LCMatrix::operator* : The LCMatrix must have the number of columns and rows = 3 (because the LCVector have 3 elements) .");
	LCVector r(MT);
	for(int ir=0; ir<3; ir++) {
		r.p[ir] = 0.;
		for (int i=0; i<3; i++) {
			//Cout << " + " << m[ir][i] * v.p[i] ;
			r.p[ir] += m[ir][i] * v.p[i];
		}
		//Cout << Endl;
	}
	return(r);
}


void LCMatrix::operator+=(const LCMatrix &a)
{
	if((nr!=a.nr)||(nc!=a.nc))
		throw std::invalid_argument("LCMatrix::operator+= : The two LCMatrix must have the same size.");
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nc; ic++)
			m[ir][ic] += a.m[ir][ic];
}


void LCMatrix::operator-=(const LCMatrix &a)
{
	if((nr!=a.nr)||(nc!=a.nc))
		throw std::invalid_argument("LCMatrix::operator-= : The two LCMatrix must have the same size.");
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nc; ic++)
			m[ir][ic] -= a.m[ir][ic];
}


void LCMatrix::operator*=(const LCMatrix &a)
{
	throw std::invalid_argument("LCMatrix::operator*= : This function is not yet implemented.");
	if((nr!=a.nr)||(nc!=a.nc))
		throw std::invalid_argument("LCMatrix::operator*= : The two LCMatrix must have the same size.");
	
	// \todo Product of matrix
	
}


void LCMatrix::operator*=(const double &v)
{
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nc; ic++)
			m[ir][ic] *= v;
}


double LCMatrix::operator()(const int &ir, const int &ic)
{
	return m[ir][ic];
}


void LCMatrix::operator()(const int &ir, const int &ic, const double &v)
{
	m[ir][ic] = v;
}


std::ostream & operator<<(std::ostream &output, const LCMatrix &a)
{
	output << Endl; 
	for(int ir=0; ir<a.nr; ir++){
		for(int ic=0; ic<a.nc; ic++){
			output.width(12);
			output.precision(6);
			output << a.m[ir][ic] << " ";
		}
		output << Endl;
	}
	output << Endl;
	return output;
}



// **************************
// * Other methods standard *
// **************************


LCMatrix LCMatrix::T()
{
	LCMatrix r(MT, nc, nr);
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nc; ic++)
			r.m[ic][ir] = m[ir][ic];
	return r;
}


LCMatrix LCMatrix::Inv(int KeyMethod)
{
	if(KeyMethod==0)
		return(InvGJ());

	if(KeyMethod==1){
//#ifdef _GSLLINKED_
		return(InvSVD());
//#else
//		throw std::invalid_argument("ERROR in LCMatrix::Inv : Sorry but we cannot use SVD decomposition because GSL is not linked !");
//#endif
	}
	
	if(KeyMethod==2){
//#ifdef _GSLLINKED_
		return(InvLURefine());
//#else
//		throw std::invalid_argument("ERROR in LCMatrix::Inv : Sorry but we cannot use LU refine because GSL is not linked !");
//#endif
	}
	
}



LCMatrix LCMatrix::InvGJ()
{
	if(nr!=nc)
		throw std::invalid_argument("LCMatrix::InvGJ : For inversion using Gauss-Jordan pivot, it must be a square matrix.");
	
	LCMatrix r(MT, nr, nr);
	
	double ** SolNullInvMat;
	SolNullInvMat = (double**) MT->AllocMemory(nr*sizeof(double*));
	for(int i=0; i<nr; i++){
		SolNullInvMat[i] = (double*) MT->AllocMemory(sizeof(double));
		SolNullInvMat[i][0] = 0.0;
		memcpy(r.m[i], m[i], nr*sizeof(double));
	}
	
	MT->Gaussj(r.m, nr, SolNullInvMat, 1);
	
	for(int i=0; i<nr; i++)
		MT->Free(SolNullInvMat[i], sizeof(double));
	MT->Free(SolNullInvMat, nr*sizeof(double*));
	
	
	return r;
}


LCMatrix LCMatrix::Extract(int ir_e, int ic_e, int nr_e, int nc_e)
{
	LCMatrix r(MT, nr_e, nc_e);
	for(int ir=0; ir<nr_e; ir++)
		for(int ic=0; ic<nc_e; ic++)
			r.m[ir][ic] = m[ir_e+ir][ic_e+ic];
	return r;
}


LCMatrix LCMatrix::Proj(int * iCom, int NCom)
{
	LCMatrix * TmpIn;
	LCMatrix * TmpOut;
	LCMatrix * TmpSwp;
	LCMatrix Res(MT, NCom, NCom);
	int NRed(nc-NCom);
	int * iRed;
	int iR;
	
	
	TmpIn = new LCMatrix(MT, nc, nr);
	TmpOut = new LCMatrix(MT, nc, nr);
	iRed = (int*) MT->AllocMemory(NRed*sizeof(int));
	
	/*! Make list of index of component to reduce */
	iR = 0;
	for(int i=0; i<nc; i++){
		bool RedC(true);
		for(int k=0; k<NCom; k++)
			if(iCom[k]==i)
				RedC = false;
		if(RedC){
			iRed[iR] = i;
			iR++;
		}
	}
	
	/*! Copy local matrix */
	for(int ic=0; ic<nc; ic++)
		for(int ir=0; ir<nr; ir++)
			TmpOut->m[ir][ic] = m[ir][ic];
	/*! Reduce the matrix */	
	for(int k=0; k<NRed; k++){
		TmpSwp = TmpOut;
		TmpOut = TmpIn;
		TmpIn = TmpSwp;
		for(int ic=0; ic<nc; ic++)
			for(int ir=0; ir<nr; ir++){
				TmpOut->m[ir][ic] = TmpIn->m[ir][ic] - TmpIn->m[ir][iRed[k]] * TmpIn->m[iRed[k]][ic] / TmpIn->m[iRed[k]][iRed[k]];
				//Cout << TmpIn->m[ir][ic] << " - " << TmpIn->m[ir][iRed[k]] << " * " << TmpIn->m[iRed[k]][ic] << " / " << TmpIn->m[iRed[k]][iRed[k]] << " = " << TmpOut->m[ir][ic] << " ==> " << *TmpOut << Endl;
			}
		//Cout << "Intermediate matrix after reducing the " << k << " component " << iRed[k] << " : " << *TmpOut << Endl;
	}
	
	Cout << "Result in full form : " << *TmpOut << Endl;
	
	/*! Put result in output LCMatrix (smaller the local one) */
	for(int ic=0; ic<NCom; ic++)
		for(int ir=0; ir<NCom; ir++)
			Res.m[ir][ic] = TmpOut->m[iCom[ir]][iCom[ic]];
	
	
	// \todo FINISH THIS FUNCTION !!!
	
	return(Res);
}


double LCMatrix::DiffWithUnit()
{
	double InvResid(0.);
	for(int ir=0; ir<nr; ir++){
		for(int ic=0; ic<nc; ic++){
			if(ic==ir)
				InvResid += fabs(m[ir][ic] - 1.0);
			else
				InvResid += fabs(m[ir][ic]);
		}
	}
	return(InvResid);
}


void LCMatrix::EigenVsym(LCMatrix & EigenVect, LCMatrix & EigenVal)
{
#ifdef _GSLLINKED_
	EigenVsymCompute(EigenVect, EigenVal);
#else
	throw std::invalid_argument("ERROR in LCMatrix::EigenVsym : Sorry but we cannot compute Eigen vetors and values because GSL is not linked !");
#endif
}



void LCMatrix::DispMathematica()
{
	for(int ir=0; ir<nr; ir++){
		if(ir==0)
			Cout << "{{";
		else 
			Cout << ",{";
		for(int ic=0; ic<nc; ic++){
			if(ic==0)
				Cout <<  m[ir][ic];
			else
				Cout << "," <<  m[ir][ic];
		}
		Cout << "}";
	}
	Cout << "}" << Endl;
}




// *******************************
// * Other methods requiring GSL *
// *******************************

//#ifdef _GSLLINKED_

LCMatrix LCMatrix::InvSVD()
{
	LCMatrix r(MT, nr, nc);
	LCMatrix u(MT,0,0), s(MT,0,0), v(MT,0,0);
	
	SVD(u,s,v);
	for(int i=0; i<nc; i++)
		s.m[i][i] = 1.0/s.m[i][i];
	r = v*s*u.T();
	
	return r;
}


void LCMatrix::SVD(LCMatrix &U, LCMatrix &S, LCMatrix &V)
{
	
	if(nc > nr )
		throw std::invalid_argument("LCMatrix::SVD : Number of rows must be larger than columns.");
	
	gsl_matrix *m1 = gsl_matrix_alloc(nr, nc);
	gsl_matrix *v = gsl_matrix_alloc(nc, nc);
	gsl_vector *s = gsl_vector_alloc(nc);
	
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nc; ic++)
			gsl_matrix_set(m1,ir,ic,m[ir][ic]);
	
	gsl_linalg_SV_decomp_jacobi(m1, v, s);
	/*
	 Cout << "m1 = " << Endl;
	 for(int ir=0; ir<nr; ir++){
	 for(int ic=0; ic<nc; ic++)
	 Cout << gsl_matrix_get(m1,ir,ic) << " \t ";
	 Cout << Endl;
	 }
	 Cout << Endl;
	 
	 Cout << "v^T = " << Endl;
	 for(int ic=0; ic<nc; ic++){
	 for(int ir=0; ir<nc; ir++)
	 Cout << gsl_matrix_get(v,ir,ic) << " \t ";
	 
	 Cout << Endl;
	 }
	 Cout << Endl;
	 
	 Cout << "s = " << Endl;
	 for(int ir=0; ir<nc; ir++){
	 Cout << gsl_vector_get(s,ir) << Endl;
	 }
	 */
	
	U.init(MT,nr,nc);
	S.init(MT,nc,nc);
	V.init(MT,nc,nc);
	
	for(int ic=0; ic<nc; ic++){
		for(int ir=0; ir<nr; ir++){
			U.m[ir][ic] = gsl_matrix_get(m1,ir,ic);
			if(ir<nc){
				V.m[ir][ic] = gsl_matrix_get(v,ir,ic);
			}
		}
		S.m[ic][ic] = gsl_vector_get(s,ic);
	}
	
	//Cout << "U.S.V^T" << U*S*V.T() << Endl;
	
	gsl_matrix_free(v);
	gsl_matrix_free(m1);
	gsl_vector_free(s);
	
}

LCMatrix LCMatrix::InvLURefine()
{
	if(nr!=nc)
		throw std::invalid_argument("LCMatrix::InvGJ : For inversion using Gauss-Jordan pivot, it must be a square matrix.");

	
	//! **** First allocation 
	LCMatrix r(MT, nr, nr);
	gsl_matrix *m0 = gsl_matrix_alloc(nr, nr);
	gsl_vector *rhs = gsl_vector_alloc(nr);
	gsl_vector *x = gsl_vector_alloc(nr);
	gsl_permutation * p = gsl_permutation_alloc(nr);
	int signum;
	gsl_matrix *m1 = gsl_matrix_alloc(nr, nr);
	gsl_vector *residual = gsl_vector_alloc(nr);
	
	//! **** Copy the matrix
	for (int ir=0; ir<nr; ir++){
		gsl_vector_set(rhs, ir, 0.0);
		for(int ic=0; ic<nc; ic++){
			gsl_matrix_set(m0, ir, ic, m[ir][ic]);
		}	
	}
	
	//! **** Sa allocation
	gsl_matrix_memcpy(m1, m0);
	
	//! **** Set LU decomposition
	int info = gsl_linalg_LU_decomp (m0, p, &signum);
	//   gsl_linalg_LU_solve(m, p, rhs, x); // third method
	//   info = gsl_linalg_LU_refine(m1, m, p, rhs, x, residual);

	//! **** 
	for (int ir=0; ir<nr; ir++){
		//    std::cout << i << "     "  <<gsl_vector_get(x, i) << "   " << gsl_vector_get(residual, i)<< std::Endl;;
		
		//! *** Set the unit 
		for (int ic=0; ic<nr; ic++){
			gsl_vector_set(rhs, ic, 0.0);
			if (ir == ic)
				gsl_vector_set(rhs, ic, 1.0);
		}
		
		//! *** Solve LU decomposition
		gsl_linalg_LU_solve(m0, p, rhs, x);
		
		//! *** Refine the result of LU
		info = gsl_linalg_LU_refine(m1, m0, p, rhs, x, residual);
		
		//! *** Put the resulting vector in the inverse
		for (int ic=0; ic<nr; ic++){
			r.m[ic][ir] = gsl_vector_get(x, ic);
		}
	}
	return r;
}


void LCMatrix::EigenVsymCompute(LCMatrix & EigenVect, LCMatrix & EigenVal)
{
	if(nc != nr )
		throw std::invalid_argument("LCMatrix::SVD : LCMatrix must be square (and symetric).");
	
	gsl_matrix *m1 = gsl_matrix_alloc(nr, nr);
	gsl_vector *eval = gsl_vector_alloc (nr);
	gsl_matrix *evec = gsl_matrix_alloc (nr, nr);
	
	for(int ir=0; ir<nr; ir++)
		for(int ic=0; ic<nr; ic++)
			gsl_matrix_set(m1,ir,ic,m[ir][ic]);
	
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (nr);
	
	gsl_eigen_symmv (m1, eval, evec, w);
	
	gsl_eigen_symmv_free (w);
	
	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
	
	
	EigenVect.init(MT, nr, nr);
	EigenVal.init(MT, nr, 1);
	
	for(int ir=0; ir<nr; ir++){
		EigenVal.m[ir][0] = gsl_vector_get(eval, ir);
		for(int ic=0; ic<nr; ic++)
			EigenVect.m[ir][ic] = gsl_matrix_get(evec,ir,ic);
	}
	
	gsl_matrix_free(m1);
	gsl_matrix_free(evec);
	gsl_vector_free(eval);
}

//#endif


// end of LISACode-Matrix.cpp
