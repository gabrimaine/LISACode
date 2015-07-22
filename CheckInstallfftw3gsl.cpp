#include <cstdlib>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <fftw3.h>


// Compilation : g++ -I/path/to/include -L/path/to/lib  -lgsl -lgslcblas -lfftw3 -lm -lstdc++ CheckInstallfftw3gsl.cpp -o CheckTest


int main(int argc, char *argv[])
{
    //! **** Test fftw
    std::cout << std::endl << "Check fftw ..." << std::endl;
    double * tf;
    tf = (double *) fftw_malloc(5*sizeof(double));
    fftw_free(tf);
    
    std::cout << std::endl << "Check fftw -> OK !" << std::endl;
    
    
    //! **** Test gsl with an svd decompostion
    std::cout << std::endl << "Check gsl ..." << std::endl;
    int nc(3);
    gsl_matrix * m1 = gsl_matrix_alloc(nc,nc);
    gsl_matrix * v  = gsl_matrix_alloc(nc,nc);
    gsl_vector * s  = gsl_vector_alloc(nc);        
    
    gsl_matrix_set(m1,0,0,3);
    gsl_matrix_set(m1,0,1,5);
    gsl_matrix_set(m1,0,2,1);
    gsl_matrix_set(m1,1,0,-1);
    gsl_matrix_set(m1,1,1,2);
    gsl_matrix_set(m1,1,2,-3);
    gsl_matrix_set(m1,2,0,4);
    gsl_matrix_set(m1,2,1,2);
    gsl_matrix_set(m1,2,2,6);
    
    gsl_linalg_SV_decomp_jacobi(m1,v,s);
    
    std::cout << std::endl << "m1 : " << std::endl;
    for(int ir=0; ir<nc; ir++){
        for(int ic=0; ic<nc; ic++)
            std::cout << " " << gsl_matrix_get(m1,ir,ic);
        std::cout << std::endl;
    }
    
    std::cout << std::endl << "v : " << std::endl;
    for(int ir=0; ir<nc; ir++){
        for(int ic=0; ic<nc; ic++)
            std::cout << " " << gsl_matrix_get(v,ir,ic);
        std::cout << std::endl;
    }
    
    std::cout << std::endl << "s : " << std::endl;
    for(int ir=0; ir<nc; ir++)
        std::cout << " " << gsl_vector_get(s,ir);
    
    /*! ** Result for checking :
     m1 :
     0.521641 -0.714815 0.465757
     -0.200707 -0.633412 -0.747332
     0.829221 0.296358 -0.473882
     
     v :
     0.575589 -0.0627363 0.815329
     0.437733 -0.818536 -0.372005
     0.690715 0.571019 -0.443679
     
     s :
     8.83011 5.18997 0.30549
     */
    
    std::cout << std::endl << "Check gsl -> OK !" << std::endl;
    
    std::cout << std::endl << "Everything seems OK !" << std::endl;
    
    return 0;
}
