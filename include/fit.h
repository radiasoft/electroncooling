#ifndef FIT_H
#define FIT_H

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <numeric>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

struct data{
  double *t;
  double *y;
  size_t n;
  //TODO: Store fit function pointer in here? then it can be called from func_f
  //TODO: Also store Jacobian function pointer for more abstraction
};


class fit{

public:
    //For gaussian fit
  static void callback(const size_t iter, void *params,
                       const gsl_multifit_nlinear_workspace *w);
  void solve_system(gsl_vector *x, gsl_multifit_nlinear_fdf *fdf,
                    gsl_multifit_nlinear_parameters *params, double *chisq);
  void gaus_fit(double *x, unsigned int n, double *amplitude, double *mean, double *sigma, double *chisq, int n_bins);
  void double_gaus_fit(double *x, unsigned int n,
                       double *amplitude1, double *mean1, double *sigma1,
                       double *amplitude2, double *mean2, double *sigma2,
                       double *chisq, int n_bins,bool print=false);


    void histogram(const double *x, int n, int n_bins, data *output);

};

//These have to be outside of the class definition
double gaussian(const double a, const double c, const double t);
double double_gaussian(const double a1, const double c1,
                       const double a2, const double c2,
                       const double t);


int func_f_gaus (const gsl_vector * x, void *params, gsl_vector * f);
int func_df_gaus (const gsl_vector * x, void *params, gsl_matrix * J);
int func_fvv_gaus (const gsl_vector * x, const gsl_vector * v, void *params, gsl_vector * fvv);

int func_f_dbl_gaus (const gsl_vector * x, void *params, gsl_vector * f);
int func_df_dbl_gaus (const gsl_vector * x, void *params, gsl_matrix * J);
int func_fvv_dbl_gaus (const gsl_vector * x, const gsl_vector * v, void *params, gsl_vector * fvv);

//Store both gaussian and double-gaussian fit results
// (extendable to other models later on)
class fit_results{

public:
    
    double amplitude_;
    double mean_;
    double sigma_;
    double chisq_;
    
    double amplitude1_;
    double amplitude2_;
    double mean1_;
    double mean2_;
    double sigma1_;
    double sigma2_;
    double chisq2_;
    
    //container for bigaussian emittance
    double emit1_ = 0.;
    double emit2_ = 0.;
    
    int n_;
    int n1_;
    int n2_;
    
    fit_results():amplitude_(0.),mean_(0.),sigma_(0.),chisq_(0.)
                    ,amplitude1_(0.),amplitude2_(0.),mean1_(0.),mean2_(0.)
                    ,sigma1_(0.),sigma2_(0.),chisq2_(0.){}
    
    
};

#endif // FIT_H
