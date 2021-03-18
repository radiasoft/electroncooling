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
                       double *chisq, int n_bins);


    void histogram(const double *x, int n, int n_bins, data *output);

};

//These have to be outside of the class definition
double gaussian(const double a, const double b, const double c, const double t);
double double_gaussian(const double a1, const double b1, const double c1,
                       const double a2, const double b2, const double c2,
                       const double t);


int func_f_gaus (const gsl_vector * x, void *params, gsl_vector * f);
int func_df_gaus (const gsl_vector * x, void *params, gsl_matrix * J);
int func_fvv_gaus (const gsl_vector * x, const gsl_vector * v, void *params, gsl_vector * fvv);

int func_f_dbl_gaus (const gsl_vector * x, void *params, gsl_vector * f);
int func_df_dbl_gaus (const gsl_vector * x, void *params, gsl_matrix * J);
int func_fvv_dbl_gaus (const gsl_vector * x, const gsl_vector * v, void *params, gsl_vector * fvv);


#endif // FIT_H
