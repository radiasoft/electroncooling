#ifndef FIT_H
#define FIT_H

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

struct data{
  double *t;
  double *y;
  size_t n;
};


class fit{

public:
    //For gaussian fit


    static void callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w);
    void solve_system(gsl_vector *x, gsl_multifit_nlinear_fdf *fdf,
             gsl_multifit_nlinear_parameters *params);
    void gaus_fit(double *x, unsigned int n, double *amplitude, double *mean, double *sigma, int n_bins);

    void histogram(const double *x, int n, int n_bins, data *output);



};

double gaussian(const double a, const double b, const double c, const double t);
int func_f (const gsl_vector * x, void *params, gsl_vector * f);
int func_df (const gsl_vector * x, void *params, gsl_matrix * J);
int func_fvv (const gsl_vector * x, const gsl_vector * v, void *params, gsl_vector * fvv);


#endif // FIT_H
