#include "fit.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <chrono>

//Use GSL Libraries to fit gaussian distributions and
// double-gaussians for emittance and IBS calculations


double gaussian(const double a, const double c, const double t)
{
    //Force the fit mean to be 0.0
    const double z = (t - 0.0) / c;
    return (abs(a) * exp(-0.5 * z * z));
}

double double_gaussian(const double a1,  const double c1,
                       const double a2, const double c2, const double t)
{
  //Force the '2' to be the component that's the widest
  double a_wide,a_narrow,c_wide,c_narrow;
  a_wide = a2;
  c_wide = c2;
  a_narrow = a1;
  c_narrow = c1;

  a_wide = abs(a_wide);
  a_narrow = abs(a_narrow);
  c_wide = abs(c_wide);
  c_narrow = abs(c_narrow);

//Force the mean to be 0.0
  const double z1 = (t - 0.0) / c_narrow;
  const double z2 = (t - 0.0) / c_wide;
    
  const double g1 = a_narrow * exp(-0.5 * z1 * z1);
  const double g2 = a_wide * exp(-0.5 * z2 * z2);
  return ( g1 + g2 );
}

//Explicity setting Gaussian function
int func_f_gaus (const gsl_vector * x, void *params, gsl_vector * f)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);
  double c = gsl_vector_get(x, 1);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double yi = d->y[i];
      double y = gaussian(a, c, ti);

      gsl_vector_set(f, i, yi - y);
    }

  return GSL_SUCCESS;
}


int func_f_dbl_gaus (const gsl_vector * x, void *params, gsl_vector * f)
{
  struct data *d = (struct data *) params;
  double a1 = gsl_vector_get(x, 0);
  double c1 = gsl_vector_get(x, 1);
  double a2 = gsl_vector_get(x, 2);
  double c2 = gsl_vector_get(x, 3);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double yi = d->y[i];
      double y = double_gaussian(a1, c1, a2, c2, ti);

      gsl_vector_set(f, i, yi - y);
    }

  return GSL_SUCCESS;
}

//df only for Gaussian fit function
int func_df_gaus (const gsl_vector * x, void *params, gsl_matrix * J)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);
  double c = gsl_vector_get(x, 1);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double zi = (ti - 0.0) / c;
      double ei = exp(-0.5 * zi * zi);

      gsl_matrix_set(J, i, 0, -ei);
      //gsl_matrix_set(J, i, 1, -(a / c) * ei * zi);
      gsl_matrix_set(J, i, 1, -(a / c) * ei * zi * zi);
    }

  return GSL_SUCCESS;
}

int func_df_dbl_gaus (const gsl_vector * x, void *params, gsl_matrix * J)
{
  struct data *d = (struct data *) params;
  double a1 = gsl_vector_get(x, 0);
  double c1 = gsl_vector_get(x, 1);
  double a2 = gsl_vector_get(x, 2);
  double c2 = gsl_vector_get(x, 3);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti   = d->t[i];
      double zi1  = (ti - 0.0) / c1;
      double zi2  = (ti - 0.0) / c2;
      double ei1  = exp(-0.5 * zi1 * zi1);
      double ei2  = exp(-0.5 * zi2 * zi2);

      gsl_matrix_set(J, i, 0, -ei1);
 //     gsl_matrix_set(J, i, 1, -(a1 / c1) * ei1 * zi1);
      gsl_matrix_set(J, i, 1, -(a1 / c1) * ei1 * zi1 * zi1);
      gsl_matrix_set(J, i, 2, -ei2);
   ///   gsl_matrix_set(J, i, 4, -(a2 / c2) * ei2 * zi2);
      gsl_matrix_set(J, i, 3, -(a2 / c2) * ei2 * zi2 * zi2);
    }

  return GSL_SUCCESS;
}

int func_fvv_gaus (const gsl_vector * x, const gsl_vector * v,
          void *params, gsl_vector * fvv)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);
  double c = gsl_vector_get(x, 1);
  double va = gsl_vector_get(v, 0);
  double vc = gsl_vector_get(v, 1);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double zi = (ti - 0.0) / c;
      double ei = exp(-0.5 * zi * zi);
      double Dab = -zi * ei / c;
      double Dac = -zi * zi * ei / c;
      double Dbb = a * ei / (c * c) * (1.0 - zi*zi);
      double Dbc = a * zi * ei / (c * c) * (2.0 - zi*zi);
      double Dcc = a * zi * zi * ei / (c * c) * (3.0 - zi*zi);
      double sum;

      sum = //2.0 * va * vb * Dab +
            2.0 * va * vc * Dac +
            //      vb * vb * Dbb +
            //2.0 * vb * vc * Dbc +
                  vc * vc * Dcc;

      gsl_vector_set(fvv, i, sum);
    }

  return GSL_SUCCESS;
}

 int func_fvv_dbl_gaus (const gsl_vector * x, const gsl_vector * v,
          void *params, gsl_vector * fvv)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);
  double c = gsl_vector_get(x, 1);
  double va = gsl_vector_get(v, 0);
  double vc = gsl_vector_get(v, 1);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double zi = (ti - 0.0) / c;
      double ei = exp(-0.5 * zi * zi);
      double Dab = -zi * ei / c;
      double Dac = -zi * zi * ei / c;
      double Dbb = a * ei / (c * c) * (1.0 - zi*zi);
      double Dbc = a * zi * ei / (c * c) * (2.0 - zi*zi);
      double Dcc = a * zi * zi * ei / (c * c) * (3.0 - zi*zi);
      double sum;

      sum = //2.0 * va * vb * Dab +
            2.0 * va * vc * Dac +
            //      vb * vb * Dbb +
            //2.0 * vb * vc * Dbc +
                  vc * vc * Dcc;

      gsl_vector_set(fvv, i, sum);
    }

  return GSL_SUCCESS;
}


void fit::callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double avratio = gsl_multifit_nlinear_avratio(w);
  double rcond;

  (void) params; /* not used */

  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);

}

void fit::solve_system(gsl_vector *x, gsl_multifit_nlinear_fdf *fdf,
             gsl_multifit_nlinear_parameters *params, double *chisq)
{
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  const size_t max_iter = 100000;
  const double xtol = 1.0e-8;
  const double gtol = 1.0e-8;
  const double ftol = 1.0e-8;
  const size_t n = fdf->n;
  const size_t p = fdf->p;
  gsl_multifit_nlinear_workspace *work =
    gsl_multifit_nlinear_alloc(T, params, n, p);
  gsl_vector * f = gsl_multifit_nlinear_residual(work);
  gsl_vector * y = gsl_multifit_nlinear_position(work);
  int info;
  double chisq0, rcond;

  /* initialize solver */
  gsl_multifit_nlinear_init(x, fdf, work);

  /* store initial cost */
  gsl_blas_ddot(f, f, &chisq0);

  /* iterate until convergence */
  gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,
                              callback, NULL, &info, work);

  /* store final cost */
  gsl_blas_ddot(f, f, chisq);

  /* store cond(J(x)) */
  gsl_multifit_nlinear_rcond(&rcond, work);

  gsl_vector_memcpy(x, y);

  gsl_multifit_nlinear_free(work);
}


void fit::histogram(const double *x, int n, int n_bins, struct data *output){

//Construct a (normalized) histogram from the macroparticle data

  gsl_histogram * h = gsl_histogram_alloc (n_bins);
  std::vector<double> vec(x,x+n);
  double x_low = *std::min_element(std::begin(vec), std::end(vec));
  double x_max = *std::max_element(std::begin(vec), std::end(vec));

  gsl_histogram_set_ranges_uniform (h, x_low, x_max);

  for(int i=0;i<n;i++){
    gsl_histogram_increment (h, x[i]);
  }

  output->t = new double[n_bins];
  output->y = new double[n_bins];
  output->n = n_bins;

  for(int i=0;i<n_bins;i++){
      double lower, upper;
      gsl_histogram_get_range(h, i, &lower, &upper);
      output->t[i] = 0.5*((lower + upper));
      output->y[i] = gsl_histogram_get(h, i) / (double)n;
  }
}

void fit::gaus_fit(double *x, unsigned int n, double *amplitude, double *mean, double *sigma, double *chisq, int n_bins=200){

  const size_t p = 2;    /* number of model parameters */

  gsl_vector *xv = gsl_vector_alloc(p); //Stores fit parameters
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

  //Histogram the distribution, store the results in our struct
  struct data fit_data;
  histogram(x,n,n_bins,&fit_data);

  //Get the moments from the histogram, to set up the fitter's initial position
  double initial_amplitude;
  double initial_mean;
  double initial_sigma;

  std::vector<double> v(x,x+n);
  initial_amplitude = *std::max_element(std::begin(v), std::end(v));
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  initial_mean = sum / v.size();
  std::vector<double> diff(v.size());
  std::transform(v.begin(), v.end(), diff.begin(), [initial_mean](double x) { return x - initial_mean; });
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  initial_sigma = std::sqrt(sq_sum / v.size());

  //Set up the fitter's struct
  fdf.f   =  &func_f_gaus;
  fdf.df  =  &func_df_gaus;
  fdf.fvv =  &func_fvv_gaus;
  fdf.n   =  n_bins;
  fdf.p   =  p;
  fdf.params = &fit_data;

  /* starting point */
  gsl_vector_set(xv, 0, initial_amplitude);
  gsl_vector_set(xv, 1, initial_sigma);

  //Use gedesic acceleration?
  fdf_params.trs = gsl_multifit_nlinear_trs_lm; //slower levinburg-marquardt
  //fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel; // with acceleration from 2nd derivative

  solve_system(xv, &fdf, &fdf_params, chisq);

  //TODO: Calculate chisquared and return FOM?
  

  //TODO: Catch NaN's?

  //Set the outputs
  *amplitude = gsl_vector_get(xv, 0);
  *mean      = 0.0; //Fixed at 0
  *sigma     = gsl_vector_get(xv, 1);

  std::cout<<"Amplitude: "<<gsl_vector_get(xv, 0)<<" sigma "<<gsl_vector_get(xv, 1)<<std::endl;

  gsl_vector_free(xv);
}

 void fit::double_gaus_fit(double *x, unsigned int n, double *amplitude1, double *mean1, double *sigma1,double *amplitude2, double *mean2, double *sigma2, double *chisq, int n_bins=200,bool print){

  const size_t p = 4;    /* number of model parameters */

  gsl_vector *xv = gsl_vector_alloc(p);
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

  if(print){
      //Output the raw data and the histogram data to separate files for QC.
      remove("Raw_distribution.txt"); //Remove a file from an earlier write
      std::ofstream myfile;
      myfile.open ("Raw_distribution.txt");

      for(int i=0;i<n;i++){
          myfile<<x[i]<<"\n";
      }
      myfile.close();
  }
     
  //Histogram the distribution
  struct data fit_data;
     
  histogram(x,n,n_bins,&fit_data);

  //Get the gaussian moments from the histogram, to set up the fitter's initial position
  double initial_amplitude;
  double initial_mean;
  double initial_sigma;

  std::vector<double> v(x,x+n);
  initial_amplitude = *std::max_element(std::begin(v), std::end(v));
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  initial_mean = sum / v.size();
  std::vector<double> diff(v.size());
  std::transform(v.begin(), v.end(), diff.begin(), [initial_mean](double x) { return x - initial_mean; });
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  initial_sigma = std::sqrt(sq_sum / v.size());

  fdf.f   =  &func_f_dbl_gaus;
  fdf.df  =  &func_df_dbl_gaus;
  fdf.n   =  n_bins;
  fdf.p   =  p;
  fdf.params = &fit_data;

  /* starting point */
  gsl_vector_set(xv, 0, initial_amplitude); //Amplitude1
  gsl_vector_set(xv, 1, initial_sigma/2);   //Width1

  //Curve 2 is defined as having the wider distribution
  gsl_vector_set(xv, 2, initial_amplitude); //Amplitude2
  gsl_vector_set(xv, 3, initial_sigma*2);   //Width2

  //Use Geodesic Acceleration?
  // If not, we don't have to deal with the 2nd derivative
  fdf_params.trs = gsl_multifit_nlinear_trs_lm;  //Slower, no acceleration

  solve_system(xv, &fdf, &fdf_params, chisq);

  double interim_amp1   = gsl_vector_get(xv, 0);
  double interim_sigma1 = gsl_vector_get(xv, 1);
  double interim_amp2   = gsl_vector_get(xv, 2);
  double interim_sigma2 = gsl_vector_get(xv, 3);

  if(print){
      //Write the histogram data out, with the top line storing the fit parameters
      //remove("histogram.txt"); //Remove a file from an earlier write
      std::ofstream myfile;
      //myfile.open ("histogram.txt");
      time_t rawtime;
      struct tm * timeinfo;
      char buffer [80];

      time (&rawtime);
      timeinfo = localtime (&rawtime);
      strftime (buffer,80,"histogram_%I%M%S.txt",timeinfo);
      myfile.open(buffer);
     //Write all fit parameters out in the first line
      for (int i=0;i<p-1;i++){
          myfile << abs(gsl_vector_get(xv,i))<< ", "; 
      }  
      myfile<< abs(gsl_vector_get(xv,3))<<"\n";
     
      for(int i=0;i<fit_data.n;i++){
          myfile << fit_data.t[i] << ", " << fit_data.y[i] << "\n";
      }
      myfile.close();
      
      
      //Also, write out the list of macroparticles that are in 
      // the central core, so we can track them over time
      std::ofstream myfile2;
//      strftime (buffer,80,"core_%I%M%S.txt",timeinfo);
      //Construct number of milliseconds since now was called for filename
      typedef std::chrono::system_clock Clock;
      auto now = Clock::now();
      auto milliseconds = std::chrono::time_point_cast<std::chrono::milliseconds>(now);
//      auto fraction = now - seconds;
      time_t cnow = Clock::to_time_t(now);
      //auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::time_point_cast<std::chrono::seconds>(now));
      sprintf(buffer,"core_%d.txt",milliseconds);
      myfile2.open(buffer);
     //Write all fit parameters out in the first line
      for (int i=0;i<n;i++){
//          if(abs(x[i])<0.0001){ //for x
          if(abs(x[i])<0.1){ //for ds
              myfile2 << i<<", "<<x[i]<< "\n"; 
          }
      }  
      myfile2.close();
      
      //Grab a spot in the tail for monitoring the core/tail ratio
      std::ofstream myfile3;
      sprintf(buffer,"tail_%d.txt",milliseconds);
      myfile3.open(buffer);
      for (int i=0;i<n;i++){
//          if(x[i]<0.0012 && x[i]>0.001){ //for x
          if(x[i]<10.2 && x[i]>10.0){ //for ds
          myfile3 << i<<", "<<x[i]<< "\n"; 
          }
      }  
      myfile3.close();
      
      
  }
     
  //Set the outputs, make sure the narrow peak comes out first
  if(abs(interim_sigma1) < abs(interim_sigma2)){
    *amplitude1 = abs(interim_amp1);
    *mean1      = 0.0;
    *sigma1     = abs(interim_sigma1);
    *amplitude2 = abs(interim_amp2);
    *mean2      = 0.0;
    *sigma2     = abs(interim_sigma2);
  }
  else{
    *amplitude1 = abs(interim_amp2);
    *mean1      = 0.0;
    *sigma1     = abs(interim_sigma2);
    *amplitude2 = abs(interim_amp1);
    *mean2      = 0.0;
    *sigma2     = abs(interim_sigma1);      
  }
     
  gsl_vector_free(xv);
}
