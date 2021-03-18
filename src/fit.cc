#include "fit.h"

//Use GSL Libraries to fit gaussian distributions and
// double-gaussians for emittance and IBS calculations

double gaussian(const double a, const double b, const double c, const double t)
{
  const double z = (t - b) / c;
  return (a * exp(-0.5 * z * z));
}

double double_gaussian(const double a1, const double b1, const double c1,
                       const double a2, const double b2, const double c2, const double t)
{
  //Force the '2' to be the component that's the widest
  double a_wide,a_narrow,b_wide,b_narrow,c_wide,c_narrow;
  a_wide = a2;
  b_wide = b2;
  c_wide = c2;
  a_narrow = a1;
  b_narrow = b1;
  c_narrow = c1;

  a_wide = abs(a_wide);
  a_narrow = abs(a_narrow);
  c_wide = abs(c_wide);
  c_narrow = abs(c_narrow);

  const double z1 = (t - b_narrow) / c_narrow;
  const double z2 = (t - b_wide) / c_wide;
  const double g1 = a_narrow * exp(-0.5 * z1 * z1);
  const double g2 = a_wide * exp(-0.5 * z2 * z2);
  return ( g1 + g2 );
}

//Explicity setting Gaussian function
int func_f_gaus (const gsl_vector * x, void *params, gsl_vector * f)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);
  double b = gsl_vector_get(x, 1);
  double c = gsl_vector_get(x, 2);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double yi = d->y[i];
      double y = gaussian(a, b, c, ti);

      gsl_vector_set(f, i, yi - y);
    }

  return GSL_SUCCESS;
}


int func_f_dbl_gaus (const gsl_vector * x, void *params, gsl_vector * f)
{
  struct data *d = (struct data *) params;
  //  double total = gsl_vector_get(x, 0);
  double a1 = gsl_vector_get(x, 0);
  double b1 = gsl_vector_get(x, 1);
  double c1 = gsl_vector_get(x, 2);
  double a2 = gsl_vector_get(x, 3);
  double b2 = gsl_vector_get(x, 4);
  double c2 = gsl_vector_get(x, 5);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double yi = d->y[i];
      double y = double_gaussian(a1, b1, c1, a2, b2, c2, ti);

      gsl_vector_set(f, i, yi - y);
    }

  return GSL_SUCCESS;
}

//df only for Gaussian fit function
int func_df_gaus (const gsl_vector * x, void *params, gsl_matrix * J)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);
  double b = gsl_vector_get(x, 1);
  double c = gsl_vector_get(x, 2);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double zi = (ti - b) / c;
      double ei = exp(-0.5 * zi * zi);

      gsl_matrix_set(J, i, 0, -ei);
      gsl_matrix_set(J, i, 1, -(a / c) * ei * zi);
      gsl_matrix_set(J, i, 2, -(a / c) * ei * zi * zi);
    }

  return GSL_SUCCESS;
}

int func_df_dbl_gaus (const gsl_vector * x, void *params, gsl_matrix * J)
{
  struct data *d = (struct data *) params;
  double a1 = gsl_vector_get(x, 0);
  double b1 = gsl_vector_get(x, 1);
  double c1 = gsl_vector_get(x, 2);
  double a2 = gsl_vector_get(x, 3);
  double b2 = gsl_vector_get(x, 4);
  double c2 = gsl_vector_get(x, 5);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti   = d->t[i];
      double zi1  = (ti - b1) / c1;
      double zi2  = (ti - b2) / c2;
      double ei1  = exp(-0.5 * zi1 * zi1);
      double ei2  = exp(-0.5 * zi2 * zi2);

      gsl_matrix_set(J, i, 0, -ei1);
      gsl_matrix_set(J, i, 1, -(a1 / c1) * ei1 * zi1);
      gsl_matrix_set(J, i, 2, -(a1 / c1) * ei1 * zi1 * zi1);
      gsl_matrix_set(J, i, 3, -ei2);
      gsl_matrix_set(J, i, 4, -(a2 / c2) * ei2 * zi2);
      gsl_matrix_set(J, i, 5, -(a2 / c2) * ei2 * zi2 * zi2);
    }

  return GSL_SUCCESS;
}

int func_fvv_gaus (const gsl_vector * x, const gsl_vector * v,
          void *params, gsl_vector * fvv)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);
  double b = gsl_vector_get(x, 1);
  double c = gsl_vector_get(x, 2);
  double va = gsl_vector_get(v, 0);
  double vb = gsl_vector_get(v, 1);
  double vc = gsl_vector_get(v, 2);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double zi = (ti - b) / c;
      double ei = exp(-0.5 * zi * zi);
      double Dab = -zi * ei / c;
      double Dac = -zi * zi * ei / c;
      double Dbb = a * ei / (c * c) * (1.0 - zi*zi);
      double Dbc = a * zi * ei / (c * c) * (2.0 - zi*zi);
      double Dcc = a * zi * zi * ei / (c * c) * (3.0 - zi*zi);
      double sum;

      sum = 2.0 * va * vb * Dab +
            2.0 * va * vc * Dac +
                  vb * vb * Dbb +
            2.0 * vb * vc * Dbc +
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
  double b = gsl_vector_get(x, 1);
  double c = gsl_vector_get(x, 2);
  double va = gsl_vector_get(v, 0);
  double vb = gsl_vector_get(v, 1);
  double vc = gsl_vector_get(v, 2);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double zi = (ti - b) / c;
      double ei = exp(-0.5 * zi * zi);
      double Dab = -zi * ei / c;
      double Dac = -zi * zi * ei / c;
      double Dbb = a * ei / (c * c) * (1.0 - zi*zi);
      double Dbc = a * zi * ei / (c * c) * (2.0 - zi*zi);
      double Dcc = a * zi * zi * ei / (c * c) * (3.0 - zi*zi);
      double sum;

      sum = 2.0 * va * vb * Dab +
            2.0 * va * vc * Dac +
                  vb * vb * Dbb +
            2.0 * vb * vc * Dbc +
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
  const size_t max_iter = 5000;
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
  //std::cout<<"chisq = "<<*chisq<<" "<<chisq0<<std::endl;
  /* store cond(J(x)) */
  gsl_multifit_nlinear_rcond(&rcond, work);

  gsl_vector_memcpy(x, y);

  gsl_multifit_nlinear_free(work);
}


void fit::histogram(const double *x, int n, int n_bins, struct data *output){

//Construct a histogram from the data

//TODO: To remove outliers, could sort by value, remove top 5% and bottom 5%
  // prior to range finding and histogram filling

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
      output->y[i] = gsl_histogram_get(h, i);
  }
}

void fit::gaus_fit(double *x, unsigned int n, double *amplitude, double *mean, double *sigma, double *chisq, int n_bins=200){

  const size_t p = 3;    /* number of model parameters */

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
  gsl_vector_set(xv, 1, initial_mean);
  gsl_vector_set(xv, 2, initial_sigma);

  //Use gedesic acceleration?
  //fdf_params.trs = gsl_multifit_nlinear_trs_lm; //slower levinburg-marquardt
  fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel; // with acceleration from 2nd derivative

  solve_system(xv, &fdf, &fdf_params, chisq);

  //TODO: Calculate chisquared and return FOM?


  //TODO: Catch NaN's?

  //Set the outputs
  *amplitude = gsl_vector_get(xv, 0);
  *mean      = gsl_vector_get(xv, 1);
  *sigma     = gsl_vector_get(xv, 2);

  //  std::cout<<"Amplitude: "<<gsl_vector_get(xv, 0)<<" mean "<<gsl_vector_get(xv, 1)<<" sigma "<<gsl_vector_get(xv, 2)<<std::endl;

  gsl_vector_free(xv);
}

 void fit::double_gaus_fit(double *x, unsigned int n, double *amplitude1, double *mean1, double *sigma1,double *amplitude2, double *mean2, double *sigma2, double *chisq, int n_bins=200){

  const size_t p = 6;    /* number of model parameters */

  gsl_vector *xv = gsl_vector_alloc(p);
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

  //Histogram the distribution
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

  //  std::cout<<initial_amplitude<<" "<<initial_mean<<" "<<initial_sigma<<std::endl;

  fdf.f   =  &func_f_dbl_gaus;
  fdf.df  =  &func_df_dbl_gaus;
  //fdf.fvv =  &func_fvv_dbl_gaus;
  fdf.n   =  n_bins;
  fdf.p   =  p;
  fdf.params = &fit_data;

  /* starting point */
  gsl_vector_set(xv, 0, 2000.0); //Amplitude1
  gsl_vector_set(xv, 1, 0.0);    //Center1
  gsl_vector_set(xv, 2, 0.01);   //Width1

  //Curve 2 is defined as having the wider distribution
  gsl_vector_set(xv, 3, 2000.0); //Amplitude2
  gsl_vector_set(xv, 4, 0.0);    //Center2
  gsl_vector_set(xv, 5, 0.5);   //Width2

  //Use Geodesic Acceleration?
  // If not, we don't have to deal with the 2nd derivative
  fdf_params.trs = gsl_multifit_nlinear_trs_lm;  //Slower, no acceleration
  //fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;

  solve_system(xv, &fdf, &fdf_params, chisq);

  double interim_amp1   = gsl_vector_get(xv, 0);
  double interim_mean1  = gsl_vector_get(xv, 1);
  double interim_sigma1 = gsl_vector_get(xv, 2);
  double interim_amp2   = gsl_vector_get(xv, 3);
  double interim_mean2  = gsl_vector_get(xv, 4);
  double interim_sigma2 = gsl_vector_get(xv, 5);
     
  //Set the outputs, make sure the narrow peak comes out first
  if(abs(interim_sigma1) < abs(interim_sigma2)){
    *amplitude1 = abs(interim_amp1);
    *mean1      = interim_mean1;
    *sigma1     = abs(interim_sigma1);
    *amplitude2 = abs(interim_amp2);
    *mean2      = interim_mean2;
    *sigma2     = abs(interim_sigma2);
  }
  else{
    *amplitude1 = abs(interim_amp2);
    *mean1      = interim_mean2;
    *sigma1     = abs(interim_sigma2);
    *amplitude2 = abs(interim_amp1);
    *mean2      = interim_mean1;
    *sigma2     = abs(interim_sigma1);      
  }
     
  gsl_vector_free(xv);
}
