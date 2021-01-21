#include "fit.h"

//Use GSL Libraries to fit gaussian distributions and
// double-gaussians for emittance and IBS calculations

double gaussian(const double a, const double b, const double c, const double t)
{
  const double z = (t - b) / c;
  return (a * exp(-0.5 * z * z));
}


int func_f (const gsl_vector * x, void *params, gsl_vector * f)
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

int func_df (const gsl_vector * x, void *params, gsl_matrix * J)
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

int func_fvv (const gsl_vector * x, const gsl_vector * v,
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

  fprintf(stderr, "iter %2zu: a = %.4f, b = %.4f, c = %.4f, |a|/|v| = %.4f cond(J) = %8.4f, |f(x)| = %.4f\n",
          iter,
          gsl_vector_get(x, 0),
          gsl_vector_get(x, 1),
          gsl_vector_get(x, 2),
          avratio,
          1.0 / rcond,
          gsl_blas_dnrm2(f));
}

void fit::solve_system(gsl_vector *x, gsl_multifit_nlinear_fdf *fdf,
             gsl_multifit_nlinear_parameters *params)
{
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  const size_t max_iter = 200;
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
  double chisq0, chisq, rcond;

  std::cout<<"Initializing solver"<<std::endl;

  /* initialize solver */

  std::cout<<x->size<<" "<<n<<" "<<p<<" "<<std::endl;
  gsl_multifit_nlinear_init(x, fdf, work);

  std::cout<<"Store initial cost"<<std::endl;

  /* store initial cost */
  gsl_blas_ddot(f, f, &chisq0);

  std::cout<<"Iterating..."<<std::endl;

  /* iterate until convergence */
  gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,
                              callback, NULL, &info, work);

  /* store final cost */
  gsl_blas_ddot(f, f, &chisq);

  /* store cond(J(x)) */
  gsl_multifit_nlinear_rcond(&rcond, work);

  gsl_vector_memcpy(x, y);

  /* print summary */

  fprintf(stderr, "NITER         = %zu\n", gsl_multifit_nlinear_niter(work));
  fprintf(stderr, "NFEV          = %zu\n", fdf->nevalf);
  fprintf(stderr, "NJEV          = %zu\n", fdf->nevaldf);
  fprintf(stderr, "NAEV          = %zu\n", fdf->nevalfvv);
  fprintf(stderr, "initial cost  = %.12e\n", chisq0);
  fprintf(stderr, "final cost    = %.12e\n", chisq);
  fprintf(stderr, "final x       = (%.12e, %.12e, %12e)\n",
          gsl_vector_get(x, 0), gsl_vector_get(x, 1), gsl_vector_get(x, 2));
  fprintf(stderr, "final cond(J) = %.12e\n", 1.0 / rcond);

  gsl_multifit_nlinear_free(work);
}


void fit::histogram(const double *x, int n, int n_bins, struct data *output){

//Construct a histogram from the data

  gsl_histogram * h = gsl_histogram_alloc (n_bins);
  std::vector<double> vec(x,x+n);
  double x_low = *std::min_element(std::begin(vec), std::end(vec));
  double x_max = *std::max_element(std::begin(vec), std::end(vec));

  gsl_histogram_set_ranges_uniform (h, x_low, x_max);

  for(int i=0;i<n;i++){
    gsl_histogram_increment (h, vec[n]);
  }

  output->t = new double[n_bins];
  output->y = new double[n_bins];
  output->n = n;

  for(int i=0;i<n_bins;i++){
      double lower, upper;
      gsl_histogram_get_range(h, i, &lower, &upper);
      output->t[i] = 0.5*((lower + upper));
      output->y[i] = gsl_histogram_get(h, i);
  }

}

void fit::gaus_fit(double *x, unsigned int n, double *amplitude, double *mean, double *sigma, int n_bins=200){


  const size_t p = 3;    /* number of model parameters */

  gsl_vector *f = gsl_vector_alloc(n);
  gsl_vector *xv = gsl_vector_alloc(p);
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

  //Histogram the distribution
  struct data fit_data;
  histogram(x,n,n_bins,&fit_data);

  for(int i=0;i<n_bins;i++){
    std::cout<<fit_data.t[i]<<" ";
      }
  std::cout<<std::endl;

  for(int i=0;i<n_bins;i++){
    if(fit_data.y[i] != 0.0){
      std::cout<<fit_data.y[i]<<" ";
    }
  }
  std::cout<<std::endl;


  fdf.f   =  &func_f;
  fdf.df  =  &func_df;
  fdf.fvv =  &func_fvv;
  fdf.n   =  n_bins;
  fdf.p   =  p;
  fdf.params = &fit_data;

  /* starting point */
  gsl_vector_set(xv, 0, 2000.0); //Amplitude
  gsl_vector_set(xv, 1, 0.0);    //Center
  gsl_vector_set(xv, 2, 0.01);   //Width

  fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;

  std::cout<<"Solving system"<<std::endl;

  solve_system(xv, &fdf, &fdf_params);

  *amplitude = gsl_vector_get(xv, 0);
  *mean      = gsl_vector_get(xv, 1);
  *sigma     = gsl_vector_get(xv, 2);

  std::cout<<amplitude<<" "<<mean<<" "<<sigma<<std::endl;

  gsl_vector_free(f);
  gsl_vector_free(xv);
}
