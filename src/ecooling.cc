#include "ecooling.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include "constants.h"
#include "cooler.h"
#include "force.h"
#include "functions.h"
#include "fit.h"

#include <fstream>
#include <vector>
#include <float.h>

//For fit test
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

std::unique_ptr<double []> x_bet, xp_bet, y_bet, yp_bet, ds, dp_p, x, y, xp, yp, ne;
std::unique_ptr<double []> force_x, force_y, force_z, v_tr, v_long;
std::unique_ptr<double []> x_spl, xp_spl, y_spl, yp_spl, ds_spl, dp_p_spl;
//double *x_bet, *xp_bet, *y_bet, *yp_bet, *ds, *dp_p, *x, *y, *xp, *yp, *ne;
//double *force_x, *force_y, *force_z, *v_tr, *v_long;
//double *x_spl, *xp_spl, *y_spl, *yp_spl, *ds_spl, *dp_p_spl;
double t_cooler;
// model_beam_count <0 : create sample ions;
//int model_beam_count = -1;
// rms_dynamic_count < 0: Initialize and Clean;  = 0: Initialize, no Clean; >=1 no Initialize, no Clean;
int rms_dynamic_count = -1;
bool dynamic_flag = false;

EcoolRateParas::EcoolRateParas(const EcoolRateParas& old_ecool){
    ion_sample_ = old_ecool.ion_sample();
    n_sample_ = old_ecool.n_sample();
    n_tr_ = old_ecool.n_tr();
    n_long_ = old_ecool.n_long();
    shift_ = old_ecool.shift();
    bunch_separate_ = old_ecool.bunch_separate();
    n_long_sample_ = old_ecool.n_long_sample();
}

//Initialize the scratch variables for electron cooling rate calculation.
int assign_ecool_scratches(unsigned int n){
    x_bet.reset(new double[n]);
    xp_bet.reset(new double[n]);
    y_bet.reset(new double[n]);
    yp_bet.reset(new double[n]);
    ds.reset(new double[n]);
    dp_p.reset(new double[n]);
    x.reset(new double[n]);
    y.reset(new double[n]);
    xp.reset(new double[n]);
    yp.reset(new double[n]);
    ne.reset(new double[n]);
    force_x.reset(new double[n]);
    force_y.reset(new double[n]);
    force_z.reset(new double[n]);
    v_tr.reset(new double[n]);
    v_long.reset(new double[n]);

    srand(time(NULL));
//    srand(3.4345878966324);
    return 0;
}

//Calculate the transverse emittance statistically
double emit(double * x, double * xp, unsigned int n){
    double emit, x_mean, xp_mean, dlt2_x, dlt2_xp, dlt_xxp;
    x_mean = 0;
    xp_mean = 0;
    for(unsigned int i=0; i<n; ++i){
        x_mean += x[i];
       xp_mean += xp[i];
    }
    x_mean /= n;
    xp_mean /= n;

    dlt2_x = 0;
    dlt2_xp = 0;
    dlt_xxp = 0;
    for(unsigned int i=0; i<n; ++i){
        double x_adj = x[i]-x_mean;
        double xp_adj = xp[i]-xp_mean;
        dlt2_x += x_adj * x_adj;
        dlt2_xp += xp_adj * xp_adj;
        dlt_xxp += x_adj * xp_adj;
    }
//    emit = sqrt(fabs(dlt2_x*dlt2_xp-dlt_xxp*dlt_xxp))/n;
    emit = sqrt( dlt2_x * dlt2_xp - dlt_xxp * dlt_xxp )/n;

    return emit;
}

void pairsort(double a[], double b[], int n)
{
    std::pair<double, double> pairt[n];

    // Storing the respective array
    // elements in pairs.
    for (int i = 0; i < n; i++)
    {
        pairt[i].first = abs(a[i]);
        pairt[i].second = b[i];
    }

    // Sorting the pair array.
    std::sort(pairt, pairt + n);

    // Modifying original arrays
    for (int i = 0; i < n; i++)
    {
        a[i] = pairt[i].first;
        b[i] = pairt[i].second;
    }
}

double emit_fit(double *x, double * xp, unsigned int n, fit_results& FR,bool print)
{
    double x_mean, xp_mean;
    double x_sigma, xp_sigma;
    double x_amplitude,xp_amplitude;
    double chisq;
    int n_bins = 400;

    fit *fitter = new fit();
    fitter->gaus_fit(x,n,&x_amplitude,&x_mean,&x_sigma,&chisq, n_bins);
    fitter->gaus_fit(xp,n,&xp_amplitude,&xp_mean,&xp_sigma,&chisq, n_bins);

    FR.amplitude_ = x_amplitude;
    FR.mean_ = x_mean;
    FR.sigma_ = x_sigma;
    FR.chisq_ = chisq;

    //Calculate the normalized emittance asssuming a gaussian fit
    double dlt2_x = x_sigma * x_sigma * n;
    double dlt2_xp = xp_sigma * xp_sigma * n;
    double dlt_xxp = 0;
    for(unsigned int i=0; i<n; ++i){
        double x_adj = x[i]-x_mean;
        double xp_adj = xp[i]-xp_mean;
        dlt_xxp += x_adj * xp_adj;
    }
    //    emit = sqrt(fabs(dlt2_x*dlt2_xp-dlt_xxp*dlt_xxp))/n;
    double emit = sqrt( dlt2_x * dlt2_xp - dlt_xxp * dlt_xxp )/n;

    //Fit a double gaussian and propogate the
    double x_mean1, x_mean2, xp_mean1, xp_mean2;
    double x_sigma1, x_sigma2, xp_sigma1, xp_sigma2;
    double x_amp1, x_amp2, xp_amp1, xp_amp2;
    double chisq1;

    //We've guaranteed that the narrowest peak is 1
    fitter->double_gaus_fit(x,n, &x_amp1, &x_mean1, &x_sigma1,
                              &x_amp2, &x_mean2, &x_sigma2, &chisq1, n_bins,print);
    fitter->double_gaus_fit(xp,n, &xp_amp1, &xp_mean1, &xp_sigma1,
                              &xp_amp2, &xp_mean2, &xp_sigma2, &chisq1, n_bins,false); //don't print

    FR.amplitude1_ = x_amp1;
    FR.amplitude2_ = x_amp2;
    FR.mean1_ = x_mean1;
    FR.mean2_ = x_mean2;
    FR.sigma1_ = x_sigma1;
    FR.sigma2_ = x_sigma2;
    FR.chisq2_ = chisq1;

    double emit1, emit2;

    if(false){
    //Calculate the normalized emittance for the inner gaussian
      double dlt2_x1 = 0.0;
      double dlt2_xp1 = 0.0;
      double dlt_xxp1 = 0.0;
      for(unsigned int i=0; i<n; ++i){
        double x1_adj = x[i]-x_mean1;
        double xp1_adj = xp[i]-xp_mean1;
        dlt_xxp1 += x1_adj * xp1_adj;
        dlt2_x1 += x1_adj * x1_adj;
        dlt2_xp1 += xp1_adj * xp1_adj;
        dlt_xxp1  += x1_adj * xp1_adj;
      }
//    emit = sqrt(fabs(dlt2_x*dlt2_xp-dlt_xxp*dlt_xxp))/n;
      emit1 = sqrt( dlt2_x1 * dlt2_xp1 - dlt_xxp1 * dlt_xxp1 )/n;

    //Calculate the normalized emittance for the outer gaussian
      double dlt2_x2  = 0.0;
      double dlt2_xp2 = 0.0;
      double dlt_xxp2 = 0.0;
      for(unsigned int i=0; i<n; ++i){
        double x2_adj = x[i]-x_mean2;
        double xp2_adj = xp[i]-xp_mean2;
        dlt_xxp2 += x2_adj * xp2_adj;
        dlt2_x2 += x2_adj * x2_adj;
        dlt2_xp2 += xp2_adj * xp2_adj;
        dlt_xxp2  += x2_adj * xp2_adj;
      }
//    emit = sqrt(fabs(dlt2_x*dlt2_xp-dlt_xxp*dlt_xxp))/n;
      emit2 = sqrt( dlt2_x2 * dlt2_xp2 - dlt_xxp2 * dlt_xxp2 )/n;

      int n1 = floor(n * (x_amp1*x_sigma1) / ( (x_amp1*x_sigma1) + (x_amp2*x_sigma2) ) );
      int n2 = n - n1;

      std::cout<<"Counts: "<<n<<" : "<<n1<<" : "<<n2<<std::endl;

    }

    if(true){
        //Split the emittance calculation up between the two distributions.
        // Integrate the two gaussian fits to determine the relative populations
        //Recall integral of gaussian = sqrt(2pi) * amplitude * sigma

        int n1 = floor(n * (x_amp1*x_sigma1) / ( (x_amp1*x_sigma1) + (x_amp2*x_sigma2) ) );
        int n2 = n - n1;

        //Sort x by centrality with a lambda. Make deep copies first
        double *x_sorted = new double[n];
        double *xp_sorted = new double[n];
        for(int i=0;i<n;i++){
            x_sorted[i] = x[i];
            xp_sorted[i] = xp[i];
        }

        //Sort ascending based on the x value, with the same ordering for xp
        pairsort(x_sorted, xp_sorted, n);

        double dlt2_x1 = 0.0;
        double dlt2_xp1 = 0.0;
        double dlt_xxp1 = 0.0;
        for(unsigned int i=0; i<n1; ++i){
            double x_adj1 = x_sorted[i]-x_mean1;
            double xp_adj1 = xp_sorted[i]-xp_mean1;
            dlt2_x1 += x_adj1 * x_adj1;
            dlt2_xp1 += xp_adj1 * xp_adj1;
            dlt_xxp1 += x_adj1 * xp_adj1;
        }

        emit1 = sqrt( dlt2_x1 * dlt2_xp1 - dlt_xxp1 * dlt_xxp1 )/n1;

        double dlt2_x2 = 0.0;
        double dlt2_xp2 = 0.0;
        double dlt_xxp2 = 0.0;
        for(unsigned int i=n1; i<n; ++i){
            double x_adj2 = x_sorted[i]-x_mean2;
            double xp_adj2 = xp_sorted[i]-xp_mean2;
            dlt2_x2 += x_adj2 * x_adj2;
            dlt2_xp2 += xp_adj2 * xp_adj2;
            dlt_xxp2 += x_adj2 * xp_adj2;
        }

        emit2 = sqrt( dlt2_x2 * dlt2_xp2 - dlt_xxp2 * dlt_xxp2 )/n2;

        std::cout<<"Counts: "<<n<<" : "<<n1<<" : "<<n2<<std::endl;
    }

    int n1 = floor(n * (x_amp1*x_sigma1) / ( (x_amp1*x_sigma1) + (x_amp2*x_sigma2) ) );
    int n2 = n - n1;

    FR.n_     = n;
    FR.n1_    = n1;
    FR.n2_    = n2;
    FR.emit1_ = emit1;
    FR.emit2_ = emit2;

    std::cout<<"Emittance: "<<emit<<" : "<<emit1<<" : "<<emit2<<std::endl;
    std::cout << "Gaussian: " << x_sigma << " :" << chisq/((double)n_bins - 3.) << std::endl;
    std::cout << " DblGaus: " << x_sigma1 << " " << x_sigma2 << " : " << chisq1/((double)n_bins - 6.) << std::endl;


    return emit;
}

//Calculate the longitudinal emittance as (dp/p)^2/n
double emit_p(double * dp_p, unsigned int n){
    double emit_p = 0;
    double dp_p_mean = 0;
    for(unsigned int i=0; i<n; ++i){
        dp_p_mean += dp_p[i];
    }
    dp_p_mean /= n;

    for(unsigned int i=0; i<n; ++i){
        double dp_p_adj = dp_p[i] - dp_p_mean;
        emit_p += dp_p_adj * dp_p_adj;
    }
    emit_p /= n;

    return emit_p;
}


double emit_p_fit(double * dp_p, unsigned int n, fit_results& FR,bool print)
{
    double emit_p = 0.0;
    double dp_p_mean, dp_p_amplitude, dp_p_sigma, chisq;
    int n_bins = 400;
    double emit1 = 0.0;
    double emit2 = 0.0;


    fit *fitter = new fit();
    fitter->gaus_fit(dp_p,n,&dp_p_amplitude,&dp_p_mean,&dp_p_sigma,&chisq, n_bins);

    FR.amplitude_ = dp_p_amplitude;
    FR.mean_ = dp_p_mean;
    FR.sigma_ = dp_p_sigma;
    FR.chisq_ = chisq;


    for(unsigned int i=0; i<n; ++i){
        double dp_p_adj = dp_p[i] - dp_p_mean;
        emit_p += dp_p_adj * dp_p_adj;
    }
    emit_p /= n;


    //Fit a double gaussian
    double dp_p_mean1, dp_p_mean2;
    double dp_p_sigma1, dp_p_sigma2;
    double dp_p_amp1, dp_p_amp2;
    double chisq1;

    //We've guaranteed that the narrowest peak is 1
    fitter->double_gaus_fit(dp_p,n, &dp_p_amp1, &dp_p_mean1, &dp_p_sigma1,
                              &dp_p_amp2, &dp_p_mean2, &dp_p_sigma2, &chisq1, n_bins,print);

    FR.amplitude1_ = dp_p_amp1;
    FR.amplitude2_ = dp_p_amp2;
    FR.mean1_ = dp_p_mean1;
    FR.mean2_ = dp_p_mean2;
    FR.sigma1_ = dp_p_sigma1;
    FR.sigma2_ = dp_p_sigma2;
    FR.chisq2_ = chisq1;

    emit1 = dp_p_sigma1;

    emit2 = dp_p_sigma2;

    int n1 = floor(n * (dp_p_amp1*dp_p_sigma1) / ( (dp_p_amp1*dp_p_sigma1) + (dp_p_amp2*dp_p_sigma2) ) );
    int n2 = n - n1;

    std::cout<<"Long. Counts: "<<n<<" : "<<n1<<" : "<<n2<<std::endl;

    FR.emit1_ = emit1;
    FR.emit2_ = emit2;

    std::cout<<"Long. Emittance: "<<emit<<" : "<<emit1<<" : "<<emit2<<std::endl;
    std::cout << "Long. Gaussian: " << dp_p_sigma << " :" << chisq/((double)n_bins - 3.) << std::endl;
    std::cout << " Long. DblGaus: " << dp_p_sigma1 << " " << dp_p_sigma2 << " : " << chisq1/((double)n_bins - 6.) << std::endl;

    return emit_p;
}

double emit_p(double * dp_p, double * ds, Ring &ring, unsigned int n){
    double emit_p = 0;
    double dp_p_mean = 0;
    for(unsigned int i=0; i<n; ++i){
        dp_p_mean += dp_p[i];
    }
    dp_p_mean /= n;

    for(unsigned int i=0; i<n; ++i){
        double dp_p_adj = dp_p[i] - dp_p_mean;
        emit_p += dp_p_adj * dp_p_adj;
    }

    double emit_s = 0;
    double ds_mean = 0;
    for(unsigned int i=0; i<n; ++i){
        ds_mean += ds[i];
    }
    ds_mean /= n;

    for(unsigned int i=0; i<n; ++i){
        double ds_adj = ds[i] - ds_mean;
        emit_s += ds_adj*ds_adj;
    }
    emit_s /= (ring.beta_s()*ring.beta_s());

    emit_p = (emit_p + emit_s)/n;

    return emit_p;
}

double emit_p_fit(double * dp_p, double * ds, Ring &ring, unsigned int n, fit_results &FRp, fit_results &FRs, bool print){
    
    //Fit the dp/p distribution and ds distribution separately

    double emit_p = 0.0;
    double dp_p_mean, dp_p_amplitude, dp_p_sigma, chisq;
    int n_bins = 400;
    double emit1 = 0.0;
    double emit2 = 0.0;

    fit *fitter = new fit();
    fitter->gaus_fit(dp_p,n,&dp_p_amplitude,&dp_p_mean,&dp_p_sigma,&chisq, n_bins);

    FRp.amplitude_ = dp_p_amplitude;
    FRp.mean_ = dp_p_mean;
    FRp.sigma_ = dp_p_sigma;
    FRp.chisq_ = chisq;

    for(unsigned int i=0; i<n; ++i){
        double dp_p_adj = dp_p[i] - dp_p_mean;
        emit_p += dp_p_adj * dp_p_adj;
    }
    emit_p /= n;


    //Fit a double gaussian
    double dp_p_mean1, dp_p_mean2;
    double dp_p_sigma1, dp_p_sigma2;
    double dp_p_amp1, dp_p_amp2;
    double chisq1;

    //We've guaranteed that the narrowest peak is 1
    fitter->double_gaus_fit(dp_p,n, &dp_p_amp1, &dp_p_mean1, &dp_p_sigma1,
                              &dp_p_amp2, &dp_p_mean2, &dp_p_sigma2, &chisq1, n_bins,print);

    FRp.amplitude1_ = dp_p_amp1;
    FRp.amplitude2_ = dp_p_amp2;
    FRp.mean1_ = dp_p_mean1;
    FRp.mean2_ = dp_p_mean2;
    FRp.sigma1_ = dp_p_sigma1;
    FRp.sigma2_ = dp_p_sigma2;
    FRp.chisq2_ = chisq1;

//Now fit ds
    double ds_mean, ds_amplitude, ds_sigma;
    double ds_mean1, ds_mean2;
    double ds_sigma1, ds_sigma2;
    double ds_amp1, ds_amp2;

    fitter->gaus_fit(ds,n,&ds_amplitude,&ds_mean,&ds_sigma,&chisq, n_bins);

    FRs.amplitude_ = ds_amplitude;
    FRs.mean_ = ds_mean;
    FRs.sigma_ = ds_sigma;
    FRs.chisq_ = chisq;

    //Fit the double gaussian
    // We've guaranteed that the narrowest peak is 1
    fitter->double_gaus_fit(ds,n, &ds_amp1, &ds_mean1, &ds_sigma1,
                              &ds_amp2, &ds_mean2, &ds_sigma2, &chisq1, n_bins,print);

    FRs.amplitude1_ = ds_amp1;
    FRs.amplitude2_ = ds_amp2;
    FRs.mean1_ = ds_mean1;
    FRs.mean2_ = ds_mean2;
    FRs.sigma1_ = ds_sigma1;
    FRs.sigma2_ = ds_sigma2;
    FRs.chisq2_ = chisq1;
        
    double emit_s1 = ds_sigma1;
    double emit_s2 = ds_sigma2;

    emit_s1 /= (ring.beta_s()*ring.beta_s());
    emit_s2 /= (ring.beta_s()*ring.beta_s());
        
    emit_p = (dp_p_sigma1 + emit_s1)/n;
//    emit_p = (dp_p_sigma2 + emit_s2)/n;
    
    return emit_p;
}



//Generate Gaussian random number in S frame with given Twiss parameters
//First, rotate to O frame where alf = 0;
//Second, Generate x and xp with Gaussian random number in O frame
//Third, rotate back to S frame
int gaussian_bet_cod(double beta_xs, double alf_xs, double emit_x, double *x_bet, double *xp_bet, unsigned int n){

    double gamma_xs = (1+alf_xs*alf_xs)/beta_xs;
    double theta = atan(2*alf_xs/(gamma_xs-beta_xs))/2;     //rotation angle between O frame and S frame

    //Transfer matrix between O and S frames
    double matrix_os[2][2], matrix_so[2][2];
    matrix_os[0][0] = cos(theta);
    matrix_os[0][1] = -sin(theta);
    matrix_os[1][0] = sin(theta);
    matrix_os[1][1] = cos(theta);
    matrix_so[0][0] = matrix_os[0][0];
    matrix_so[0][1] = -matrix_os[0][1] ;
    matrix_so[1][0] = -matrix_os[1][0];
    matrix_so[1][1] = matrix_os[1][1];

    //Calculate beta and sigma in O frame
    double beta_xo = matrix_so[0][0]*matrix_so[0][0] * beta_xs-2*matrix_so[0][0]*matrix_so[0][1]*alf_xs+
                     matrix_so[0][1]*matrix_so[0][1]*gamma_xs;
    double sigma_xo = sqrt(emit_x*beta_xo);
    double sigma_xpo = sqrt(emit_x/beta_xo);
    //Generate x and xp in O frame
    gaussian_random(n, x_bet, sigma_xo);
    gaussian_random(n, xp_bet, sigma_xpo);
    gaussian_random_adjust(n, x_bet, sigma_xo);
    gaussian_random_adjust(n, xp_bet, sigma_xpo);

    //Rotate back to S frame
    for(unsigned int i=0; i<n;++i){
        double x = matrix_os[0][0]*x_bet[i]+matrix_os[0][1]*xp_bet[i];
        double xp = matrix_os[1][0]*x_bet[i]+matrix_os[1][1]*xp_bet[i];
        x_bet[i] = x;
        xp_bet[i] = xp;
    }
    return 0;
}

int adjust_disp(double dx, double *x_bet, double *dp_p, double *x, unsigned int n){
    for(unsigned int i=0; i<n; ++i) x[i] = x_bet[i]+dx*dp_p[i];
    return 0;
}

int assign_single_particle_scratches(EcoolRateParas &ecool_paras, Beam &ion) {
    unsigned int n_tr = ecool_paras.n_tr();
    if (x_spl.get() == nullptr) x_spl.reset(new double[n_tr]);
    if (y_spl.get() == nullptr) y_spl.reset(new double[n_tr]);
    if (xp_spl.get() == nullptr) xp_spl.reset(new double[n_tr]);
    if (yp_spl.get() == nullptr) yp_spl.reset(new double[n_tr]);
//    x_spl = new double[n_tr];
//    y_spl = new double[n_tr];
//    xp_spl = new double[n_tr];
//    yp_spl = new double[n_tr];
    if(ion.bunched()) {
        unsigned int n_long = ecool_paras.n_long();
        if (ds_spl.get() == nullptr) ds_spl.reset(new double[n_long]);
        if (dp_p_spl.get() == nullptr) dp_p_spl.reset(new double[n_long]);
//        ds_spl = new double[n_long];
//        dp_p_spl = new double[n_long];
    }
    return 0;
}

int config_ecooling(EcoolRateParas &ecool_paras, Beam &ion) {
    switch (ecool_paras.ion_sample()) {
        case IonSample::MONTE_CARLO: {
            assign_ecool_scratches(ecool_paras.n_sample());
            break;
        }
        case IonSample::SINGLE_PARTICLE: {
            if (!ion.bunched()) {
                ecool_paras.set_n_long(2);
                std::cout<<"Continuous ion beam! Only 2 possible velocities"<<std::endl;
            }
            assign_ecool_scratches(ecool_paras.n_sample());
            assign_single_particle_scratches(ecool_paras, ion);
            break;
        }
        default: {
            assert(false);
        }
    }
    return 0;
}

int single_particle_grid(unsigned int n_tr, unsigned int n_long, Beam &ion, Cooler &cooler){
    double alf_x = cooler.alpha_h();
    double alf_y = cooler.alpha_v();
    double dphi = 2.0*k_pi/n_tr;
    double phi = 0;
    for(unsigned int i=0; i<n_tr; ++i){
        x_spl[i] = sin(phi);
        y_spl[i] = sin(phi);
        xp_spl[i] = cos(phi)-alf_x*sin(phi);
        yp_spl[i] = cos(phi)-alf_y*sin(phi);
        phi += dphi;
    }

    if(ion.bunched()){
        phi = 0;
        dphi = 2.0*k_pi/n_long;
        for(unsigned int i=0; i<n_long; ++i){
            ds_spl[i] = sin(phi);
            dp_p_spl[i] = cos(phi);
            phi += dphi;
        }
    }
    return 0;
}

int ion_beam_model_single_particle(unsigned int n_tr, unsigned int n_long, Beam &ion, Cooler &cooler){
    double emit_x = ion.emit_x();
    double emit_y = ion.emit_y();
    double sigma_p = ion.dp_p();
    double beta_x = cooler.beta_h();
    double beta_y = cooler.beta_v();
    double dx = cooler.disp_h();
    double dy = cooler.disp_v();
    double dpx = cooler.der_disp_h();
    double dpy = cooler.der_disp_v();

    double y_amp = sqrt(2.0*emit_y*beta_y);
    double yp_amp = sqrt(2.0*emit_y/beta_y);
    double x_amp = sqrt(2.0*emit_x*beta_x);
    double xp_amp = sqrt(2.0*emit_x/beta_x);

    double ds_amp, dp_amp;
    if(ion.bunched()){  //bunched beam
        double sigma_s = ion.sigma_s();
        ds_amp = sqrt(2.0)*sigma_s;
        dp_amp = sqrt(2.0)*sigma_p;
    }

    unsigned int cnt = 0;
    for(unsigned int i=0; i<n_tr; ++i){
        double y_spl_tmp = y_amp*y_spl[i];
        double yp_spl_tmp = yp_amp*yp_spl[i];
        for(unsigned int j=0; j<n_tr; ++j){
            double x_spl_tmp = x_amp*x_spl[j];
            double xp_spl_tmp = xp_amp*xp_spl[j];
            if(ion.bunched()){  //bunched beam
                for(unsigned int k=0; k<n_long; ++k){
                    double ds_spl_tmp = ds_amp*ds_spl[k];
                    double dp_spl_tmp = dp_amp*dp_p_spl[k];
                    x_bet[cnt] = x_spl_tmp;
                    xp_bet[cnt] = xp_spl_tmp;
                    y_bet[cnt] = y_spl_tmp;
                    yp_bet[cnt] = yp_spl_tmp;
                    ds[cnt] = ds_spl_tmp;
                    dp_p[cnt] = dp_spl_tmp;
                    ++cnt;
                }
            }
            else{   //coasting beam, ds=s-s0 is set to be zero!
                for(int k=-1; k<3; k +=2){
                    x_bet[cnt] = x_spl_tmp;
                    xp_bet[cnt] = xp_spl_tmp;
                    y_bet[cnt] = y_spl_tmp;
                    yp_bet[cnt] = yp_spl_tmp;
                    dp_p[cnt] = k*sigma_p;
                    ++cnt;
                }
            }
        }
    }
    for(unsigned int i=0; i<cnt; ++i){
        x[i] = x_bet[i] + dp_p[i]*dx;
        xp[i] = xp_bet[i] + dp_p[i]*dpx;
        y[i] = y_bet[i] + dp_p[i]*dy;
        yp[i] = yp_bet[i] + dp_p[i]*dpy;
    }
    return 0;
}

int ion_beam_model_MonteCarlo_Gaussian(unsigned int n_sample, Beam &ion, EBeam &ebeam, Cooler &cooler){

    double beta_xs = cooler.beta_h();
    double beta_ys = cooler.beta_v();
    double alf_xs = cooler.alpha_h();
    double alf_ys = cooler.alpha_v();
    double emit_x = ion.emit_x();
    double emit_y = ion.emit_y();

//    gaussian_bet_cod(beta_xs, alf_xs, emit_x, x_bet, xp_bet, n_sample);
//    gaussian_bet_cod(beta_ys, alf_ys, emit_y, y_bet, yp_bet, n_sample);
//
//    double sigma_p = ion.dp_p();
//    gaussian_random(n_sample, dp_p, sigma_p);
//    gaussian_random_adjust(n_sample, dp_p, sigma_p);

    gaussian_bet_cod(beta_xs, alf_xs, emit_x, x_bet.get(), xp_bet.get(), n_sample);
    gaussian_bet_cod(beta_ys, alf_ys, emit_y, y_bet.get(), yp_bet.get(), n_sample);

    double sigma_p = ion.dp_p();
    gaussian_random(n_sample, dp_p.get(), sigma_p);
    gaussian_random_adjust(n_sample, dp_p.get(), sigma_p);

    //longitudinal sampling
    if(ion.bunched()) {
        double sigma_s = ion.sigma_s();
//        gaussian_random(n_sample, ds, sigma_s);
//        gaussian_random_adjust(n_sample, ds, sigma_s);
//
        gaussian_random(n_sample, ds.get(), sigma_s);
        gaussian_random_adjust(n_sample, ds.get(), sigma_s);
    }
//    else if(ebeam.bunched()) {      //Bunched electron beam to cool coasting ion beam
//        double length = ebeam.shape_->length();
//        uniform_random(n_sample, ds, -0.5*length, 0.5*length);
//        uniform_random_adjust(n_sample, ds);
//    }

    double dx = cooler.disp_h();
    double dy = cooler.disp_v();
    double dpx = cooler.der_disp_h();
    double dpy = cooler.der_disp_v();
//    adjust_disp(dx, x_bet, dp_p, x, n_sample);
//    adjust_disp(dy, y_bet, dp_p, y, n_sample);
//    adjust_disp(dpx, xp_bet, dp_p, xp, n_sample);
//    adjust_disp(dpy, yp_bet, dp_p, yp, n_sample);

    adjust_disp(dx, x_bet.get(), dp_p.get(), x.get(), n_sample);
    adjust_disp(dy, y_bet.get(), dp_p.get(), y.get(), n_sample);
    adjust_disp(dpx, xp_bet.get(), dp_p.get(), xp.get(), n_sample);
    adjust_disp(dpy, yp_bet.get(), dp_p.get(), yp.get(), n_sample);


    return 0;
}

int ion_beam_model_MonteCarlo_Gaussian(unsigned int n_sample, Beam &ion, Twiss &twiss){

    double beta_xs = twiss.bet_x;
    double beta_ys = twiss.bet_y;
    double alf_xs = twiss.alf_x;
    double alf_ys = twiss.alf_y;
    double emit_x = ion.emit_x();
    double emit_y = ion.emit_y();

//    gaussian_bet_cod(beta_xs, alf_xs, emit_x, x_bet, xp_bet, n_sample);
//    gaussian_bet_cod(beta_ys, alf_ys, emit_y, y_bet, yp_bet, n_sample);
//
//    double sigma_p = ion.dp_p();
//    gaussian_random(n_sample, dp_p, sigma_p);
//    gaussian_random_adjust(n_sample, dp_p, sigma_p);
//
//    //longitudinal sampling
//    if(ion.bunched()) {
//        double sigma_s = ion.sigma_s();
//        gaussian_random(n_sample, ds, sigma_s);
//        gaussian_random_adjust(n_sample, ds, sigma_s);
//    }

    gaussian_bet_cod(beta_xs, alf_xs, emit_x, x_bet.get(), xp_bet.get(), n_sample);
    gaussian_bet_cod(beta_ys, alf_ys, emit_y, y_bet.get(), yp_bet.get(), n_sample);

    double sigma_p = ion.dp_p();
    gaussian_random(n_sample, dp_p.get(), sigma_p);
    gaussian_random_adjust(n_sample, dp_p.get(), sigma_p);

    //longitudinal sampling
    if(ion.bunched()) {
        double sigma_s = ion.sigma_s();
        gaussian_random(n_sample, ds.get(), sigma_s);
        gaussian_random_adjust(n_sample, ds.get(), sigma_s);
    }

//    else if(ebeam.bunched()) {      //Bunched electron beam to cool coasting ion beam
//        double length = ebeam.shape_->length();
//        uniform_random(n_sample, ds, -0.5*length, 0.5*length);
//        uniform_random_adjust(n_sample, ds);
//    }

    double dx = twiss.disp_x;
    double dy = twiss.disp_y;
    double dpx = twiss.disp_dx;
    double dpy = twiss.disp_dy;
    adjust_disp(dx, x_bet.get(), dp_p.get(), x.get(), n_sample);
    adjust_disp(dy, y_bet.get(), dp_p.get(), y.get(), n_sample);
    adjust_disp(dpx, xp_bet.get(), dp_p.get(), xp.get(), n_sample);
    adjust_disp(dpy, yp_bet.get(), dp_p.get(), yp.get(), n_sample);

//    adjust_disp(dx, x_bet, dp_p, x, n_sample);
//    adjust_disp(dy, y_bet, dp_p, y, n_sample);
//    adjust_disp(dpx, xp_bet, dp_p, xp, n_sample);
//    adjust_disp(dpy, yp_bet, dp_p, yp, n_sample);

    return 0;
}

int ion_beam_model_MonteCarlo_Gaussian(unsigned int n_sample, Beam &ion, Cooler &cooler){

    double beta_xs = cooler.beta_h();
    double beta_ys = cooler.beta_v();
    double alf_xs = cooler.alpha_h();
    double alf_ys = cooler.alpha_v();
    double emit_x = ion.emit_x();
    double emit_y = ion.emit_y();

    gaussian_bet_cod(beta_xs, alf_xs, emit_x, x_bet.get(), xp_bet.get(), n_sample);
    gaussian_bet_cod(beta_ys, alf_ys, emit_y, y_bet.get(), yp_bet.get(), n_sample);

    double sigma_p = ion.dp_p();
    gaussian_random(n_sample, dp_p.get(), sigma_p);
    gaussian_random_adjust(n_sample, dp_p.get(), sigma_p);

    //longitudinal sampling
    if(ion.bunched()) {
        double sigma_s = ion.sigma_s();
        gaussian_random(n_sample, ds.get(), sigma_s);
        gaussian_random_adjust(n_sample, ds.get(), sigma_s);
    }

    double dx = cooler.disp_h();
    double dy = cooler.disp_v();
    double dpx = cooler.der_disp_h();
    double dpy = cooler.der_disp_v();
    adjust_disp(dx, x_bet.get(), dp_p.get(), x.get(), n_sample);
    adjust_disp(dy, y_bet.get(), dp_p.get(), y.get(), n_sample);
    adjust_disp(dpx, xp_bet.get(), dp_p.get(), xp.get(), n_sample);
    adjust_disp(dpy, yp_bet.get(), dp_p.get(), yp.get(), n_sample);

    return 0;
}

int ion_sample(EcoolRateParas &ecool_paras, Beam &ion, Ring &ring, Cooler &cooler, EBeam &ebeam) {
    switch (ecool_paras.ion_sample()) {
        case IonSample::SINGLE_PARTICLE: {
            if (rms_dynamic_count<1) single_particle_grid(ecool_paras.n_tr(), ecool_paras.n_long(), ion, cooler);
            ion_beam_model_single_particle(ecool_paras.n_tr(), ecool_paras.n_long(), ion, cooler);
            break;
        }
        case IonSample::MONTE_CARLO: {
            ion_beam_model_MonteCarlo_Gaussian(ecool_paras.n_sample(), ion, ebeam, cooler);
            break;
        }
        default: {
            perror("Error in ion beam model definition!");
            break;
        }
    }
    return 0;
}

int ion_sample(EcoolRateParas &ecool_paras, Beam &ion, Ring &ring, Cooler &cooler) {
    switch (ecool_paras.ion_sample()) {
        case IonSample::SINGLE_PARTICLE: {
            if (rms_dynamic_count<1) single_particle_grid(ecool_paras.n_tr(), ecool_paras.n_long(), ion, cooler);
            ion_beam_model_single_particle(ecool_paras.n_tr(), ecool_paras.n_long(), ion, cooler);
            break;
        }
        case IonSample::MONTE_CARLO: {
            ion_beam_model_MonteCarlo_Gaussian(ecool_paras.n_sample(), ion, cooler);
            break;
        }
        default: {
            perror("Error in ion beam model definition!");
            break;
        }
    }
    return 0;
}

int electron_density(EcoolRateParas &ecool_paras, Beam &ion, EBeam &ebeam) {
    if(ecool_paras.shift()){
        double cx, cy, cz;
        ion.center(cx, cy, cz);
//        ebeam.shape_->density(x, y, ds, ebeam, ne, ecool_paras.n_sample(), cx, cy, cz);
        ebeam.shape_->density(x.get(), y.get(), ds.get(), ebeam, ne.get(), ecool_paras.n_sample(), cx, cy, cz);
    }
    else{
//        ebeam.shape_->density(x, y, ds, ebeam, ne, ecool_paras.n_sample());
         ebeam.shape_->density(x.get(), y.get(), ds.get(), ebeam, ne.get(), ecool_paras.n_sample());
    }
    return 0;
}

int space_to_dynamic(unsigned int n_sample, Beam &ion) {
    double v = ion.beta()*k_c;
    for(unsigned int i=0; i<n_sample; ++i) {
//        v_long[i] = dp_p[i]*v/(ion.gamma()*ion.gamma());  //Convert from dp/p to dv/v
//        v_long[i] /= (1-(v_long[i]+v)*ion.beta()/k_c);    //Convert to beam frame, when v_long<<v, canceled with the above line.
        v_long[i] = dp_p[i]*v;
        v_tr[i] = sqrt(xp[i]*xp[i]+yp[i]*yp[i])*v;
    }
    return 0;
}

int space_to_dynamic(unsigned int n_sample, Beam &ion, EBeam &ebeam) {
    double v = ion.beta()*k_c;
    ParticleBunch* prtl_bunch = nullptr;
    if(ebeam.shape_->shape() == Shape::PARTICLE_BUNCH)
        prtl_bunch = dynamic_cast<ParticleBunch*>(ebeam.shape_);
    for(unsigned int i=0; i<n_sample; ++i) {
        v_tr[i] = sqrt(xp[i]*xp[i]+yp[i]*yp[i])*v;
        v_long[i] = dp_p[i]*v;
    }

    if(prtl_bunch->corr())
        for(unsigned int i=0; i<n_sample; ++i)
            v_long[i] -= prtl_bunch->v_avg_z.at(i);
    return 0;
}

int restore_velocity(unsigned int n_sample, EBeam &ebeam) {
    ParticleBunch* prtl_bunch = nullptr;
    if(ebeam.shape_->shape() == Shape::PARTICLE_BUNCH) {
        prtl_bunch = dynamic_cast<ParticleBunch*>(ebeam.shape_);    //Electron beam is defined by particles from file.
        if(prtl_bunch->corr())
            for(unsigned int i=0; i<n_sample; ++i)
                v_long[i] += prtl_bunch->v_avg_z.at(i);         //Restore the original longitudinal ion velocity.
    }
    return 0;
}

//TODO: Add a check to make sure we never convert to beam_frame from beam_frame
int beam_frame(unsigned int n_sample, double gamma_e) {
    double gamma_e_inv = 1./gamma_e;
    for(unsigned int i=0; i<n_sample; ++i){
        v_tr[i] *= gamma_e;
        ne[i] *= gamma_e_inv;
    }
    t_cooler /= gamma_e;
    return 0;
}

int lab_frame(unsigned int n_sample, double gamma_e) {
    double gamma_e_inv = 1/gamma_e;
    for(unsigned int i=0; i<n_sample; ++i){
            force_x[i] *= gamma_e_inv;
            v_tr[i] *= gamma_e_inv;
    }
    t_cooler *= gamma_e;
    return 0;
}

//Calculate friction force
int force(unsigned int n_sample, Beam &ion, EBeam &ebeam, Cooler &cooler, ForceParas &force_paras) {
    //set parameters for friction force calculation
    force_paras.set_magnetic_field(cooler.magnetic_field());
    force_paras.set_time_cooler(t_cooler);

    if(ebeam.shape_->shape() == Shape::PARTICLE_BUNCH) {
        ParticleBunch* prtl_bunch = dynamic_cast<ParticleBunch*>(ebeam.shape_);
        force_paras.set_ptr_d_paral_e(prtl_bunch->v_rms_l);
        force_paras.set_ptr_d_perp_e(prtl_bunch->v_rms_t);
    }
    else {
        force_paras.set_d_perp_e(ebeam.v_rms_tr());
        force_paras.set_d_paral_e(ebeam.v_rms_long());
    }

    //Calculate friction force
//    friction_force(ion.charge_number(),n_sample,v_tr, v_long, ne, force_paras, force_x, force_z);
    friction_force(ion.charge_number(),n_sample,v_tr.get(), v_long.get(), ne.get(), force_paras, force_x.get(), force_z.get());

    return 0;
}

//Distribute to x and y direction (in lab frame)
int force_distribute(unsigned int n_sample, Beam &ion) {
    double v0 = ion.beta()*k_c;
    for(unsigned int i=0; i<n_sample; ++i){
        //Dividing by v_tr here leads to a problem when v_tr = 0
        if(v_tr[i] == 0.0){
            force_y[i] = 0.0;
            force_x[i] = 0.0;
        }
        else{
            force_y[i] = yp[i]!=0?force_x[i]*yp[i]*v0/v_tr[i]:0;
            force_x[i] = xp[i]!=0?force_x[i]*xp[i]*v0/v_tr[i]:0;
        }
    }
    return 0;
}

int original_emittance(EcoolRateParas &ecool_paras, Beam &ion, double &emit_x0, double &emit_y0, double &emit_z0) {
    switch(ecool_paras.ion_sample()){
        case IonSample::SINGLE_PARTICLE: {
            emit_x0 = ion.emit_x();
            emit_y0 = ion.emit_y();
            if(ion.bunched()) emit_z0 = (2*ion.dp_p()*ion.dp_p());
            else emit_z0 = (ion.dp_p()*ion.dp_p());
            break;
        }
        case IonSample::MONTE_CARLO: {
            unsigned int n_sample = ecool_paras.n_sample();
            emit_x0 = emit(x_bet.get(), xp_bet.get(), n_sample);
            emit_y0 = emit(y_bet.get(), yp_bet.get(), n_sample);
            emit_z0 = emit_p(dp_p.get(), n_sample);
            break;
        }
        default:{
        perror("Error in defining ion beam model!");
        break;
        }
    }
    return 0;
}


int original_emittance(EcoolRateParas &ecool_paras, Ring &ring, Beam &ion, double &emit_x0, double &emit_y0, double &emit_z0) {
    switch(ecool_paras.ion_sample()){
        case IonSample::SINGLE_PARTICLE: {
            emit_x0 = ion.emit_x();
            emit_y0 = ion.emit_y();
            if(ion.bunched()) emit_z0 = (2*ion.dp_p()*ion.dp_p());
            else emit_z0 = (ion.dp_p()*ion.dp_p());
            break;
        }
        case IonSample::MONTE_CARLO: {
            unsigned int n_sample = ecool_paras.n_sample();
            emit_x0 = emit(x_bet.get(), xp_bet.get(), n_sample);
            emit_y0 = emit(y_bet.get(), yp_bet.get(), n_sample);
            if (ion.bunched()) emit_z0 = emit_p(dp_p.get(), ds.get(), ring, n_sample);
            else emit_z0 = emit_p(dp_p.get(), n_sample);
            break;
        }
        default:{
        perror("Error in defining ion beam model!");
        break;
        }
    }
    return 0;
}

int apply_kick(unsigned int n_sample, Beam &ion) {
    double p0 = ion.gamma()*ion.mass()*1e6*k_e*ion.beta()/k_c;
    for(unsigned int i=0; i<n_sample; ++i){
        xp[i] = !iszero(xp[i])?xp[i]*exp(force_x[i]*t_cooler/(p0*xp[i])):xp[i];
        yp[i] = !iszero(yp[i])?yp[i]*exp(force_y[i]*t_cooler/(p0*yp[i])):yp[i];
        dp_p[i] = !iszero(dp_p[i])?dp_p[i]*exp(force_z[i]*t_cooler/(p0*dp_p[i])):dp_p[i];
    }
    return 0;
}

int new_emittance(EcoolRateParas ecool_paras, Beam ion, Ring &ring, Cooler &cooler, double &emit_x, double &emit_y,
                  double &emit_z) {
    double dx = cooler.disp_h();
    double dy = cooler.disp_v();
    double dpx = cooler.der_disp_h();
    double dpy = cooler.der_disp_v();

    unsigned int n_sample = ecool_paras.n_sample();
    for(unsigned int i=0; i<n_sample; ++i) {
        x_bet[i] = x[i] - dx*dp_p[i];
        xp_bet[i] = xp[i] - dpx*dp_p[i];
        y_bet[i] = y[i] - dy*dp_p[i];
        yp_bet[i] = yp[i] - dpy*dp_p[i];
    }

    switch (ecool_paras.ion_sample()) {
        case IonSample::SINGLE_PARTICLE: {
            double alf_x = cooler.alpha_h();
            double alf_y = cooler.alpha_v();
            double beta_x = cooler.beta_h();
            double beta_y = cooler.beta_v();
            double gamma_x = (1+alf_x*alf_x)/beta_x;
            double gamma_y = (1+alf_y*alf_y)/beta_y;
            emit_x = 0;
            emit_y = 0;
            emit_z = 0;
            for(unsigned int i=0; i<n_sample; ++i) {
                emit_x += beta_x*xp_bet[i]*xp_bet[i]+2*alf_x*x_bet[i]*xp_bet[i]+gamma_x*x_bet[i]*x_bet[i];
                emit_y += beta_y*yp_bet[i]*yp_bet[i]+2*alf_y*y_bet[i]*yp_bet[i]+gamma_y*y_bet[i]*y_bet[i];
                emit_z += dp_p[i]*dp_p[i];
            }
            if(ion.bunched()) {
                double bs2 = ring.beta_s();
                bs2 *= bs2;
                double bs2_inv = 1/bs2;
                for(unsigned int i=0; i<n_sample; ++i) emit_z+= ds[i]*ds[i]*bs2_inv;
            }
            emit_x /= 2*n_sample;
            emit_y /= 2*n_sample;
            emit_z /= n_sample;
            break;
        }
        case IonSample::MONTE_CARLO: {
            emit_x = emit(x_bet.get(), xp_bet.get(), n_sample);
            emit_y = emit(y_bet.get(), yp_bet.get(), n_sample);
            if(ion.bunched()) emit_z = emit_p(dp_p.get(), ds.get(), ring, n_sample);
            else emit_z = emit_p(dp_p.get(), n_sample);
            break;
        }
        default: {
            perror("Error in defining ion beam model!");
            break;
        }
    }
    return 0;
}


int new_emittance(EcoolRateParas ecool_paras, Beam ion, Cooler &cooler, double &emit_x, double &emit_y,
                  double &emit_z) {
    double dx = cooler.disp_h();
    double dy = cooler.disp_v();
    double dpx = cooler.der_disp_h();
    double dpy = cooler.der_disp_v();

    unsigned int n_sample = ecool_paras.n_sample();
    for(unsigned int i=0; i<n_sample; ++i) {
        x_bet[i] = x[i] - dx*dp_p[i];
        xp_bet[i] = xp[i] - dpx*dp_p[i];
        y_bet[i] = y[i] - dy*dp_p[i];
        yp_bet[i] = yp[i] - dpy*dp_p[i];
    }

    switch (ecool_paras.ion_sample()) {
        case IonSample::SINGLE_PARTICLE: {
            double alf_x = cooler.alpha_h();
            double alf_y = cooler.alpha_v();
            double beta_x = cooler.beta_h();
            double beta_y = cooler.beta_v();
            double gamma_x = (1+alf_x*alf_x)/beta_x;
            double gamma_y = (1+alf_y*alf_y)/beta_y;
            emit_x = 0;
            emit_y = 0;
            emit_z = 0;
            for(unsigned int i=0; i<n_sample; ++i) {
                emit_x += beta_x*xp_bet[i]*xp_bet[i]+2*alf_x*x_bet[i]*xp_bet[i]+gamma_x*x_bet[i]*x_bet[i];
                emit_y += beta_y*yp_bet[i]*yp_bet[i]+2*alf_y*y_bet[i]*yp_bet[i]+gamma_y*y_bet[i]*y_bet[i];
                emit_z += dp_p[i]*dp_p[i];
            }
//            if(ion.bunched()) {
//                double bs2 = ring.beta_s();
//                bs2 *= bs2;
//                double bs2_inv = 1/bs2;
//                for(unsigned int i=0; i<n_sample; ++i) emit_z+= ds[i]*ds[i]*bs2_inv;
//            }
            emit_x /= 2*n_sample;
            emit_y /= 2*n_sample;
            emit_z /= n_sample;
            break;
        }
        case IonSample::MONTE_CARLO: {
            emit_x = emit(x_bet.get(), xp_bet.get(), n_sample);
            emit_y = emit(y_bet.get(), yp_bet.get(), n_sample);
            emit_z = emit_p(dp_p.get(), n_sample);
            break;
        }
        default: {
            perror("Error in defining ion beam model!");
            break;
        }
    }
    return 0;
}

//For bunched electron beam to cool coasting ion beam, adjust the cooling rate wit the duty factor.
int adjust_rate(EcoolRateParas &ecool_paras, Beam &ion, Ring &ring, Cooler &cooler, EBeam &ebeam, double &rate_x,
         double &rate_y, double &rate_s) {
    if(ebeam.bunched()&&(!ion.bunched())) {
        double sample_length = ebeam.length();
        double bunch_separate = ecool_paras.bunch_separate();
        if(sample_length<0) perror("electron bunch length must be positive!");
        if(bunch_separate>=sample_length) {
            double duty_factor = sample_length/bunch_separate;
            rate_x *= duty_factor;
            rate_y *= duty_factor;
            rate_s *= duty_factor;
        }
        else {
            perror("Electron bunch length is larger than the distance between electron bunches");
        }
    }
    return 0;
}
int bunched_to_coasting(EcoolRateParas &ecool_paras, Beam &ion, EBeam &ebeam, Cooler &cooler,
                        ForceParas &force_paras){
    unsigned int n_sample = ecool_paras.n_sample();
    int count = 1;
    double *force_tr_rcd = new double[n_sample]();
    double *force_long_rcd = new double[n_sample]();
//    std::copy(force_x, force_x+n_sample, force_tr_rcd);
//    std::copy(force_z, force_z+n_sample, force_long_rcd);
    double cz_rcd = ion.center(2);

    ecool_paras.set_shift(true);
    int n_long = ecool_paras.n_long_sample();
    double length = ebeam.length();
    double step = length/n_long;
    double gamma_e_inv = 1./ebeam.gamma();
    for(double cz = cz_rcd-0.5*length; cz <= cz_rcd+0.5*length; cz += step) {
        ion.set_center(2,cz);
        electron_density(ecool_paras, ion, ebeam);
        for(unsigned int i=0; i<n_sample; ++i) ne[i] *= gamma_e_inv;
        force(n_sample, ion, ebeam, cooler, force_paras);
        for(unsigned int i=0; i<n_sample; ++i) {
            force_tr_rcd[i] += force_x[i];
            force_long_rcd[i] += force_z[i];
        }
        ++count;
    }
//    for(double cz = cz_rcd-0.5*length; cz <= cz_rcd+0.5*length; cz += step) {
//        if(cz!=0) {
//            ion.set_center(2,cz);
//            electron_density(ecool_paras, ion, ebeam);
//            for(unsigned int i=0; i<n_sample; ++i) ne[i] *= gamma_e_inv;
//            force(n_sample, ion, ebeam, cooler, force_paras);
//            for(unsigned int i=0; i<n_sample; ++i) {
//                force_tr_rcd[i] += force_x[i];
//                force_long_rcd[i] += force_z[i];
//            }
//            ++count;
//        }
//    }
    ion.set_center(2, cz_rcd);
    ecool_paras.set_shift(false);
    double count_inv = 1.0/static_cast<double>(count);
    for(unsigned int i=0; i<n_sample; ++i) {
        force_x[i] = force_tr_rcd[i]*count_inv;
        force_z[i] = force_long_rcd[i]*count_inv;
    }
    delete[] force_tr_rcd;
    delete[] force_long_rcd;
    return 0;
}


int ecooling_rate(EcoolRateParas &ecool_paras, ForceParas &force_paras, Beam &ion, Cooler &cooler, EBeam &ebeam,
                  Ring &ring, double &rate_x, double &rate_y, double &rate_s) {
    //Initialize the scratch variables
//    if(rms_dynamic_count<1) config_ecooling(ecool_paras, ion);
    if(!dynamic_flag) config_ecooling(ecool_paras, ion);
    //Create the ion samples
//    if(model_beam_count<0) ion_sample(ecool_paras, ion, ring, cooler, ebeam);
//    if(model_beam_count<0) ion_sample(ecool_paras, ion, ring, cooler);
    if(!dynamic_flag) ion_sample(ecool_paras, ion, ring, cooler);
    //Calculate the electron density for each ion
    electron_density(ecool_paras, ion, ebeam);

    unsigned int n_sample = ecool_paras.n_sample();
    //Phase space variables to dynamic variables
    space_to_dynamic(n_sample, ion);
    //Time through the cooler
    t_cooler = cooler.length()/(ion.beta()*k_c);
    //Transfer into e- beam frame
//    beam_frame(ecool_paras.n_sample(), ebeam, ion);
    beam_frame(n_sample, ebeam.gamma());
    //Calculate friction force
    force(n_sample, ion, ebeam, cooler, force_paras);
    //Restore the longitudinal velocity if it has been changed due to electron velocity gradient
    restore_velocity(n_sample, ebeam);

    //Special treatment for bunched electron beam to cool coasting ion beam
    if(!ion.bunched()&&ebeam.bunched()) {
        //bunched_to_coasting is calling force() within it (again?)
        std::cout<<"Bunched to coasting!"<<std::endl;
        bunched_to_coasting(ecool_paras, ion, ebeam, cooler, force_paras);
        //TODO: Should this really be calculating force for a 3rd time?
        force(n_sample, ion, ebeam, cooler, force_paras);
    }

    //Transfer back to lab frame
//    lab_frame(ecool_paras.n_sample(), ebeam, ion);
    lab_frame(n_sample, ebeam.gamma());
    //Distribute the friction force into x,y direction.
    force_distribute(n_sample,ion);
    //Original emittance
    double emit_x0, emit_y0, emit_z0;

    //Archive a copy of the parameters
    double *xp_archive = new double[n_sample];
    double *yp_archive = new double[n_sample];
    double *dp_p_archive = new double[n_sample];
    for(unsigned int i=0;i<n_sample;i++){
        xp_archive[i] = xp[i];
        yp_archive[i] = yp[i];
        dp_p_archive[i] = dp_p[i];
    }

//    original_emittance(ecool_paras, ion, emit_x0, emit_y0, emit_z0);
    original_emittance(ecool_paras, ring, ion, emit_x0, emit_y0, emit_z0);
    //Apply kick
    apply_kick(n_sample, ion); //Note: This mutates the global xp,yp,dp_p values
    //New emittance
    double emit_x, emit_y, emit_z;
//    new_emittance(ecool_paras, ion, cooler, emit_x, emit_y, emit_z);
    new_emittance(ecool_paras, ion, ring, cooler, emit_x, emit_y, emit_z);
    //rate
    rate_x = emit_x/emit_x0-1;
    rate_y = emit_y/emit_y0-1;
    rate_s = emit_z/emit_z0-1;
    double freq = k_c*ion.beta()/ring.circ()*cooler.section_number();
    rate_x *= freq;
    rate_y *= freq;
    rate_s *= freq;
    adjust_rate(ecool_paras, ion, ring, cooler, ebeam, rate_x, rate_y, rate_s);

    //Replace with archived values so cooling kick is not double-counted
    for(unsigned int i=0;i<n_sample;i++){
        xp[i] = xp_archive[i];
        yp[i] = yp_archive[i];
        dp_p[i] = dp_p_archive[i];
    }

    delete xp_archive;
    delete yp_archive;
    delete dp_p_archive;

    return 0;
}


//Fix V_tr or V_long and calculate the force at different values.
// Using the built-in functions within ecooling to initialize all the
// variables as in an ecooling calculation, but then swapping out the
// velocities to calculate the forces needed for plotting.
int CalculateForce(EcoolRateParas &ecool_paras, ForceParas &force_paras, Beam &ion, Cooler &cooler, EBeam &ebeam,
                  Ring &ring){


    //Initialize the scratch variables based on the beam configuration
    // defined by all the function arguments
    config_ecooling(ecool_paras, ion);

    //Create the ion samples
    ion_sample(ecool_paras, ion, ring, cooler);

    unsigned int n_sample = ecool_paras.n_sample();

    //Calculate the electron density for each ion
    electron_density(ecool_paras, ion, ebeam);

    //Phase space variables to dynamic variables
    space_to_dynamic(n_sample, ion);

    //Time through the cooler
    t_cooler = cooler.length()/(ion.beta()*k_c);

    //Transfer into e- beam frame
    beam_frame(n_sample, ebeam.gamma());

    //In the beam frame, determine the max/min velocities in v_long
    // and v_tr in the specified beam configuration
    double v_long_max, v_tr_max, ne_max = - DBL_MAX;
    double v_long_min, v_tr_min, ne_min = DBL_MAX;

    for(int i=0;i<n_sample;i++){
        if(v_long[i] > v_long_max) {
            if(v_long[i] < k_c) v_long_max = v_long[i];
        }
        if(v_tr[i]   > v_tr_max) {
            if(v_tr[i] < k_c)   v_tr_max = v_tr[i];
        }

        if(v_long[i] < v_long_min){
            if(v_long[i] > -k_c) v_long_min = v_long[i];
        }
        if(v_tr[i]   < v_tr_min) {
            if(v_tr[i] > -k_c)   v_tr_min = v_tr[i];
        }

        if(ne[i] > ne_max) ne_max = ne[i];
        if(ne[i] < ne_min) ne_min = ne[i];
    }

    //Protect against some extraneous large values
    // that wreck our v_interval. This was observed
    //only once and wasn't reproducible.
    //Enforce speed limits.
    if(v_tr_min < -k_c) v_tr_min = -k_c;
    if(v_tr_max > k_c) v_tr_max = k_c;

    //TODO: This is a temporary (permenant?) hack to fix the velocities that are calculated
    v_tr_max = 6e5;
    v_tr_min = 0;
    v_long_max = 6e5;
    v_long_min = -6e5;


    //Define the values of the velocities for reporting
    // on the dependence of v_long. We set v_long = 0 for
    // calculating f_tr, and v_tr = 0 for calculating f_long.
    double v_interval = (v_tr_max - v_tr_min)/n_sample;

#ifndef NDEBUG
    std::cout<<"Calculate Force: "<<std::endl;
    std::cout<<"V_tr min = "<<v_tr_min<<" V_tr max = "<<v_tr_max<<std::endl;
    std::cout<<"V_long min = "<<v_long_min<<" V_long max = "<<v_long_max<<std::endl;
    std::cout<<"ne min = "<<ne_min<<" ne max = "<<ne_max<<std::endl;
    std::cout<<"n_sample = "<<n_sample<<" v trans interval = "<<v_interval<<std::endl;
#endif

    //Our v_tr values are already set. Zero out v_long
    for(int i=0;i<n_sample;i++){
        v_tr[i] = v_tr_min + (double)i*v_interval;
        v_long[i] = 0.0;
        //yp[i]  = 0.0 ; // Focus on transverse only in x direction
        ne[i] = ne_max;  //We need a constant n_e for a smooth plot
    }

    //Calculate friction force f_tr when v_long = 0.
    force(n_sample, ion, ebeam, cooler, force_paras);

    //Store the values so they don't get overwritten when we
    // calculate f_long
    std::vector< double > v_tmp;
    std::vector< double > f_tmp;

    for(int i=0;i<n_sample;i++){
        v_tmp.push_back(v_tr[i]);
        f_tmp.push_back(force_x[i]);
    }

    //Prepare for calculating force_z (still in the beam frame)
    v_interval = (v_long_max - v_long_min)/n_sample;
    for(int i=0;i<n_sample;i++){
        v_long[i]  = v_long_min + (double)i*v_interval;
        v_tr[i]    = 0.0;
        force_x[i] = 0.0;
        force_y[i] = 0.0;
        force_z[i] = 0.0;
        ne[i]      = ne_max;  //We need a constant n_e for a smooth plot
    }
#ifndef NDEBUG
    std::cout<<"n_sample = "<<n_sample<<" v long interval = "<<v_interval<<std::endl;
#endif
    //Calculate friction force f_long when v_tr = 0
    force(n_sample, ion, ebeam, cooler, force_paras);
    //This writes force_z and v_long

    //Replace the values calculated previously. Sort them by velocity value
    std::vector< std::pair<double,double> > zip;
    for(size_t i=0; i<v_tmp.size(); ++i){
        zip.push_back(std::make_pair(v_tmp[i], f_tmp[i]));
    }

    // Sort the vector of pairs
    std::sort(std::begin(zip), std::end(zip),
              [&](std::pair<double,double> a, std::pair<double,double> b){
                  return a.first < b.first;
              });
    //Over-write the global variables with sorted values
    for(size_t i=0; i<zip.size(); i++){
        v_tmp[i] = zip[i].first;
        f_tmp[i] = zip[i].second;
    }
    //Now write to the global variables
    for(size_t i=0;i<zip.size();i++){
        v_tr[i]    = v_tmp[i];
        force_x[i] = f_tmp[i];
    }

    //We don't want to transfer the velocities back to lab frame
    //lab_frame(n_sample, ebeam.gamma());
    //end_ecooling(ecool_paras, ion);

    //TODO: Betacool shows this as a 2d contour plot. Should we instead calculate
    // on a grid to prepare for this kind of plot in the future?


    return 0;
}
