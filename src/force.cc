#include "force.h"

#include <cstdio>
#include <cassert>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>


//These functions calculate the force from the derrived class, and output the test file
// from a filename stored in the derived class header.
int ForceParas::ApplyForce(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                            double temperature, double magnetic_field, double *d_perp_e, double *d_paral_e, double time_cooler,
                            double *force_tr, double *force_long) {

//Compiler should ignore #pragma it doesn't understand
    #pragma omp parallel for  
    for(unsigned long int i=0; i<ion_number; ++i){
        double result_trans,result_long;
        force(v_tr[i],v_long[i],d_perp_e[i],d_paral_e[i],temperature,charge_number,
              density_e[i],time_cooler,magnetic_field, result_trans, result_long);
        //the force in Newtons
        force_tr[i] = result_trans; 
        force_long[i] = result_long; 
    } 

    return 0;
}

int ForceParas::ApplyForce(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                            double temperature, double magnetic_field, double d_perp_e, double d_paral_e, double time_cooler,
                            double *force_tr, double *force_long, bool do_test) {
    
//Compiler should ignore #pragma it doesn't understand
    #pragma omp parallel for //ordered schedule(dynamic,10)
    for(unsigned long int i=0; i<ion_number; ++i){
        double result_trans,result_long;
        force(v_tr[i],v_long[i],d_perp_e,d_paral_e,temperature,charge_number,
              density_e[i],time_cooler,magnetic_field, result_trans, result_long);

        //the force in newtons
        force_tr[i] = result_trans;
        force_long[i] = result_long;
    }
    
    //Do this separately to ensure that ofstream is threadsafe
    if(do_test){
        //Open up an output file if we're in the testing phase
        std::ofstream outfile;
        outfile.open(test_filename_);
        
        for(unsigned long int i=0;i<ion_number;++i){
            outfile<< f_const_ <<", "<<v_tr[i]<<", "<<v_long[i]
                   <<", "<<density_e[i]<<", "<<force_tr[i]<<", "<<force_long[i]<<"\n";
        }
        outfile.close();
    }
    
    return 0;
}

void ForceParas::EvalIntegral(double (*func)(double*, size_t, void*), int_info &params,double *xl,double *xu,size_t dim,
                                  double &result,double &error){

    //This is a wrapper for a number of triple-integrals that appear in these force calculations
    
    const gsl_rng_type *T;
    gsl_rng *r;
    
    gsl_monte_function G;
    G.f = func;
    G.dim = dim;
    G.params = &params; 

    //Initialize the random number generator, which we'll need regardless of the
    // Monte Carlo scheme.
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    //Hard-code this for now while we're testing
    bool plain_ = true;
    bool miser_ = false;
    bool vegas_ = false;

    //This uses the PLAIN algorithm, which samples a fixed number of points and should be fast
    if(plain_){
        gsl_monte_plain_state *s = gsl_monte_plain_alloc (dim);
        gsl_monte_plain_integrate (&G, xl, xu, dim, calls_, r, s,
                                   &result, &error);
        gsl_monte_plain_free (s);
    }
    
    //This uses the MISER algorithm, which is a balance between speed and accuracy.
    else if(miser_){
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(dim);
        gsl_monte_miser_integrate (&G, xl, xu, dim, calls_, r, s,
                               &result, &error);
        gsl_monte_miser_free (s);
    }
    
    //This uses the VEGAS algorithm for Monte Carlo evaluation of the integral. 
    // This method is a trade-off of speed for accuracy.
    else if(vegas_){
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);
        //Some burn-in for importance-weighting
        gsl_monte_vegas_integrate(&G,xl,xu,dim,calls_/5,r,s, &result,&error);

        //Sample until the result is 'good enough' and no more
        do
        {
            gsl_monte_vegas_integrate (&G, xl, xu, dim, calls_/5, r, s,
                                       &result, &error);
        }
        while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

        gsl_monte_vegas_free (s);
    }
    else{
        assert(false&&"Need to set a MC Integration method!");
    }
    
    gsl_rng_free (r);   
}


double ForceParas::max_impact_factor(double v_dlt, int charge_number,double density_e,double time_cooler){
    
    //double wp_const = 4 * k_pi * k_c*k_c * k_e * k_ke / (k_me*1e6);
    double wp = sqrt(4 * k_pi * k_c*k_c * k_re * density_e);
    double rho_max = v_dlt / wp; //The shelding radius, rho_sh
    double rho_max_2 = pow(3 * charge_number / density_e, 1.0/3);
    if(rho_max<rho_max_2) rho_max = rho_max_2;
    double rho_max_3 = v_dlt * time_cooler;
    if(rho_max>rho_max_3) rho_max = rho_max_3;
    
    return rho_max;
}

//This function is used for calculating the forces on single particles or arrays of particles
void Force_Parkhomchuk::force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,
                              int charge_number, double density_e, double time_cooler, double magnetic_field, 
                              double &force_result_trans, double &force_result_long){

    double v2 = v_tr*v_tr + v_long*v_long;    
    double force = 0.0;
    
    if(v2>0){
        //In SI units, this is m*v_tr /(qB). In CGS units, there's an extra k_c in the numerator.
        double rho_lamor = k_me_kg * d_perp_e / ( magnetic_field * k_e ); //in units of m
        double v2_eff_e = temperature * k_c*k_c / (k_me*1e6);
        
        //Use the effective delta v specific to the Parkhomchuk model
        double dlt2_eff_e = d_paral_e*d_paral_e + v2_eff_e;
        double dlt = v2 + dlt2_eff_e;

        //double rho_min_const = charge_number * k_e * k_ke * k_c*k_c / (k_me*1e6);
        double rho_min_const = charge_number * k_re * k_c*k_c; //SI, in units of m
        double rho_min = rho_min_const/dlt;
        dlt = sqrt(dlt);

        double rho_max = max_impact_factor(dlt,charge_number,density_e,time_cooler); //dlt already sqrt

        double lc = log( ( rho_max + rho_min + rho_lamor ) / ( rho_min + rho_lamor ) );   //Coulomb Logarithm =~6

        //Calculate friction force
        force = f_const_ * charge_number*charge_number * density_e * lc / ( dlt*dlt*dlt );
    }

    //This factor of the force must be multiplied by either v_tr or v_long to truly be
    // equal to the force in newtons
    force_result_trans = force * v_tr;
    force_result_long = force * v_long;

}    

double Force_DS::trans_integrand(double alpha,void *params){
    //Using Pestrikov's trick to change the integration variable

    int_info *p = (struct int_info *)params;    
    double y = p->V_long / p->width;
    double z = p->V_trans / p->width;
    double x = y + abs(z)*tan(alpha);
    
    double integrand = tan(alpha)*(y*cos(alpha) + abs(z)*sin(alpha));
    //The F_const already has a minus sign, so the integrand term
    // should carry the sign of z
    integrand *= exp(-0.5*pow(y + abs(z)*tan(alpha),2));
    integrand *= copysign(1,z);
    
    return integrand;    
}

double Force_DS::long_integrand(double alpha, void *params){
    //Using Pestrikov's trick to change the integration variable

    struct int_info *p = (struct int_info *)params;    
    
    double y = p->V_long / p->width;
    double z = p->V_trans / p->width;
    double x = y + abs(z)*tan(alpha);
    
    double integrand = y*cos(alpha) + abs(z)*sin(alpha);
    integrand *= exp(-0.5*pow(y+abs(z)*tan(alpha),2));
    
   return integrand; 
}

void Force_DS::EvalIntegral(double (*func)(double, void*), int_info &params,double &result,double &error){
    unsigned int space_size = 100;
    gsl_set_error_handler_off();
    gsl_function F;

    F.params = &params;
    F.function = func;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(space_size);
    double result_trans,result_long, error_trans,error_long;

    //In D&S/Pestrikov, all integrals have the same limits, so we can hard-code it here
    int status = gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-6,space_size,1,w,&result,&error);
  
    if(status == GSL_EDIVERGE){ //If we didn't get divergene quickly, try again in a more time-consuming way
        status = gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-10,space_size,1,w,&result,&error);
    }

    gsl_integration_workspace_free(w);
}

//This function is used for calculating the forces on single particles or arrays of particles
void Force_DS::force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,
                     int charge_number,double density_e,double time_cooler,double magnetic_field, 
                     double &force_result_tr, double &force_result_long){
    
  double rho_L      = k_me_kg * d_perp_e / ( magnetic_field * k_e ); //SI units, in meters
    
  //Calculate rho_max as in the Parkhomchuk model
  double v_mag      = sqrt(v_tr*v_tr + v_long*v_long);
  double rho_max = max_impact_factor(v_mag,charge_number,density_e,time_cooler);
  double lm = log(rho_max / rho_L); //Coulomb Log ~ 5

  //Save ourselves some cycles
  if(lm < 0. || density_e <= 0.0){
      //lm = 0.0;
      force_result_long = 0.0;
      force_result_tr = 0.0;
      return;
  }
    
  int_info params;
  params.V_trans = v_tr; //Ion velocity components
  params.V_long  = v_long;
  params.width = d_paral_e; //The electron bunch width RMS
      
  double result_trans,result_long, error_trans,error_long;

  //Solve indefinite integrals with GSL
  EvalIntegral(&Force_DS::trans_integrand,params,result_trans,error_trans);
  EvalIntegral(&Force_DS::long_integrand,params,result_long,error_long);
     
  force_result_tr = f_const_ * density_e * charge_number*charge_number * lm * result_trans;
  force_result_tr /= (d_paral_e*d_paral_e); 

  force_result_long = f_const_ * density_e * charge_number*charge_number * lm * result_long;
  force_result_long /= (d_paral_e*d_paral_e);        
}

void Force_Meshkov::force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,
                          int charge_number,double density_e,double time_cooler,double magnetic_field,
                          double &force_result_trans, double &force_result_long){
    
  double rho_L    = k_me_kg * d_perp_e / ( magnetic_field * k_e ); //SI units, in m
  double wp_const = 4 * k_pi * k_c*k_c * k_re;
  double wp       = sqrt(wp_const * density_e);
    
  double rho_min_const = charge_number * k_re * k_c*k_c; //SI units
    
  double v2      = v_tr*v_tr + v_long*v_long;
  double v3      = pow(sqrt(v2),3);

  //The Number of multiple adiabatic collisions with a single electron
  double N_col = 1 + d_perp_e/(k_pi*sqrt(v2 + d_paral_e*d_paral_e));
      
  //dynamic shielding radius
  double rho_sh = sqrt(v2 + d_paral_e*d_paral_e) / wp;
  //minimum impact parameter  
  double rho_min = rho_min_const / (v2 + d_paral_e*d_paral_e);  //in units of m
  //intermediate impact parameter
  double rho_F = rho_L * sqrt(v2 + d_paral_e*d_paral_e) / d_perp_e;
  
  double rho_max = max_impact_factor(sqrt(v2 + d_paral_e*d_paral_e),charge_number,density_e,time_cooler);

  //k_ is a fudge factor to smooth the friction force shape. Using the
  // "Classical" default definition here from the BETACOOL documentation

  //Coulomb Logarithms
  double L_M = log(rho_max / (k_*rho_L)); //there's a typo in the BETACOOL guide, the numerator here is rho_max
  double L_A = log((k_*rho_L) / rho_F);
  double L_F = log(rho_F / rho_min);

  //std::cout<<"Meshkov: L_M="<<L_M<<" L_A="<<L_A<<" L_F="<<L_F<<std::endl;

  //If it's < 0, the interaction is absent
  if( L_M < 0. ) L_M = 0.0; //This controls the low ion V_long rise, ~10-20
  if( L_A < 0. ) L_A = 0.0; 
  if( L_F < 0. ) L_F = 0.0;  //This one is responsible for discontinuity at high ion V_long ~5.5
    
  //Define the regions of the ion velocity domains
  double result_trans,result_long;
    
  double ellipse = (v_tr*v_tr)/(d_perp_e*d_perp_e) + (v_long*v_long)/(d_paral_e*d_paral_e);
     
  if(sqrt(v2) > d_perp_e) { //Region I
       result_trans = 2*L_F + L_M*(v_tr*v_tr - 2*v_long*v_long)/v2;
       result_trans /= v3;
         
        result_long = 2 + 2*L_F + L_M*(3*v_tr*v_tr/v2);
        result_long /= v3;
  }   
  else if(sqrt(v2) < d_paral_e) { //Region III
        result_trans = 2*(L_F + N_col*L_A)/pow(d_perp_e,3) + L_M/pow(d_paral_e,3);
        result_long = 2*(L_F + N_col*L_A)/(d_perp_e*d_perp_e * d_paral_e) + L_M/pow(d_paral_e,3);
      }  

  //Region II (a or b)
  else{ //This constrains to a donut shape (d_paral_e < sqrt(v2) < d_perp_e).
        result_trans = ((v_tr*v_tr - 2*v_long*v_long)/v2) * (L_M/v3);
        result_trans += (2/pow(d_perp_e,3)) * (L_F + N_col*L_A);
          
        if( ellipse <= 1.0 ){ // Region IIb
              //This is the same as result_long in Region 3 
              result_long = 2*(L_F + N_col*L_A)/(d_perp_e*d_perp_e * d_paral_e) + L_M/pow(d_paral_e,3);
        }
        else{ //It must be Region IIa
              result_long = (((3*v_tr*v_tr*L_M)/v2) + 2) / v3;
              result_long += 2 * (L_F + N_col*L_A)/(d_perp_e*d_perp_e * v_long);
        }
  } 
        
  force_result_trans = f_const_ * charge_number*charge_number * density_e * result_trans * v_tr;
  force_result_long = f_const_ * charge_number*charge_number * density_e * result_long * v_long;
}

double Force_Unmagnetized::normalization_factor(double *k, size_t dim, void *params){

    (void) (dim); //avoid unused parameter warnings
    
    double v_t = k[0]; //= integrating variable v_trans (non-negative)
    double v_l = k[1];  //= integrating variable v_long (all values)
    double phi = k[2];   //= integrating angle phi (radians)
    struct int_info *p = (struct int_info *)params;    
    double d_perp_e = p->d_perp_e;
    double d_paral_e = p->d_paral_e;
    
    double U = exp((-v_t*v_t)/(2*d_perp_e*d_perp_e)-(v_l*v_l)/(2*d_paral_e*d_paral_e));
    U *= v_t;
   
    return U;
}

double Force_Unmagnetized::trans_integrand(double *k, size_t dim, void *params){
   
    (void) (dim); //avoid unused parameter warnings
    
    double v_t = k[0]; //= integrating variable v_trans (non-negative)
    double v_l = k[1];  //= integrating variable v_long (all values)
    double phi = k[2];   //= integrating angle phi (radians)
    struct int_info *p = (struct int_info *)params;    
    double y = p->V_long - v_l;
    double z = p->V_trans - v_t*cos(phi);
    double d_perp_e = p->d_perp_e;
    double d_paral_e = p->d_paral_e;
    double rho_max = p->rho_max;
    double rho_min_const = p->Z * k_re * k_c*k_c; //SI units, m
    double rho_min = rho_min_const / ( y*y + z*z + pow(v_t*sin(phi),2));
    
    double L_C = log(rho_max/rho_min);
    if(L_C<0.0) L_C = 0.0;
    
    double U = L_C * z;
    U *= exp(-((v_t*v_t)/(2*d_perp_e*d_perp_e)) - ((v_l*v_l)/(2*d_paral_e*d_paral_e)));
    U /= pow(y*y + z*z + pow(v_t*sin(phi),2),3/2);
    U *= v_t;
    return U;

}

double Force_Unmagnetized::long_integrand(double *k, size_t dim, void *params){
   
    (void) (dim); //avoid unused parameter warnings
    
    double v_t = k[0]; //= integrating variable v_trans (non-negative)
    double v_l = k[1];  //= integrating variable v_long (all values)
    double phi = k[2];   //= integrating angle phi (radians)
    struct int_info *p = (struct int_info *)params;    
    double y = p->V_long - v_l;
    double z = p->V_trans - v_t*cos(phi);
    double d_perp_e = p->d_perp_e;
    double d_paral_e = p->d_paral_e;
    double rho_max = p->rho_max;
    double rho_min_const = p->Z * k_re * k_c*k_c; //SI units, m
    double rho_min = rho_min_const / ( y*y + z*z + pow(v_t*sin(phi),2));

    double L_C = log(rho_max/rho_min);
    if(L_C<0.0) L_C = 0.0;
    
    double U = L_C * y;
    U *= exp(-((v_t*v_t)/(2*d_perp_e*d_perp_e)) - ((v_l*v_l)/(2*d_paral_e*d_paral_e)));
    U /= pow(y*y + z*z + pow(v_t*sin(phi),2),3/2);
    U *= v_t;
    return U;

}


void Force_Unmagnetized::force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,
                         int charge_number,double density_e,double time_cooler, double magnetic_field, 
                         double &force_result_trans, double &force_result_long){

    double v2 = v_tr*v_tr + v_long*v_long;
    double v_mag = sqrt(v2);

    //Save some CPU cycles
    if(density_e>0.0 && v2 >0.0){
        double rho_max = max_impact_factor(v_mag,charge_number,density_e,time_cooler);
        
        if(approximate_){
            //this gives the value to within 2%
            
            //First, calculate the normalization factor.
            int_info params;
            params.V_trans = v_tr; //Ion velocity components
            params.V_long  = v_long;
            params.d_paral_e = d_paral_e; //The electron bunch width RMS
            params.d_perp_e = d_perp_e;
            params.Z = charge_number;
            params.rho_max = rho_max;
            
            //Lower limits to the integrals 
            double xl[3] = {0.,-3*d_paral_e,0.}; //in order: v_trans,v_long,phi
            //Upper limits to the variables
            double xu[3] = {3*d_perp_e,3*d_paral_e,k_pi};
            double norm_factor,error_norm;
            
            EvalIntegral(&normalization_factor,params,xl,xu,(size_t) 3,norm_factor,error_norm);

            //Now calculate the actual integrals (limits are the same as above)           
            double trans_result,trans_error;
            double long_result,long_error;
            
            EvalIntegral(&trans_integrand,params,xl,xu,(size_t) 3,trans_result,trans_error);
            EvalIntegral(&long_integrand,params,xl,xu,(size_t) 3,long_result,long_error);

            double constants = f_const_ * charge_number*charge_number * density_e;
            //std::cout<<" norm="<<norm_factor<<" trans_result="<<trans_result<<" long_result="<<long_result<<std::endl;
            force_result_trans = constants * trans_result / norm_factor;
            force_result_long  = constants * long_result / norm_factor;
        }
        else{
            //Full integrals, will take longer
            int_info params;
            params.V_trans = v_tr; //Ion velocity components
            params.V_long  = v_long;
            params.d_paral_e = d_paral_e; //The electron bunch width RMS
            params.d_perp_e = d_perp_e;
            params.Z = charge_number;
            params.rho_max = rho_max;
            
            //Lower limits to the integrals 
            double xl[3] = {0.,-cutoff_,0.}; //in order: v_trans,v_long,phi
            //Upper limits to the variables
            double xu[3] = {cutoff_,cutoff_,2*k_pi};
            double norm_factor,error_norm;
            
            //Now calculate the actual integrals (limits are the same)           
            double trans_result,trans_error;
            double long_result,long_error;
            
            EvalIntegral(&trans_integrand,params,xl,xu,(size_t) 3,trans_result,trans_error);
            EvalIntegral(&long_integrand,params,xl,xu,(size_t) 3,long_result,long_error);

            double constants = f_const_ * charge_number*charge_number * density_e;
            force_result_trans = constants * trans_result / (d_perp_e*d_perp_e * d_paral_e);
            force_result_long  = constants * long_result / (d_perp_e*d_perp_e * d_paral_e);           
        }
    }
    else{
        force_result_trans = 0.0;
        force_result_long = 0.0;
    }
}

void Force_Budker::force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,
                         int charge_number,double density_e,double time_cooler, double magnetic_field, 
                         double &force_result_trans, double &force_result_long){
    
    //Here we are assuming a maxwellian distribution, d_paral_e = d_perp_e
    // and ignoring magnetic field
        
    double delta_e = d_perp_e; 
    double v2 = v_tr*v_tr + v_long*v_long;
    double v_mag = sqrt(v2);

    //Save some CPU cycles
    if(density_e>0.0 && v2 >0.0){
        double rho_max = max_impact_factor(v_mag,charge_number,density_e,time_cooler);
        double rho_min_const = charge_number * k_re * k_c*k_c; //Gives LC ~ 9. SI units, in m
        double rho_min = rho_min_const / (v2 + delta_e*delta_e);

        double lc = log( rho_max / rho_min );
        if(lc<0.) lc = 0.;
    
        double arg = v_mag / delta_e;
        //double phi = sqrt( 2 / k_pi ) * ( erf(arg/sqrt(2)) - arg * exp((-arg*arg)/2) );
        double phi = sqrt( 2 / k_pi ) * ( erf(arg) - arg * exp((-arg*arg)/2) );
        double result = phi * pow(v_mag,-3);
 
        force_result_trans = f_const_ * charge_number*charge_number * lc * density_e * result * v_tr;
        force_result_long = f_const_ * charge_number*charge_number * lc * density_e * result * v_long;
 
    }
    else{
        force_result_trans = 0.0;
        force_result_long = 0.0;
    }
}

double Force_Erlangen::fast_trans(double *k, size_t dim, void *params){

    (void) (dim); //avoid unused parameter warnings
    
    double v_t = k[0]; //= integrating variable v_trans (non-negative)
    double v_l = k[1];  //= integrating variable v_long (all values)
    double phi = k[2];   //= integrating angle phi (radians)
    double L_C;
    struct int_info *p = (struct int_info *)params;    
    double y = p->V_long - v_l;
    double z = p->V_trans - v_t*cos(phi);
    double d_perp_e = p->d_perp_e;
    double d_paral_e = p->d_paral_e;
   
    //The Coulomb log
    double alt_rho_max = k_me_kg * v_t / (k_e * p->B_field); //SI units, m
    double rho_min_const = p->Z * k_re * k_c*k_c; //SI units, m
//    double rho_min_const = p->Z * k_e*k_e * k_ke / k_me; //CGS units
    double rho_min = rho_min_const / ( y*y + z*z + pow(v_t*sin(phi),2));
    
    //for some regions of integration, use the smaller rho_max
    if( alt_rho_max < p->rho_max ) L_C = log( alt_rho_max / rho_min );
    else L_C = log( p->rho_max / rho_min );
    if(L_C < 0.0) L_C = 0.0;
    
    double U = z * v_t * L_C * exp(-((v_t*v_t)/(2*d_perp_e*d_perp_e)) - ((v_l*v_l)/(2*d_paral_e*d_paral_e)));
    U /= pow(y*y + z*z + pow(v_t*sin(phi),2),3/2);
   
    return U;
}

double Force_Erlangen::fast_long(double *k, size_t dim, void *params){
    
    (void) (dim); //avoid unused parameter warnings
    
    double v_t = k[0]; //= integrating variable v_trans (non-negative)
    double v_l = k[1];  //= integrating variable v_long (all values)
    double phi = k[2];   //= integrating angle phi (radians)
    double L_C; 
    struct int_info *p = (struct int_info *)params;  
    double y = p->V_long - v_l;
    double z = p->V_trans - v_t*cos(phi);
    double d_perp_e = p->d_perp_e;
    double d_paral_e = p->d_paral_e;
    
    //The Coulomb log
    double alt_rho_max = k_me_kg * v_t / (k_e * p->B_field); //SI units, m
    double rho_min_const = p->Z * k_re * k_c*k_c; //SI units, m
//    double rho_min_const = p->Z * k_e*k_e * k_ke / k_me;
    double rho_min = rho_min_const / ( y*y + z*z + pow(v_t*sin(phi),2));
    
    //for some regions of integration, use the smaller rho_max
    if( alt_rho_max < p->rho_max ) L_C = log( alt_rho_max / rho_min );
    else L_C = log( p->rho_max / rho_min );
    if(L_C < 0.0) L_C = 0.0;
    
    double U = y * v_t * L_C * exp(-((v_t*v_t)/(2*d_perp_e*d_perp_e)) - ((v_l*v_l)/(2*d_paral_e*d_paral_e)));
    U /= pow(y*y + z*z + pow(v_t*sin(phi),2),3/2);
    
    return U;
}

double Force_Erlangen::stretched_trans(double *k, size_t dim, void *params){

    (void) (dim); //avoid unused parameter warnings
        
    double v_t = k[0]; //= integrating variable v_trans
    double v_l = k[1];  //= integrating variable v_long
    double L_C;
    struct int_info *p = (struct int_info *)params;     
    double y = p->V_long - v_l;
    double d_perp_e = p->d_perp_e;
    double d_paral_e = p->d_paral_e;
    double vt2 = p->V_trans*p->V_trans;
    
    //The Coulomb log
    double delta = k_me_kg * sqrt(vt2 + y*y) / (k_e * p->B_field);
    double alt_rho =  k_me_kg * v_t / (k_e * p->B_field);

    double numerator = p->rho_max;
    double denominator = p->rho_max;
    if( alt_rho < p->rho_max ) denominator = alt_rho;
    if( delta < p->rho_max ) numerator = delta;
    if( numerator < denominator ) return 0.0;
    else L_C = log( numerator / denominator);
    
    double U = v_t * p->V_trans * pow(vt2 + y*y,-3/2) * exp(-v_l*v_l/(2*d_paral_e*d_paral_e));
    U *= L_C * exp(-v_t*v_t/(2*d_perp_e*d_perp_e));
    
    return U;
}

double Force_Erlangen::stretched_long(double *k, size_t dim, void *params){
    
    (void) (dim); //avoid unused parameter warnings
    
    double v_t = k[0]; //= integrating variable v_trans 
    double v_l = k[1]; //= integrating variable v_long
    double L_C;
    struct int_info *p = (struct int_info *)params;  
    double y = p->V_long - v_l;
    double d_perp_e = p->d_perp_e;
    double d_paral_e = p->d_paral_e;
    double vt2 = p->V_trans*p->V_trans;
    
        //The Coulomb log
    double delta = k_me_kg * sqrt(vt2 + y*y) / (k_e * p->B_field);
    double alt_rho =  k_me_kg * v_t / (k_e * p->B_field);

    double numerator = p->rho_max;
    double denominator = p->rho_max;
    if( alt_rho < p->rho_max ) denominator = alt_rho;
    if( delta < p->rho_max ) numerator = delta;
    if( numerator < denominator ) return 0.0;
    else L_C = log( numerator / denominator);
    
    double U = v_t * y * pow(vt2 + y*y,-3/2) * exp(-v_l*v_l/(2*d_paral_e*d_paral_e));
    U *= L_C * exp(-v_t*v_t/(2*d_perp_e*d_perp_e));
    
    return U;
}

double Force_Erlangen::tight_trans(double *k, size_t dim, void *params){

    (void) (dim); //avoid unused parameter warnings
    
    double v_t = k[0]; //= integrating variable v_trans
    double v_l = k[1];  //= integrating variable v_long
   
    struct int_info *p = (struct int_info *)params;     
    double y = p->V_long - v_l;
    double d_perp_e = p->d_perp_e;
    double d_paral_e = p->d_paral_e;
    double vt2 = p->V_trans*p->V_trans;
    
    //The Coulomb log
    double L_C;
    double delta = k_me_kg * sqrt(vt2 + y*y) / (k_e * p->B_field);
    double rho_min = k_me_kg * v_t / (k_e * p->B_field);
    
    
    if( rho_min < delta ) rho_min = delta;
    if( p->rho_max < rho_min ) return 0.0;
    else L_C = log(p->rho_max / rho_min);
    
    double U = p->V_trans * (vt2 - y*y) * 0.5 * pow(vt2 + y*y,-5/2) * exp(-v_l*v_l/(2*d_paral_e*d_paral_e));
    U *= v_t * L_C * exp(-v_t*v_t/(2*d_perp_e*d_perp_e));
    
    return U;
}

double Force_Erlangen::tight_long(double *k, size_t dim, void *params){
    
    (void) (dim); //avoid unused parameter warnings
        
    double v_t = k[0]; //=integrating variable v_trans
    double v_l = k[1];  //=integrating variable v_long
    
    struct int_info *p = (struct int_info *)params;
    double y = p->V_long - v_l;
    double d_perp_e = p->d_perp_e;
    double d_paral_e = p->d_paral_e;
    double vt2 = p->V_trans*p->V_trans;
    
    //The Coulomb log
    double L_C;
    double delta = k_me_kg * sqrt(vt2 + y*y) / (k_e * p->B_field);
    double rho_min = k_me_kg * v_t / (k_e * p->B_field);  
    
    if( rho_min < delta ) rho_min = delta;
    if( p->rho_max < rho_min ) return 0.0;
    else L_C = log( p->rho_max / rho_min);
    
    double U = v_t * vt2 * y * pow(vt2 + y*y,-5/2) * exp(-v_l*v_l/(2*d_paral_e*d_paral_e));
    U *= L_C * exp(-v_t*v_t/(2*d_perp_e*d_perp_e));

    return U;
}



void Force_Erlangen::force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,
                           int charge_number,double density_e,double time_cooler, double magnetic_field, 
                           double &force_result_trans, double &force_result_long){

    double v_mag = sqrt(v_tr*v_tr + v_long*v_long);

    double final_result_trans,final_result_long = 0.0;
    double result_trans,result_long, error_trans,error_long;

    if(v_mag > 0. && density_e > 0.){
        int_info params;
        params.V_trans = v_tr; //Ion velocity components
        params.V_long  = v_long;
        params.d_paral_e = d_paral_e; //The electron bunch width RMS
        params.d_perp_e = d_perp_e;
        params.B_field = magnetic_field;
        params.Z = charge_number;
        
        double temp = k_N_eVm * f_const_ * charge_number*charge_number * density_e / (d_perp_e*d_perp_e * d_paral_e);

        double v_eff = sqrt(v_tr*v_tr + v_long*v_long + d_paral_e*d_paral_e);
        double rho_max = max_impact_factor(v_eff,charge_number,density_e,time_cooler);
        params.rho_max = rho_max;
        
        if(fast_){
            //Lower limits to the integrals 
            double xl[3] = {0.,-cutoff_,0.}; //in order: v_trans,v_long,phi
            //Upper limits to the variables
            double xu[3] = {cutoff_,cutoff_,2*k_pi};

            EvalIntegral(&fast_trans,params,xl,xu,(size_t) 3,result_trans,error_trans);
            EvalIntegral(&fast_long, params,xl,xu,(size_t) 3,result_long, error_long);
            //std::cout<<"Fast trans = "<<result_trans * temp <<" long = "<<result_long * temp<<std::endl;
                      
            final_result_trans += sqrt(2/k_pi) * result_trans;
            final_result_long  += sqrt(2/k_pi) * result_long;
        }
        
        if(stretched_){
            //Integration limits
            double xl[2] = {0.,-cutoff_}; //in order: v_trans,v_long
            double xu[2] = {cutoff_,cutoff_};

            //There is an approximation when V>>d_paral_e

            //A guess at the limit of approximation
            if(v_mag >(50*d_paral_e)){

                //Calculate the coulomb log
                double delta = k_me_kg * v_mag / (k_e * magnetic_field);
                double w_p = sqrt(4*k_pi * density_e * k_e*k_e / k_me_kg);  //plasma frequency
                double w_b = k_e * magnetic_field/k_me_kg; //cyclotron frequency    
                
                //Several factors cancel out in the coulomb log
                double L_C = log( rho_max / std::max(v_tr,delta) ) + log(w_p/w_b);              
                
                result_trans = v_tr * L_C * pow(v_mag*v_mag + d_paral_e*d_paral_e,-3/2);
                result_long = v_long * L_C * pow(v_mag*v_mag + d_paral_e*d_paral_e,-3/2);
          
            }
            else{
                EvalIntegral(&stretched_trans,params,xl,xu,(size_t) 2,result_trans,error_trans);
                EvalIntegral(&stretched_long, params,xl,xu,(size_t) 2,result_long, error_long);
                //std::cout<<"stretched trans = "<<result_trans * temp<<" long = "<<result_long * temp<<std::endl;
                
            }
            final_result_trans += 4 * k_pi * result_trans / sqrt(2*k_pi);
            final_result_long  += 4 * k_pi * result_long / sqrt(2*k_pi);
        }
        
        if(tight_){
            //Integration limits
            double xl[2] = {0.,-cutoff_};
            double xu[2] = {cutoff_,cutoff_};
            double delta = k_me_kg * v_mag / (k_e * magnetic_field);
    
            //An approximation for speed
            if(v_mag >(500*d_paral_e)){
                 //Calculate the coulomb log
                double delta = k_me_kg * v_mag / (k_e * magnetic_field);
                //Several factors cancel out in the coulomb log
                double L_C = log( v_tr / std::max(v_tr,delta) );

                result_trans = v_long * L_C * v_tr*v_tr / pow(v_mag,5);
                result_long = v_long * L_C * (v_long*v_long - v_tr*v_tr) / pow(v_mag,5);
            }
            else{
                EvalIntegral(&tight_trans,params,xl,xu,(size_t) 2,result_trans,error_trans);
                EvalIntegral(&tight_long, params,xl,xu,(size_t) 2,result_long, error_long);
                //std::cout<<"Tight trans = "<<result_trans * temp<<" long = "<<result_long * temp<<std::endl;                
            }
            final_result_trans += 4 * k_pi * result_trans / sqrt(2*k_pi);
            final_result_long  += 4 * k_pi * result_long / sqrt(2*k_pi);
        }
    
        force_result_trans = f_const_ * charge_number*charge_number * density_e * final_result_trans;
        force_result_trans /= d_perp_e*d_perp_e * d_paral_e;

        force_result_long = f_const_ * charge_number*charge_number * density_e * final_result_long;
        force_result_long /= d_perp_e*d_perp_e * d_paral_e;    
    }
    else{
        force_result_trans = 0.0;
        force_result_long = 0.0;
    }
}
    
int friction_force(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                   ForceParas &force_paras, double *force_tr, double *force_long){

    double temperature = force_paras.park_temperature_eff();
    double time_cooler = force_paras.time_cooler();
    double magnetic_field = force_paras.magnetic_field();

    if(force_paras.ptr_d_perp_e()||force_paras.ptr_d_paral_e()) {
        double *d_perp_e = force_paras.ptr_d_perp_e();
        double *d_paral_e = force_paras.ptr_d_paral_e();
        force_paras.ApplyForce(charge_number, ion_number, v_tr, v_long, density_e, temperature, magnetic_field, d_perp_e,
                            d_paral_e, time_cooler, force_tr, force_long);

    }
    else {
        double d_perp_e = force_paras.d_perp_e();   //From ebeam.v_rms_tr
        double d_paral_e = force_paras.d_paral_e(); //From ebeam.v_rms_long
              
        force_paras.ApplyForce(charge_number, ion_number, v_tr, v_long, density_e, temperature, magnetic_field, d_perp_e,
                            d_paral_e, time_cooler, force_tr, force_long, force_paras.do_test());
    }
    
    return 0;
}

//Define our switch on the derived class in only one place
ForceParas* ChooseForce(ForceFormula force_formula){

    ForceParas *force_paras;
    
    switch(force_formula){
        case (ForceFormula::PARKHOMCHUK):
            force_paras = new Force_Parkhomchuk();
            break;
        
        case(ForceFormula::DERBENEVSKRINSKY):
            force_paras = new Force_DS();
            break;

        case(ForceFormula::MESHKOV):
            force_paras = new Force_Meshkov();
            break;
        
        case(ForceFormula::UNMAGNETIZED):
            force_paras = new Force_Unmagnetized();
            break;
            
        case(ForceFormula::BUDKER):
            force_paras = new Force_Budker();
            break;
            
        case(ForceFormula::ERLANGEN):
            force_paras = new Force_Erlangen();
            break;
            
    }
    return force_paras;

}