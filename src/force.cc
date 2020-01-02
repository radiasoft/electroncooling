#include "force.h"

#include <cstdio>
#include <cassert>

#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte.h>
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
    
    //Open up an output file if we're in the testing phase
    std::ofstream outfile;
    if(do_test){
      outfile.open(test_filename_);
    }   

//Compiler should ignore #pragma it doesn't understand
//TODO: Check the effect of this on do_test output
    #pragma omp parallel for ordered schedule(dynamic,10)

    for(unsigned long int i=0; i<ion_number; ++i){
        double result_trans,result_long;
        force(v_tr[i],v_long[i],d_perp_e,d_paral_e,temperature,charge_number,
              density_e[i],time_cooler,magnetic_field, result_trans, result_long);

        //the force in newtons
        force_tr[i] = result_trans;
        force_long[i] = result_long;

        if(do_test){
            #pragma omp ordered
            outfile<<f_const_ * charge_number*charge_number <<", "<<v_tr[i]<<", "<<v_long[i]
                   <<", "<<density_e[i]<<", "<<force_tr[i]<<", "<<force_long[i]<<"\n";
        }
    }

    if(do_test){
        outfile.close();
    }
    
    return 0;
}

//This function is used for calculating the forces on single particles or arrays of particles
void Force_Parkhomchuk::force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,int charge_number,
                              double density_e,double time_cooler,double magnetic_field, double &force_result_trans, double &force_result_long){

    double v2 = v_tr*v_tr + v_long*v_long;    
    double force = 0.0;
    
    if(v2>0){
        double rho_lamor = k_me_kg * d_perp_e / ( magnetic_field * k_e );
        double v2_eff_e = temperature*k_c*k_c/(k_me*1e6);
        double dlt2_eff_e = d_paral_e*d_paral_e+v2_eff_e;

        double dlt = v2 + dlt2_eff_e;

        double rho_min_const = charge_number*k_e*k_ke*k_c*k_c/(k_me*1e6);
        double rho_min = rho_min_const/dlt;
        dlt = sqrt(dlt);
        double wp_const = 4*k_pi*k_c*k_c*k_e*k_ke/(k_me*1e6);
        double wp = sqrt(wp_const*density_e);

        double rho_max = dlt/wp;
        double rho_max_2 = pow(3*charge_number/density_e, 1.0/3);
        if(rho_max<rho_max_2) rho_max = rho_max_2;
        double rho_max_3 = dlt*time_cooler;
        if(rho_max>rho_max_3) rho_max = rho_max_3;

        double lc = log((rho_max+rho_min+rho_lamor)/(rho_min+rho_lamor));   //Coulomb Logarithm

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

    //struct int_info *p = (struct int_info *)params;    
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

//This function is used for calculating the forces on single particles or arrays of particles
void Force_DS::force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,int charge_number,
                   double density_e,double time_cooler,double magnetic_field, double &force_result_tr, double &force_result_long){
    
  double v2_eff_e   = temperature * k_c*k_c / (k_me*1e6);
  double dlt2_eff_e = d_paral_e*d_paral_e + v2_eff_e;
  double rho_L      = k_me_kg * d_perp_e / ( magnetic_field * k_e );
  double wp_const   = 4*k_pi * k_c*k_c * k_e * k_ke/(k_me*1e6);
    
  //Calculate rho_max as in the Parkhomchuk model
  double v2      = v_tr*v_tr + v_long*v_long;
  double dlt     = sqrt(v2 + dlt2_eff_e);
  double wp      = sqrt(wp_const*density_e);

  double rho_max   = dlt/wp;
  double rho_max_2 = pow(3*charge_number/density_e, 1.0/3);
  if (rho_max < rho_max_2) rho_max = rho_max_2;
  double rho_max_3 = dlt*time_cooler;
  if(rho_max > rho_max_3) rho_max = rho_max_3;
        
  //Maximum impact parameter log ratio, BETACOOL eq. (3.47)
  //  This has been brought outside of the integrand
  double lm = log(rho_max/rho_L); //value is ~11      
            
  // If rho_max < rho_L the coulomb log is < 1, the
  // interaction doesn't happen
  if( lm < 1. ) lm = 0.0;
        
  int_info params;
  params.V_trans = v_tr; //Ion velocity components
  params.V_long  = v_long;
  params.width = d_paral_e; //The electron bunch width RMS
      
  unsigned int space_size = 100;
  gsl_set_error_handler_off();
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(space_size);
  double result_trans,result_long, error_trans,error_long;
  gsl_function F;

  F.params = &params;
  F.function = &Force_DS::trans_integrand;

  int status = gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-6,space_size,1,w,&result_trans,&error_trans);
  
  if(status == GSL_EDIVERGE){
    status = gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-10,space_size,1,w,&result_trans,&error_trans);
  }
       
  F.function = &Force_DS::long_integrand;

  status = gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-6,space_size,1,w,&result_long,&error_long);
  if(status == GSL_EDIVERGE){
    status = gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-10,space_size,1,w,&result_long,&error_long);
  }
     
  force_result_tr = f_const_ * density_e * charge_number*charge_number * lm * result_trans;
  force_result_tr /= (d_paral_e*d_paral_e); 

  force_result_long = f_const_ * density_e * charge_number*charge_number * lm * result_long;
  force_result_long /= (d_paral_e*d_paral_e);        
          
  gsl_integration_workspace_free(w);
   
}

void Force_Meshkov::force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,int charge_number,
                          double density_e,double time_cooler,double magnetic_field, double &force_result_trans, double &force_result_long){
    
  double v2_eff_e   = temperature * k_c*k_c / (k_me*1e6);
  double dlt2_eff_e = d_paral_e*d_paral_e + v2_eff_e;
  double rho_L      = (k_me*1e6) * d_perp_e / (k_c*k_c * magnetic_field);
  double wp_const   = 4*k_pi * k_c*k_c * k_e * k_ke/(k_me*1e6);
    
  //A fudge factor to smooth the friction force shape. Using the
  // "Classical" definition here from the BETACOOL documentation
  double k = 2;
  double rho_min_const = charge_number * k_e*k_e * k_c*k_c / (k_me * 1e6);
    
  double v2      = v_tr*v_tr + v_long*v_long;
  double v3      = pow(sqrt(v2),3);
  double dlt     = sqrt(v2 + dlt2_eff_e);
  double wp      = sqrt(wp_const*density_e);

  //The Number of multiple adiabatic collisions with a single electron
  double N_col = 1 + d_perp_e/(k_pi*sqrt(v2 + d_paral_e*d_paral_e));
      
  //dynamic shielding radius
  double rho_sh = sqrt(v2 + d_paral_e*d_paral_e) / wp;
  //minimum impact parameter  
  double rho_min = rho_min_const / (v2 + d_paral_e*d_paral_e); 
  //intermediate impact parameter
  double rho_F = rho_L * sqrt(v2 + d_paral_e*d_paral_e) / d_perp_e;
        
  double rho_max   = dlt/wp;
  double rho_max_2 = pow(3*charge_number/density_e, 1.0/3);
  if (rho_max < rho_max_2) rho_max = rho_max_2;
  double rho_max_3 = dlt*time_cooler;
  if(rho_max > rho_max_3) rho_max = rho_max_3;
        
  //Coulomb Logarithms
  double L_M = log(rho_max / (k*rho_L)); //there's a typo in the BETACOOL guide, the numerator here is rho_max
  double L_A = log((k*rho_L) / rho_F);
  double L_F = log(rho_F / rho_min);
    
  //If it's < 0, the interaction is absent
  if( L_M < 1. ) L_M = 0.0;
  if( L_A < 1. ) L_A = 0.0;
  if( L_F < 1. ) L_F = 0.0;  

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
        
  //The factor of pi/2 comes from the difference between the constants
  // used in Parkhomchuk with the constants used in the Meshkov representation
  force_result_trans = f_const_ * charge_number*charge_number * density_e * result_trans * v_tr;
  force_result_long = f_const_ * charge_number*charge_number * density_e * result_long * v_long;
}


void Force_Budker::force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,int charge_number,
                  double density_e,double time_cooler, double magnetic_field, double &force_result_trans, double &force_result_long){

    //Here we are assuming a maxwellian distribution, d_paral_e = d_perp_e
    // and ignoring magnetic field
    
    double delta_e = d_paral_e; 

    double v_mag = sqrt(v_tr*v_tr + v_long*v_long);
        
    double rho_max = (v_mag*v_mag + delta_e*delta_e)*k_me;
    rho_max /= sqrt(4 * k_pi * density_e * k_e*k_e);
    
    double rho_max_2 = pow(3*charge_number / density_e,1/3);
    if(rho_max < rho_max_2) rho_max = rho_max_2;
    double rho_max_3 = v_mag*time_cooler;
    if(rho_max > rho_max_3) rho_max = rho_max_3;
    
    double rho_min = charge_number * k_e*k_e / (k_me*(v_mag*v_mag + delta_e*delta_e));
    double lc = log(rho_max/rho_min);
    if(lc<0.) lc = 0.;
    
    double arg = v_mag/delta_e;
    double phi = sqrt(2/k_pi) * (erf(arg/sqrt(2)) - arg *exp((-arg*arg)/2));
    double result = phi * pow(v_mag,-3);
    
    force_result_trans = f_const_ * charge_number*charge_number * lc * density_e * result * v_tr;
    force_result_long = f_const_ * charge_number*charge_number * lc * density_e * result * v_long;
    
}

double Force_Erlangen::fast_trans(double *k, size_t dim, void *params){

    (void) (dim); //avoid unused parameter warnings
    
    //alpha[0] = integrating variable v_trans 
    //alpha [1] = integrating variable v_long
    //phi = integrating angle phi
    
    double alpha = k[0];
    double beta = k[1];
    double phi = k[2];
   
    struct int_info *p = (struct int_info *)params;    
    
    double y = p->V_long - alpha;
    double z = p->V_trans - beta*cos(phi);
    double d_perp_e = p->d_perp_e;
    double d_paral_e = p->d_paral_e;
   
    double U = z*exp(-alpha*alpha/(2*d_perp_e*d_perp_e) - beta*beta/(2*d_paral_e*d_paral_e));
    U /= pow(y*y + z*z + pow(beta*sin(phi),2),3/2);
    U *= alpha;
    
    return U;
}

double Force_Erlangen::fast_long(double *k, size_t dim, void *params){
    
    (void) (dim); //avoid unused parameter warnings
    
    //alpha[0] = integrating variable v_trans 
    //alpha [1] = integrating variable v_long
    //phi = integrating angle phi
    
    double alpha = k[0];
    double beta = k[1];
    double phi = k[2];
    
    struct int_info *p = (struct int_info *)params;  
    
//    std::cout<<alpha<<" "<<beta<<" "<<phi<<" "<<p->V_long<<" "<<p->V_trans<<" "<<p->d_perp_e<<" "<<p->d_paral_e<<std::endl;
    
    double y = p->V_long - alpha;
    double z = p->V_trans - beta*cos(phi);
    double d_perp_e = p->d_perp_e;
    double d_paral_e = p->d_paral_e;
    
    double U = y*exp(-alpha*alpha/(2*d_perp_e*d_perp_e) - beta*beta/(2*d_paral_e*d_paral_e));
    U /= pow(y*y + z*z + pow(beta*sin(phi),2),3/2);
    U *= alpha;
    
    return U;
}


void Force_Erlangen::force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,int charge_number,
                  double density_e,double time_cooler, double magnetic_field, double &force_result_trans, double &force_result_long){

    double v_mag = sqrt(v_tr*v_tr + v_long*v_long);
    double result_trans,result_long;
    double delta_e = d_paral_e; 

    if(v_mag > 0.){
        double rho_max = (v_mag*v_mag + delta_e*delta_e)*k_me / sqrt(4*k_pi * density_e * k_e*k_e);
        double rho_max_2 = pow(3*charge_number / density_e,1/3);
        double rho_max_3 = sqrt(v_mag*v_mag + delta_e*delta_e)*time_cooler;
        if(rho_max < rho_max_2) rho_max = rho_max_2;
        if(rho_max > rho_max_3) rho_max = rho_max_3;

        double result_trans,result_long, error_trans,error_long;
        
        if(fast_){
            
            int_info params;
            params.V_trans = v_tr; //Ion velocity components
            params.V_long  = v_long;
            params.d_paral_e = d_paral_e; //The electron bunch width RMS
            params.d_perp_e = d_perp_e;
            
            //Lower limits to the variables 
            double xl[3] = {0.,-cutoff_,0.};
            //Upper limits to the variables
            double xu[3] = {cutoff_,cutoff_,2*k_pi};

            const gsl_rng_type *T;
            gsl_rng *r;
            
            gsl_monte_function G = {&fast_trans, 3, &params}; 
            
            size_t calls = 5e5;
            
            gsl_rng_env_setup();
            
            T = gsl_rng_default;
            r = gsl_rng_alloc(T);
            
            
            gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(3);

            gsl_monte_vegas_integrate(&G,xl,xu,3,1e5,r,s, &result_trans,&error_trans);

            do
              {
                gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
                                           &result_trans, &error_trans);
//                printf ("result = % .6f sigma = % .6f "
//                        "chisq/dof = %.1f\n", result_trans, error_trans, gsl_monte_vegas_chisq (s));
              }
            while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

            G = {&fast_long, 3, &params}; //what are these?
            
            do
              {
                gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
                                           &result_long, &error_long);
//                printf ("result = % .6f sigma = % .6f "
//                        "chisq/dof = %.1f\n", result_long, error_long, gsl_monte_vegas_chisq (s));
              }
            while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
            
            gsl_monte_vegas_free (s);
            gsl_rng_free (r);
            
        }
        
//        std::cout<<result_trans<<" "<<result_long<<std::endl;
        
        if(stretched_){}
        
        if(tight_){}
        
    }
    
    force_result_trans = f_const_ * charge_number*charge_number * density_e * result_trans;
    force_result_trans /= d_perp_e*d_perp_e * d_paral_e;
    force_result_long = f_const_ * charge_number*charge_number * density_e * result_long;
    force_result_long /= d_perp_e*d_perp_e * d_paral_e;
    
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

        case(ForceFormula::BUDKER):
            force_paras = new Force_Budker();
            break;
            
        case(ForceFormula::ERLANGEN):
            force_paras = new Force_Erlangen();
            break;
            
    }
    return force_paras;

}