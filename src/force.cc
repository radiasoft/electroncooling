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
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

//These functions calculate the force from the derrived class, and output the test file
// from a filename stored in the derived class header.
int ForceParas::ApplyForce(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                            double temperature, double magnetic_field, double *d_perp_e, double *d_paral_e, double time_cooler,
                            double *force_tr, double *force_long) {

//This function definition differs from the next only in the definition
// of d_perp_e or d_paral_e as double or double*.  The other implementation
// is the only one currently called by our testing suite
    
    
//Compiler should ignore #pragma it doesn't understand, but we'll place a condition
// so OpenMP can be easily turned off for debugging
#ifdef _OPENMP
    #pragma omp parallel for  
#endif
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
    
#ifdef _OPENMP
    #pragma omp parallel for 
#endif
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


double max_impact_factor(double v_dlt, int charge_number,double density_e,double time_cooler){
    
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


//Abstracted functions for the GSL integration
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
  
    if(status == GSL_EDIVERGE){ //Try again a more time-consuming way
        status = gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-10,space_size,1,w,&result,&error);
    }

    gsl_integration_workspace_free(w);
  
}

void Force_DS::force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,
                     int charge_number,double density_e,double time_cooler,double magnetic_field, 
                     double &force_result_tr, double &force_result_long){
    
  double rho_L      = k_me_kg * d_perp_e / ( magnetic_field * k_e ); //SI units, in meters
    
  //Calculate rho_max as in the Parkhomchuk model
  double v_eff      = sqrt(v_tr*v_tr + v_long*v_long);
  double rho_max = max_impact_factor(v_eff,charge_number,density_e,time_cooler);
  double lm = log(rho_max / rho_L); //Coulomb Log ~ 5

  if(lm < 0.) lm = 0.0;
    
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
    
  //A fudge factor to smooth the friction force shape. Using the
  // "Classical" default definition here from the BETACOOL documentation
  double k = 2;

  //double rho_min_const = charge_number * k_e*k_e * k_c*k_c / (k_me * 1e6); //CGS units
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

  //Coulomb Logarithms
  double L_M = log(rho_max / (k*rho_L)); //there's a typo in the BETACOOL guide, the numerator here is rho_max
  double L_A = log((k*rho_L) / rho_F);
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
        
  //The factor of pi/2 comes from the difference between the constants
  // used in Parkhomchuk with the constants used in the Meshkov representation
  force_result_trans = f_const_ * charge_number*charge_number * density_e * result_trans * v_tr;
  force_result_long = f_const_ * charge_number*charge_number * density_e * result_long * v_long;
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


    }
    return force_paras;

}