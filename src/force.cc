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

/*  
    double f_con = -4*k_c*k_c*k_ke*k_ke*k_e*k_e*k_e/(k_me*1e6);
    double v2_eff_e = temperature*k_c*k_c/(k_me*1e6);
    double wp_const = 4*k_pi*k_c*k_c*k_e*k_ke/(k_me*1e6);
    double rho_min_const = charge_number*k_e*k_ke*k_c*k_c/(k_me*1e6);
    
    double dlt2_eff_e = d_paral_e*d_paral_e+v2_eff_e;
    double rho_lamor = k_me*1e6*d_perp_e/(magnetic_field*k_c*k_c);
    
    double v2 = v_tr*v_tr+v_long*v_long;
    if(v2>0){
        double dlt = v2+dlt2_eff_e;
        //Calculate rho_min
        double rho_min = rho_min_const/dlt;
        dlt = sqrt(dlt);
        double wp = sqrt(wp_const*density_e);

        //Calculate rho_max
        double rho_max = dlt/wp;
        double rho_max_2 = pow(3*charge_number/density_e, 1.0/3);
        if(rho_max<rho_max_2) rho_max = rho_max_2;
        double rho_max_3 = dlt*time_cooler;
        if(rho_max>rho_max_3) rho_max = rho_max_3;

        double lc = log((rho_max+rho_min+rho_lamor)/(rho_min+rho_lamor));   //Coulomb Logarithm =~5
        std::cout<<lc<<std::endl;
        //Calculate friction force
        double f = f_con*density_e*lc/(dlt*dlt*dlt);
        force_result_trans = f*v_tr;
        force_result_long = f*v_long;
    }
    else{
        force_result_trans = 0;
        force_result_long = 0;
    }
*/
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

//This function is used for calculating the forces on single particles or arrays of particles
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

    }
    return force_paras;

}