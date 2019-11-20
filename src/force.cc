#include "force.h"

#include <cmath>
#include <cstdio>

#include <iostream>
#include <fstream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

//This function is used for calculating the forces on single particles or arrays of particles
double parkhomchuk_force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,int charge_number,
                         double density_e,double time_cooler,double magnetic_field){

    double v2 = v_tr*v_tr + v_long*v_long;    
    double force = 0.0;
    
    if(v2>0){
        double f_const = -4 * charge_number*charge_number * k_me_kg * pow(k_re*k_c*k_c,2);
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
        force = f_const * density_e * lc / ( dlt*dlt*dlt );
    }

    //This factor of the force must be multiplied by either v_tr or v_long to truly be
    // equal to the force in newtons
    return force;
}    

int parkhomchuk(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                double temperature, double magnetic_field, double *d_perp_e, double *d_paral_e, double time_cooler,
                double *force_tr, double *force_long) {

    for(unsigned long int i=0; i<ion_number; ++i){

        double f = parkhomchuk_force(v_tr[i],v_long[i],d_perp_e[i],d_paral_e[i],temperature,
                                     charge_number,density_e[i],time_cooler,magnetic_field);
        //the force in Newtons
        force_tr[i] = f*v_tr[i]; 
        force_long[i] = f*v_long[i]; 
    }

    return 0;
}

int parkhomchuk(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                double temperature, double magnetic_field, double d_perp_e, double d_paral_e, double time_cooler,
                double *force_tr, double *force_long, bool do_test) {
    
    //TODO: Remove the expectation of f_const from the test output file, then remove this definition   
    double f_const = -4 * charge_number*charge_number * k_me_kg * pow(k_re*k_c*k_c,2);

    //Open up an output file if we're in the testing phase
    std::ofstream outfile;
    if(do_test){
      outfile.open("Parkhomchuk.txt");
    }   
    
    for(unsigned long int i=0; i<ion_number; ++i){

        double f = parkhomchuk_force(v_tr[i],v_long[i],d_perp_e,d_paral_e,temperature,
                                     charge_number,density_e[i],time_cooler,magnetic_field);

        //the force in newtons
        force_tr[i] = f*v_tr[i];
        force_long[i] = f*v_long[i];

        if(do_test){
            outfile<<f_const<<", "<<v_tr[i]<<", "<<v_long[i]<<", "<<density_e[i]<<", "<<force_tr[i]<<", "<<force_long[i]<<"\n";
        }
    }
    
    outfile.close();

    return 0;
}


//Convenience struct to pass all of the information
// needed to calculate the integral
struct int_info {
	double V_trans;
	double V_long;
	double width; 
};

double bunch_shape(double v_e,void *params){

    struct int_info *p = (struct int_info *)params;
    double pdf = 1.0 / p->width;
	pdf *= exp(-(v_e * v_e) / (2.0 * p->width*p->width));
    
    return pdf;
}

double DS_trans_integrand(double alpha,void *params){
    //Using Pestrikov's trick to change the integration variable

    struct int_info *p = (struct int_info *)params;    
    
    double y = p->V_long / p->width;
    double z = p->V_trans / p->width;
    double x = y + abs(z)*tan(alpha);
    
    double integrand = tan(alpha)*(y*cos(alpha) + abs(z)*sin(alpha));
    integrand *= exp(-0.5*pow(y + abs(z)*tan(alpha),2));
    integrand *= copysign(1,z);
    
    return integrand;    
}

double DS_long_integrand(double alpha, void *params){
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
int DS_force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,int charge_number,
                         double density_e,double time_cooler,double magnetic_field, double &force_tr, double &force_long){
    
  double f_const    = -2 * charge_number*charge_number * density_e * k_me_kg * pow(k_re*k_c*k_c,2);
    
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
       
  //In the middle range where we must evaluate the integral
        
  int_info params;
  params.V_trans = v_tr; //Ion velocity components
  params.V_long  = v_long;
  params.width = d_paral_e; //The electron bunch width RMS
      
  unsigned int space_size = 100;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(space_size);
  double result_trans,result_long, error_trans,error_long;
  gsl_function F;

  F.params = &params;
  F.function = &DS_trans_integrand;

  gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-7,space_size,1,w,&result_trans,&error_trans);
        
  F.function = &DS_long_integrand;

  gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-7,space_size,1,w,&result_long,&error_long);
     
  //The factor of pi/(2sqrt(2pi)) comes from the difference between the constants
  // used in Parkhomchuk with the constants used in the Pestrikov D&S integrals
  force_tr = f_const * (k_pi/(2*sqrt(2*k_pi))) * lm * result_trans / (d_paral_e*d_paral_e); 
  force_long = f_const * (k_pi/(2*sqrt(2*k_pi))) * lm * result_long / (d_paral_e*d_paral_e);        
          
  gsl_integration_workspace_free(w);
   
  return 0;
}

int DerbenevSkrinsky(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
        double temperature, double magnetic_field, double *d_perp_e, double *d_paral_e, double time_cooler,
        double *force_tr, double *force_long) {    

  for(unsigned long int i=0;i<ion_number; ++i){
      double result_trans,result_long;
      
      DS_force(v_tr[i], v_long[i], d_perp_e[i], d_paral_e[i], temperature, charge_number,
               density_e[i], time_cooler, magnetic_field, result_trans,result_long);

      //The force in newtons
      force_tr[i] = result_trans;
      force_long[i] = result_long;
      
    }
    
    return 0;
}

int DerbenevSkrinsky(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
        double temperature, double magnetic_field, double d_perp_e, double d_paral_e, double time_cooler,
        double *force_tr, double *force_long, bool do_test) {    

  double f_const = -2 * charge_number*charge_number * k_me_kg * pow(k_re*k_c*k_c,2);
    
  //Open up an output file if we're in the testing phase
  std::ofstream outfile;
  if(do_test){
      outfile.open("DerbenevSkrinsky.txt");
  }   
    
  for(unsigned long int i=0;i<ion_number; ++i){
      double result_trans,result_long;
      
      DS_force(v_tr[i], v_long[i], d_perp_e, d_paral_e, temperature, charge_number,
               density_e[i], time_cooler, magnetic_field, result_trans,result_long);

      force_tr[i] = result_trans;
      force_long[i] = result_long;
      
      if(do_test){
          outfile<<f_const<<", "<<v_tr[i]<<", "<<v_long[i]<<", "<<density_e[i]<<", "<<force_tr[i]<<", "<<force_long[i]<<"\n";
      }
    }
    
    outfile.close();
       
    return 0;
}


int friction_force(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                   ForceParas &force_paras, double *force_tr, double *force_long){
    switch (force_paras.formula()) {
        case ForceFormula::PARKHOMCHUK: {
            double temperature = force_paras.park_temperature_eff();
            double time_cooler = force_paras.time_cooler();
            double magnetic_field = force_paras.magnetic_field();
            //Use this overloaded function if the ebeam is bunched
            if(force_paras.ptr_d_perp_e()||force_paras.ptr_d_paral_e()) {
                double *d_perp_e = force_paras.ptr_d_perp_e();
                double *d_paral_e = force_paras.ptr_d_paral_e();
                parkhomchuk(charge_number, ion_number, v_tr, v_long, density_e, temperature, magnetic_field, d_perp_e,
                            d_paral_e, time_cooler, force_tr, force_long);
            }
            else {
                double d_perp_e = force_paras.d_perp_e();
                double d_paral_e = force_paras.d_paral_e();
                parkhomchuk(charge_number, ion_number, v_tr, v_long, density_e, temperature,  magnetic_field, d_perp_e,
                            d_paral_e, time_cooler, force_tr, force_long, force_paras.do_test());
            }
            break;
        }
        case ForceFormula::DERBENEVSKRINSKY: {
            double temperature = force_paras.park_temperature_eff();
            double time_cooler = force_paras.time_cooler();
            double magnetic_field = force_paras.magnetic_field();
        	//Use this overloaded function if the ebeam is bunched
            if(force_paras.ptr_d_perp_e()||force_paras.ptr_d_paral_e()) {
                double *d_perp_e = force_paras.ptr_d_perp_e();
                double *d_paral_e = force_paras.ptr_d_paral_e();
                DerbenevSkrinsky(charge_number, ion_number, v_tr, v_long, density_e, temperature, magnetic_field, d_perp_e,
                         d_paral_e, time_cooler, force_tr, force_long);
        		break;
            }
            else {
                double d_perp_e = force_paras.d_perp_e();   //From ebeam.v_rms_tr
                double d_paral_e = force_paras.d_paral_e(); //From ebeam.v_rms_long

                DerbenevSkrinsky(charge_number, ion_number, v_tr, v_long, density_e, temperature,  magnetic_field, d_perp_e,
                            	d_paral_e, time_cooler, force_tr, force_long, force_paras.do_test());
            }
            break;
        }


        default: {
            perror("Choose your formula for friction force calculation!");
            break;
        }
    }
    return 0;
}
