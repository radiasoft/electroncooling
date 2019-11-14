#include "force.h"

#include <cmath>
#include <cstdio>

#include <iostream>
#include <fstream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>


int parkhomchuk(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                double temperature, double magnetic_field, double *d_perp_e, double *d_paral_e, double time_cooler,
                double *force_tr, double *force_long, bool do_test) {
    double f_const = -4*k_c*k_c*k_ke*k_ke*k_e*k_e*k_e/(k_me*1e6);
    double v2_eff_e = temperature*k_c*k_c/(k_me*1e6);
//    double dlt2_eff_e = dlt_paral_e*dlt_paral_e+V2_eff_e;
//    double rho_Lamor = k_me * 1e6 * (*d_perp_e)/(magnetic_field*k_c*k_c);
    double wp_const = 4*k_pi*k_c*k_c*k_e*k_ke/(k_me*1e6);
    double rho_min_const = charge_number*k_e*k_ke*k_c*k_c/(k_me*1e6);

    //Open up an output file if we're in the testing phase
    std::ofstream outfile;
    if(do_test){
      outfile.open("Parkhomchuk.txt");
    }   
    
    double *dlt2_eff_e = new double[ion_number];
    for(unsigned long int i=0; i<ion_number; ++i){
        dlt2_eff_e[i] = d_paral_e[i]*d_paral_e[i]+v2_eff_e;
    }
    double *rho_lamor = new double[ion_number];
    for(unsigned long int i=0; i<ion_number; ++i) {
        rho_lamor[i] = k_me*1e6*d_perp_e[i]/(magnetic_field*k_c*k_c);
    }

    for(unsigned long int i=0; i<ion_number; ++i){
        double v2 = v_tr[i]*v_tr[i]+v_long[i]*v_long[i];
        if(v2>0){
            double dlt = v2+dlt2_eff_e[i];
            //Calculate rho_min
            double rho_min = rho_min_const/dlt;
            dlt = sqrt(dlt);
            double wp = sqrt(wp_const*density_e[i]);

            //Calculate rho_max
            double rho_max = dlt/wp;
            double rho_max_2 = pow(3*charge_number/density_e[i], 1.0/3);
            if(rho_max<rho_max_2) rho_max = rho_max_2;
            double rho_max_3 = dlt*time_cooler;
            if(rho_max>rho_max_3) rho_max = rho_max_3;

            double lc = log((rho_max+rho_min+rho_lamor[i])/(rho_min+rho_lamor[i]));   //Coulomb Logarithm
            //Calculate friction force
            double f = f_const*density_e[i]*lc/(dlt*dlt*dlt);
            force_tr[i] = f*v_tr[i];
            force_long[i] = f*v_long[i];
        }
        else{
            force_tr[i] = 0;
            force_long[i] = 0;
        }
        if(do_test){
            outfile<<f_const<<", "<<v_tr[i]<<", "<<v_long[i]<<", "<<density_e[i]<<", "<<force_tr[i]<<", "<<force_long[i]<<"\n";
        }
    }
    outfile.close();
    return 0;
}

int parkhomchuk(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                double temperature, double magnetic_field, double d_perp_e, double d_paral_e, double time_cooler,
                double *force_tr, double *force_long, bool do_test) {
    double f_const = -4 * charge_number*charge_number * k_c*k_c * k_ke*k_ke * k_e*k_e*k_e /(k_me*1e6);
    double v2_eff_e = temperature*k_c*k_c/(k_me*1e6);
    double dlt2_eff_e = d_paral_e*d_paral_e + v2_eff_e;
    double rho_lamor = (k_me*1e6)*d_perp_e/(magnetic_field*k_c*k_c);
    double wp_const = 4*k_pi * k_c*k_c * k_e * k_ke/(k_me*1e6);
    double rho_min_const = charge_number * k_e * k_ke * k_c*k_c/(k_me*1e6);

    //Open up an output file if we're in the testing phase
    std::ofstream outfile;
    if(do_test){
      outfile.open("Parkhomchuk.txt");
    }   
    
    for(unsigned long int i=0; i<ion_number; ++i){
        double v2 = v_tr[i]*v_tr[i]+v_long[i]*v_long[i];
        if(v2>0){
            double dlt = v2+dlt2_eff_e;
            //Calculate rho_min
            double rho_min = rho_min_const/dlt;
            dlt = sqrt(dlt);
            double wp = sqrt(wp_const*density_e[i]);

            //Calculate rho_max
            double rho_max = dlt/wp;
            double rho_max_2 = pow(3*charge_number/density_e[i], 1.0/3);
            if(rho_max<rho_max_2) rho_max = rho_max_2;
            double rho_max_3 = dlt*time_cooler;
            if(rho_max>rho_max_3) rho_max = rho_max_3;

            double lc = log((rho_max+rho_min+rho_lamor)/(rho_min+rho_lamor));   //Coulomb Logarithm
            //Calculate friction force
            double f = f_const*density_e[i]*lc/(dlt*dlt*dlt);
            force_tr[i] = f*v_tr[i];
            force_long[i] = f*v_long[i];
        }
        else{
            force_tr[i] = 0;
            force_long[i] = 0;
        }
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

double DS_trans_integrand(double alpha,void *params){
    //Using Pestrikov's trick to change the integration variable

    struct int_info *p = (struct int_info *)params;    
    
    double y = p->V_long / p->width;
    double z = p->V_trans / p->width;
    double x = y + abs(z)*tan(alpha);
    
    double integrand = tan(alpha)*(y*cos(alpha) + abs(z)*sin(alpha));
    integrand *= exp(-0.5*pow(y+abs(z)*tan(alpha),2));
    //The F_const already has a minus sign, so the integrand term
    // should carry the sign of z
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

int DerbenevSkrinsky(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
        double temperature, double magnetic_field, double d_perp_e, double d_paral_e, double time_cooler,
        double *force_tr, double *force_long, bool do_test) {    
    
  //Constant term for Parkhomchuk functions, above
  double f_const    = -4 * charge_number*charge_number * k_c*k_c * k_ke*k_ke * k_e*k_e*k_e /(k_me*1e6);
  //double f_const = 2 * k_pi*charge_number*charge_number * pow(k_e,4)/(k_me * 1e6);
  double v2_eff_e   = temperature * k_c*k_c / (k_me*1e6);
  double dlt2_eff_e = d_paral_e*d_paral_e + v2_eff_e;
  double rho_L      = (k_me*1e6) * d_perp_e / (k_c*k_c * magnetic_field);
  double wp_const   = 4*k_pi * k_c*k_c * k_e * k_ke/(k_me*1e6);
  
  //Open up an output file if we're in the testing phase
  std::ofstream outfile;
  if(do_test){
      outfile.open("DerbenevSkrinsky.txt");
  }   
    
  for(unsigned long int i=0;i<ion_number; ++i){

      //Calculate rho_max as in the Parkhomchuk model
      double v2      = v_tr[i]*v_tr[i] + v_long[i]*v_long[i];
      double dlt     = sqrt(v2 + dlt2_eff_e);
      double wp      = sqrt(wp_const*density_e[i]);

      double rho_max   = dlt/wp;
      double rho_max_2 = pow(3*charge_number/density_e[i], 1.0/3);
      if (rho_max < rho_max_2) rho_max = rho_max_2;
      double rho_max_3 = dlt*time_cooler;
      if(rho_max > rho_max_3) rho_max = rho_max_3;
        
      //Maximum impact parameter log ratio, BETACOOL eq. (3.47)
      //  This has been brought outside of the integrand
      double lm = log(rho_max/rho_L); //value is ~11      

      // If rho_max < rho_L the coulomb log is < 1, the
      // interaction doesn't happen
      if( lm < 1. ) lm = 0.0;
        
      //In the middle range where we must evaluate the integral
         
      int_info params;
      params.V_trans = v_tr[i]; //Ion velocity components
      params.V_long  = v_long[i];
      params.width = d_paral_e; //The electron bunch width RMS
        
      unsigned int space_size = 100;
      gsl_integration_workspace *w = gsl_integration_workspace_alloc(space_size);
      double result_trans,result_long,error_trans,error_long;
      gsl_function F;

      F.params = &params;
      F.function = &DS_trans_integrand;

      gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-7,space_size,1,w,&result_trans,&error_trans);
          
      F.function = &DS_long_integrand;
          
      gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-7,space_size,1,w,&result_long,&error_long);
          
      //The factor of pi/(2sqrt(2pi)) comes from the difference between the constants
      // used in Parkhomchuk with the constants used in the Pestrikov D&S integrals      
      force_tr[i] = f_const * (k_pi/(2*sqrt(2*k_pi))) * lm * density_e[i] * result_trans / (d_paral_e*d_paral_e); 
      force_long[i] = f_const * (k_pi/(2*sqrt(2*k_pi))) * lm * density_e[i] * result_long / (d_paral_e*d_paral_e);        
          
      gsl_integration_workspace_free(w);

      if(do_test){
          outfile<<f_const<<", "<<v_tr[i]<<", "<<v_long[i]<<", "<<density_e[i]<<", "<<force_tr[i]<<", "<<force_long[i]<<"\n";
      }
    }
    
    outfile.close();
       
    return 0;
}

//This version for when d_perp_e and d_paral_e are doubles
int Meshkov(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
        double temperature, double magnetic_field, double d_perp_e, double d_paral_e, 
        double time_cooler, double *force_tr, double *force_long,bool do_test) {
    
    //Open up an output file if we're in the testing phase
    std::ofstream outfile;
    if(do_test){
      outfile.open("Meshkov.txt");
    }      
    
  //Constant term for Parkhomchuk functions, above
  double f_const    = -4 * charge_number*charge_number * k_c*k_c * k_ke*k_ke * k_e*k_e*k_e /(k_me*1e6);
  double v2_eff_e   = temperature * k_c*k_c / (k_me*1e6);
  double dlt2_eff_e = d_paral_e*d_paral_e + v2_eff_e;
  double rho_L      = (k_me*1e6) * d_perp_e / (k_c*k_c * magnetic_field);
  double wp_const   = 4*k_pi * k_c*k_c * k_e * k_ke/(k_me*1e6);
    
  //A fudge factor to smooth the friction force shape. Using the
  // "Classical" definition here from the BETACOOL documentation
  double k = 2;
  double rho_min_const = charge_number * k_e*k_e * k_c*k_c / (k_me * 1e6);
    
  for(unsigned long int i=0;i<ion_number; ++i){

      double v2      = v_tr[i]*v_tr[i] + v_long[i]*v_long[i];
      double v3      = pow(sqrt(v2),3);
      double dlt     = sqrt(v2 + dlt2_eff_e);
      double wp      = sqrt(wp_const*density_e[i]);

      //The Number of multiple adiabatic collisions with a single electron
      double N_col = 1 + d_perp_e/(k_pi*sqrt(v2 + d_paral_e*d_paral_e));
        
      //dynamic shielding radius
      double rho_sh = sqrt(v2 + d_paral_e*d_paral_e) / wp;
      //minimum impact parameter  
      double rho_min = rho_min_const / (v2 + d_paral_e*d_paral_e); 
      //intermediate impact parameter
      double rho_F = rho_L * sqrt(v2 + d_paral_e*d_paral_e) / d_perp_e;
        
      double rho_max   = dlt/wp;
      double rho_max_2 = pow(3*charge_number/density_e[i], 1.0/3);
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
    
      double ellipse = (v_tr[i]*v_tr[i])/(d_perp_e*d_perp_e) + (v_long[i]*v_long[i])/(d_paral_e*d_paral_e);
      
      if(sqrt(v2) > d_perp_e) { //Region I
              result_trans = 2*L_F + L_M*(v_tr[i]*v_tr[i] - 2*v_long[i]*v_long[i])/v2;
              result_trans /= v3;
          
              result_long = 2 + 2*L_F + L_M*(3*v_tr[i]*v_tr[i]/v2);
              result_long /= v3;
      }   
      else if(sqrt(v2) < d_paral_e) { //Region III
              result_trans = 2*(L_F + N_col*L_A)/pow(d_perp_e,3) + L_M/pow(d_paral_e,3);
              result_long = 2*(L_F + N_col*L_A)/(d_perp_e*d_perp_e * d_paral_e) + L_M/pow(d_paral_e,3);
      }  

      //Region II (a or b)
      else{ //This constrains to a donut shape (d_paral_e < sqrt(v2) < d_perp_e).
          result_trans = ((v_tr[i]*v_tr[i] - 2*v_long[i]*v_long[i])/v2) * (L_M/v3);
          result_trans += (2/pow(d_perp_e,3)) * (L_F + N_col*L_A);
          
          if( ellipse <= 1.0 ){ // Region IIb
              //This is the same as result_long in Region 3 
              result_long = 2*(L_F + N_col*L_A)/(d_perp_e*d_perp_e * d_paral_e) + L_M/pow(d_paral_e,3);
          }
          else{ //It must be Region IIa
              result_long = (((3*v_tr[i]*v_tr[i]*L_M)/v2) + 2) / v3;
              result_long += 2 * (L_F + N_col*L_A)/(d_perp_e*d_perp_e * v_long[i]);
          }
      } 
        
      //The factor of pi/2 comes from the difference between the constants
      // used in Parkhomchuk with the constants used in the Meshkov representation
      force_tr[i] = f_const * (k_pi/2) * density_e[i] * result_trans * v_tr[i];
      force_long[i] = f_const * (k_pi/2) * density_e[i] * result_long * v_long[i];
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
                            d_paral_e, time_cooler, force_tr, force_long, force_paras.do_test());
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
        	if(force_paras.ptr_d_perp_e()||force_paras.ptr_d_paral_e()) {
        		perror("This kind of parameterization not supported yet!");
//        		double *d_perp_e = force_paras.ptr_d_perp_e();
//              double *d_paral_e = force_paras.ptr_d_paral_e();
//              DerbenevSkrinsky(charge_number, ion_number, v_tr, v_long, density_e, temperature, magnetic_field, d_perp_e,
 //                         d_paral_e, time_cooler, force_tr, force_long, force_paras.do_test());
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
        case ForceFormula::MESHKOV: {
            double temperature = force_paras.park_temperature_eff();
            double time_cooler = force_paras.time_cooler();
            double magnetic_field = force_paras.magnetic_field();
        	if(force_paras.ptr_d_perp_e()||force_paras.ptr_d_paral_e()) {
                perror("Must use constant V_RMS for electrons for Meshkov model");
//                double *d_perp_e = force_paras.ptr_d_perp_e();
//                double *d_paral_e = force_paras.ptr_d_paral_e();
//                Meshkov(charge_number, ion_number, v_tr, v_long, density_e, temperature, magnetic_field, d_perp_e,
//                              d_paral_e, time_cooler, force_tr, force_long);
        		break;
            }
            else {
                double d_perp_e = force_paras.d_perp_e();   //From ebeam.v_rms_tr
                double d_paral_e = force_paras.d_paral_e(); //From ebeam.v_rms_long
                
                Meshkov(charge_number, ion_number, v_tr, v_long, density_e, temperature,  magnetic_field, d_perp_e,
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
