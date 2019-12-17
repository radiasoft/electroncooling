#include "force.h"

#include <cstdio>

#include <iostream>
#include <fstream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

//This function is used for calculating the forces on single particles or arrays of particles
double parkhomchuk_force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,
                         int charge_number, double density_e,double time_cooler,double magnetic_field){

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
        force = f_const * charge_number*charge_number * density_e * lc / ( dlt*dlt*dlt );
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
            outfile<<f_const * charge_number*charge_number <<", "<<v_tr[i]<<", "<<v_long[i]<<", "<<density_e[i]<<", "<<force_tr[i]<<", "<<force_long[i]<<"\n";
        }
    }

    if(do_test){
        outfile.close();
    }
    
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
    //The F_const already has a minus sign, so the integrand term
    // should carry the sign of z
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
                         double density_e,double time_cooler,double magnetic_field, double &DS_force_tr, double &DS_force_long){
    
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
  F.function = &DS_trans_integrand;

  int status = gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-6,space_size,1,w,&result_trans,&error_trans);
  
  if(status == GSL_EDIVERGE){
    status = gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-10,space_size,1,w,&result_trans,&error_trans);
  }
       
  F.function = &DS_long_integrand;

  status = gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-6,space_size,1,w,&result_long,&error_long);
  if(status == GSL_EDIVERGE){
    status = gsl_integration_qag(&F, -k_pi/2, k_pi/2, 1, 1e-10,space_size,1,w,&result_long,&error_long);
  }
     
  //The factor of (0.5*pi/(2sqrt(2pi))) comes from the difference between the constants
  // used in Parkhomchuk with the constants used in the Pestrikov D&S integrals
  DS_force_tr = f_const * density_e * charge_number*charge_number * 0.5*(k_pi/(2*sqrt(2*k_pi))) * lm * result_trans;
  DS_force_tr /= (d_paral_e*d_paral_e); 
  DS_force_long = f_const * density_e * charge_number*charge_number * 0.5*(k_pi/(2*sqrt(2*k_pi))) * lm * result_long;
  DS_force_long /= (d_paral_e*d_paral_e);        
          
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
          //TODO: Remove f_const as a column
          outfile<<f_const<<", "<<v_tr[i]<<", "<<v_long[i]<<", "<<density_e[i]<<", "<<force_tr[i]<<", "<<force_long[i]<<"\n";
      }
    }

    if(do_test){
        outfile.close();
    }
    
    return 0;
}

int Meshkov_force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,int charge_number,
                  double density_e,double time_cooler,double magnetic_field, double &Mesh_force_tr, double &Mesh_force_long){
    
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
  Mesh_force_tr = f_const * charge_number*charge_number * (k_pi/2) * density_e * result_trans * v_tr;
  Mesh_force_long = f_const * charge_number*charge_number * (k_pi/2) * density_e * result_long * v_long;

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

    for(unsigned long int i=0;i<ion_number; ++i){
        double result_trans,result_long;
      
        Meshkov_force(v_tr[i], v_long[i], d_perp_e, d_paral_e, temperature, charge_number,
                   density_e[i], time_cooler, magnetic_field, result_trans,result_long);

        force_tr[i] = result_trans;
        force_long[i] = result_long;
      
        if(do_test){
          outfile<<f_const<<", "<<v_tr[i]<<", "<<v_long[i]<<", "<<density_e[i]<<", "<<force_tr[i]<<", "<<force_long[i]<<"\n";
      }
    }

    if(do_test){
        outfile.close();
    }
    
    return 0;
}

//This version for when d_perp_e and d_paral_e are arrays
int Meshkov(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
        double temperature, double magnetic_field, double *d_perp_e, double *d_paral_e, 
        double time_cooler, double *force_tr, double *force_long) {

    for(unsigned long int i=0;i<ion_number; ++i){
        double result_trans,result_long;
      
        Meshkov_force(v_tr[i], v_long[i], d_perp_e[i], d_paral_e[i], temperature, charge_number,
                   density_e[i], time_cooler, magnetic_field, result_trans,result_long);

        force_tr[i] = result_trans;
        force_long[i] = result_long;
      
    }
    
    return 0;
}


int Budker_force(double v_tr, double v_long, double d_perp_e, double d_paral_e, double temperature,int charge_number,
                  double density_e,double time_cooler, double &Bud_force_tr, double &Bud_force_long){

    //Here we are assuming a maxwellian distribution, d_paral_e = d_perp_e
    // and ignoring magnetic field
    
    double delta_e = d_paral_e; 

    double v_mag = sqrt(v_tr*v_tr + v_long*v_long);
    
/*
    //The mean minimal impact parameter
    double rho_min = charge_number * k_e*k_e / k_me_kg;
    rho_min /= v_tr*v_tr + v_long*v_long + 2*delta_e*delta_e;

    double v2_eff_e = temperature*k_c*k_c/(k_me*1e6);
    double dlt2_eff_e = delta_e*delta_e+v2_eff_e;

    double dlt = v_mag*v_mag + dlt2_eff_e;
    dlt = sqrt(dlt);
    double wp_const = 4*k_pi*k_c*k_c*k_e*k_ke/(k_me*1e6);
    double wp = sqrt(wp_const*density_e);
    
    double rho_max = dlt/wp; 
    double rho_max_2 = pow(3*charge_number/density_e, 1.0/3);
    if(rho_max<rho_max_2) rho_max = rho_max_2;
    double rho_max_3 = dlt*time_cooler;
    if(rho_max>rho_max_3) rho_max = rho_max_3;
    
    double lc = log(rho_max/rho_min);
*/
    
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
    
    Bud_force_tr = f_const * charge_number*charge_number * (k_pi/2) * lc * density_e * result * v_tr;
    Bud_force_long = f_const * charge_number*charge_number * (k_pi/2) * lc * density_e * result * v_long;
    
    return 0;
}

//This version for when d_perp_e and d_paral_e are doubles
int Budker(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
        double temperature, double d_perp_e, double d_paral_e, 
        double time_cooler, double *force_tr, double *force_long,bool do_test) {
    
    //Open up an output file if we're in the testing phase
    std::ofstream outfile;

    if(do_test){
      outfile.open("Budker.txt");
    }          

    for(unsigned long int i=0;i<ion_number; ++i){
        double result_trans,result_long;
      
        Budker_force(v_tr[i], v_long[i], d_perp_e, d_paral_e, temperature, charge_number,
                   density_e[i], time_cooler, result_trans,result_long);

        force_tr[i] = result_trans;
        force_long[i] = result_long;
      
        if(do_test){
          outfile<<f_const<<", "<<v_tr[i]<<", "<<v_long[i]<<", "<<density_e[i]<<", "<<force_tr[i]<<", "<<force_long[i]<<"\n";
      }
    }
    if(do_test){
        outfile.close();
    }
    
    return 0;
}

//This version for when d_perp_e and d_paral_e are arrays
int Budker(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
        double temperature, double *d_perp_e, double *d_paral_e, 
        double time_cooler, double *force_tr, double *force_long) {

    for(unsigned long int i=0;i<ion_number; ++i){
        double result_trans,result_long;
      
        Budker_force(v_tr[i], v_long[i], d_perp_e[i], d_paral_e[i], temperature, charge_number,
                   density_e[i], time_cooler, result_trans,result_long);

        force_tr[i] = result_trans;
        force_long[i] = result_long;
      
    }
    
    return 0;
}

    
int friction_force(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                   ForceParas &force_paras, double *force_tr, double *force_long){

    double temperature = force_paras.park_temperature_eff();
    double time_cooler = force_paras.time_cooler();
    double magnetic_field = force_paras.magnetic_field();

    switch (force_paras.formula()) {
        case ForceFormula::PARKHOMCHUK: {
            //Use this overloaded function if the ebeam is user-defined
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
        	//Use this overloaded function if the ebeam is user-defined
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
        case ForceFormula::MESHKOV: {       	
            if(force_paras.ptr_d_perp_e()||force_paras.ptr_d_paral_e()) {
                double *d_perp_e = force_paras.ptr_d_perp_e();
                double *d_paral_e = force_paras.ptr_d_paral_e();
                Meshkov(charge_number, ion_number, v_tr, v_long, density_e, temperature, magnetic_field, d_perp_e,
                              d_paral_e, time_cooler, force_tr, force_long);
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
        case ForceFormula::BUDKER: {        	
            if(force_paras.ptr_d_perp_e()||force_paras.ptr_d_paral_e()) {
                double *d_perp_e = force_paras.ptr_d_perp_e();
                double *d_paral_e = force_paras.ptr_d_paral_e();
                Budker(charge_number, ion_number, v_tr, v_long, density_e, temperature, d_perp_e,
                              d_paral_e, time_cooler, force_tr, force_long);
        		break;
            }
            else {
                double d_perp_e = force_paras.d_perp_e();   //From ebeam.v_rms_tr
                double d_paral_e = force_paras.d_paral_e(); //From ebeam.v_rms_long
                
                Budker(charge_number, ion_number, v_tr, v_long, density_e, temperature, d_perp_e,
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
