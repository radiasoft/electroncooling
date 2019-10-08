#include "force.h"
#include <cmath>
#include <cstdio>

#include <iostream>
#include <fstream>
#include <gsl/gsl_integration.h>

int parkhomchuk(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                double temperature, double magnetic_field, double *d_perp_e, double *d_paral_e, double time_cooler,
                double *force_tr, double *force_long) {
    double f_const = -4*k_c*k_c*k_ke*k_ke*k_e*k_e*k_e/(k_me*1e6);
    double v2_eff_e = temperature*k_c*k_c/(k_me*1e6);
//    double dlt2_eff_e = dlt_paral_e*dlt_paral_e+V2_eff_e;
    double rho_Lamor = k_me * 1e6 * (*d_perp_e)/(magnetic_field*k_c*k_c);
    double wp_const = 4*k_pi*k_c*k_c*k_e*k_ke/(k_me*1e6);
    double rho_min_const = charge_number*k_e*k_ke*k_c*k_c/(k_me*1e6);

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
    }
    return 0;
}

int parkhomchuk(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                double temperature, double magnetic_field, double d_perp_e, double d_paral_e, double time_cooler,
                double *force_tr, double *force_long) {
    double f_const = -4*charge_number*charge_number*k_c*k_c*k_ke*k_ke*k_e*k_e*k_e/(k_me*1e6);
    double v2_eff_e = temperature*k_c*k_c/(k_me*1e6);
    double dlt2_eff_e = d_paral_e*d_paral_e+v2_eff_e;
    double rho_lamor = k_me*1e6*d_perp_e/(magnetic_field*k_c*k_c);
    double wp_const = 4*k_pi*k_c*k_c*k_e*k_ke/(k_me*1e6);
    double rho_min_const = charge_number*k_e*k_ke*k_c*k_c/(k_me*1e6);

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
    }
    return 0;
}


//Convenience struct to pass all of the information
// needed to calculate the integral
struct int_info {
	double V_trans;
	double V_long;
	double v_e;     //The integrating variable
	double Delta_e;
	double width;   //A priori 1e5
};

double DS_trans_integrand(double v_e, void * params){
	//Integrating over v_e

	struct int_info *p = (struct int_info *)params;
	double U = sqrt(p->V_trans*p->V_trans + (p->V_long - p->v_e)*(p->V_long - p->v_e));
	double integrand = (p->V_trans * (p->V_trans*p->V_trans - 2*(p->V_long - p->v_e)*(p->V_long - p->v_e))) / pow(U,5);

	//Now calculate gaussian PDF of the electron bunch (centered on 0 with standard deviation width)
	double pdf = 1.0 / (sqrt(2.0 * k_pi));
	pdf *= exp(-(p->v_e * p->v_e) / (2.0 * p->width));

	return integrand * pdf;
}

double DS_long_integrand(double v_e, void * params){
	//Integrating over v_e
	struct int_info *p = (struct int_info *)params;

	double U = sqrt(p->V_trans*p->V_trans + (p->V_long-p->v_e)*(p->V_long-p->v_e));
	double integrand;
	double V_ion = sqrt(p->V_trans*p->V_trans + p->V_long*p->V_long);
	if (p->V_trans < 1.0){
		integrand = -2 * p->V_long * exp((-p->V_long*p->V_long) / (2*p->Delta_e*p->Delta_e));
		integrand /= sqrt(2*k_pi) * pow(p->Delta_e,3);
	}

	else if (V_ion > (100*p->v_e)){
		integrand = 3*(p->V_trans*p->V_trans)*(p->V_long * p->v_e)/pow(U,5);
		integrand += 2*(p->V_long - p->v_e) / pow(U,3);
	}
	else{
		integrand = 3*(p->V_trans*p->V_trans)*(p->V_long * p->v_e)/pow(U,5);
	}

	//Now calculate gaussian PDF of the electron bunch (centered on 0 with standard deviation width)
	double pdf = 1.0 / (sqrt(2.0 * k_pi));
	pdf *= exp(-(p->v_e * p->v_e) / (2.0 * p->width));

	return integrand * pdf;
}

int DerbenevSkrinsky(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
        double temperature, double magnetic_field, double d_perp_e, double d_paral_e, double time_cooler,
        double *force_tr, double *force_long) {

	double f_const = -2 * k_pi /(k_me*1e6);
    double v2_eff_e = temperature*k_c*k_c / (k_me*1e6);
    double dlt2_eff_e = d_paral_e*d_paral_e + v2_eff_e;
    double rho_min_const = charge_number*k_e*k_ke*k_c*k_c/(k_me*1e6);

	for(unsigned long int i=0;i<ion_number; ++i){

		double rho_L  = k_c * (k_me*1e6) * d_perp_e / (k_e * magnetic_field);
		double v2     = v_tr[i]*v_tr[i] + v_long[i]*v_long[i];
		double dlt    = v2 + dlt2_eff_e;
        double rho_min = rho_min_const/dlt;
		double omega_p = sqrt(4 * k_pi * density_e[i] * k_e*k_e / (k_me*1e6));
		double rho_sh = sqrt(v2)/omega_p;

		double rho_max_2 = pow(3*charge_number/density_e[i], 1.0/3);
		double rho_max = rho_sh;
		if (rho_sh < rho_max_2) rho_max = rho_max_2;
		double rho_max_3 = dlt*time_cooler;
        if(rho_max > rho_max_3) rho_max = rho_max_3;

        double lc = log((rho_max+rho_min+rho_L)/(rho_min+rho_L));   //Coulomb Logarithm

        int_info params;
        params.V_trans = v_tr[i];
        params.V_long  = v_long[i];
        params.Delta_e = sqrt(dlt2_eff_e);
        params.width   = 1E5;

        gsl_integration_workspace *w = gsl_integration_workspace_alloc(100);
        double result_trans,result_long,error;
        gsl_function F;

        F.function = &DS_trans_integrand;
        F.params = &params;
        //Integrate over the infinite interval
        gsl_integration_qagi(&F, 1, 1e-7,100,w,&result_trans,&error);

        F.function = &DS_long_integrand;
        gsl_integration_qagi(&F, 1, 1e-7,100,w,&result_long,&error);


		force_tr[i] = f_const * charge_number * charge_number * result_trans;
		force_long[i] = f_const * charge_number * charge_number * result_long;
	}
	return 0;
}

int friction_force(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                   ForceParas &force_paras, double *force_tr, double *force_long){
    switch (force_paras.formula()) {
        case ForceFormula::PARKHOMCHUK: {
            double temperature = force_paras.park_temperature_eff();
            double time_cooler = force_paras.time_cooler();
            double magnetic_field = force_paras.magnetic_field();
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
                            d_paral_e, time_cooler, force_tr, force_long);
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
 //                         d_paral_e, time_cooler, force_tr, force_long);
        		break;
            }
            else {
                double d_perp_e = force_paras.d_perp_e();
                double d_paral_e = force_paras.d_paral_e();
                DerbenevSkrinsky(charge_number, ion_number, v_tr, v_long, density_e, temperature,  magnetic_field, d_perp_e,
                            	d_paral_e, time_cooler, force_tr, force_long);
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
