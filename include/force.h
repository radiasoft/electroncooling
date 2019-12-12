#ifndef FORCE_H
#define FORCE_H

#include <vector>
#include <cmath>
#include "constants.h"

enum class ForceFormula {PARKHOMCHUK,DERBENEVSKRINSKY,MESHKOV,BUDKER};

class ForceParas{
    ForceFormula formula_;
    double park_temperature_eff_ = 0;
    double magnetic_field_ = 0;
    double d_perp_e_;
    double d_paral_e_;
    double *ptr_d_perp_e_ = nullptr;
    double *ptr_d_paral_e_ = nullptr;
    double time_cooler_;
    bool do_test_ = false;
 public:
    ForceFormula formula(){return formula_;}
    double park_temperature_eff(){return park_temperature_eff_;}
    double magnetic_field(){return magnetic_field_;}
    double d_perp_e(){return d_perp_e_;}
    double d_paral_e(){return d_paral_e_;}
    double *ptr_d_perp_e(){return ptr_d_perp_e_;}
    double *ptr_d_paral_e(){return ptr_d_paral_e_;}
    double time_cooler(){return time_cooler_;}
    bool do_test(){return do_test_;}
    int set_park_temperature_eff(double x){park_temperature_eff_ = x; return 0;}
    int set_magnetic_field(double x){magnetic_field_=x; return 0;}
    int set_d_perp_e(double x){d_perp_e_ = x; return 0;}
    int set_d_paral_e(double x){d_paral_e_ = x; return 0;}
    int set_ptr_d_perp_e(double* x){ptr_d_perp_e_ = x; return 0;}
    int set_ptr_d_paral_e(double* x){ptr_d_paral_e_ = x; return 0;}
    int set_ptr_d_perp_e(std::vector<double>& x){ptr_d_perp_e_ = &*x.begin(); return 0;}
    int set_ptr_d_paral_e(std::vector<double>& x){ptr_d_paral_e_ = &*x.begin(); return 0;}
    int set_time_cooler(double x){time_cooler_ = x; return 0;}
    int set_do_test(bool b){do_test_ = b; return 0;}
    ForceParas(ForceFormula formula):formula_(formula){};
};

//Force-dependent constants
const double f_const = -4 * k_me_kg * pow(k_re*k_c*k_c,2); //for Parkhomchuk

int friction_force(int charge_number, unsigned long int ion_number, double *v_tr, double *v_z, double *density_e,
                  ForceParas &force_paras, double *force_tr, double *force_long);

#endif // FORCE_H
