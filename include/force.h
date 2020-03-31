#ifndef FORCE_H
#define FORCE_H

#include <vector>
#include <cmath>
#include "constants.h"

enum class ForceFormula {PARKHOMCHUK,DERBENEVSKRINSKY,MESHKOV,UNMAGNETIZED,BUDKER};

class ForceParas{
 protected:
    ForceFormula formula_;
    double park_temperature_eff_ = 0;
    double magnetic_field_       = 0;
    double d_perp_e_;
    double d_paral_e_;
    double *ptr_d_perp_e_  = nullptr;
    double *ptr_d_paral_e_ = nullptr;
    double time_cooler_;
    bool do_test_ = false;
    const char* test_filename_;
    double f_const_;

    //Variables related to the integration required for some force models
    double cutoff_ = 1e5; //cutoff value for indefinite integrals
    size_t calls_  = 5e5; //The number of MC samples when calculating integrals  
    struct int_info {
        double V_trans;
        double V_long;
        double width; 
        double d_perp_e;
        double d_paral_e;
        double rho_max;
        double B_field;
        int Z;
    };
 
    //For Erlangen model
    bool fast_      = true;
    bool stretched_ = false;
    bool tight_     = false;

    //For the Un=magnetized model
    bool approximate_ = false; //Use a simplification in the integration
    bool binney_ = true; //A further simplification
    
 public:
    ForceFormula formula()        const {return formula_;}
    double park_temperature_eff() const {return park_temperature_eff_;}
    double magnetic_field()       const {return magnetic_field_;}
    double d_perp_e()             const {return d_perp_e_;}
    double d_paral_e()            const {return d_paral_e_;}
    double *ptr_d_perp_e()        const {return ptr_d_perp_e_;}
    double *ptr_d_paral_e()       const {return ptr_d_paral_e_;}
    double time_cooler()          const {return time_cooler_;}
    bool do_test()                const {return do_test_;} 
    
    int set_park_temperature_eff(double x) {park_temperature_eff_ = x; return 0;}
    int set_magnetic_field(double x)       {magnetic_field_=x; return 0;}
    int set_d_perp_e(double x)             {d_perp_e_ = x; return 0;}
    int set_d_paral_e(double x)            {d_paral_e_ = x; return 0;}
    int set_ptr_d_perp_e(double* x)        {ptr_d_perp_e_ = x; return 0;}
    int set_ptr_d_paral_e(double* x)       {ptr_d_paral_e_ = x; return 0;}
    int set_ptr_d_perp_e(std::vector<double>& x) {ptr_d_perp_e_ = &*x.begin(); return 0;}
    int set_ptr_d_paral_e(std::vector<double>& x){ptr_d_paral_e_ = &*x.begin(); return 0;}
    int set_time_cooler(double x)          {time_cooler_ = x; return 0;}
    int set_do_test(bool b)                {do_test_ = b; return 0;}
    int set_filename(const char* f)        {test_filename_=f; return 0;}
    int set_cutoff(double c)               {cutoff_ = c; return 0;}
    int set_calls(size_t c)                {calls_ = c; return 0;} 
    //For Erlangen model    
    int set_fast(bool k)                   {fast_ = k; return 0;}
    int set_tight(bool k)                  {tight_ = k; return 0;}
    int set_stretched(bool k)              {stretched_ = k; return 0;}
    //For the un-magnetized model
    int set_approximate(bool k)            {approximate_ = k; return 0;}
    int set_binney(bool k)                 {binney_ = k; return 0;}

    //Constructors
    ForceParas(ForceFormula formula):formula_(formula){};
    ForceParas(ForceFormula formula,double f_const,const char* test_filename):
                formula_(formula),
                f_const_(f_const),
                test_filename_(test_filename){};
    
    //Thse are the wrappers for OpenMP parallel loops
    int ApplyForce(int charge_number, unsigned long int ion_number, 
                   double *v_tr, double *v_long, double *density_e,
                   double temperature, double magnetic_field, 
                   double *d_perp_e, double *d_paral_e, double time_cooler,
                   double *force_tr, double *force_long);
    
    int ApplyForce(int charge_number, unsigned long int ion_number, 
                   double *v_tr, double *v_long, double *density_e,
                   double temperature, double magnetic_field, 
                   double d_perp_e, double d_paral_e, double time_cooler,
                   double *force_tr, double *force_long, bool do_test);
<<<<<<< HEAD
    
    //This is a wrapper for the multidimensional integrals that show up 
    // in these force calculations
    void EvalIntegral(double (*func)(double*, size_t, void*), int_info &params,
                          double *xl, double *xu, size_t dim, 
                          double &result, double &error);
    
=======
    
    //This is a wrapper for the multidimensional integrals that show up 
    // in these force calculations
    void EvalIntegral(double (*func)(double*, size_t, void*), int_info &params,
                          double *xl, double *xu, size_t dim, 
                          double &result, double &error);
    
>>>>>>> 8947a02860fdcba41633e419a1ef2462ae623315
    //A 1d version of the eval integral function
    void EvalIntegral(double (*func)(double, void*), int_info &params,
                          double xl, double xu, double &result, double &error);

    //A calculation common to most (all?) friction force models
    double max_impact_factor(double v_dlt, int charge_number,
                             double density_e,double time_cooler);
    
    //This virtual function will be overloaded by each individual force model
    virtual void force(double v_tr, double v_long, double d_perp_e, 
                       double d_paral_e, double temperature, int charge_number,
                       double density_e,double time_cooler,double magnetic_field,
                       double &result_trans, double &result_long) = 0;      
};

class Force_Parkhomchuk : public ForceParas{
    private:
    //Force-dependent constants
        //double f_const = -4 * k_me_kg * pow(k_re*k_c*k_c,2); //for Parkhomchuk
        //const char* test_filename = "Parkhomchuk.txt";
    public:
        Force_Parkhomchuk():ForceParas(ForceFormula::PARKHOMCHUK,-4 * k_me_kg * pow(k_re*k_c*k_c,2),"Parkhomchuk.txt"){};

        virtual void force(double v_tr, double v_long, double d_perp_e, 
                           double d_paral_e, double temperature, int charge_number,
                           double density_e,double time_cooler,double magnetic_field,
                           double &force_result_trans, double &force_result_long);
};

class Force_DS : public ForceParas{
    private:
        //The factor of (0.5*pi/(2sqrt(2pi))) comes from the difference between the constants
        // used in Parkhomchuk with the constants used in the Pestrikov D&S integrals
        //const double f_const = -2 * (k_pi/(2*sqrt(2*k_pi))) * k_me_kg * pow(k_re*k_c*k_c,2);
        //const char* test_filename = "DerbenevSkrinsky.txt";
    
        //These must be static so they can be passed to GSL integration
        static double trans_integrand(double alpha,void *params);
        static double long_integrand(double alpha, void *params);
        
    public:
        Force_DS():ForceParas(ForceFormula::DERBENEVSKRINSKY,-sqrt(2*k_pi) * k_me_kg * pow(k_re*k_c*k_c,2),"DerbenevSkrinsky.txt"){};
        virtual void force(double v_tr, double v_long, double d_perp_e, 
                           double d_paral_e, double temperature, int charge_number,
                           double density_e,double time_cooler,double magnetic_field,
                           double &force_result_trans, double &force_result_long);
};

class Force_Meshkov : public ForceParas{
    private:
        double k_ = 2; //A smoothing factor between regions. 2 is the default in betacool
    public:
        Force_Meshkov():ForceParas(ForceFormula::MESHKOV,-2 * k_pi * k_me_kg * pow(k_re*k_c*k_c,2),"Meshkov.txt"){};
    
        int set_k(double k){k_ = k; return 0;}
        virtual void force(double v_tr, double v_long, double d_perp_e, 
                           double d_paral_e, double temperature, int charge_number,
                           double density_e,double time_cooler,double magnetic_field, 
                           double &force_result_trans, double &force_result_long);
};

class Force_Unmagnetized : public ForceParas{    
    public:    
        Force_Unmagnetized():ForceParas(ForceFormula::UNMAGNETIZED, -k_me_kg * pow(k_re*k_c*k_c,2), "Unmagnetized.txt"){};
        static double normalization_factor(double *k, size_t dim, void *params);    
        static double trans_integrand(double *k, size_t dim, void *params);
        static double long_integrand(double *k, size_t dim, void *params);
        static double full_trans_integrand(double *k, size_t dim, void *params);
        static double full_long_integrand(double *k, size_t dim, void *params);
        static double Binney_trans(double alpha, void *params);
        static double Binney_long(double alpha, void *params);
    
    virtual void force(double v_tr, double v_long, double d_perp_e,
                       double d_paral_e, double temperature, int charge_number,
                       double density_e,double time_cooler,double magnetic_field,
                       double &force_result_trans, double &force_result_long); 
};

class Force_Budker : public ForceParas{
    private:
    //Force-dependent constants
        //The factor of pi comes from the difference between the constants
        // used in Parkhomchuk with the constants used in the Meshkov representation
  
    public:
        Force_Budker():ForceParas(ForceFormula::BUDKER,-4 * k_pi * k_me_kg * pow(k_re*k_c*k_c,2),"Budker.txt"){};
    
        virtual void force(double v_tr, double v_long, double d_perp_e,
                           double d_paral_e, double temperature, int charge_number,
                           double density_e,double time_cooler,double magnetic_field,
                           double &force_result_trans, double &force_result_long);
};

<<<<<<< HEAD

=======
class Force_Erlangen : public ForceParas{
    private:
        //Definitions of integrals that need to be evaluated through Monte Carlo:
        static double fast_trans(double *k, size_t dim, void *params);
        static double fast_long(double *k, size_t dim, void *params);
        static double Alt_fast_trans(double *k, size_t dim, void *params);
        static double Alt_fast_long(double *k, size_t dim, void *params);
        static double tight_trans(double *k, size_t dim, void *params);
        static double tight_long(double *k, size_t dim, void *params);
        static double stretched_trans(double *k, size_t dim, void *params);
        static double stretched_long(double *k, size_t dim, void *params); 
 
    public:
        Force_Erlangen():ForceParas(ForceFormula::ERLANGEN, -k_me_kg * pow(k_re*k_c*k_c,2),"Erlangen.txt"){};

        int set_fast(bool k)      {fast_ = k; return 0;}
        int set_tight(bool k)     {tight_ = k; return 0;}
        int set_stretched(bool k) {stretched_ = k; return 0;}    
    
         virtual void force(double v_tr, double v_long, double d_perp_e,
                            double d_paral_e, double temperature, int charge_number,
                            double density_e, double time_cooler, double magnetic_field,
                            double &force_result_trans, double &force_result_long);
};
>>>>>>> 8947a02860fdcba41633e419a1ef2462ae623315

//A frequently-used case structure for switching
int friction_force(int charge_number, unsigned long int ion_number, double *v_tr, double *v_z, double *density_e,
                  ForceParas &force_paras, double *force_tr, double *force_long);

ForceParas* ChooseForce(ForceFormula force_formula);

#endif // FORCE_H
