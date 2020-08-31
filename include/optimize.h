#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <algorithm>
#include <chrono>
#include <fstream>
#include <cmath>
#include "dynamic.h"
#include "ecooling.h"
#include "ibs.h"
#include "ring.h"
#include "gsl/gsl_multimin.h"
#include <random>
#include <string>
#include <map>
#include <iomanip>

//This object stores all optimizable variables
// as doubles. Seems like a duplication of 
// storage in this class, but I don't know a better
// way at this moment.

//We also need to store the step size and RMS for the 
// Gaussian search

class Optimize{
    
    std::vector<std::string> FitVariables;
    std::map<const std::string, double> FitStepSize;
    std::map<const std::string, double> BestFits;
    std::map<const std::string, double> InitialValues;
    double best_fval = DBL_MAX; //Any fit will decrease this

    static double fit_fcn(const gsl_vector *v, void *params);
    void Randomize();
        
    int n_trials = 15;
         
    
  protected:
     //Store some parameters we don't scan over,
     // including ion beam info.
     // It's easier to define these things in a struct
     // because the whole thing can be passed into the 
     // static function for the fitter.
     struct opt_info{
//         int Z_                 = 1;
//         double m0_             = 938.272;
         double magnetic_field_ = 1.0;    // in Tesla
         double length_         = 130.0;  // in meter
         double section_number_ = 1;
         double beta_h_         = 10.0;   // 10.0 in meter
         double beta_v_         = 10.0;   // in meter
         double disp_h_         = 0.0;    // 0.1, in meter
         double disp_v_         = 0.0;
         double alpha_h_        = 0.0;
         double alpha_v_        = 0.0;
         double disp_der_h_     = 0.0;
         double disp_der_v_     = 0.0;
         double sigma_x_        = 1e-4;
         double sigma_y_        = 1e-4;
         double sigma_s_        = 0.15;
         double temp_tr_        = 0.01;
         double temp_long_      = 0.01;
         double n_electron_     = 1.5; // in 1e10
    
         double cool_target_    = 20.; //The time in minutes to aim for
         
         int n_sample           = 1e6;
         
         
         std::string lattice_filename = "eRHIC.tfs";    
         //ForceFormula ff = ForceFormula::PARKHOMCHUK;
        ForceFormula ff = ForceFormula::BUDKER;
//        ForceFormula ff = ForceFormula::UNMAGNETIZED;
         
         Lattice *lattice;
         ForceParas *force_paras;
         EcoolRateParas *ecool_paras;
         Beam *beam;
         
         std::vector<std::string> FitVariables_working;
         std::map<const std::string, double> FitStepSize_working;
         std::map<const std::string, double> BestFits_working;
         std::map<const std::string, double> InitialValues_working;
         
        } fitter_values;
    
    //opt_info fitter_values;
    
    public:

        void InitializeFitter(std::vector<std::string>, 
                              std::vector<double>, 
                              Lattice*, 
                              Beam*, 
                              ForceFormula ff);
    
        
        void OptimizeTrial();
        void ManyTrials();
        //A function to access the optimization from the UI 
        int Optimize_From_UI(std::vector<std::string> Params, 
                             std::vector<double> InitialValues, 
                             Beam &ion, 
                             Cooler &cooler, 
                             EBeam &ebeam, 
                             Ring &ring, 
                             ForceFormula ff);
    
        Optimize(){};
    
    
    
};
#endif // OPTIMIZE_H