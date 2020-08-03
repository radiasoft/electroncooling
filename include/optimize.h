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

//This object stores all optimizable variables
// as doubles. Seems like a duplication of 
// storage in this class, but I don't know a better
// way at this moment.

//We also need to store the step size and RMS for the 
// Gaussian search

class Optimize{
    double length_;      // in meter
    double section_number_;
    double magnetic_field_; // 1.0, in Tesla
    double beta_h_;      // 10.0 in meter
    double beta_v_;      // in meter
    double disp_h_;      // 0.1, in meter
    double disp_v_;
    double alpha_h_;
    double alpha_v_;
    double der_disp_h_;
    double der_disp_v_;
    double sigma_x_; //1e-4
    double sigma_y_; //1e-4
    double sigma_s_; //0.15
    double temp_tr_; //0.01
    double temp_long_;
    double n_electron_; //1.5, in 1e10
    
    std::vector<std::string> FitVariables;
    std::map<const std::string, double> FitStepSize;
    std::map<const std::string, double> BestFits;
    double best_fval = DBL_MAX; //Any fit will decrease this

    static double fit_fcn(const gsl_vector *v, void *params);
    std::map<std::string, double> Randomize();
    
    
     protected:
        //Store some parameters we don't scan over,
        // including ion beam info
        struct opt_info{
            int Z = 1;
            double m0 = 938.272;
            std::string lattice_filename;    
        };
    
    
    public:
        double length()const {return length_;}
        double section_number()const {return section_number_;}
        double magnetic_field()const {return magnetic_field_;}
        double beta_h()const {return beta_h_;}
        double beta_v()const {return beta_v_;}
        double alpha_h()const {return alpha_h_;}
        double alpha_v()const {return alpha_v_;}
        double disp_h()const {return disp_h_;}
        double disp_v()const {return disp_v_;}
        double der_disp_h()const {return der_disp_h_;}
        double der_disp_v()const {return der_disp_v_;}
        double sigma_x()const {return sigma_x_;}
        double sigma_y()const {return sigma_y_;}
        double sigma_s()const {return sigma_s_;}
        double temp_tr()const {return temp_tr_;}
        double temp_long()const {return temp_long_;}
        
        void OptimizeTrial();
        void ManyTrials();
    
        Optimize(){};
    
    
    
};
#endif // OPTIMIZE_H