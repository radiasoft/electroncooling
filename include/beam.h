#ifndef BEAM_H
#define BEAM_H

#include "constants.h"
#include <cstdio>
#include <memory>
#include <string>
#include <vector>
#include "beamshape.h"



enum class Velocity {CONST, USER_DEFINE, SPACE_CHARGE}  ;
enum class Temperature {CONST, USER_DEFINE, SPACE_CHARGE}  ;

class BeamShape;

class Beam{
    int charge_number_;   //Number of charges
    double mass_number_;       //mass = A*u [MeV/c^2]
    double mass_;    //unit in MeV/c^2
    double r_;       //classical radius, in m
    double kinetic_energy_;      //kinetic energy, in MeV
    double beta_;    //Lorentz factors
    double gamma_;   //Lorentz factors
    double emit_nx_; //normalized horizontal emittance, in m
    double emit_ny_; //normalized vertical emittance, in m
    double emit_x_;  //geometrical horizontal emittance, in m
    double emit_y_;  //geometrical vertical emittance, in m
    double dp_p_;     //momentum spread dp/p
    double energy_spread_;       // dE/E
    double sigma_s_; //RMS bunch length. set it to -1 for coasting beam, in m
    double particle_number_; //number of particles
    bool bunched_;   //Return true if beam is bunched.
    double center_[3] = {0,0,0};

public:
    int set_emit_nx(double x){emit_nx_ = x; emit_x_ = emit_nx_/(beta_*gamma_); return 0;}
    int set_emit_ny(double x){emit_ny_ = x; emit_y_ = emit_ny_/(beta_*gamma_); return 0;}
    int set_emit_x(double x){emit_x_ = x; emit_nx_ = beta_*gamma_*emit_x_; return 0;}
    int set_emit_y(double x){emit_y_ = x; emit_ny_ = beta_*gamma_*emit_y_; return 0;}
    int set_dp_p(double x){dp_p_ = x; energy_spread_ = beta_*beta_*dp_p_; return 0;}
    int set_sigma_s(double x){sigma_s_ = x; return 0;}
    int set_center(double cx, double cy, double cz){center_[0] = cx; center_[1] = cy; center_[2] = cz; return 0;}
    int set_center(int i, double x);

    int set_bunched(bool s){bunched_ = s; return 0;}
    int charge_number() const {return charge_number_;}
    double mass() const {return mass_;}
    double kinetic_energy() const {return kinetic_energy_;}
    double beta() const {return beta_;}
    double gamma() const {return gamma_;}
    double emit_nx() const {return emit_nx_;}
    double emit_ny() const {return emit_ny_;}
    double emit_x() const {return emit_x_;}
    double emit_y() const {return emit_y_;}
    double dp_p() const {return dp_p_;}
    double energy_spread() const {return energy_spread_;}
    double sigma_s() const {return sigma_s_;}
    double r() const {return r_;}
    double particle_number() const {return particle_number_;}
    double mass_number() const {return mass_number_;}
    double mass_J() const {return mass_*1e6*k_e;}
    bool bunched()const {return bunched_;}

    int center(double &cx, double &cy, double &cz){cx = center_[0]; cy = center_[1]; cz = center_[2]; return 0;}
    double center(int i){ if (i<3) return center_[i]; else perror("Error index for electron beam center!"); return 1.0;}
    Beam(int charge_number, double mass_number, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
        double sigma_s, double n_particle);
    Beam(int charge_number, double mass_number, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
        double n_particle);
};



class EBeam:public Beam{
    double tmp_tr_;            //Transverse temperature, in eV
    double tmp_long_;          //Longitudinal temperature, in eV
    double v_rms_tr_;        //Transverse RMS velocity, in m/s
    double v_rms_long_;      //Longitudinal RMS velocity, in m/s
    Velocity velocity_ = Velocity::CONST;
    Temperature temperature_ = Temperature::CONST;

 public:
    BeamShape *shape_;            //Shape of the electron beam
    double tmp_tr(){return tmp_tr_;}
    double tmp_long(){return tmp_long_;}
    double v_rms_tr(){return v_rms_tr_;}
    double v_rms_long(){return v_rms_long_;}
    void set_velocity(Velocity velocity){velocity_ = velocity;}
    void set_temperature(Temperature temp){temperature_ = temp;}
    Velocity velocity(){return velocity_;}
    Temperature temperature(){return temperature_;}

    int emit_nx(){perror("This function is not defined for cooling electron beam"); return 1;}
    int emit_ny(){perror("This function is not defined for cooling electron beam"); return 1;}
    int emit_x(){perror("This function is not defined for cooling electron beam"); return 1;}
    int emit_y(){perror("This function is not defined for cooling electron beam"); return 1;}
    int dp_p(){perror("This function is not defined for cooling electron beam"); return 1;}
    int sigma_s(){perror("This function is not defined for cooling electron beam"); return 1;}
    int n_particle(){perror("This function is not defined for cooling electron beam"); return 1;}
    bool bunched();
    double length();


    EBeam(double gamma, double tmp_tr, double tmp_long, BeamShape &shape_defined);
    EBeam(double gamma, BeamShape &shape_defined);
};
#endif // BEAM_H
