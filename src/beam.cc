
#include "beam.h"
#include <cmath>
#include <cstring>
#include "arbitrary_electron_beam.h"

#include <fstream>

Beam::Beam(int charge_number, double mass_number, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
           double sigma_s, double n_particle): charge_number_(charge_number), mass_number_(mass_number),
           kinetic_energy_(kinetic_energy), emit_nx_(emit_nx), emit_ny_(emit_ny), dp_p_(dp_p), sigma_s_(sigma_s),
           particle_number_(n_particle) {
    mass_ = mass_number*k_u;
    gamma_ = 1+kinetic_energy_/mass_;
    beta_ = sqrt(gamma_*gamma_-1)/gamma_;
    r_ = k_ke*charge_number_*charge_number_*k_e*1e-6/mass_;
    bunched_ = (sigma_s_>0)?true:false;
    emit_x_ = emit_nx_/(beta_*gamma_);
    emit_y_ = emit_ny_/(beta_*gamma_);
    energy_spread_ = beta_*beta_*dp_p_;
}

Beam::Beam(int charge_number, double mass_number, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
           double n_particle): charge_number_(charge_number), mass_number_(mass_number), kinetic_energy_(kinetic_energy),
           emit_nx_(emit_nx), emit_ny_(emit_ny), dp_p_(dp_p), particle_number_(n_particle) {
    mass_ = mass_number_*k_u;
    gamma_ = 1+kinetic_energy_/mass_;
    beta_ = sqrt(gamma_*gamma_-1)/gamma_;
    r_ = k_ke*charge_number_*charge_number_*k_e*1e-6/mass_;
    bunched_ = false;
    sigma_s_ = -1;
    emit_x_ = emit_nx_/(beta_*gamma_);
    emit_y_ = emit_ny_/(beta_*gamma_);
    energy_spread_ = beta_*beta_*dp_p_;
}

int Beam::set_center(int i, double x) {
    if(i<3) {
        center_[i] = x;
        return 0;
    }
    else {
        perror("Error index for electron beam center!");
        return 1;
    }
}


EBeam::EBeam(double gamma, double tmp_tr, double tmp_long, BeamShape &shape_defined):
    Beam(-1, k_me/k_u, (gamma-1)*k_me, 0, 0, 0, 0),tmp_tr_(tmp_tr),tmp_long_(tmp_long){
    shape_ = &shape_defined;
    v_rms_long_ = sqrt(tmp_long_/this->mass())*0.001*k_c;
    v_rms_tr_ = sqrt(tmp_tr_/this->mass())*0.001*k_c;
    if(shape_defined.shape()==Shape::PARTICLE_BUNCH) {
        velocity_ = Velocity::USER_DEFINE;
        temperature_ = Temperature::USER_DEFINE;
    }
}

EBeam::EBeam(double gamma, BeamShape &shape_defined):Beam(-1, k_me/k_u, (gamma-1)*k_me, 0, 0, 0, 0) {
    shape_ = &shape_defined;
    tmp_tr_ = 0;
    tmp_long_ = 0;
    v_rms_long_ = sqrt(tmp_long_/this->mass())*0.001*k_c;
    v_rms_tr_ = sqrt(tmp_tr_/this->mass())*0.001*k_c;
    if(shape_defined.shape()==Shape::PARTICLE_BUNCH) {
        velocity_ = Velocity::USER_DEFINE;
        temperature_ = Temperature::USER_DEFINE;
    }
}

bool EBeam::bunched(){return shape_->bunched();}
double EBeam::length(){return shape_->length();}

