#ifndef PARTICLE_MODEL_HPP
#define PARTICLE_MODEL_HPP

#include <vector>
#include "beam.h"
#include "cooler.h"
#include "ring.h"
#include "fit.h"

void initialize_particle_model(Beam &ion);
void referece_twiss(Cooler &cooler);
void particle_model_ibs_scratches();
void ibs_kick(unsigned int n_sample, double rate, double twiss, double dt, double emit, double* p);
void ibs_kick(Beam &ion, double rate1, double rate2, double twiss, double dt, double *x, double* p, fit_results &FR);
void restore_cord(double t_cooler);
void adjust_freq(double &freq, EBeam ebeam);
void apply_cooling_kick(double t_cooler, double freq, double dt);
void apply_ibs_kick(double dt, Beam &ion, std::vector<double> &r_ibs);
void apply_ibs_kick_DG(double dt, Beam &ion, std::vector<double> &r_ibs);
void move_particles(Beam &ion, Ring &ring);
void update_beam_parameters(Beam &ion,Ring &ring);

#endif // IBS_HPP
