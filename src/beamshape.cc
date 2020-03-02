#include "beamshape.h"
#include <cmath>
#include <cstring>
#include "arbitrary_electron_beam.h"

#include <fstream>

int UniformCylinder::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle){
    int nq = ebeam.charge_number();
    if (nq<0) nq *= -1;
    double r2 = radius_*radius_;
    double density = current_/(k_pi*r2*nq*k_e*ebeam.beta()*k_c);
    memset(ne, 0, n_particle*sizeof(double));
    
    #pragma omp parallel for
    for(unsigned int i=0; i<n_particle; ++i){
        if(x[i]*x[i]+y[i]*y[i]<=r2) ne[i] = density;
    }
    return 0;
}

int UniformCylinder::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx,
                             double cy, double cz){
    int nq = ebeam.charge_number();
    if (nq<0) nq *= -1;
    double r2 = radius_*radius_;
    double density = current_/(k_pi*r2*nq*k_e*ebeam.beta()*k_c);
    memset(ne, 0, n_particle*sizeof(double));
    //ion_center - electron_center
    cx -= ebeam.center(0);
    cy -= ebeam.center(1);
    cz -= ebeam.center(2);
    
    #pragma omp parallel for
    for(unsigned int i=0; i<n_particle; ++i){
        if((x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy)<=r2) ne[i] = density;
    }
    return 0;
}

int UniformHollow::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle) {
    int nq = ebeam.charge_number();
    double out_r2 = out_radius_*out_radius_;
    double in_r2 = in_radius_*in_radius_;
    double area = k_pi*(out_r2-in_r2);
    double density;
    if (area!=0)
        density = current_/(area*nq*k_e*ebeam.beta()*k_c);
    else
        density = 0;
    if (density<0) density *= -1;
    memset(ne, 0, n_particle*sizeof(double));

    #pragma omp parallel for
    for(unsigned int i=0; i<n_particle; ++i){
        double r2 = x[i]*x[i]+y[i]*y[i];
        if(r2<=out_r2 && r2>=in_r2) ne[i] = density;
    }
    return 0;
}

int UniformHollow::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx,
                             double cy, double cz) {
    int nq = ebeam.charge_number();
    double out_r2 = out_radius_*out_radius_;
    double in_r2 = in_radius_*in_radius_;
    double area = k_pi*(out_r2-in_r2);
    double density;
    if (area!=0)
        density = current_/(area*nq*k_e*ebeam.beta()*k_c);
    else
        density = 0;
    if (density<0) density *= -1;
    memset(ne, 0, n_particle*sizeof(double));
     //ion_center - electron_center
    cx -= ebeam.center(0);
    cy -= ebeam.center(1);
    cz -= ebeam.center(2);
    
    #pragma omp parallel for
    for(unsigned int i=0; i<n_particle; ++i){
        double r2 = (x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy);
        if(r2<=out_r2 && r2>=in_r2) ne[i] = density;
    }
    return 0;
}



int UniformBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle){

    int nq = ebeam.charge_number();
    if (nq<0) nq *= -1;
    double r2 = radius_*radius_;
    double density = current_/(k_pi*r2*nq*k_e*ebeam.beta()*k_c);
    memset(ne, 0, n_particle*sizeof(double));
    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    
    #pragma omp parallel for
    for(unsigned int i=0; i<n_particle; ++i){
        if(z[i]<=right_end && z[i]>= left_end && x[i]*x[i]+y[i]*y[i]<=r2)
            ne[i] = density;
    }
    return 0;
}

int UniformBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx,
                          double cy, double cz){
    int nq = ebeam.charge_number();
    if (nq<0) nq *= -1;
    double r2 = radius_*radius_;
    double density = current_/(k_pi*r2*nq*k_e*ebeam.beta()*k_c);

    //ion_center - electron_center
    cx -= ebeam.center(0);
    cy -= ebeam.center(1);
    cz -= ebeam.center(2);
    memset(ne, 0, n_particle*sizeof(double));
    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    
    #pragma omp parallel for
    for(unsigned int i=0; i<n_particle; ++i){
        if((z[i]+cz)<=right_end && (z[i]+cz)>=left_end && (x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy)<=r2)
            ne[i] = density;
    }
    return 0;
}

int UniformHollowBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle){

    int nq = ebeam.charge_number();
    double out_r2 = out_radius_*out_radius_;
    double in_r2 = in_radius_*in_radius_;
    double area = k_pi*(out_r2-in_r2);
    double density;
    if (area!=0)
        density = current_/(area*nq*k_e*ebeam.beta()*k_c);
    else
        density = 0;
    if (density<0) density *= -1;
    memset(ne, 0, n_particle*sizeof(double));

    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    
    #pragma omp parallel for
    for(unsigned int i=0; i<n_particle; ++i){
        double r2 = x[i]*x[i]+y[i]*y[i];
        if(z[i]<=right_end && z[i]>=left_end &&r2<=out_r2 && r2>=in_r2) ne[i] = density;
    }
    return 0;
}

int UniformHollowBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx,
                          double cy, double cz){

    int nq = ebeam.charge_number();
    double out_r2 = out_radius_*out_radius_;
    double in_r2 = in_radius_*in_radius_;
    double area = k_pi*(out_r2-in_r2);
    double density;
    if (area!=0)
        density = current_/(area*nq*k_e*ebeam.beta()*k_c);
    else
        density = 0;
    if (density<0) density *= -1;
    //ion_center - electron_center
    cx -= ebeam.center(0);
    cy -= ebeam.center(1);
    cz -= ebeam.center(2);
    memset(ne, 0, n_particle*sizeof(double));

    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    
    #pragma omp parallel for
    for(unsigned int i=0; i<n_particle; ++i){
        double r2 = (x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy);
        double z_shifted = z[i]+cz;
        if(z_shifted<=right_end && z_shifted>=left_end &&r2<=out_r2 && r2>=in_r2) ne[i] = density;
    }
    return 0;
}



int EllipticUniformBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle){

    int nq = ebeam.charge_number();
    if (nq<0) nq *= -1;
    double density = current_/(k_pi*rh_*rv_*nq*k_e*ebeam.beta()*k_c);
    memset(ne, 0, n_particle*sizeof(double));
    double inv_rh2 = 1.0/(rh_*rh_);
    double inv_rv2 = 1.0/(rv_*rv_);
    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    
    #pragma omp parallel for
    for(unsigned int i=0; i<n_particle; ++i){
        if(z[i]<=right_end && z[i]>=left_end && inv_rh2*x[i]*x[i]+inv_rv2*y[i]*y[i]<=1)
            ne[i] = density;
    }
    return 0;
}

int EllipticUniformBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle,
                                  double cx, double cy, double cz){
    int nq = ebeam.charge_number();
    if (nq<0) nq *= -1;
    double density = current_/(k_pi*rh_*rv_*nq*k_e*ebeam.beta()*k_c);

    //ion_center - electron_center
    cx -= ebeam.center(0);
    cy -= ebeam.center(1);
    cz -= ebeam.center(2);
    memset(ne, 0, n_particle*sizeof(double));
    double inv_rh2 = 1.0/(rh_*rh_);
    double inv_rv2 = 1.0/(rv_*rv_);
    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    
    #pragma omp parallel for
    for(unsigned int i=0; i<n_particle; ++i){
        if((z[i]+cz)<=right_end && (z[i]+cz)>=left_end &&
           inv_rh2*(x[i]+cx)*(x[i]+cx)+inv_rv2*(y[i]+cy)*(y[i]+cy)<=1)
            ne[i] = density;
    }
    return 0;
}

int GaussianBunch::density(double *x, double *y, double *z, Beam &beam, double *ne, unsigned int n_particle){
    double amp = n_electron_/(sqrt(8*k_pi*k_pi*k_pi)*sigma_x_*sigma_y_*sigma_s_);
    double sigma_x2 = -1/(2*sigma_x_*sigma_y_);
    double sigma_y2 = -1/(2*sigma_y_*sigma_y_);
    double sigma_s2 = -1/(2*sigma_s_*sigma_s_);
    
    #pragma omp parallel for
    for(unsigned int i=0; i<n_particle; ++i){
//        ne[i] = amp*exp(-0.5*(x[i]*x[i]/sigma_x2+y[i]*y[i]/sigma_y2+z[i]*z[i]/sigma_s2));
        ne[i] = amp*exp(x[i]*x[i]*sigma_x2+y[i]*y[i]*sigma_y2+z[i]*z[i]*sigma_s2);
    }
    return 0;
}

int GaussianBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx,
                           double cy, double cz){
    double amp = n_electron_/(sqrt(8*k_pi*k_pi*k_pi)*sigma_x_*sigma_y_*sigma_s_);
    double sigma_x2 = -1/(2*sigma_x_*sigma_y_);
    double sigma_y2 = -1/(2*sigma_y_*sigma_y_);
    double sigma_s2 = -1/(2*sigma_s_*sigma_s_);
    //ion_center - electron_center
    cx -= ebeam.center(0);
    cy -= ebeam.center(1);
    cz -= ebeam.center(2);
    
    #pragma omp parallel for
    for(unsigned int i=0; i<n_particle; ++i){
        ne[i] = amp*exp((x[i]+cx)*(x[i]+cx)*sigma_x2+(y[i]+cy)*(y[i]+cy)*sigma_y2+(z[i]+cz)*(z[i]+cz)*sigma_s2);
//        ne[i] = amp*exp(-0.5*((x[i]+cx)*(x[i]+cx)/sigma_x2+(y[i]+cy)*(y[i]+cy)/sigma_y2+(z[i]+cz)*(z[i]+cz)/sigma_s2));
    }
    return 0;
}


void ParticleBunch::load_particle(long int n) {
    n_ = load_electrons(x, y, z, vx, vy, vz, filename_, n, line_skip_, binary_, buffer_);
    create_e_tree(x, y, z, n_, s_, tree_, list_e_);
    if(length_==0) {
        auto itr = z.begin();
        double z_max = *itr;
        double z_min = *itr;
        ++itr;
        for(; itr!=z.end(); ++itr) {
            if(*itr>z_max) z_max = *itr;
            if(*itr<z_min) z_min = *itr;
        }
        length_ = z_max - z_min;
    }
}

void ParticleBunch::load_particle() {
    load_particle(0);
}

int ParticleBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n) {
    double rate = n_electron_/n_;
    std::vector<unsigned int> list_i;
    unsigned int idx_out;
    create_ion_tree(x, y, z, n, tree_, list_i, idx_out);
    if (v_x_corr_)
        {::density(tree_, list_e_, vx, vy, vz, n_, list_i, idx_out, n, ne, v_avg_z, v_rms_t, v_rms_l);}
    else
        {::density(tree_, list_e_, vx, vy, vz, n_, list_i, idx_out, n, ne, v_rms_t, v_rms_l);}
//    tpr_l.clear();
//    tpr_t.clear();
//    double c2_inv = 1.0/(k_c*k_c);
//    for(auto& v: v_rms_l) tpr_l.push_back(k_me*v*v*1e6*c2_inv);
//    for(auto& v: v_rms_t) tpr_t.push_back(k_me*v*v*1e6*c2_inv);

    for(unsigned int i=0; i<n; ++i) ne[i] *= rate;

    return 0;
}

int ParticleBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx, double cy,
                       double cz) {
    double rate = n_electron_/n_;
    for(unsigned int i=0; i<n; ++i) {
        x[i] -= cx;
        y[i] -= cy;
        z[i] -= cz;
    }
    std::vector<unsigned int> list_i;
    unsigned int idx_out;
    create_ion_tree(x, y, z, n, tree_, list_i, idx_out);
    if (v_x_corr_)
        {::density(tree_, list_e_, vx, vy, vz, n_, list_i, idx_out, n, ne, v_avg_z, v_rms_t, v_rms_l);}
    else
        {::density(tree_, list_e_, vx, vy, vz, n_, list_i, idx_out, n, ne, v_rms_t, v_rms_l);}
    for(unsigned int i=0; i<n; ++i) {
        x[i] += cx;
        y[i] += cy;
        z[i] += cz;
    }
    double c2_inv = 1.0/(k_c*k_c);
    for(auto& v: v_rms_l) tpr_l.push_back(k_me*v*v*1e6*c2_inv);
    for(auto& v: v_rms_t) tpr_t.push_back(k_me*v*v*1e6*c2_inv);
    for(unsigned int i=0; i<n; ++i) ne[i] *= rate;
    return 0;
}
