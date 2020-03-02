#ifndef BEAMSHAPE_H
#define BEAMSHAPE_H

#include "constants.h"
#include <cstdio>
#include <memory>
#include <string>
#include <vector>
#include "arbitrary_electron_beam.h"
#include "beam.h"

class Beam;

enum class Shape {UNIFORM_CYLINDER, GAUSSIAN_BUNCH, UNIFORM_BUNCH, GAUSSIAN_CYLINDER, ELLIPTIC_UNIFORM_BUNCH,
    UNIFORM_HOLLOW, UNIFORM_HOLLOW_BUNCH, PARTICLE_BUNCH};

class BeamShape{
 public:
    //Calculate the charge density for a given position (x,y,z) in Lab frame.
    virtual int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n)=0;
    //(cx, cy, cz) is the relative position of the ion beam center to the electron beam center.
    virtual int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx,
                            double cy, double cz)=0;
    virtual Shape shape()=0;
    virtual bool bunched() = 0;
    virtual double length()=0; //For bunched electron beam, return full length of the electron bunch.
    virtual double neutralisation()=0;
    BeamShape(){};
};

class UniformCylinder: public BeamShape{
    double current_;                   //Current of the beam in A
    double radius_;              //Radius of the beam in meter
    double neutralisation_;
 public:
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle);
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx, double cy,
                double cz);
    double current(){return current_;}
    double radius(){return radius_;}
    double neutralisation(){return neutralisation_;}
    Shape shape(){return Shape::UNIFORM_CYLINDER;}
    double length(){perror("length() not defined for UniformCylinder, which is coasting"); return 0;}
    bool bunched(){return false;}
    UniformCylinder(double current, double radius, double neutralisation=2):current_(current),radius_(radius),
                    neutralisation_(neutralisation){};
};

class UniformHollow: public BeamShape {
    double current_;    //Peak current, the current as if the beam is coasting.
    double in_radius_;
    double out_radius_;
    double neutralisation_;
 public:
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle);
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx, double cy,
                double cz);
    double current(){return current_;}
    double out_radius(){return out_radius_;}
    double in_radius(){return in_radius_;}
    double neutralisation(){return neutralisation_;}
    Shape shape(){return Shape::UNIFORM_HOLLOW;}
    double length(){perror("length() not defined for UniformHollow, which is coasting"); return 0;}
    bool bunched(){return false;}
    UniformHollow(double current, double in_radius, double out_radius, double neutralisation=2):current_(current),
        in_radius_(in_radius), out_radius_(out_radius),neutralisation_(neutralisation){};
};

class UniformHollowBunch: public BeamShape {
    double current_;
    double in_radius_;
    double out_radius_;
    double neutralisation_;
    double length_;
 public:
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle);
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx, double cy,
                double cz);
    double current(){return current_;}
    double out_radius(){return out_radius_;}
    double in_radius(){return in_radius_;}
    double neutralisation(){return neutralisation_;}
    Shape shape(){return Shape::UNIFORM_HOLLOW_BUNCH;}
    double length(){return length_;}
    bool bunched(){return true;}
    UniformHollowBunch(double current, double in_radius, double out_radius, double length, double neutralisation=2):current_(current),
        in_radius_(in_radius), out_radius_(out_radius), neutralisation_(neutralisation), length_(length) {}
};

class GaussianBunch: public BeamShape{
    double n_electron_;
    double sigma_x_;
    double sigma_y_;
    double sigma_s_;
    double neutralisation_;
 public:
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle);
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx, double cy,
                double cz);
    Shape shape(){return Shape::GAUSSIAN_BUNCH;}
    double length(){return 6*sigma_s_;}
    bool bunched(){return true;}
    double neutralisation(){return neutralisation_;}
    GaussianBunch(double n_electron, double sigma_x, double sigma_y, double sigma_s):n_electron_(n_electron),
                sigma_x_(sigma_x),sigma_y_(sigma_y),sigma_s_(sigma_s){};

};

class UniformBunch: public BeamShape{
    double current_;                   //Current of the beam in A, assuming the beam is DC.
    double radius_;              //Radius of the beam in meter
    double length_;
    double neutralisation_;
public:
    //Calculate the charge density for a given position (x,y,z) in Lab frame.
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n);
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx, double cy, double cz);
    Shape shape(){return Shape::UNIFORM_BUNCH;}
    double length(){return length_;}
    bool bunched(){return true;}
    double current(){return current_;}
    double radius(){return radius_;}
    double neutralisation(){return neutralisation_;}
//    UniformCylinder(double I, double radius, double neutralisation):Shape(ShapeList::uniformCylinder),I(I),radius(radius),neutralisation(neutralisation){};
    UniformBunch(double current, double radius, double length, double neutralisation=2):current_(current),radius_(radius),
            length_(length), neutralisation_(neutralisation){};

};

class EllipticUniformBunch: public BeamShape{
    double current_;
    double rh_;         //half horizontal axis
    double rv_;         //half vertical axis
    double length_;     //bunch length
    double neutralisation_;
public:
    //Calculate the charge density for a given position (x,y,z) in Lab frame.
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n);
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx, double cy, double cz);
    Shape shape(){return Shape::ELLIPTIC_UNIFORM_BUNCH;}
    double length(){return length_;}
    bool bunched(){return true;}
    double neutralisation(){return neutralisation_;}
    EllipticUniformBunch(double current, double rh, double rv, double length, double neutralisation=2):current_(current),
            rh_(rh),rv_(rv),length_(length),neutralisation_(neutralisation){};
};

class ParticleBunch: public BeamShape {
    double n_electron_;
    std::string filename_;
    unsigned long int n_ = 0;
    double length_ = 0;
    bool v_x_corr_ = false;    //Velocity position correlation
    double neutralisation_;
    int line_skip_ = 0;
    vector<Box> tree_;
    vector<unsigned long int> list_e_;
    int s_ = 100;
    bool binary_ = false;
    int buffer_ = 1000;
public:
    std::vector<double> x, y, z, vx, vy, vz;  //Electron phase space coordinates
    std::vector<double> v_avg_z, v_rms_l, v_rms_t, tpr_l, tpr_t;  //Velocity and temperate w.r.t. ions.
    //Calculate the charge density for a given position (x,y,z) in Lab frame.
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n);
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx, double cy, double cz);
    Shape shape(){return Shape::PARTICLE_BUNCH;}
    double length(){return length_;}
    bool bunched(){return true;}
    bool corr(){return v_x_corr_;}
    void set_corr(bool corr = true){v_x_corr_ = corr;}
    void set_buffer(int n) {buffer_ = n;}
    void set_s(int s) {s_ = s;}
    void set_binary(bool b) {binary_ = b;}
    void set_skip(int n) {line_skip_ = n;}
    void set_neutralisation(double x) { neutralisation_ = x;}
    double neutralisation(){return neutralisation_;}

    ParticleBunch(double n_electron, std::string filename, double length):n_electron_(n_electron),
        filename_(filename),length_(length){};
    ParticleBunch(double n_electron, std::string filename):n_electron_(n_electron),
        filename_(filename){};
    void load_particle(long int n);
    void load_particle();

};

#endif // BEAMSHAPE_H
