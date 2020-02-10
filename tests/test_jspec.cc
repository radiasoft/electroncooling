#include <math.h>
#include <Base.hh>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "force.h"
#include "ecooling.h"
#include "ibs.h"
#include "dynamic.h"
#include "cooler.h"

using std::string;

extern DynamicParas * dynamic_paras;
extern EcoolRateParas * ecool_paras;
extern ForceParas * force_paras;
extern IBSSolver * ibs_solver;

int both(){

    JSPEC_TEST_BEGIN("Electron cooling + IBS:");
    // define proton beam;
    double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
    int Z;
    Z = 1;
    m0 = 938.272;
    KE = 100e3;
    emit_nx0 = 1.2e-6;
    emit_ny0 = 0.6e-6;
    dp_p0 = 5e-4;
    sigma_s0 = 0.84e-2;

    N_ptcl = 6.56E9;
    Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl);

    // define the lattice of the proton ring
    //std::string filename = "MEICColliderRedesign1IP.tfs";
    std::string filename = CMAKE_SOURCE_DIR + std::string("/data/Booster.tfs");
    Lattice lattice(filename);

    //Define the ring
    Ring ring(lattice, p_beam);

    //Set IBS parameters.
    int nu = 100;
    int nv = 100;
    int nz = 40;
    double log_c = 39.9/2;
    double k = 0.5;
    ibs_solver = new IBSSolver_Martini(nu, nv, nz, log_c, k);

//  //Calculate IBS rate.
//
//  double rx_ibs, ry_ibs, rz_ibs;
//  ibs_rate(lattice, p_beam, *ibs_solver, rx_ibs, ry_ibs, rz_ibs);
//  std::cout<<"ibs rate: "<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;


    //define the cooler
    double cooler_length = 60;
    double n_section = 1;
    double magnetic_field = 1;
    double beta_h = 100;
    double beta_v = 100;
    double dis_h = 0;
//  double dis_h = 2.5;
    double dis_v = 0;
    Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);

////
    //define electron beam
//  double n_electron = 2.62E9;
//  double sigma_x = 0.035E-2;
//  double sigma_y = 0.035E-2;
//  double sigma_s = 0.84E-2;
//  GaussianBunch gaussian_bunch(n_electron, sigma_x, sigma_y, sigma_s);
//  double gamma_e = p_beam.gamma();
//  double tmp_tr = 0.1;
//  double tmp_long = 0.5;
//  EBeam e_beam(gamma_e, tmp_tr, tmp_long, gaussian_bunch);

    std::string electron_file = "electrons.dat";
    double ns = 1e6;
    double n_electron = 2.62E9;
    int s = 200;
    int line_skip = 0;
    ParticleBunch particle_bunch(n_electron, electron_file);
    particle_bunch.set_s(s);
    particle_bunch.set_skip(line_skip);
    particle_bunch.set_binary(true);
    particle_bunch.load_particle(ns);

//                ParticleBunch particle_bunch(n_electron, electron_file, ns, line_skip, false, 1000, s);
    double gamma_e = p_beam.gamma();
    EBeam e_beam(gamma_e, particle_bunch);

    unsigned int n_sample = 10000;
    ecool_paras = new EcoolRateParas(n_sample);
    //define friction force formula
    force_paras = ChooseForce(ForceFormula::PARKHOMCHUK);

    double rate_x, rate_y, rate_s;
    ecooling_rate(*ecool_paras, *force_paras, p_beam, cooler, e_beam, ring, rate_x, rate_y, rate_s);
    std::cout<<std::endl;
    std::cout<<"rate_x = "<<rate_x<<" rate_y = "<<rate_y<<" rate_s = "<<rate_s<<std::endl;

//                rate_x += rx_ibs;
//                rate_y += ry_ibs;
//                rate_s += rz_ibs;
//                std::cout<<"rate_x = "<<rate_x<<" rate_y = "<<rate_y<<" rate_s = "<<rate_s<<std::endl;



    //            double t = 7200;
    //            int n_step = 720;
    //            bool ibs = true;
    //            bool ecool = true;
    //            dynamic_paras = new DynamicParas(t, n_step, ibs, ecool);
    //            dynamic_paras->set_model(DynamicModel::MODEL_BEAM);
    //
    //            char file[100] = "tr-0.5eV.txt";
    //            std::ofstream outfile;
    //            outfile.open(file);
    //            dynamic(p_beam, cooler, e_beam, ring, outfile);
    ////            dynamic(p_beam, cooler, e_beam, ring, outfile);
    //            outfile.close();

    JSPEC_TEST_END();

    return 0;
}


int ecool(ForceFormula ff, double *check_rates){
    JSPEC_TEST_BEGIN("Electron cooling:");

    //define coasting 12C6+ beam
    int n_charge = 6;
    double n_mass = 12;
    double kinetic_energy = 30*n_mass;
    double gamma = 1+kinetic_energy/(n_mass*k_u);
    double beta = sqrt(1-1/(gamma*gamma));
    double emit_nx0 = beta*gamma*5e-6;
    double emit_ny0 = beta*gamma*5e-6;
    double dp_p0 = 0.0004;
    double n_ptcl = 5E8;
    Beam c_beam(n_charge, n_mass, kinetic_energy, emit_nx0, emit_ny0, dp_p0, n_ptcl);

    //define the ring
//    std::string filename = "csrm.tfs";
    std::string filename = CMAKE_SOURCE_DIR + std::string("/data/Booster.tfs");
    Lattice lattice(filename);
    Ring ring(lattice, c_beam);

    //define the cooler
    double cooler_length = 3.4;
    double n_section = 1;
    double magnetic_field = 0.039;
    double beta_h = 10;
    double beta_v = 17;
    double dis_h = 0;
    double dis_v = 0;
    Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);

    //define electron beam
    double current = 0.03;
    double radius = 0.025;
    double neutralisation = 0;
    UniformCylinder uniform_cylinder(current, radius, neutralisation);
    double gamma_e = c_beam.gamma();
    double tmp_tr = 0.05;
    double tmp_long = 0.1;
    EBeam e_beam(gamma_e, tmp_tr, tmp_long, uniform_cylinder);


    //define cooling model
//            unsigned int n_sample = 5000;
//            EcoolRateParas ecool_rate_paras(n_sample);
    unsigned int n_tr = 100;
    unsigned int n_long = 100;
    EcoolRateParas ecool_rate_paras(n_tr, n_long);

    ForceParas *force_paras = ChooseForce(ff);
    double rate_x, rate_y, rate_s;
    ecooling_rate(ecool_rate_paras, *force_paras, c_beam, cooler, e_beam, ring, rate_x, rate_y, rate_s);
    std::cout<<std::endl;
    std::cout<<"rate_x = "<<rate_x<<" rate_y = "<<rate_y<<" rate_s = "<<rate_s<<std::endl;


//This set of test values was satisfied with the original F_const definitions
// in terms of k_ke in Parkhomchuk
//    JSPEC_ASSERT_THROW(abs(rate_x + 0.00865669) < 1e-5);
//    JSPEC_ASSERT_THROW(abs(rate_y + 0.00900383) < 1e-5);
//    JSPEC_ASSERT_THROW(abs(rate_s + 0.0190261) < 1e-5);

//This set of test values is satisfied with the F_const defined
// in terms of r_e in Parkhomchuk
    JSPEC_ASSERT_THROW(abs(rate_x - check_rates[0]) < 1e-5);
    JSPEC_ASSERT_THROW(abs(rate_y - check_rates[1]) < 1e-5);
    JSPEC_ASSERT_THROW(abs(rate_s - check_rates[2]) < 1e-5);


    std::cout<<"Force calculation"<<std::endl;
    //electrons and ions must be bunched for this to work
    double sigma_x = 1.5e-2;
    double sigma_y = 1.5e-2;
    double sigma_s = 2e-2;
    double n_particle = 1e7;
    GaussianBunch gb(n_particle,sigma_x,sigma_y,sigma_s);
    EBeam e_beam2(gamma_e, tmp_tr, tmp_long, gb);
    Beam c_beam2(n_charge, n_mass, kinetic_energy, emit_nx0, emit_ny0, dp_p0,
                 sigma_s, //This makes sure c_beam is bunched and allows that v_long can have > 2 values
                 n_ptcl);

    n_long = 2;
    EcoolRateParas ecool_rate_paras2(n_tr, n_long);
    //ForceParas force_paras2(ForceFormula::DERBENEVSKRINSKY);

    ForceParas *force_paras2 = ChooseForce(ff);
    CalculateForce(ecool_rate_paras2, *force_paras2, c_beam2, cooler, e_beam2, ring);

    JSPEC_TEST_END();

    return 0;
}

int ibs(){
    JSPEC_TEST_BEGIN("IBS:");
    double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
    int Z;
//                Z = 1;
//                m0 = 938.272;
//                KE = 100e3;
//                emit_nx0 = 6.00332e-006;
//                emit_ny0 = 3.01154e-007;
//                dp_p0 = 0.000997401;
//                sigma_s0 = 0.0284972;
//                N_ptcl = 6.56E9;

    Z = 1;
    m0 = 938.272;
    KE = 8000;
    emit_nx0 = 2.2e-6;
    emit_ny0 = 2.2e-6;
    dp_p0 = 0.001*0.6;
    N_ptcl = 6.58E11;
    sigma_s0 = 7;
    Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl);

    // define the lattice of the proton ring
//    std::string filename = "MEICColliderRedesign1IP.tfs";
    std::string filename = CMAKE_SOURCE_DIR + std::string("/data/Booster.tfs");
    Lattice lattice(filename);

    //Define the ring
    Ring ring(lattice, p_beam);

    //Set IBS parameters.
    int nu = 100;
    int nv = 100;
    int nz = 40;
    double log_c = 39.2/2;
    double k = 0.0;
    IBSSolver_Martini ibs_solver(nu, nv, nz, log_c, k);

//Calculate IBS rate.
//  std::chrono::steady_clock::time_point start, end;
//  start = std::chrono::steady_clock::now();
//
    double rx_ibs, ry_ibs, rz_ibs;
//  ibs_rate(lattice, p_beam, ibs_solver, rx_ibs, ry_ibs, rz_ibs);
//
//  end = std::chrono::steady_clock::now();
//  auto t1 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//  std::cout<<"IBS 3D integral: "<<t1<<std::endl;
//
//  std::cout<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;

    ibs_solver.rate(lattice, p_beam, rx_ibs, ry_ibs, rz_ibs);
    std::cout<<std::endl;
    std::cout<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;

    JSPEC_ASSERT_THROW(abs(rx_ibs - 0.000737069) < 1e-4);
    JSPEC_ASSERT_THROW(abs(ry_ibs + 8.61496e-6) < 1e-6);
    JSPEC_ASSERT_THROW(abs(rz_ibs - 0.000651936) < 1e-5);

    // Test a different coefficient
    ibs_solver.set_k(0.2);
    ibs_solver.rate(lattice, p_beam, rx_ibs, ry_ibs, rz_ibs);
    std::cout<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;

    JSPEC_ASSERT_THROW( abs(rx_ibs - 0.0006625) < 1e-4 );
    JSPEC_ASSERT_THROW( abs(ry_ibs - 6.59534e-05) < 1e-6 );
    JSPEC_ASSERT_THROW( abs(rz_ibs - 0.000651936) < 1e-5 );

    JSPEC_TEST_END();
    return 0;
}

int dynamicibsbunched(){
    JSPEC_TEST_BEGIN("Dynamic IBS Bunched:")

    // define proton beam;
    double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
    int Z;
    Z = 1;
    m0 = 938.272;
    KE = 3e4;
    emit_nx0 = 0.4962695094e-6;
    emit_ny0 = 0.4962695094e-6;
    dp_p0 = 4e-4;
    sigma_s0 = 1.994525702e-2;
    N_ptcl = 6.56E9;
    Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl);

    // define the lattice of the proton ring
    //std::string filename = "MEICColliderRedesign1IP.tfs";
    std::string filename = CMAKE_SOURCE_DIR + std::string("/data/Booster.tfs");
    Lattice lattice(filename);

    //Define the ring
    Ring ring(lattice, p_beam);

    //Set IBS parameters.
    int nu = 200;
    int nv = 200;
    int nz = 40;
    double log_c = 39.9/2;
    double k = 0.5;
    ibs_solver = new IBSSolver_Martini(nu, nv, nz, log_c, k);

    dynamic_paras = new DynamicParas(3600, 360, true, false);

    Cooler *cooler=nullptr;
    EBeam *e_beam=nullptr;

    //Skip the Force visualization calculation
    dynamic_paras->set_test(true);
    dynamic(p_beam, *cooler, *e_beam, ring);

    JSPEC_TEST_END();
    return 0;
}

int dynamicibs(){
    JSPEC_TEST_BEGIN("Dynamic IBS:")

    double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
    int Z;
    Z = 1;
    m0 = 938.272;
    KE = 800;
    emit_nx0 = 1.039757508e-6;
    emit_ny0 = 1.039757508e-6;
    dp_p0 = 2e-3;
    N_ptcl = 3.6E11;
    Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, N_ptcl);

    // define the lattice of the proton ring
    //std::string filename = "MEICBoosterRedesign.tfs";
    std::string filename = CMAKE_SOURCE_DIR + std::string("/data/Booster.tfs");
    Lattice lattice(filename);

    //Define the ring
    Ring ring(lattice, p_beam);

    //Set IBS parameters.
    int nu = 200;
    int nv = 200;
    int nz = 40;
    double log_c = 44.8/2;
    ibs_solver = new IBSSolver_Martini(nu, nv, nz, log_c, 1.0);

    dynamic_paras = new DynamicParas(3600, 360, true, false);

//  char file[100] = "test_dynamic_ibs.txt";
//  std::ofstream outfile;
//  outfile.open(file);
    Cooler *cooler=nullptr;
    EBeam *e_beam=nullptr;
//  dynamic(p_beam, *cooler, *e_beam, ring, outfile);
//  outfile.close();

    dynamic_paras->set_output_file("test_dynamic_ibs.txt");
    //Skip the Force visualization calculation
    dynamic_paras->set_test(true);
    dynamic(p_beam, *cooler, *e_beam, ring);


    JSPEC_TEST_END();
    return 0;
}

int dynamicecool(){

    JSPEC_TEST_BEGIN("Dynamic Ecool:");

    // define proton beam;
    double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
    int Z;
    Z = 1;
    m0 = 938.272;
    KE = 800;
    emit_nx0 = 1.039757508e-6;
    emit_ny0 = 1.039757508e-6;
    dp_p0 = 0.002012615391;
    N_ptcl = 3.6E11;
    Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, N_ptcl);

    // define the lattice of the proton ring
//    std::string filename = "MEICBoosterRedesign.tfs";
    std::string filename = CMAKE_SOURCE_DIR + std::string("/data/Booster.tfs");
    Lattice lattice(filename);

    //Define the ring
    Ring ring(lattice, p_beam);

    //define the cooler
    double cooler_length = 10;
    double n_section = 1;
    double magnetic_field = 0.1;
    double beta_h = 10;
    double beta_v = 10;
    double dis_h = 0;
    double dis_v = 0;
    Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);

    //define electron beam
    double current = 2;
    double radius = 0.008;
    double neutralisation = 0;
    UniformCylinder uniform_cylinder(current, radius, neutralisation);
    double gamma_e = p_beam.gamma();
    double tmp_tr = 0.1;
    double tmp_long = 0.1;
    EBeam e_beam(gamma_e, tmp_tr, tmp_long, uniform_cylinder);

    ////define cooling model: single particle
    // unsigned int n_tr = 100;
    // unsigned int n_long = 100;
    // ecool_paras = new EcoolRateParas(n_tr, n_long);
    //define cooling model: monte carlo
    unsigned int n_sample = 40000;
    ecool_paras = new EcoolRateParas(n_sample);
    //define friction force formula
    force_paras = ChooseForce(ForceFormula::PARKHOMCHUK);
    //define dynamic simulation
    dynamic_paras = new DynamicParas(60, 120, false, true);
    dynamic_paras->set_model(DynamicModel::MODEL_BEAM);

//  char file[100] = "test_dynamic_ecool_DC_model_beam.txt";
//  std::ofstream outfile;
//  outfile.open(file);
//  dynamic(p_beam, cooler, e_beam, ring, outfile);
//  outfile.close();

    dynamic_paras->set_output_file("test_dynamic_ecool_DC_model_beam.txt");
    //Skip the Force visualization calculation
    dynamic_paras->set_test(true);
    dynamic(p_beam, cooler, e_beam, ring);


    JSPEC_TEST_END();
    return 0;
}

int dynamicboth(){

    JSPEC_TEST_BEGIN("Dynamic IBS + ECOOL:");

    srand(time(NULL));
    //            srand(0);

    // define proton beam;
    double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
    int Z;
    Z = 1;
    m0 = 938.272;
    KE = 100e3;
//            // CM energy 44.7 GeV
//            emit_nx0 = 0.5e-6;
//            emit_ny0 = 0.12e-6;
//            dp_p0 = 0.0008;
//            N_ptcl = 0.98E10*0.59;
//            sigma_s0 = 2.5E-2;

    double factor = 2.25;
    // CM energy 44.7 GeV
    emit_nx0 = 0.5e-6*factor;
    emit_ny0 = 0.15e-6*factor;
    dp_p0 = 0.0008;
    N_ptcl = 0.98E10*0.93;
    sigma_s0 = 1.5E-2;

//            // CM energy 63.3 GeV
//            emit_nx0 = 1.25e-6;
//            emit_ny0 = 0.38e-6;
//            dp_p0 = 0.0008;
//            N_ptcl = 3.9E10*0.25;
//            sigma_s0 = 2.5E-2;
    Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl);
    std::cout<<"Normalized emittance: "<<p_beam.emit_nx()<<' '<<p_beam.emit_ny()<<std::endl;
    std::cout<<"Geometric emittance: "<<p_beam.emit_x()<<' '<<p_beam.emit_y()<<std::endl;

    // define the lattice of the proton ring
//            std::string filename = "MEICBoosterRedesign.tfs";
//    std::string filename = "MEICColliderRedesign1IP.tfs";
    std::string filename = CMAKE_SOURCE_DIR + std::string("/data/Booster.tfs");
    Lattice lattice(filename);

    //Define the ring
    Ring ring(lattice, p_beam);

//            //define the cooler
    double cooler_length = 60;
    double n_section = 1;
    double magnetic_field = 1;
    double beta_h = 60;
    double beta_v = 200;
    double dis_h = 2.0;
    double dis_v = 0.6;
    Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);
    std::cout<<"Ion beam size at cooler: "<<sqrt(cooler.beta_h()*p_beam.emit_x())
            <<' '<<sqrt(cooler.beta_v()*p_beam.emit_y())<<std::endl<<std::endl;

//            //define electron beam
    double length = 0.02;
    double radius = 0.000528*sqrt(factor);
    std::cout<<"Electron beam radius: "<<radius<<std::endl;
    double q_e = 2.0e-9;
    double current = q_e*p_beam.beta()*k_c/length;
    UniformBunch uniform_bunch(current, radius, length);
    double gamma_e = p_beam.gamma();
    double tmp_tr = 0.246;
    double tmp_long = 0.184;
    EBeam e_beam(gamma_e, tmp_tr, tmp_long, uniform_bunch);

//            double n_e = 1.248E10;
//            double sigma_e_x = 0.035E-2;
//            double sigma_e_y = 0.035E-2;
//            double sigma_e_s = 0.84E-2;
//            GaussianBunch gaussian_bunch(n_e, sigma_e_x, sigma_e_y, sigma_e_s);
//            double gamma_e = p_beam.gamma();
//            double tmp_tr = 0.5;
//            double tmp_long = 0.1;
//            EBeam e_beam(gamma_e, tmp_tr, tmp_long, gaussian_bunch);
//
    //define cooling model: monte carlo
    unsigned int n_sample = 100000;
    ecool_paras = new EcoolRateParas(n_sample);
//            //define friction force formula
    force_paras = ChooseForce(ForceFormula::PARKHOMCHUK);

    double rate_x, rate_y, rate_s;
    ecooling_rate(*ecool_paras, *force_paras, p_beam, cooler, e_beam, ring, rate_x, rate_y, rate_s);
    std::cout<<"cooling rate: [1/s] "<<rate_x<<" "<<rate_y<<" "<<rate_s<<std::endl;

    //Set IBS parameters.
    int nu = 100;
    int nv = 100;
    int nz = 40;
    double log_c = 40.4/2;      //100 GeV, 63.3 GeV CM Energy
    double k = 0.4;
//            double log_c = 39.2/2;    //100 GeV
    ibs_solver = new IBSSolver_Martini(nu, nv, nz, log_c, k);
//            ibs_solver = new IBSSolver(nu, nv, nz);

    double rx_ibs, ry_ibs, rz_ibs;
    ibs_solver->rate(lattice, p_beam, rx_ibs, ry_ibs, rz_ibs);
    std::cout<<"IBS rate: [1/s] "<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;
    std::cout<<"Total rate: [1/s] "<<rx_ibs+rate_x<<' '<<ry_ibs+rate_y<<' '<<rz_ibs+rate_s<<std::endl<<std::endl;

//                return 0;

    //define dynamic simulation
    double t = 1200*2;
    int n_step = 600*2;
    bool ibs = true;
    bool ecool = true;
    dynamic_paras = new DynamicParas(t, n_step, ibs, ecool);
    dynamic_paras->set_model(DynamicModel::MODEL_BEAM);
//                dynamic_paras->set_model(DynamicModel::RMS);

//                char file[100] = "Collider_100GeV_strong_cooling_baseline_44Ecm_2nC_rms_02_long_cool_compensated.txt";
//                std::ofstream outfile;
//                outfile.open(file);
//                Cooler *cooler;
//                EBeam *e_beam;
//                if(twiss_ref.get()==nullptr) twiss_ref.reset(new Twiss());
//                twiss_ref->bet_x = 10;
//                twiss_ref->bet_y = 10;
//                twiss_ref->disp_x = 5;
    dynamic_paras->twiss_ref.bet_x = 10;
    dynamic_paras->twiss_ref.bet_y = 10;
    dynamic_paras->twiss_ref.disp_x = 5;
//                dynamic_paras->set_n_sample(5000);
//                dynamic(p_beam, *cooler, *e_beam, ring, outfile);
//                dynamic(p_beam, cooler, e_beam, ring, outfile);
//                outfile.close();

    dynamic_paras->set_output_file("Collider_100GeV_strong_cooling_baseline_44Ecm_2nC_rms_02_long_cool_compensated.txt");
    //Skip the Force visualization calculation
    dynamic_paras->set_test(true);
    dynamic(p_beam, cooler, e_beam, ring);


    JSPEC_TEST_END();
    return 0;
}


int main(int, char**)
{

    double check_rates[3];
    //  check_rates[0] = -0.0087216;
   //check_rates[1] = -0.00907129;
    //check_rates[2] = -0.0191686;
    //ecool(ForceFormula::PARKHOMCHUK,check_rates);

    //    check_rates[0] = -0.00503489;
    //check_rates[1] = -0.00544656;
    //check_rates[2] = -0.0122515;
    //ecool(ForceFormula::DERBENEVSKRINSKY,check_rates);


//    check_rates[0] = -1.72854;
//    check_rates[1] = -1.72863;
//    check_rates[2] = -2.61073;
//    ecool(ForceFormula::MESHKOV,check_rates);

//    dynamicibsbunched();
//    dynamicibs();
//    dynamicecool();
//    ibs();
    //These both take a long time
    //dynamicboth();
   // both();

}
