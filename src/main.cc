#include <algorithm>
#include <chrono>
#include <fstream>
#include <map>
#include "dynamic.h"
#include "ecooling.h"
#include "ibs.h"
#include "math_parser.h"
#include "muParserDLL.h"
#include "ring.h"
#include "ui.h"


using std::string;

extern DynamicParas * dynamic_paras;
extern IBSParas * ibs_paras;
extern EcoolRateParas * ecool_paras;
extern ForceParas * force_paras;
extern Luminosity *luminosity_paras;
//extern std::unique_ptr<Twiss> twiss_ref;
//extern int n_ion_model;

extern std::map<std::string, Section> sections;
extern muParserHandle_t math_parser;

//extern std::vector<std::string> sections;
//extern std::vector<std::string> ion_pars;


enum class Test {IBS, ECOOL, BOTH, DYNAMICIBS, DYNAMICECOOL, DYNAMICBOTH, DYNAMICIBSBUNCHED};

int main(int argc, char** argv) {

  Test test = Test::ECOOL;

  //********************************
  // Test Electron Cooling Rate
  //********************************
    
  //define coasting gold ion beam
  int n_charge = 79;
  double n_mass = 197;
  double kinetic_energy = 100000 * n_mass; //120*1000;
  double gamma = 1+kinetic_energy/(n_mass*k_u);
  double beta = sqrt(1-1/(gamma*gamma));
  std::cout<<"gamma = "<<gamma<<" beta = "<<beta<<std::endl;
  double emit_nx0 = beta*gamma*5e-6;
  double emit_ny0 = beta*gamma*5e-6;
  double dp_p0 = 0.002; //0.0004;
  double n_ptcl = 5E8;

    //This constructor does not define a sigma s, so it's got a flag set for bunched=false
  //Beam c_beam(n_charge, n_mass, kinetic_energy, emit_nx0, emit_ny0, dp_p0, n_ptcl);

  double sigma_s_ion = 2e-2;  //What's a reasonable value here?
  Beam c_beam(n_charge, n_mass, kinetic_energy, emit_nx0, emit_ny0, dp_p0, sigma_s_ion, n_ptcl);
    
  //define the ring
  std::string filename = "Booster.tfs";
  Lattice lattice(filename);
  Ring ring(lattice, c_beam);

  //define the cooler
  double cooler_length = 10;//3.4;
  double n_section = 1;
  double magnetic_field = 5;//0.039;
  double beta_h = 10;
  double beta_v = 10;//17;
  double dis_h = 0; //dispersion - definitions not required
  double dis_v = 0;
  Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);

 //unused - vesigial from original main.cc
  double current = 0.03;
  double radius = 0.025;
  double neutralisation = 0;
  UniformCylinder uniform_cylinder(current, radius, neutralisation);

  //We want to test a gaussian bunch for Derbenev & Skrinsky
  //made up numbers, just trying to figure out what gaussian bunch does
  double sigma_x = 1.5e-2; //RMS bunch size in meters
  double sigma_y = 1.5e-2;
  double sigma_s = 2e-2;   //RMS bunch length (longitudinal) in meters
  double n_particle = 1e7;

  GaussianBunch gb(n_particle,sigma_x,sigma_y,sigma_s);

  double gamma_e = c_beam.gamma(); //ebeam and ion beam are co-travelling
  double tmp_tr = 0.1; //0.05; //Transverse temperature
  double tmp_long = 0.1; //longitudinal temperature
  EBeam e_beam(gamma_e, tmp_tr, tmp_long, gb);

  //define cooling model

  //Here, the constructor signifies the monte carlo method for rate calculation
  //unsigned int n_sample = 5000;
  //EcoolRateParas ecool_rate_paras(n_sample);

  //Here, the constructor signifies the single-particle method of rate calculation.
  // The n's are the number of different values for single particles sampling the phase space
  unsigned int n_tr = 4; //number of phi for transverse direction
  unsigned int n_long = 400; //number of phi for longitudinal direction
  EcoolRateParas ecool_rate_paras(n_tr, n_long);

  double rate_x, rate_y, rate_s;
  if ( argc > 1){
      ForceParas force_paras(ForceFormula::PARKHOMCHUK);
      ecooling_rate(ecool_rate_paras, force_paras, c_beam, cooler, e_beam, ring, rate_x, rate_y, rate_s);
      std::cout<<"Parkhomchuk model:"<<std::endl;
      std::cout<<"rate_x = "<<rate_x<<" rate_y = "<<rate_y<<" rate_s = "<<rate_s<<std::endl;
  }
   else{
      ForceParas force_parasDS(ForceFormula::DERBENEVSKRINSKY);
      ecooling_rate(ecool_rate_paras, force_parasDS, c_beam, cooler, e_beam, ring, rate_x, rate_y, rate_s);
      std::cout<<"Derbenev-Skrinsky model:"<<std::endl;
      std::cout<<"rate_x = "<<rate_x<<" rate_y = "<<rate_y<<" rate_s = "<<rate_s<<std::endl;
   }


    //Pause the system
//    std::cout<<std::endl<<"Press the Enter key to close the window."<<std::endl;
//    std::cin.ignore( std::numeric_limits< std::streamsize >::max( ), '\n' );

    return 0;
}
