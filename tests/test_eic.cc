#include <math.h>
#include <Base.hh>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdio.h>

#include "force.h"
#include "ecooling.h"

#include <gsl/gsl_fit.h>

using std::string;

extern EcoolRateParas * ecool_paras;
extern ForceParas * force_paras;

//Run a cooling rate calculation relevant to EIC and compare 
// it to a stored result. 

//To perform this test, run with a specific setup, selecting
// a specific force function, and output results to a file.
// Then, perform a regression test on those results compared to a
// 'golden' dataset stored in the data subdirectory.


void SetupModel(ForceFormula ff)
{
  //define coasting gold ion beam
  int n_charge = 79;
  double n_mass = 197;
  double kinetic_energy = 1e5*n_mass;
  double gamma = 1+kinetic_energy/(n_mass*k_u);
  double beta = sqrt(1-1/(gamma*gamma));
  double emit_nx0 = beta*gamma*5e-6;
  double emit_ny0 = beta*gamma*5e-6;
  std::cout<<std::endl;
  std::cout<<"gamma ="<<gamma<<" beta ="<<beta<<" emit_nx0="<<emit_nx0<<std::endl;
  double dp_p0 = 0.004; //0.0004
  double n_ptcl = 5E8;
  double sigma_s_ion = 2e-2;
  Beam ion_beam(n_charge, n_mass, kinetic_energy, emit_nx0, emit_ny0, dp_p0, sigma_s_ion, n_ptcl);

  //define the ring (not really used for friction force calculation)
  string filename = CMAKE_SOURCE_DIR + std::string("/data/Booster.tfs");
  
  Lattice lattice(filename);
  Ring ring(lattice, ion_beam);

  //define the cooler
  double cooler_length = 1;
  double n_section = 1;
  double magnetic_field = 5; 
  double beta_h = 10;
  double beta_v = 17;
  double dis_h = 0;
  double dis_v = 0;
  Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);

  //define electron bunch
  double sigma_x = 1.5e-2;
  double sigma_y = 1.5e-2;
  double sigma_s = 2e-2;
  double n_particle = 1e7;
  GaussianBunch gb(n_particle,sigma_x,sigma_y,sigma_s);
            
  double gamma_e = ion_beam.gamma();
  double tmp_tr = 0.1; //0.05
  double tmp_long = 0.1;
  EBeam e_beam(gamma_e, tmp_tr, tmp_long, gb);
    
  //This choice tells ecool_rate_paras to select
  // 4 transverse velocities spanning the whole space,
  // from which we can whittle down to the one we are
  // interested in. As yet, there is no way to select
  // a single velocity in a straightforward way.
  unsigned int n_tr = 4;
  unsigned int n_long = 400;
  EcoolRateParas ecool_rate_paras(n_tr, n_long); //This sets IonSample::SingleParticle

  double rate_x, rate_y, rate_s;
            
  ForceParas *force_paras = ChooseForce(ff);
    
  force_paras->set_do_test(true);
  ecooling_rate(ecool_rate_paras, *force_paras, ion_beam, cooler, e_beam, ring, rate_x, rate_y, rate_s);
  std::cout<<"rate_x = "<<rate_x<<" rate_y = "<<rate_y<<" rate_s = "<<rate_s<<std::endl;  
    
  return;
}


void SetupBoosterModel(ForceFormula ff)
{
  //Set up the same model as the Booster example on Sirepo
  int n_charge = 1;
  double n_mass = 1;
  double kinetic_energy = 7942.25; //gamma = 9.526
  double gamma = 1+kinetic_energy/(n_mass*k_u);
  double beta = sqrt(1-1/(gamma*gamma));
  double emit_nx0 = 2.2e-6;
  double emit_ny0 = 2.2e-6;
  std::cout<<std::endl;
  std::cout<<"gamma ="<<gamma<<" beta ="<<beta<<" emit_nx0="<<emit_nx0<<std::endl;
  double dp_p0 = 6e-4;
  double n_ptcl = 6.58e13;
  double sigma_s_ion = 6e-4;
  Beam ion_beam(n_charge, n_mass, kinetic_energy, emit_nx0, emit_ny0, dp_p0, sigma_s_ion, n_ptcl);

  //define the ring (not really used for friction force calculation)
  string filename = CMAKE_SOURCE_DIR + std::string("/data/Booster.tfs");
  
  Lattice lattice(filename);
  Ring ring(lattice, ion_beam);

  //define the cooler
  double cooler_length = 3.4;
  double n_section = 1;
  double magnetic_field = 0.039;
  double beta_h = 10;
  double beta_v = 10;
  double dis_h = 0;
  double dis_v = 0;
  Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);

  //define electron beam
  double current = 2; //Amps 
  double radius = 0.004; //m
  double neutralisation = 2;
  UniformCylinder uc(current,radius,neutralisation);

  double gamma_e = ion_beam.gamma();
  double tmp_tr = 0.1;
  double tmp_long = 0.01;
  EBeam e_beam(gamma_e, tmp_tr, tmp_long, uc);
    
  //This choice tells ecool_rate_paras to select
  // 4 transverse velocities spanning the whole space,
  // from which we can whittle down to the one we are
  // interested in. As yet, there is no way to select
  // a single velocity in a straightforward way.
  unsigned int n_tr = 100;
  unsigned int n_long = 400;
  EcoolRateParas ecool_rate_paras(n_tr, n_long); //This sets IonSample::SingleParticle

  double rate_x, rate_y, rate_s;
            
  ForceParas *force_paras = ChooseForce(ff);
    
  force_paras->set_do_test(true);
  ecooling_rate(ecool_rate_paras, *force_paras, ion_beam, cooler, e_beam, ring, rate_x, rate_y, rate_s);
  std::cout<<"rate_x = "<<rate_x<<" rate_y = "<<rate_y<<" rate_s = "<<rate_s<<std::endl;  
    
  return;
}
void BetacoolBoosterModel(ForceFormula ff)
{
  //Set up the benchmarking model provided by He Zhang
  int n_charge = 1;
  double n_mass = 1;
  double kinetic_energy = 14097.; //MeV, Betacool uses GeV/u, this fixes momentum @ 15 GeV/c
  double gamma = 1+kinetic_energy/(n_mass*k_u); // = 16.13405
  double beta = sqrt(1-1/(gamma*gamma));
  double emit_nx0 = 1.039757/1e6; //normalized emittance in m, converted from 1.039 pi-mm-mrad (normalized), in colliding beam window
  double emit_ny0 = emit_nx0;
  std::cout<<std::endl;
  std::cout<<"gamma ="<<gamma<<" beta ="<<beta<<" emit_nx0="<<emit_nx0<<std::endl;
  double dp_p0 = 0.00291917972; //GeV/c, corresponding to 0.002 in unitless longitudinal scale (?)
  double n_ptcl = 3.6e11; //from coliding beam window
  double sigma_s_ion = 0.3; //m, from 30cm, setting this defines the beam state as bunched

  //This constructor defines a bunched ion beam with rms width sigma_s_ion
  //Beam ion_beam(n_charge, n_mass, kinetic_energy, emit_nx0, emit_ny0, dp_p0, sigma_s_ion, n_ptcl);

  //This constructor defines a continuous ion beam
  Beam ion_beam(n_charge, n_mass, kinetic_energy, emit_nx0, emit_ny0, dp_p0, n_ptcl); 

  //define the ring (not really used for friction force calculation)
  string filename = CMAKE_SOURCE_DIR + std::string("/data/Booster.tfs");
  
  Lattice lattice(filename);
  Ring ring(lattice, ion_beam);

  //define the cooler
  double cooler_length = 10; //in m
  double n_section = 1;
  double magnetic_field = 0.1; //Tesla, from 1kG 
  double beta_h = 10; //m
  double beta_v = 10; //m
  double dis_h = 0;
  double dis_v = 0;
  Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);

  //define electron beam
  double current = 2; //Amps 
  double radius = 0.008; //m, from 0.8cm
  double neutralisation = 2; //%from 200%
  UniformCylinder uc(current,radius,neutralisation);

  double gamma_e = ion_beam.gamma();
  double tmp_tr = 0.1; //for parkhomchuck 
  double tmp_long = 0.1;
  EBeam e_beam(gamma_e, tmp_tr, tmp_long, uc);
    
  //This choice tells ecool_rate_paras to select
  // 4 transverse velocities spanning the whole space,
  // from which we can whittle down to the one we are
  // interested in. As yet, there is no way to select
  // a single velocity in a straightforward way.
  unsigned int n_tr = 100;
  unsigned int n_long = 400;
  EcoolRateParas ecool_rate_paras(n_tr, n_long); //This sets IonSample::SingleParticle

  double rate_x, rate_y, rate_s;
            
  ForceParas *force_paras = ChooseForce(ff);

    /*
  if(force_paras->formula() == ForceFormula::ERLANGEN){
      std::string suffix;
      //Start with a clean slate
      dynamic_cast<Force_Erlangen *>(force_paras)->set_fast(false);
      dynamic_cast<Force_Erlangen *>(force_paras)->set_tight(false);
      dynamic_cast<Force_Erlangen *>(force_paras)->set_stretched(false);
                              
      if(Erlangen_fast){
          dynamic_cast<Force_Erlangen *>(force_paras)->set_fast(true);
          suffix += "F";
      }
      if(Erlangen_tight) {
          dynamic_cast<Force_Erlangen *>(force_paras)->set_tight(true);
          suffix += "T";
      }
      if(Erlangen_stretched){
          dynamic_cast<Force_Erlangen *>(force_paras)->set_stretched(true);
          suffix += "S";
      }
      force_paras->set_filename(((std::string)"Erlangen_" + suffix + (std::string)".txt").c_str());
      
      //Speed this up
      dynamic_cast<Force_Erlangen *>(force_paras)->set_calls(50000);
      
  }
  */
    
  force_paras->set_do_test(true);
  ecooling_rate(ecool_rate_paras, *force_paras, ion_beam, cooler, e_beam, ring, rate_x, rate_y, rate_s);
  std::cout<<"rate_x = "<<rate_x<<" rate_y = "<<rate_y<<" rate_s = "<<rate_s<<std::endl;  
    
  return;
}


//A method to sort rows of a matrix based on a specific column.
// This modifies the original matrix
void sortrows(std::vector<std::vector<double>>& matrix, int col) {    
    std::sort(matrix.begin(),
              matrix.end(),
              [col](const std::vector<double>& lhs, const std::vector<double>& rhs) {
                  return lhs[col] > rhs[col];
              });
}

//A standard C++ CSV parser to handle the input files we'll need
// for these comparisons
vector< vector<double> > ReadCSV(string filename){
    
    std::ifstream infile(filename);

    string line;
    vector<double> vec;
    vector< vector<double> > matrix;

    while( std::getline(infile,line) ){
        std::stringstream ss(line);
        vector<double> row;
        std::string data;
        while ( getline(ss,data,',') ){

            //Handle lines that don't match the format
            try{
                row.push_back( stod(data) );
            }
            catch(const std::exception& e){
             //Something other than the double we expect has snuck in.   
             // Use the testing framework to point us to this line # for guidance
              JSPEC_ASSERT_THROW( 1 == 0);
            }
        }
        //The variables are F_const, v_trans, v_long, electron density, force_trans, force_long
        JSPEC_ASSERT_THROW( row.size() == 6 );
        matrix.push_back(row);
    }    

    //now that we've ingested the whole file,
    // make cuts to get the range we want to test.

    vector< vector<double> > out_matrix;

    for(int i=0;i<matrix.size();i++){
        vector<double> vec = matrix[i];
        if(vec.size()>0){
            if (vec[1] < 1){ //Pick a specific low v_transverse
                out_matrix.push_back(vec);
            }
        }
    }
    
    //Now we need to sort by v_long and return it to the
    // comparison function
    sortrows(out_matrix,2);
    
    return out_matrix;
}

//Parse the file generated in SetupModel and compare it to the 
// 'golden' dataset that is stored in the data directory
double CompareOutput(string filename_golden,string filename_test){

    vector< vector<double> > golden = ReadCSV(filename_golden);
    vector< vector<double> > test = ReadCSV(filename_test);

    //Calculate the difference between these trends, so there's just
    // one trend to fit    
    vector<double> x_vec;
    vector<double> y_vec;
    vector<double> y_vec_g;
    int n = 0;

    for(int i = 0; i<golden.size(); i++){
        
        double v_long_g = golden[i][2];
        double f_long_g = golden[i][5];
     
        //Find the matching point for the test dataset
        // & add it to the fit vector. Our data is already
        // sorted by v_long.
        double min = DBL_MAX; //The closest point
        int min_index = 0;
        
        for(int j = 0; j<test.size();j++){
            double v_long_t = test[j][2];
            //can't expect the sampled velocities to be exactly
            // the same, but check that they're close

            if(abs(v_long_t - v_long_g) < min){
                min = abs(v_long_t - v_long_g);
                min_index = j;
            }
        }

        JSPEC_ASSERT_THROW( min < DBL_MAX );
        double f_long_t = test[min_index][5];
        x_vec.push_back(test[min_index][2]);
        y_vec.push_back(f_long_t - f_long_g);
        y_vec_g.push_back(f_long_g);
        n++;               
    }

    //Use the GSL library we've already introduced as a depencency
    // to compare these trends. Initialize the inputs and outputs:
    double cov00,cov01,cov11;
    double c0,c1,sumsq;
    
    //Strides refer to the separation within
    // the array of consecutive points
    int xstride = 1;
    int ystride = 1; 
    
    //The difference between the reference trend and
    // the test trend should fit to a linear fit with 
    // slope ~ 0 and intercept ~ 0. 
    assert(x_vec.size() == y_vec.size());
    gsl_fit_linear(
      x_vec.data(),
      xstride,
      y_vec.data(),
      ystride,
      x_vec.size(),
      &c0,
      &c1,
      &cov00,
      &cov01,
      &cov11,
      &sumsq
    );
    
    std::cout<<std::endl;
    std::cout<<n<<" points in comparison"<<std::endl;
    std::cout<<"Fit: intercept = "<<c0<<" slope = "<<c1<<" sumsq = "<<sumsq<<std::endl;
    
    return c1;    
}

void testForce(){

  JSPEC_TEST_BEGIN("Magnetized Electron Cooling:");

  BetacoolBoosterModel(ForceFormula::PARKHOMCHUK);
  //Run the quick simulation for the model
  SetupModel(ForceFormula::DERBENEVSKRINSKY);
  //Get the output and compare via regression
  string data_path = CMAKE_SOURCE_DIR + std::string("/data/dumpDS.txt");
  string test_path = CMAKE_SOURCE_DIR + std::string("/build/tests/DerbenevSkrinsky.txt");
              
  double slope = CompareOutput(data_path,test_path);
  //JSPEC_ASSERT_THROW( abs(slope) < 1e-40 );
  
  //A negative control: compare D&S forces to Parkhomchuk forces
  SetupModel(ForceFormula::DERBENEVSKRINSKY);  
  data_path = CMAKE_SOURCE_DIR + std::string("/data/dumpPark.txt");
  test_path = CMAKE_SOURCE_DIR + std::string("/build/tests/DerbenevSkrinsky.txt");  
  slope = CompareOutput(data_path,test_path);
  JSPEC_ASSERT_THROW( abs(slope) > 1e-40 );
  
  SetupModel(ForceFormula::PARKHOMCHUK);
  data_path = CMAKE_SOURCE_DIR + std::string("/data/dumpPark.txt");
  test_path = CMAKE_SOURCE_DIR + std::string("/build/tests/Parkhomchuk.txt");
  slope = CompareOutput(data_path,test_path);
  JSPEC_ASSERT_THROW( abs(slope) < 1e-28 );
  
    
  SetupModel(ForceFormula::MESHKOV);
  data_path = CMAKE_SOURCE_DIR + std::string("/data/dumpMesh.txt");
  test_path = CMAKE_SOURCE_DIR + std::string("/build/tests/Meshkov.txt");
  slope = CompareOutput(data_path,test_path);
  JSPEC_ASSERT_THROW( abs(slope) < 1e-27 );
    
  JSPEC_TEST_END();

  //TODO: Clean up after our test, delete the test files  
    
}

int main(int, char**)
{
 testForce();
}
