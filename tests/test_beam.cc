#include <Base.hh>
#include "beam.h"
#include <numeric>
#include <math.h>
#include <vector>
#include <algorithm>

void test_beam(){
    JSPEC_TEST_BEGIN( "Beam" );

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

    //Nothing really to test for a constructor...
    
    

    JSPEC_TEST_END();
}

void test_uniform_cylinder(){
    JSPEC_TEST_BEGIN( "Uniform Cylinder" );

    

    JSPEC_TEST_END();
}

int main(int,char**){
    test_beam();
    test_uniform_cylinder();
    //test_uniform_hollow();
    //test_uniform_bunch();
    //test_elliptic_uniform_bunch();
    //test_EBeam();
    //test_GaussianBunch();
    //test_ParticleBunch();
    
    
    
}
