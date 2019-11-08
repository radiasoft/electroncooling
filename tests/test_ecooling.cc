#include <Base.hh>
#include <ecooling.h>
#include <ring.h>
#include <numeric>
#include <math.h>
#include <vector>
#include <algorithm>

using std::string;

void test_emit(){
    JSPEC_TEST_BEGIN( "Emittance" );

    //emit() and emit_p() are just calculating the standard
    // deviation, lets skip unit testing those. This emittance
    // function is the only one that depends on a JSPEC class.
    
    int n = 10;
    std::vector<double> dp_p = {1,2,3,4,5,6,7,8,9,10};
    std::vector<double> ds = {1,2,3,4,5,6,7,8,9,10};
    
    //define coasting gold ion beam
    // (only the value of beta passes from ion_beam to ring to emit_p)
    int n_charge = 79;
    double n_mass = 197;
    double kinetic_energy = 1e5*n_mass;
    double gamma = 1+kinetic_energy/(n_mass*k_u);
    double beta = sqrt(1-1/(gamma*gamma));
    double emit_nx0 = beta*gamma*5e-6;
    double emit_ny0 = beta*gamma*5e-6;
    double dp_p0 = 0.004; //0.0004
    double n_ptcl = 5E8;
    double sigma_s_ion = 2e-2;
    Beam ion_beam(n_charge, n_mass, kinetic_energy, emit_nx0, emit_ny0, dp_p0, sigma_s_ion, n_ptcl);

    string filename = "/home/vagrant/jupyter/github/electroncooling/data/Booster.tfs";
    Lattice lattice(filename);
    Ring ring(lattice, ion_beam);
      
    double emit_p_test = emit_p(&dp_p[0], &ds[0], ring,n);

    JSPEC_ASSERT_EQUAL(emit_p_test,8.58);
    
    JSPEC_TEST_END();
}


int main(int,char**){
    test_emit();
}