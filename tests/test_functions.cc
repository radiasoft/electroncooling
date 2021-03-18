#include <Base.hh>
#include <functions.h>
#include <numeric>
#include <math.h>
#include <vector>
#include <algorithm>

//Testing against a random number generator makes for interesting comparisons.
// Mostly, we want to make sure the random seed is working and the
// shifting functions don't violate the initial constraints of the distributions.

void test_gaussian_random(){
    JSPEC_TEST_BEGIN( "Gaussian Random" );

    int n = 50000;
    double random_num[n];
    double sigma = 1.5;
    double avg = 0.;

    gaussian_random(n,random_num,sigma,avg);

    std::vector<double> v(random_num, random_num + n);
    double sum = std::accumulate(std::begin(v), std::end(v), 0.0);
    double m =  sum / v.size();

    double accum = 0.0;
    std::for_each (std::begin(v), std::end(v), [&](const double d) {
        accum += (d - m) * (d - m);
    });

    double stdev = sqrt(accum / (v.size()-1.0));

    std::cout<<std::endl;
    std::cout<<"Require delta mean < 0.01, delta_sigma < 2%"<<std::endl;
    std::cout<<"delta mean = "<<abs(m-avg)<<" delta sigma = "<<abs(stdev-sigma)/stdev<<"%"<<std::endl;

    JSPEC_ASSERT_THROW( abs(m - avg) < 0.01 );
    JSPEC_ASSERT_THROW(abs(stdev - sigma)/stdev < 0.02 );

    //If the seed isn't really random, a second draw
    // will be exactly the same with the same moments

    double random_num2[n];
    gaussian_random(n,random_num2,sigma,avg);

    std::vector<double> v2(random_num2, random_num2 + n);
    double sum2 = std::accumulate(std::begin(v2), std::end(v2), 0.0);
    double m2 =  sum2 / v2.size();

    double accum2 = 0.0;
    std::for_each (std::begin(v2), std::end(v2), [&](const double d) {
        accum2 += (d - m2) * (d - m2);
    });

    double stdev2 = sqrt(accum2 / (v2.size()-1.0));

    std::cout<<"Require delta mean != 0, delta sigma != 0"<<std::endl;
    std::cout<<"delta mean = "<<abs(m - m2)<<" delta sigma = "<<abs(stdev - stdev2)<<std::endl;

    JSPEC_ASSERT_THROW(m != m2);
    JSPEC_ASSERT_THROW(stdev != stdev2);

    JSPEC_TEST_END();
}

void test_uniform_random(){
    JSPEC_TEST_BEGIN( "Uniform Random" );

    int n = 50000;
    double random_num[n];
    double r_min = -10.;
    double r_max = 10.;

    uniform_random(n, random_num, r_min,r_max);

    std::vector<double> v(random_num, random_num + n);

    double max = *std::max_element(v.begin(),v.end());
    double min = *std::min_element(v.begin(),v.end());

    std::cout<<std::endl;
    std::cout<<"Require max < r_max, min > r_min"<<std::endl;
    std::cout<<"max = "<<max<<" min = "<<min<<" r_max = "<<r_max<<" r_min = "<<r_min<<std::endl;

    JSPEC_ASSERT_THROW( max < r_max);
    JSPEC_ASSERT_THROW( min > r_min);

    //If the seed isn't really random, a second draw
    // will be exactly the same with the same moments

    double random_num2[n];
    uniform_random(n,random_num2,r_min,r_max);

    double sum = std::accumulate(std::begin(v), std::end(v), 0.0);
    double m =  sum / v.size();

    double accum = 0.0;
    std::for_each (std::begin(v), std::end(v), [&](const double d) {
        accum += (d - m) * (d - m);
    });

    double stdev = sqrt(accum / (v.size()-1.0));

    std::vector<double> v2(random_num2, random_num2 + n);
    double sum2 = std::accumulate(std::begin(v2), std::end(v2), 0.0);
    double m2 =  sum2 / v2.size();

    double accum2 = 0.0;
    std::for_each (std::begin(v2), std::end(v2), [&](const double d) {
        accum2 += (d - m2) * (d - m2);
    });

    double stdev2 = sqrt(accum2 / (v2.size()-1.0));

    std::cout<<"Require delta mean != 0, delta sigma != 0"<<std::endl;
    std::cout<<"delta mean = "<<abs(m - m2)<<" delta sigma = "<<abs(stdev - stdev2)<<std::endl;

    JSPEC_ASSERT_THROW(m != m2);
    JSPEC_ASSERT_THROW(stdev != stdev2);

    JSPEC_TEST_END();
}

void test_gaussian_random_adjust(){
    JSPEC_TEST_BEGIN( "Gaussian Adjust" );
    int n = 50000;
    double random_num[n];
    double sigma = 1.5;
    double avg = 0.;

    gaussian_random(n,random_num,sigma,avg);

    sigma = 3.5;
    avg = 1.0;

    gaussian_random_adjust(n,random_num,sigma,avg);

    std::vector<double> v(random_num, random_num + n);
    double sum = std::accumulate(std::begin(v), std::end(v), 0.0);
    double m =  sum / v.size();

    double accum = 0.0;
    std::for_each (std::begin(v), std::end(v), [&](const double d) {
        accum += (d - m) * (d - m);
    });

    double stdev = sqrt(accum / (v.size()-1.0));

    std::cout<<std::endl;
    std::cout<<" Require delta mean < 0.01, delta sigma < 2%"<<std::endl;
    std::cout<<"delta mean = "<<abs(m-avg)<<" delta sigma = "<<abs(stdev-sigma)/stdev<<"%"<<std::endl;

    JSPEC_ASSERT_THROW( abs(m - avg) < 0.01 );
    JSPEC_ASSERT_THROW( abs(stdev - sigma)/stdev < 0.02 );

    JSPEC_TEST_END();
}

void test_uniform_random_adjust(){
    JSPEC_TEST_BEGIN( "Uniform Adjust" );

    int n = 50000;
    double random_num[n];
    double r_min = -10.;
    double r_max = 10.;

    uniform_random(n, random_num, r_min,r_max);

    double avg_shift = 2.5;

    uniform_random_adjust(n, random_num, avg_shift);

    std::vector<double> v(random_num, random_num + n);

    double max = *std::max_element(v.begin(),v.end());
    double min = *std::min_element(v.begin(),v.end());

    std::cout<<std::endl;
    std::cout<<"Require min + 2.5 > shifted_min, max + 2.5 < shifted_max"<<std::endl;
    std::cout<<"Min = "<<min<<" Max = "<<max<<std::endl;
    std::cout<<"Shifted min = "<<(r_min + avg_shift)<<" Shifted max = "<<(r_max + avg_shift)<<std::endl;

    JSPEC_ASSERT_THROW( max < (r_max + avg_shift) );
    JSPEC_ASSERT_THROW( min > (r_min + avg_shift) );

    JSPEC_TEST_END();
}

void test_rd(){
    JSPEC_TEST_BEGIN( "Carlson Integral (untested)" );

    //TODO: Choose some values from a real use-case

    double test_val = rd(4,5,6);

    JSPEC_ASSERT_THROW((test_val>0.0807146) && (test_val < 0.0807148));

    // The integral carries some specific constraints
    JSPEC_EXPECT_EXCEPTION(rd(4,5,1e-50),std::invalid_argument);
    JSPEC_EXPECT_EXCEPTION(rd(-4,5,6),std::invalid_argument);
    JSPEC_EXPECT_EXCEPTION(rd(1e-50,1e-50,6),std::invalid_argument);

    JSPEC_TEST_END();
}


int main(int,char**){
    test_gaussian_random();
    test_uniform_random();
    test_gaussian_random_adjust();
    test_uniform_random_adjust(); //Need to fix the main function
    test_rd(); //Undefined, carlson elliptical integral
}
