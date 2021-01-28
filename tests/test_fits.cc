#include <Base.hh>
#include <fit.h>
#include <functions.h>
#include <numeric>
#include <math.h>
#include <vector>
#include <algorithm>


//Generate fake datasets to test the gaussian and double-gaussian fits

void test_gaussian_fit(){
  JSPEC_TEST_BEGIN("Gaussian Fit\n");

  unsigned int n = 500000;
  double random_num[n];
  double sigma = 1.5;
  double avg = 0.;

  gaussian_random(n,random_num,sigma,avg);

  double fit_amplitude, fit_mean, fit_sigma;
  double n_bins = 100;

  fit *fitter = new fit();
  fitter->gaus_fit(&(random_num[0]), n, &fit_amplitude, &fit_mean, &fit_sigma,n_bins);

  std::cout << "amplitude = " << fit_amplitude << " mean " << fit_mean << " sigma " << fit_sigma << std::endl;

  JSPEC_ASSERT_THROW(abs( avg - fit_mean) < 0.01);
  JSPEC_ASSERT_THROW(abs( sigma - fit_sigma) < 0.01);


  JSPEC_TEST_END();
}


void test_double_gaussian_fit(){
  JSPEC_TEST_BEGIN("Double Gaussian Fit\n");

  unsigned int n = 100000;
  double random_num1[n];
  double sigma1 = 0.5;
  double avg1 = 0.;
  double random_num2[n];
  double sigma2 = 2.5;
  double avg2 = 0.;

  //Generate two gaussian distributions with different widths
  gaussian_random(n,random_num1,sigma1,avg1);
  gaussian_random(n,random_num2,sigma2,avg2);

  //Add them together
  double sum_distribution[n*2];
  for(int i=0;i<n;i++){
    sum_distribution[i] = random_num1[i];
    sum_distribution[i+n] = random_num2[i];
  }


  double total_amplitude, fit_amplitude1, fit_mean1, fit_sigma1, fit_amplitude2, fit_mean2, fit_sigma2;
  double n_bins = 100;

  fit *fitter = new fit();
  fitter->double_gaus_fit(&(sum_distribution[0]), n, &fit_amplitude1, &fit_mean1, &fit_sigma1,
                          &fit_amplitude2, &fit_mean2, &fit_sigma2, n_bins);

  std::cout << "amplitude = " << fit_amplitude1 << " mean " << fit_mean1 << " sigma " << fit_sigma1 << std::endl;
  std::cout << "amplitude = " << fit_amplitude2 << " mean " << fit_mean2 << " sigma " << fit_sigma2 << std::endl;

  JSPEC_TEST_END();
}

int main(int,char**){
  //  test_gaussian_fit();
  test_double_gaussian_fit();
}
