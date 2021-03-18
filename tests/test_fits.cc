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
  double chisq;
    
  gaussian_random(n,random_num,sigma,avg);

  double fit_amplitude, fit_mean, fit_sigma;
  double n_bins = 100;

  fit *fitter = new fit();
  fitter->gaus_fit(&(random_num[0]), n, &fit_amplitude, &fit_mean, &fit_sigma, &chisq, n_bins);

  std::cout << "amplitude = " << fit_amplitude << " mean " << fit_mean 
      << " sigma " << fit_sigma << " chisq "<< chisq <<std::endl;

  JSPEC_ASSERT_THROW(abs( avg - fit_mean) < 0.01);
  JSPEC_ASSERT_THROW(abs( sigma - fit_sigma) < 0.01);


  JSPEC_TEST_END();
}


void test_double_gaussian_fit(){
  JSPEC_TEST_BEGIN("Double Gaussian Fit\n");

  unsigned int n = 100000;
  double random_num1[n];
  double sigma1 = 0.05;
  double avg1 = 0.;
  double random_num2[n];
  double sigma2 = 1.5;
  double avg2 = 0.5;
  double chisq1;
    
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

  //Beacuse there is an element of randomness here, we expect to get the occasional failure here 
  // even under normal operation. Thus, run this a cuple of times to see if it still fails
  for(int i=0;i<10;i++){
      fitter->double_gaus_fit(&(sum_distribution[0]), 2*n, &fit_amplitude1, &fit_mean1, &fit_sigma1,
                              &fit_amplitude2, &fit_mean2, &fit_sigma2, &chisq1, n_bins);
      
      std::cout <<"Attempt #"<<i<<std::endl;
      std::cout << "amplitude1 = " << abs(fit_amplitude1) 
          << " mean1 " << fit_mean1 << " sigma1 " << abs(fit_sigma1) 
          << " chisq "<<chisq1<< std::endl;
      std::cout << "amplitude2 = " << abs(fit_amplitude2) 
          << " mean2 " << fit_mean2 << " sigma2 " << abs(fit_sigma2) 
          << std::endl;
      
      //we got a successful attempt, get out
      if( (abs(fit_mean1-avg1) < 0.05) & (abs(fit_mean2-avg2) < 0.05) &
          (abs(fit_sigma1-sigma1) < 0.05) & abs(fit_sigma2-sigma2) < 0.05)
      {
          break;
      }
  }
    
  //If we got here, none of the attempts were successfull. Throw the error.
  JSPEC_ASSERT_THROW(abs(fit_mean1-avg1) < 0.05);
  JSPEC_ASSERT_THROW(abs(fit_mean2-avg2) < 0.05);
  JSPEC_ASSERT_THROW(abs(fit_sigma1-sigma1) < 0.05);
  JSPEC_ASSERT_THROW(abs(fit_sigma2-sigma2) < 0.05);

  JSPEC_TEST_END();
}

void test_cross_gaussian_fit(){
  JSPEC_TEST_BEGIN("Cross Gaussian Fit\n");

  unsigned int n = 100000;
  double random_num1[n];
  double sigma1 = 0.05;
  double avg1 = 0.;
  double random_num2[n];
  double sigma2 = 1.5;
  double avg2 = 0.5;
  double chisq1, chisq2;
    
  //Generate two gaussian distributions with different widths
  gaussian_random(n,random_num1,sigma1,avg1);
  gaussian_random(n,random_num2,sigma2,avg2);

  //Add them together
  double sum_distribution[n*2];
  for(int i=0;i<n;i++){
    sum_distribution[i] = random_num1[i];
    sum_distribution[i+n] = random_num2[i];
  }


  double total_amplitude, fit_amplitude1, fit_mean1, fit_sigma1, fit_amplitude2, fit_mean2, fit_sigma2, fit_amplitude,fit_mean,fit_sigma;
  double n_bins = 100;

  fit *fitter = new fit();

  //Beacuse there is an element of randomness here, we expect to get the occasional failure here 
  // even under normal operation. Thus, run this a cuple of times to see if it still fails
  for(int i=0;i<10;i++){
      fitter->double_gaus_fit(&(sum_distribution[0]), 2*n, &fit_amplitude1, &fit_mean1, &fit_sigma1,
                              &fit_amplitude2, &fit_mean2, &fit_sigma2, &chisq1, n_bins);
      
        fit *fitter = new fit();
          fitter->gaus_fit(&(sum_distribution[0]), n, &fit_amplitude, &fit_mean, &fit_sigma, &chisq2, n_bins);

      
      std::cout <<"Attempt #"<<i<<std::endl;
      std::cout << "amplitude1 = " << abs(fit_amplitude1) << " mean1 " << fit_mean1 << " sigma1 " << abs(fit_sigma1) << " chisq "<< chisq1 << std::endl;
      std::cout << "amplitude2 = " << abs(fit_amplitude2) << " mean2 " << fit_mean2 << " sigma2 " << abs(fit_sigma2) << std::endl;
      std::cout << "amplitude = " << abs(fit_amplitude) << " mean " << fit_mean << " sigma " << abs(fit_sigma) << " chisq "<<chisq2 << std::endl;
      
      //we got a successful attempt, get out
      if( (abs(fit_mean1-avg1) < 0.05) & (abs(fit_mean2-avg2) < 0.05) &
          (abs(fit_sigma1-sigma1) < 0.05) & abs(fit_sigma2-sigma2) < 0.05)
      {
          break;
      }
  }
}


int main(int,char**){
  test_gaussian_fit();
  test_double_gaussian_fit();
  test_cross_gaussian_fit();
}
