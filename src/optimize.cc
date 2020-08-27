#include "optimize.h"

extern DynamicParas * dynamic_paras;
extern IBSSolver * ibs_solver;
//extern EcoolRateParas * ecool_paras;
//extern ForceParas * force_paras;
extern Luminosity *luminosity_paras;

void Optimize::InitializeFitter(std::vector<std::string> Params, std::vector<double> IV, Lattice *lattice, Beam *ion){

    assert(Params.size() != IV.size() && "Mismatch in optimization parameters!");
    
    FitVariables = Params; //a deep copy, for storage
    
    for(int i=0;i<Params.size();i++) InitialValues[FitVariables[i]] = IV[i];
    
    //Hard-coding this here, temporarily(ish)
    // Tuning these parameters takes some care. If the step size
    // is too small, the optimizer will slow-walk that parameter.
    // If it's too big, it might step over the minimum.
    FitStepSize.clear();
    FitStepSize["bfield"]     = 0.1;
    FitStepSize["sigma_x"]    = 5e-5;
    FitStepSize["sigma_y"]    = 5e-5;
    FitStepSize["sigma_s"]    = 5e-3;
    FitStepSize["beta_h"]     = 0.5;
    FitStepSize["beta_v"]     = 0.5;
    FitStepSize["alpha_h"]    = 0.5;
    FitStepSize["alpha_v"]    = 0.5;
    FitStepSize["disp_h"]     = 0.1;
    FitStepSize["disp_v"]     = 0.1;
    FitStepSize["disp_der_h"] = 0.1;
    FitStepSize["disp_der_v"] = 0.1;
    FitStepSize["temp_tr"]    = 0.001;
    FitStepSize["temp_long"]  = 0.001;
    FitStepSize["n_electron"] = 0.1;
    
    //Parse the lattice file once, store it for each 
    // call to fit_fcn
    fitter_values.lattice = lattice;

    //Parse the ion beam that's pre-defined
    fitter_values.beam = ion;
    
    //Initialize other classes that don't change
    fitter_values.ecool_paras = new EcoolRateParas(fitter_values.n_sample);
    fitter_values.force_paras = ChooseForce(fitter_values.ff);//ForceFormula::PARKHOMCHUK);
    
    //Now that we've stored the user's input values,
    // copy to a 'working' version that the optimizing function will use
    fitter_values.FitVariables_working = FitVariables;
    fitter_values.FitStepSize_working  = FitStepSize;
    fitter_values.InitialValues_working = InitialValues;
}

double Optimize::fit_fcn(const gsl_vector *v, void *params){

    //Load in the struct holding our run parameters
    struct opt_info *p = (struct opt_info *)params;
    
    //Depending on the values that are passed to the 
    // optimizer, overwrite our set values
    
    //Fixed values for params come from the object...
    double magnetic_field = p->magnetic_field_;
    double sigma_x  = p->sigma_x_; 
    double sigma_y  = p->sigma_y_; 
    double sigma_s  = p->sigma_s_; 
    double beta_h   = p->beta_h_; 
    double beta_v   = p->beta_v_; 
    double alpha_h  = p->alpha_h_; 
    double alpha_v  = p->alpha_v_; 
    double dis_h    = p->disp_h_; 
    double dis_v    = p->disp_v_; 
    double dis_der_h= p->disp_der_h_; 
    double dis_der_v= p->disp_der_v_; 
    double tmp_tr   = p->temp_tr_;
    double tmp_long = p->temp_long_;
    double n_electron = p->n_electron_;

    //Reminder: Some other parameters are set directly from the pointer below (like m0 and Z)
    
    //Grab our list of parameters that we allow to vary....
    std::map<const std::string, double> IV = p->InitialValues_working;
    std::vector<std::string> FV = p->FitVariables_working;
    
    // Then overwrite the set values with the values the fitter is trying:
    if(IV.find("bfield")     != IV.end() ) magnetic_field = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "bfield")));
    if(IV.find("sigma_x")    != IV.end() )        sigma_x = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "sigma_x"))); 
    if(IV.find("sigma_y")    != IV.end() )        sigma_y = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "sigma_y")));
    if(IV.find("sigma_s")    != IV.end() )        sigma_s = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "sigma_s"))); 
    if(IV.find("beta_h")     != IV.end() )         beta_h = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "beta_h")));
    if(IV.find("beta_v")     != IV.end() )         beta_v = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "beta_v")));
    if(IV.find("alpha_h")    != IV.end() )        alpha_h = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "alpha_h")));
    if(IV.find("alpha_v")    != IV.end() )        alpha_v = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "alpha_v")));
    if(IV.find("disp_h")     != IV.end() )          dis_h = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "disp_h"))); 
    if(IV.find("disp_v")     != IV.end() )          dis_v = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "disp_v"))); 
    if(IV.find("disp_der_h") != IV.end() )      dis_der_h = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "disp_der_h"))); 
    if(IV.find("disp_der_v") != IV.end() )      dis_der_v = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "disp_der_v"))); 
    if(IV.find("temp_tr")    != IV.end() )         tmp_tr = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "temp_tr"))); 
    if(IV.find("temp_long")  != IV.end() )       tmp_long = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "temp_long"))); 
    if(IV.find("n_electron") != IV.end() )     n_electron = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "n_electron")));
       
    
    //Punish unphysical parameters
    if (magnetic_field < 0.0 ||
        sigma_x < 0.0 || sigma_y < 0.0 || sigma_s < 0.0 ||
        tmp_tr < 1e-6 || tmp_long < 1e-6 || n_electron > 5.0){
        std::cout<<"Unphysical"<<std::endl;
        return 100000.0;
    }
    
    
    //We're really only varying the scale factor
    // TODO: Convert this to nanocoulombs
    n_electron = n_electron * 1e10; 
    
    // define ion beam;
    /*
    double KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
    KE = 2.5e4;
    emit_nx0 = 2.5e-6;
    emit_ny0 = 2.5e-6;
    dp_p0= 1e-4;
    sigma_s0 = 0.7;
    N_ptcl = 1.34e12;
    Beam p_beam(p->Z_,
                p->m0_/k_u, 
                KE, 
                emit_nx0, emit_ny0, 
                dp_p0, sigma_s0, N_ptcl);
*/
    
    // The lattice is defined once in initialization & we just pass the reference around
    Lattice lattice = *(p->lattice);

    //Define the ring
//    Ring ring(lattice, p_beam);
    Ring ring(lattice,*(p->beam));
    
    //Set IBS parameters.
    int nu = 100;
    int nv = 100;
    int nz = 40;
    double log_c = 20.6;
    ibs_solver = new IBSSolver_Martini(nu, nv, nz, log_c, 0);

    //define the cooler
    double cooler_length = p->length_;
    double n_section = p->section_number_;
    //Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);
    Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v,alpha_h,alpha_v,dis_der_h,dis_der_v);
    
    //define electron beam
    GaussianBunch gaussian_bunch(n_electron, sigma_x, sigma_y, sigma_s);
    double gamma_e = p->beam->gamma();
    EBeam e_beam(gamma_e, tmp_tr, tmp_long, gaussian_bunch);

    double rate_x, rate_y, rate_s;
    ecooling_rate(*(p->ecool_paras), *(p->force_paras),
                  *(p->beam), cooler, e_beam, ring,
                  rate_x, rate_y, rate_s);
     
    //NOTE: We're only cooling in the y direction here
    // Chuck it if we get a heating rate (rate_y>0)
    if(rate_y > 0.0 || isnan(rate_y)){
        return 100000.0;
    }
    else{
        return abs(p->cool_target_ - ((-1./(rate_y))/60.));
    }
}


void Optimize::Randomize()
{           
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);

    //TODO: Stash these RMS values in a Map too
    std::normal_distribution<double> bfield_distribution  (InitialValues["bfield"],3.0);
    std::normal_distribution<double> sigma_distribution   (InitialValues["sigma_x"],1e-4);
    std::normal_distribution<double> sigma_s_distribution (InitialValues["sigma_s"],0.2);
    std::normal_distribution<double> beta_distribution    (InitialValues["beta_h"],10.);
    std::normal_distribution<double> alpha_distribution   (InitialValues["alpha_h"],10.);
    std::normal_distribution<double> disp_distribution    (InitialValues["disp_h"],0.3);
    std::normal_distribution<double> disp_der_distribution(InitialValues["disp_der_h"],0.3);
    std::normal_distribution<double> temp_distribution    (InitialValues["temp_tr"],0.01);
    std::normal_distribution<double> e_distribution       (InitialValues["n_electron"],0.5);
    
    //Burn-in the random number generators
    for(int i=0;i<500;i++){
        bfield_distribution(generator);
        sigma_distribution(generator);
        sigma_s_distribution(generator);
        beta_distribution(generator);
        alpha_distribution(generator);
        disp_distribution(generator);
        temp_distribution(generator);
        e_distribution(generator);
    }

    //Select random values, but add some constraints
    if(std::find(FitVariables.begin(),FitVariables.end(),"bfield") != FitVariables.end() ){
        double bfield_tmp = bfield_distribution(generator);
        while(bfield_tmp <= 0.) bfield_tmp = bfield_distribution(generator);
        fitter_values.magnetic_field_ = bfield_tmp;
    }
    
    if(std::find(FitVariables.begin(),FitVariables.end(),"sigma_x") != FitVariables.end() ){
        double sigma_x_tmp = sigma_distribution(generator);
        while(sigma_x_tmp <= 1e-5) sigma_x_tmp = sigma_distribution(generator);
        fitter_values.sigma_x_ = sigma_x_tmp;
    }

    if(std::find(FitVariables.begin(),FitVariables.end(),"sigma_y") != FitVariables.end() ){
        double sigma_y_tmp = sigma_distribution(generator);
        while(sigma_y_tmp <= 1e-5) sigma_y_tmp = sigma_distribution(generator);
        fitter_values.sigma_y_ = sigma_y_tmp;
    }

    if(std::find(FitVariables.begin(),FitVariables.end(),"sigma_s") != FitVariables.end() ){
        double sigma_s_tmp = sigma_s_distribution(generator);
        while(sigma_s_tmp <= 0.) sigma_s_tmp = sigma_s_distribution(generator);
        fitter_values.sigma_s_ = sigma_s_tmp;
    }
    
    if(std::find(FitVariables.begin(),FitVariables.end(),"beta_h") != FitVariables.end() ){
        double beta_h_tmp = beta_distribution(generator);
        while(beta_h_tmp <= 0.) beta_h_tmp = beta_distribution(generator);
        fitter_values.beta_h_ = beta_h_tmp;
    }
    
    if(std::find(FitVariables.begin(),FitVariables.end(),"beta_v") != FitVariables.end() ){
        double beta_v_tmp = beta_distribution(generator);
        while(beta_v_tmp <= 0.) beta_v_tmp = beta_distribution(generator);
        fitter_values.beta_v_ = beta_v_tmp;
    }

    if(std::find(FitVariables.begin(),FitVariables.end(),"alpha_h") != FitVariables.end() ){
        double alpha_h_tmp = alpha_distribution(generator);
        while(alpha_h_tmp <= 0.) alpha_h_tmp = alpha_distribution(generator);
        fitter_values.alpha_h_ = alpha_h_tmp;
    }
    
    if(std::find(FitVariables.begin(),FitVariables.end(),"alpha_v") != FitVariables.end() ){
        double alpha_v_tmp = alpha_distribution(generator);
        while(alpha_v_tmp <= 0.) alpha_v_tmp = alpha_distribution(generator);
        fitter_values.alpha_v_ = alpha_v_tmp;
    }

    
    if(std::find(FitVariables.begin(),FitVariables.end(),"disp_h") != FitVariables.end() ){
        double disp_h_tmp = disp_distribution(generator);
        while(disp_h_tmp <= 1e-5) disp_h_tmp = disp_distribution(generator);
        fitter_values.disp_h_ = disp_h_tmp;
    }
    
    if(std::find(FitVariables.begin(),FitVariables.end(),"disp_v") != FitVariables.end() ){
        double disp_v_tmp = disp_distribution(generator);
        while(disp_v_tmp <= 1e-5) disp_v_tmp = disp_distribution(generator);
        fitter_values.disp_v_ = disp_v_tmp;
    }
    
    if(std::find(FitVariables.begin(),FitVariables.end(),"disp_der_h") != FitVariables.end() ){
        double disp_der_h_tmp = disp_der_distribution(generator);
        while(disp_der_h_tmp <= 1e-5) disp_der_h_tmp = disp_der_distribution(generator);
        fitter_values.disp_der_h_ = disp_der_h_tmp;
    }
    
    if(std::find(FitVariables.begin(),FitVariables.end(),"disp_der_v") != FitVariables.end() ){
        double disp_der_v_tmp = disp_der_distribution(generator);
        while(disp_der_v_tmp <= 1e-5) disp_der_v_tmp = disp_der_distribution(generator);
        fitter_values.disp_der_v_ = disp_der_v_tmp;
    }

    if(std::find(FitVariables.begin(),FitVariables.end(),"temp_tr") != FitVariables.end() ){
        double temp_tr_tmp = temp_distribution(generator);
        while(temp_tr_tmp < 1e-5) temp_tr_tmp = temp_distribution(generator);
        fitter_values.temp_tr_ = temp_tr_tmp;
    }

    if(std::find(FitVariables.begin(),FitVariables.end(),"temp_long") != FitVariables.end() ){
        double temp_long_tmp = temp_distribution(generator);
        while(temp_long_tmp < 1e-5) temp_long_tmp = temp_distribution(generator);
        fitter_values.temp_long_ = temp_long_tmp;
    }
    
    if(std::find(FitVariables.begin(),FitVariables.end(),"n_electron") != FitVariables.end() ){
        double n_electron_tmp = e_distribution(generator);
        while( n_electron_tmp <= 0) n_electron_tmp = e_distribution(generator);
        fitter_values.n_electron_ = n_electron_tmp; 
    }
    /*
    #ifndef NDEBUG
    std::cout<<"Starting values:"<<std::endl;
    printf ("Bfield = %.5f sigma_x = %.5f sigma_y = %.5f sigma_s = %.5f beta_h = %.5f beta_v = %.5f \n",
            bfield_tmp, sigma_x_tmp,sigma_y_tmp,sigma_z_tmp,beta_h_tmp,beta_v_tmp);       
    printf( "disp_h = %.5f disp_v = %.5f temp_tr = %.5f temp_long = %.5f n_electrons = %.5 e10 \n", 
           disp_h_tmp,disp_v_tmp,temp_tr_tmp,temp_long_tmp,n_electron_tmp);
    #endif
    */

}

void Optimize::OptimizeTrial(){
                
    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;

    gsl_vector *x;
    gsl_vector *step_size;

    gsl_multimin_function my_func;

    //TODO: Implement some checks here so GSL doesn't coredump
    
    int n_pars = FitVariables.size();
    my_func.n = n_pars;
    my_func.params = &fitter_values;
    my_func.f = &Optimize::fit_fcn;

    T = gsl_multimin_fminimizer_nmsimplex2rand;
    s = gsl_multimin_fminimizer_alloc (T, n_pars);
    
    try{
        //Set the best guess
        x = gsl_vector_alloc (n_pars);
        step_size = gsl_vector_alloc(n_pars);

        for(int i = 0; i < n_pars; i++){
            gsl_vector_set (x, i, fitter_values.InitialValues_working[FitVariables[i]] );
            gsl_vector_set (step_size, i, FitStepSize[FitVariables[i]] );
        }
                       
        gsl_multimin_fminimizer_set (s, &my_func, x, step_size);
       
        int iter = 0;
        int status;
                
        do{
            iter++;
            status = gsl_multimin_fminimizer_iterate (s);

            if(iter % 20 == 0 ){
                std::cout<<iter;

                for(int i=0; i < n_pars; i++){
                    std::cout<<" "<<FitVariables[i]<<": "<<std::setprecision(5)<<gsl_vector_get (s->x,i);
                }
                std::cout<<" Score: "<<s->fval<<std::endl;
            }
        } while(iter<125);
//                while (status == GSL_CONTINUE && iter < 100);

        if(s->fval < best_fval){
            for(int i = 0;i < n_pars; i++) BestFits[FitVariables[i]] = gsl_vector_get(s->x, i);            
            best_fval = s->fval;
        }
    }
    catch (...) { 
        gsl_multimin_fminimizer_free (s);
        gsl_vector_free (x);
        return;
    }

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    return;
}
    
    
void Optimize::ManyTrials()
{                   
    for(int i = 0;i<n_trials;i++){
        std::cout<<"Random Iteration "<<i<<std::endl;
        Randomize();
        OptimizeTrial();
    }                 
    
    std::cout<<"Best values:"<<std::endl;
    for(int i=0; i < BestFits.size(); i++){
       std::cout<<" "<<FitVariables[i]<<": "<<std::setprecision(5)<<BestFits[FitVariables[i]];
    }
    std::cout<<" Score: "<<best_fval<<std::endl;
                
    //std::ofstream myfile;
    //myfile.open ("bests.txt");
    //myfile << best_b<<" "<< best_sigma_x<<" "<<best_sigma_y<<" "<<best_sigma_z<<" "<<
    //        best_beta_h<<" "<<best_beta_v<<" "<< best_disp_h<<" "<<best_disp_v<<" "<<
    //        best_temp_tr<<" "<<best_temp_long<<" "<< best_n_electron<<" "<< best_fval<<"\n";
    //myfile.close();    
}


int Optimize::Optimize_From_UI(std::vector<std::string> Params, std::vector<double>InitialValues, Beam &ion, Cooler &cooler, EBeam &ebeam, Ring &ring){
    
    
    Lattice *lattice = ring.lattice_; //"eRHIC.tfs"
    
    //Params is a vector of string ID's for parameters
    //InitialValues is a vector of doubles, matched 1:1 with the params

    
    this->InitializeFitter(Params, InitialValues, lattice, &ion);
    this->ManyTrials();
        
    return 1;
}