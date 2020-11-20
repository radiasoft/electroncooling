#include "optimize.h"

extern DynamicParas * dynamic_paras;
extern IBSSolver * ibs_solver;
//extern EcoolRateParas * ecool_paras;
//extern ForceParas * force_paras;
extern Luminosity *luminosity_paras;

void Optimize::InitializeFitter(std::vector<std::string> Params, std::vector<double> IV, Lattice *lattice, Beam *ion, Cooler *cooler, EBeam *ebeam, ForceFormula ff){

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
    FitStepSize["log_c"]      = 0.2;

    //Store the values that were pre-set in the input file, before
    // the optimization section was defined.
    // If they were re-defined in the optimization section they'll
    // be over-written
    fitter_values.magnetic_field_ = cooler->magnetic_field();
    fitter_values.beta_h_         = cooler->beta_h();
    fitter_values.beta_v_         = cooler->beta_v();
    fitter_values.alpha_h_        = cooler->alpha_h();
    fitter_values.alpha_v_        = cooler->alpha_v();
    fitter_values.disp_h_         = cooler->disp_h();
    fitter_values.disp_v_         = cooler->disp_v();
    fitter_values.disp_der_h_     = cooler->der_disp_h();
    fitter_values.disp_der_v_     = cooler->der_disp_v();
    fitter_values.temp_tr_        = ebeam->tmp_tr();
    fitter_values.temp_long_      = ebeam->tmp_long();
    fitter_values.sigma_x_        = ebeam->sigma_x();
    fitter_values.sigma_y_        = ebeam->sigma_y();
    fitter_values.sigma_s_        = ebeam->length(); //Assuming it's bunched!
    //Assuming Gaussian bunch!
    // Also, to conform the n_electron pulled from the object (no scaling)
    // with the n_electron we typically get from n_electron in the optimization
    // section (scaled 1/1e10), we scale down the value here.
    fitter_values.n_electron_     = ebeam->n_electron() / 1e10; 
    
    //Parse the lattice file once, store it in memory for each
    // subsequent call to fit_fcn
    fitter_values.lattice = lattice;

    //Parse the ion beam that's pre-defined
    fitter_values.beam = ion;
    fitter_values.ff = ff;    

    //Initialize other classes that don't (yet) have floating parameters
    fitter_values.ecool_paras = new EcoolRateParas(fitter_values.n_sample);
    fitter_values.force_paras = ChooseForce(fitter_values.ff);

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
    double n_electron = p->n_electron_; //scaled by 1/1e10
    double log_c    = p->log_c_;

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
    if(IV.find("log_c") != IV.end() )               log_c = gsl_vector_get(v,std::distance(FV.begin(),std::find(FV.begin(), FV.end(), "log_c")));
       
    
    //Punish unphysical parameters
    if (magnetic_field < 0.0 ||
        sigma_x < 0.0 || sigma_y < 0.0 || sigma_s < 0.0 ||
        tmp_tr < 1e-6 || tmp_long < 1e-6 || n_electron > 5.0){
        std::cout<<"Unphysical"<<std::endl;
        return DBL_MAX;
    }
    
    
    //We're really only varying the scale factor
    // TODO: Convert this to nanocoulombs? (6.242e9 electrons/nC)
    n_electron = n_electron * 1e10; 
    //n_electron = n_electron * 6.242e9; 
    
    //Punish unphysical parameters
    if (magnetic_field < 0.0 ||
        sigma_x < 0.0 || sigma_y < 0.0 || sigma_s < 0.0 ||
        tmp_tr < 1e-6 || tmp_long < 1e-6 || 
        n_electron < 0.0 || n_electron > 10.0 ||
        log_c < 0){
        std::cout<<"Unphysical"<<std::endl;
        return DBL_MAX;
    }

    //We're really only varying the scale factor
    // TODO: Convert this to nanocoulombs? (6.242e9 electrons/nC)
    n_electron = n_electron * 1e10;

    // The lattice is defined once in initialization & we just pass the reference around
    Lattice lattice = *(p->lattice);

    //Define the ring
    Ring ring(lattice,*(p->beam));
    
    //Set IBS parameters.
    int nu = 100;
    int nv = 100;
    int nz = 40;
    ibs_solver = new IBSSolver_Martini(nu, nv, nz, log_c, 0);

    //define the cooler
    double cooler_length = p->length_;
    double n_section = p->section_number_;
    Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v,alpha_h,alpha_v,dis_der_h,dis_der_v);

    //define electron beam
    GaussianBunch gaussian_bunch(n_electron, sigma_x, sigma_y, sigma_s);
    double gamma_e = p->beam->gamma();
    EBeam e_beam(gamma_e, tmp_tr, tmp_long, gaussian_bunch);
    
    //Calculate the delta rate from IBS
    double rx_ibs, ry_ibs, rz_ibs;
    ibs_solver->rate(lattice, *(p->beam), rx_ibs, ry_ibs, rz_ibs);
    //std::cout<<"IBS:" << rx_ibs<<" "<<ry_ibs<<" "<<rz_ibs<<std::endl;
    
    double rate_x, rate_y, rate_s;
    ecooling_rate(*(p->ecool_paras), *(p->force_paras),
                  *(p->beam), cooler, e_beam, ring,
                  rate_x, rate_y, rate_s);

    //std::cout<<rate_x<<" "<<rate_y<<" "<<rate_s<<std::endl;

    double total_rate_x = rate_x + rx_ibs;
    double total_rate_y = rate_y + ry_ibs;
    double total_rate_z = rate_s + rz_ibs;
        
    //NOTE: We're only optimizing cooling in the y direction here (so far)
    // Chuck the trial if we get a heating rate (rate_y>0)
    if(total_rate_y > 0.0 || isnan(total_rate_y)){
        return DBL_MAX;
    }
    else{
        return abs(p->cool_target_ - ((-1./(total_rate_y))/60.));
    }
}

//Acts like the np.linspace function in Python/Numpy
template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num);

  for(int i=0; i < num; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}


std::map<int, vector<double> >Optimize::ParameterScan(string scan_par, double par_min, double par_max, int n_steps, opt_info params)
{
    //Fix all but one values, and scan the rates returned from varying that parameter
    
    //Fixed values for params come from the object...
    double magnetic_field = params.magnetic_field_;
    double sigma_x  = params.sigma_x_;
    double sigma_y  = params.sigma_y_;
    double sigma_s  = params.sigma_s_;
    double beta_h   = params.beta_h_;
    double beta_v   = params.beta_v_;
    double alpha_h  = params.alpha_h_;
    double alpha_v  = params.alpha_v_;
    double dis_h    = params.disp_h_;
    double dis_v    = params.disp_v_;
    double dis_der_h= params.disp_der_h_;
    double dis_der_v= params.disp_der_v_;
    double tmp_tr   = params.temp_tr_;
    double tmp_long = params.temp_long_;
    double n_electron = params.n_electron_*1e10; //Was already scaled by 1/1e10 (see InitializeFitter)
    double log_c    = params.log_c_;

    //Reminder: Some other parameters are set directly from the pointer below (like m0 and Z)

    //get the range of our parameter scan once 
    std::vector<double> scan_samples = linspace(par_min,par_max,n_steps);
    
    //for(int i=0;i<scan_samples.size();i++){
    //    std::cout<<scan_samples[i]<<" ";
    //}
    //std::cout<<std::endl;
    
    std::map< int,std::vector<double> > results;

    for(unsigned int i=0;i<n_steps;++i){
        // Then overwrite the set value of the scanned variable
        if     (scan_par == "bfield" )      magnetic_field = scan_samples[i];
        else if(scan_par == "sigma_x")      sigma_x        = scan_samples[i];
        else if(scan_par == "sigma_y")      sigma_y        = scan_samples[i];
        else if(scan_par == "sigma_s")      sigma_s        = scan_samples[i];
        else if(scan_par == "beta_h")       beta_h         = scan_samples[i];
        else if(scan_par == "beta_v")       beta_v         = scan_samples[i];
        else if(scan_par == "alpha_h")      alpha_h        = scan_samples[i];
        else if(scan_par == "alpha_v")      alpha_v        = scan_samples[i];
        else if(scan_par == "disp_h")       dis_h          = scan_samples[i];
        else if(scan_par == "disp_v")       dis_v          = scan_samples[i];
        else if(scan_par == "disp_der_h")   dis_der_h      = scan_samples[i];
        else if(scan_par == "disp_der_v")   dis_der_v      = scan_samples[i];
        else if(scan_par == "temp_tr")      tmp_tr         = scan_samples[i];
        else if(scan_par == "temp_long")    tmp_long       = scan_samples[i];
        else if(scan_par == "n_electron")   n_electron     = scan_samples[i]*1e10; // We scale it by 1e10 
        else if(scan_par == "log_c")        log_c          = scan_samples[i];
        else{
            std::cout<<"I don't recognize variable "<<scan_par<<std::endl;
        }
        
        std::cout<<"Testing. "<<n_electron<<std::endl;

        // The lattice is defined once in initialization & we just pass the reference around
        Lattice lattice = *(params.lattice);

        //Define the ring
        Ring ring(lattice,*(params.beam));

        //Set IBS parameters.
        int nu = 100;
        int nv = 100;
        int nz = 40;
        ibs_solver = new IBSSolver_Martini(nu, nv, nz, log_c, 0);

        //define the cooler
        double cooler_length = params.length_;
        double n_section = params.section_number_;
        Cooler cooler(cooler_length,n_section,magnetic_field,
                      beta_h,beta_v,
                      dis_h, dis_v,
                      alpha_h,alpha_v,
                      dis_der_h,dis_der_v);

        //define electron beam
        GaussianBunch gaussian_bunch(n_electron, sigma_x, sigma_y, sigma_s);
        double gamma_e = params.beam->gamma();
        EBeam e_beam(gamma_e, tmp_tr, tmp_long, gaussian_bunch);

        //We need the sample to be big to eliminate statistical errors
        params.ecool_paras->set_n_sample(1e7);
            
        double rx_ibs, ry_ibs, rz_ibs = 0.0;
        if(do_ibs_){    
            //Calculate the delta rate from IBS
            ibs_solver->rate(lattice, *(params.beam), rx_ibs, ry_ibs, rz_ibs);
            //std::cout<<"IBS:" << rx_ibs<<" "<<ry_ibs<<" "<<rz_ibs<<std::endl;
        }
        
        double rate_x, rate_y, rate_s;
        ecooling_rate(*(params.ecool_paras), *(params.force_paras),
                      *(params.beam), cooler, e_beam, ring,
                      rate_x, rate_y, rate_s);

        std::cout<<rate_x<<" "<<rate_y<<" "<<rate_s<<std::endl;
        
        double total_rate_x = rate_x + rx_ibs;
        double total_rate_y = rate_y + ry_ibs;
        double total_rate_s = rate_s + rz_ibs;

        results[0].push_back(total_rate_x);
        results[1].push_back(total_rate_y);
        results[2].push_back(total_rate_s);
    }

    std::ofstream myfile;
    myfile.open ("scan.txt");
    for(int i=0; i < results[0].size(); i++){
        myfile<<" "<<scan_samples[i]<<", "<<std::setprecision(5)<<(results[0])[i]<<", "<<(results[1])[i]<<", "<<(results[2])[i]<<"\n";
    }
    myfile.close();
    
    return results;
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
    std::normal_distribution<double> log_c_distribution   (InitialValues["log_c"],0.5);

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
        log_c_distribution(generator);
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
    if(std::find(FitVariables.begin(),FitVariables.end(),"log_c") != FitVariables.end() ){
        double log_c_tmp = log_c_distribution(generator);
        while( log_c_tmp <= 0) log_c_tmp = log_c_distribution(generator);
        fitter_values.log_c_ = log_c_tmp;
    }
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

            /*
            if(iter % 20 == 0 ){
                std::cout<<iter;

                for(int i=0; i < n_pars; i++){
                    std::cout<<" "<<FitVariables[i]<<": "<<std::setprecision(5)<<gsl_vector_get (s->x,i);
                }
                std::cout<<" Score: "<<s->fval<<std::endl;
            }
            
                std::cout<<" Score: "<<s->fval<<" minutes from target cooling time"<<std::endl;
            */
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
        std::cout<<"Random Iteration "<<i+1<<"/"<<n_trials<<std::endl;
        Randomize();
        OptimizeTrial();
    }                 
    
    std::cout<<"Best values:"<<std::endl;
    for(int i=0; i < BestFits.size(); i++){
       std::cout<<" "<<FitVariables[i]<<": "<<std::setprecision(5)<<BestFits[FitVariables[i]];
    }
    std::cout<<" Score: "<<best_fval<<" minutes from target cooling time"<<std::endl;

    std::ofstream myfile;
    myfile.open ("bests.txt");
    for(int i=0; i < BestFits.size(); i++){
        myfile<<" "<<FitVariables[i]<<": "<<std::setprecision(5)<<BestFits[FitVariables[i]]<<"\n";
    }
    //myfile << best_b<<" "<< best_sigma_x<<" "<<best_sigma_y<<" "<<best_sigma_z<<" "<<
    //        best_beta_h<<" "<<best_beta_v<<" "<< best_disp_h<<" "<<best_disp_v<<" "<<
    //        best_temp_tr<<" "<<best_temp_long<<" "<< best_n_electron<<" "<< best_fval<<"\n";
    myfile.close();
}


int Optimize::Optimize_From_UI(std::vector<std::string> Params, std::vector<double>InitialValues, Beam &ion, Cooler &cooler, EBeam &ebeam, Ring &ring, ForceFormula ff){
    Lattice *lattice = ring.lattice_; //"eRHIC.tfs"

    //Params is a vector of string ID's for parameters
    //InitialValues is a vector of doubles, matched 1:1 with the params


    this->InitializeFitter(Params, InitialValues, lattice, &ion, &cooler,&ebeam, ff);
    this->ManyTrials();

    return 1;
}


int Optimize::ParameterScan_From_UI(std::vector<std::string> Params, std::vector<double>min_max, int n_steps, Beam &ion, Cooler &cooler, EBeam &ebeam, Ring &ring, ForceFormula ff){


    Lattice *lattice = ring.lattice_; //"eRHIC.tfs"

    //Params is a vector of string ID's for parameters
    //InitialValues is a vector of doubles, matched 1:1 with the params

    std::vector<double> dummy_d;
    std::vector<string> dummy_s;
    this->InitializeFitter(dummy_s, dummy_d, lattice, &ion, &cooler,&ebeam, ff);
    this->ParameterScan(Params[0], 
                        min_max[0], 
                        min_max[1],
                        n_steps,
                        this->fitter_values); //defined with InitializeFitter
    
    return 1;
}
