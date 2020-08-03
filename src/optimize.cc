#include "optimize.h"

extern DynamicParas * dynamic_paras;
extern IBSSolver * ibs_solver;
extern EcoolRateParas * ecool_paras;
extern ForceParas * force_paras;
extern Luminosity *luminosity_paras;

double Optimize::fit_fcn(const gsl_vector *v, void *params){

    //Load in the struct holding our run parameters
    struct opt_info *p = (struct opt_info *)params;
    
    double magnetic_field = gsl_vector_get(v, 0);
    double sigma_x = gsl_vector_get(v, 1);
    double sigma_y = gsl_vector_get(v, 2);
    double sigma_s = gsl_vector_get(v, 3);
    double beta_h  = gsl_vector_get(v, 4);
    double beta_v  = gsl_vector_get(v, 5);
    double dis_h  = gsl_vector_get(v, 6);
    double dis_v  = gsl_vector_get(v, 7);
    double tmp_tr = gsl_vector_get(v, 8);
    double tmp_long = gsl_vector_get(v, 9);
    double n_electron = gsl_vector_get(v, 10);
    
    //Punish unphysical parameters
    if (magnetic_field < 0.0 || 
        sigma_x < 0.0 || sigma_y < 0.0 || sigma_x < 0.0 ||
        tmp_tr < 1e-5 || tmp_long < 1e-5 ||
        n_electron > 2.0){
        return 100000.0;
    }
    
    //We're really only varying the scale factor
    n_electron = n_electron * 1e10; 
    
    // define proton beam;
    double KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
    //m0 = 938.272;
    KE = 2.5e4;
    emit_nx0 = 2.5e-6;
    emit_ny0 = 2.5e-6;
    dp_p0= 1e-4;
    sigma_s0 = 0.7;
    N_ptcl = 1.34e12;
    Beam p_beam(p->Z,p->m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl);

    // define the lattice 
    Lattice lattice(p->lattice_filename);

    //Define the ring
    Ring ring(lattice, p_beam);

    //Set IBS parameters.
    int nu = 100;
    int nv = 100;
    int nz = 40;
    double log_c = 20.6;
    ibs_solver = new IBSSolver_Martini(nu, nv, nz, log_c, 0);

    //define the cooler
    double cooler_length = 130.0;
    double n_section = 1;
    Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);

    //define electron beam
    GaussianBunch gaussian_bunch(n_electron, sigma_x, sigma_y, sigma_s);
    double gamma_e = p_beam.gamma();
    EBeam e_beam(gamma_e, tmp_tr, tmp_long, gaussian_bunch);

    unsigned int n_sample = 5000000;
    ecool_paras = new EcoolRateParas(n_sample);
    //define friction force formula
    force_paras = ChooseForce(ForceFormula::PARKHOMCHUK);

    double rate_x, rate_y, rate_s;
    ecooling_rate(*ecool_paras, *force_paras, p_beam, cooler, e_beam, ring, rate_x, rate_y, rate_s);
        
    if(isnan(rate_y)){
        return 100000.0;
    }
    else{
        return abs(20. - ((-1./(rate_y))/60.));
    }
}


std::map<std::string, double> Optimize::Randomize()
{           
    //Select a new starting point in the parameter space. For efficiency we
    // throw all variables, even if one is frozen in the implementation
        
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);

    //TODO: Stash these RMS values in a Map too
    std::normal_distribution<double> bfield_distribution (magnetic_field_,3.0);
    std::normal_distribution<double> sigma_distribution (sigma_x_,1e-4);
    std::normal_distribution<double> sigma_s_distribution (sigma_s_,0.2);
    std::normal_distribution<double> beta_distribution (beta_h_,10.);
    std::normal_distribution<double> disp_distribution (disp_h_,.3);
    std::normal_distribution<double> temp_distribution (temp_tr_,0.01);
    std::normal_distribution<double> e_distribution (n_electron_,0.5);

    std::map<std::string, double> output;
    
    //Select random values, but add some constraints
    double bfield_tmp = bfield_distribution(generator);
    while(bfield_tmp <= 0.) bfield_tmp = bfield_distribution(generator);
    output["bfield"] = bfield_tmp;
    
    double sigma_x_tmp = sigma_distribution(generator);
    while(sigma_x_tmp <= 1e-4) sigma_x_tmp = sigma_distribution(generator);
    output["sigma_x"] = sigma_x_tmp;
    
    double sigma_y_tmp = sigma_distribution(generator);
    while(sigma_y_tmp <= 1e-4) sigma_y_tmp = sigma_distribution(generator);
    output["sigma_y"] = sigma_y_tmp;
    
    double sigma_s_tmp = sigma_s_distribution(generator);
    while(sigma_s_tmp <= 0.) sigma_s_tmp = sigma_s_distribution(generator);
    output["sigma_s"] = sigma_s_tmp;
    
    double beta_h_tmp = beta_distribution(generator);
    while(beta_h_tmp <= 0.) beta_h_tmp = beta_distribution(generator);
    output["beta_h"] = beta_h_tmp;
    
    double beta_v_tmp = beta_distribution(generator);
    while(beta_v_tmp <= 0.) beta_v_tmp = beta_distribution(generator);
    output["beta_v_tmp"] = beta_v_tmp;
    
    double disp_h_tmp = disp_distribution(generator);
    while(disp_h_tmp <= 1e-5) disp_h_tmp = disp_distribution(generator);
    output["disp_h"] = disp_h_tmp;
    
    double disp_v_tmp = disp_distribution(generator);
    while(disp_v_tmp <= 1e-5) disp_v_tmp = disp_distribution(generator);
    output["disp_v"] = disp_h_tmp;
    
    double temp_tr_tmp = temp_distribution(generator);
    while(temp_tr_tmp <= 1e-4) temp_tr_tmp = temp_distribution(generator);
    output["temp_tr"] = temp_tr_tmp;
    
    double temp_long_tmp = temp_distribution(generator);
    while(temp_long_tmp <= 1e-4) temp_long_tmp = temp_distribution(generator);
    output["temp_long"] = temp_long_tmp;
    
    double n_electron_tmp = e_distribution(generator);
    while( n_electron_tmp <= 1e-4) n_electron_tmp = e_distribution(generator);
    output["n_electron"] = n_electron_tmp;                
    
    //Hard-coding this here, temporarily
    FitStepSize.clear();
    FitStepSize["bfield"] = 0.1;
    FitStepSize["sigma_x"] = 5e-5;
    FitStepSize["sigma_y"] = 5e-5;
    FitStepSize["sigma_s"] = 5e-3;
    FitStepSize["beta_h"] = 0.5;
    FitStepSize["beta_v"] = 0.5;
    FitStepSize["disp_h"] = 0.1;
    FitStepSize["disp_v"] = 0.1;
    FitStepSize["temp_tr"] = 0.01;
    FitStepSize["temp_long"] = 0.01;

    #ifndef NDEBUG
    std::cout<<"Starting values:"<<std::endl;
    printf ("Bfield = %.5f sigma_x = %.5f sigma_y = %.5f sigma_s = %.5f beta_h = %.5f beta_v = %.5f \n",
            bfield_tmp, sigma_x_tmp,sigma_y_tmp,sigma_z_tmp,beta_h_tmp,beta_v_tmp);       
    printf( "disp_h = %.5f disp_v = %.5f temp_tr = %.5f temp_long = %.5f n_electrons = %.5 e10 \n", 
           disp_h_tmp,disp_v_tmp,temp_tr_tmp,temp_long_tmp,n_electron_tmp);
    #endif
    return output;
}

void Optimize::OptimizeTrial(){
                
    std::map<std::string, double> InitialPositions = Randomize();
    
    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;

    gsl_vector *x;
    gsl_vector *step_size;

    gsl_multimin_function my_func;

    int n_pars = FitVariables.size();
    my_func.n = n_pars;
    
    my_func.f = &Optimize::fit_fcn;
                    
    try{
        //Set the best guess
        x = gsl_vector_alloc (n_pars);
        step_size = gsl_vector_alloc(n_pars);

        for(int i = 0; i < n_pars; i++){
            gsl_vector_set (x, i, InitialPositions[FitVariables[i]] );
            gsl_vector_set (x, i, FitStepSize[FitVariables[i]] );
        }
        
        T = gsl_multimin_fminimizer_nmsimplex2rand;
        s = gsl_multimin_fminimizer_alloc (T, n_pars);
               
        gsl_multimin_fminimizer_set (s, &my_func, x, step_size);
       
        int iter = 0;
        int status;
                
        do{
            iter++;
            status = gsl_multimin_fminimizer_iterate (s);

            if (status)
                break;    
            if(iter % 10 == 0 ){
                printf ("%5d ",iter);
                for(int i=0; i < n_pars; i++) printf("%.5f ",gsl_vector_get (s->x,i));
                printf("%10.5f\n", s->fval);
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
}
    
    
void Optimize::ManyTrials()
{                   
    for(int i = 0;i<150;i++){
        std::cout<<"Random Iteration "<<i<<std::endl;

        OptimizeTrial();
    }                 
    
    std::cout<<"Best values:"<<std::endl;
    printf ("%5d ",BestFits.size());
    for(int i=0; i < BestFits.size(); i++) printf("%.5f ",BestFits[FitVariables[i]]);
    printf("%10.5f\n", best_fval);
    
            
    //std::ofstream myfile;
    //myfile.open ("bests.txt");
    //myfile << best_b<<" "<< best_sigma_x<<" "<<best_sigma_y<<" "<<best_sigma_z<<" "<<
    //        best_beta_h<<" "<<best_beta_v<<" "<< best_disp_h<<" "<<best_disp_v<<" "<<
    //        best_temp_tr<<" "<<best_temp_long<<" "<< best_n_electron<<" "<< best_fval<<"\n";
    //myfile.close();    
}