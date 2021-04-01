#include "general_libraries.h"
#include "cell_abm.h"
#include "main_abm_oxy.h"

//###########################################################################
//                    0 - Dead Cells
//                    1 - Tumor Cells
//                    2 - Proliferative Tumor Cells
//                    3 - Hypoxic Tumor Cells
//                    4 - Dying Tumor Cells
//                    5 - G1 Tumor Cells
//                    6 - Normoxic Cells
//###########################################################################

double scalar_prod(const Point& e1,const Point& e2)
{
  return e1(0)*e2(0)+e1(1)*e2(1);
}

void main_code(LibMeshInit &init, vector<double>& LIVE_CONFL, vector<double>& DEAD_CONFL, vector<double> Parameters, const unsigned int rand_seed)
{

  list <Cell> Cells_local;
  
  //********** Read arguments **********
  //Read parameters from options file; if not in file, values are assigned to default values shown in next ~25ish lines
  GetPot input_file("options.in");
  const unsigned int n_timesteps = input_file("n_timesteps",100);
  const unsigned int initial_tum = input_file("initial_tum",20);
  const unsigned int print_inter = input_file("print_inter",100);
  const unsigned int verbose     = input_file("verbose",0);
  const unsigned int print_sa    = input_file("print_sa",0);
  const double cluster_scale     = input_file("cluster_scale",1.0);
  const double domain_diameter   = input_file("domain_diameter",1.0);
  const double time_step         = input_file("time_step",0.05);
  const double nucleus_radius    = input_file("nucleus_radius",5.295);
  const double cell_radius       = input_file("cell_radius",9.953);
  const double action_prop       = input_file("action_prop",1.214);
  const double c_ccr             = input_file("c_ccr",10.0);
  const double c_tta             = input_file("c_tta",0.488836);
  const double c_hha             = input_file("c_hha",0.588836);
  
  //Seed the RNG for stochastic processes
  Ran ran(rand_seed);
  int outside_cells = 0;
  int total_tumor = initial_tum;
  
  //Declare output file to store cell data in
  ofstream out_data;
  if(verbose || print_sa){
    char bufferf[30];
    sprintf(bufferf, "cells_data%d.txt",1);
    out_data.open(bufferf);
  }
  
  //********** Create a uniform mesh **********
  SerialMesh mesh(init.comm());
  MeshTools::Generation::build_cube (mesh,
                                     20,
                                     20,
                                     0,
                                     0., domain_diameter,
                                     0., domain_diameter,
                                     0., 0.,
                                     QUAD4);
  //********** Create an equation systems object and set a simulation-specific parameter **********
  EquationSystems eq_sys(mesh);
  eq_sys.parameters.set<Real> ("nucleus_radius")      = cluster_scale*nucleus_radius;
  eq_sys.parameters.set<Real> ("cell_radius")         = cluster_scale*cell_radius;
  eq_sys.parameters.set<Real> ("action_prop")         = action_prop;
  eq_sys.parameters.set<Real> ("action_radius")       = action_prop*cell_radius*cluster_scale;
  eq_sys.parameters.set<Real> ("domain_diameter")     = domain_diameter;
  eq_sys.parameters.set<Real> ("lambda_cell")         = 0.1;
  eq_sys.parameters.set<Real> ("initial_con_live")    = Parameters[2]; // Must be decimal
  eq_sys.parameters.set<Real> ("initial_con_dead")    = Parameters[3]; // Must be decimal
  eq_sys.parameters.set<Real> ("cellc_time")          = 18.0;
  eq_sys.parameters.set<Real> ("g1_time")             = 9.0;
  eq_sys.parameters.set<Real> ("prol_intes")          = Parameters[0];
  eq_sys.parameters.set<Real> ("apop_intes")          = Parameters[1];
  eq_sys.parameters.set<Real> ("apop_time")           = Parameters[4];
  eq_sys.parameters.set<Real> ("delta_tt")            = 1.0;
  eq_sys.parameters.set<Real> ("c_ccr")               = c_ccr;
  eq_sys.parameters.set<Real> ("c_tta")               = c_tta;
  eq_sys.parameters.set<Real> ("c_hha")               = c_hha;

  //********** Initial conditions for cells **********
  init_cond_cells(Cells_local, eq_sys.parameters, ran);
  compute_initial_forces(Cells_local, eq_sys, domain_diameter, outside_cells, total_tumor);
  total_tumor = Cells_local.size();
  
  //********** Initialize the data structures for the equation system **********
  double confluence = 0.;
  double dead_conf = 0.;
  int num_dead = 0;
  int qtd_c4,qtd_c1,qtd_cp,qtd_cm,qtd_gp;
  double mean_f,std_f;
  qtd_c4 = qtd_c1 = qtd_cp = qtd_cm = qtd_gp = 0;
  mean_f = std_f = 0.0;
  std::list<Cell>::iterator cg;
  for(cg = Cells_local.begin(); cg != Cells_local.end(); cg++){
    if((*cg).state==4){
	  dead_conf += std::pow((*cg).C_radius, 2) / std::pow(0.5 * domain_diameter, 2);
	  qtd_c4++;
    }
    else if((*cg).state==1){
	  confluence += std::pow((*cg).C_radius, 2) / std::pow(0.5 * domain_diameter, 2);
	  mean_f += (*cg).abs_nf;
	  qtd_c1++;
    }
    else if((*cg).state==8){
	  confluence += std::pow((*cg).C_radius, 2) / std::pow(0.5 * domain_diameter, 2);
	  qtd_cm++;
    }
    else if((*cg).state==2){
	  confluence += std::pow((*cg).C_radius, 2) / std::pow(0.5 * domain_diameter, 2);
	  mean_f += (*cg).abs_nf;
	  qtd_cp++;
    }
    else{
	  confluence += std::pow((*cg).C_radius, 2) / std::pow(0.5 * domain_diameter, 2);
	  qtd_gp++;
    }
  }
  mean_f = mean_f/(qtd_c1+qtd_cp);
  for(cg = Cells_local.begin(); cg != Cells_local.end(); cg++){
    if((*cg).state==1 || (*cg).state==2){
	  std_f += std::pow(mean_f-(*cg).abs_nf,2);
    }
  }
  std_f = sqrt(std_f/(qtd_c1+qtd_cp));
  num_dead = qtd_c4;
  if(confluence+dead_conf>1){
    double val = confluence+dead_conf;
    confluence=confluence/val;
    dead_conf = dead_conf/val;
  }
  LIVE_CONFL.push_back(confluence);
  DEAD_CONFL.push_back(dead_conf);
  //cout << dead_conf << endl;
  if(verbose)
    save_cells(Cells_local,domain_diameter,"saida",1,0);
  if(verbose || print_sa)
    out_data << 0.0 << " " << confluence << " " << dead_conf << " " << qtd_c1 << " " << qtd_cp << " " << qtd_gp << " " << qtd_c4 << " " << qtd_cm << " " << mean_f << " " << std_f << endl;
    
  //********** Loop over time **********
  unsigned int t_step = 0;
  
  do{
    //********** Increase time_step counter **********
    t_step++;

    //********** Solve ABM system **********
    update_states(Cells_local, t_step, ran, eq_sys, outside_cells, total_tumor, num_dead);
    compute_forces(Cells_local, eq_sys, domain_diameter, outside_cells, total_tumor, ran);

    // Confluence
    confluence = 0.;
    dead_conf = 0.;
    qtd_c4 = qtd_c1 = qtd_cp = qtd_cm = qtd_gp = 0;
    mean_f = std_f = 0.0;
    for(cg = Cells_local.begin(); cg != Cells_local.end(); cg++){
      if((*cg).state==4){
	    dead_conf += std::pow((*cg).C_radius, 2) / std::pow(0.5 * domain_diameter, 2);
	    qtd_c4++;
      }
      else if((*cg).state==1){
	    confluence += std::pow((*cg).C_radius, 2) / std::pow(0.5 * domain_diameter, 2);
	    mean_f += (*cg).abs_nf;
  	    qtd_c1++;
      }
      else if((*cg).state==8){
	    confluence += std::pow((*cg).C_radius, 2) / std::pow(0.5 * domain_diameter, 2);
	    qtd_cm++;
      }
      else if((*cg).state==2){
	    confluence += std::pow((*cg).C_radius, 2) / std::pow(0.5 * domain_diameter, 2);
	    mean_f += (*cg).abs_nf;
	    qtd_cp++;
      }
      else{
	    confluence += std::pow((*cg).C_radius, 2) / std::pow(0.5 * domain_diameter, 2);
	    qtd_gp++;
      }
    }
    mean_f = mean_f/(qtd_c1+qtd_cp);
    for(cg = Cells_local.begin(); cg != Cells_local.end(); cg++){
      if((*cg).state==1 || (*cg).state==2){
	    std_f += std::pow(mean_f-(*cg).abs_nf,2);
      }
    }
    std_f = sqrt(std_f/(qtd_c1+qtd_cp));
    if(confluence+dead_conf>1){
      double val = confluence+dead_conf;
      confluence=confluence/val;
      dead_conf = dead_conf/val;
    }
    // Output live confluence info if verbose
    if(verbose){
	  //if(t_step%print_inter==0 || t_step>283)
	  if(t_step%print_inter==0)
	    save_cells(Cells_local,domain_diameter,"saida",1,t_step);
	  //cout << "Time            = " << t_step*time_step << endl;
	    /*
	  cout << "**************************************************" << endl;
  	  cout << "T Confluence    = " << confluence + dead_conf << endl;
	  cout << "L Confluence    = " << confluence << endl;
	  cout << "D Confluence    = " << dead_conf << endl;
	  cout << "Number of cells = " << Cells_local.size() << endl;
	  cout << "Tumor cells     = " << total_tumor << endl;
	  cout << "Outside cells   = " << outside_cells << endl;
	  cout << "Time iteration  = " << t_step << endl;
	  cout << "Time            = " << t_step*time_step << endl;
	  cout << "**************************************************" << endl;
	  */
    }
    if(verbose || print_sa)
      out_data << t_step << " " << confluence << " " << dead_conf << " " << qtd_c1 << " " << qtd_cp << " " << qtd_gp << " " << qtd_c4 << " " << qtd_cm << " " << mean_f << " " << std_f << endl;
    if(t_step%3==0){
	  LIVE_CONFL.push_back(confluence);
	  DEAD_CONFL.push_back(dead_conf);
    }
  }while(t_step<n_timesteps);
  if(verbose)
    save_cells(Cells_local, domain_diameter, "saida", 1, t_step);
  if(verbose || print_sa)
    out_data.close();
}
