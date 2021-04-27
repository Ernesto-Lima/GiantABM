#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <functional>
#include <fstream>
#include <string>
#include "main_abm_oxy.h"
#include "general_libraries.h"
#include "libapg.h"

struct likelihoodRoutine_DataType{
  std::vector< std::vector< std::vector<double> > > Full_Data;
  LibMeshInit *init;
};

using namespace std;

double likelihoodRoutine(const std::vector<double> paramValues,
                         const void* functionDataPtr);

inline void read_full_data(std::vector< std::vector< std::vector<double> > >& v_data,
                           const std::string data_file);
                           
int main(int argc, char* argv[]){
  MPI_Init(&argc,&argv);
  {
    Parameter_Space mesh(argc,argv);
    mesh.set_maxlvl(21);
    mesh.set_inter(11);
    mesh.add_parameter("proliferation rate",0,1);
    mesh.add_parameter("standard deviation",0.01,50);
    if(!mesh.world_rank){
      cout << "Number of intervals........" << mesh.get_inter() << endl;
      cout << "Maximum number of levels..." << mesh.get_maxlvl() << endl;
      cout << "Number of parameters......." << mesh.get_dim() << endl;
    }
    mesh.print_parameters();
    mesh.grid_init();
    mesh.grid_save("grid_model_calibration");
    //###################################################################
    //--------------- Initializing libmesh and likelihood ---------------
    LibMeshInit init (argc, argv, MPI_COMM_SELF);
    likelihoodRoutine_DataType likelihoodRoutine_Data;
    likelihoodRoutine_Data.init = &init;
    read_full_data(likelihoodRoutine_Data.Full_Data,"../Data/small_cells.dat");
    if(!mesh.world_rank){
      cout << "Number of datasets........." << likelihoodRoutine_Data.Full_Data.size() << endl;
      for (unsigned int i = 0; i<likelihoodRoutine_Data.Full_Data.size(); i++)
        cout << "Number of time points......" << likelihoodRoutine_Data.Full_Data[i][0].size() << endl;
    }
    //###################################################################
    mesh.solve(&likelihoodRoutine,&likelihoodRoutine_Data); // real 0m25.792s - MLE: P(0.195883|11.3736) = -138.985
    //mesh.solve_serial(&likelihoodRoutine,&likelihoodRoutine_Data); // real 1m29.590s - MLE: P(0.195883|11.3736) = -138.985
    if(!mesh.world_rank){
      mesh.grid_save("solved_model_calibration");
      vector<double> mle = mesh.get_mle();
      std::cout << "MLE: P(" << mle[0];
      for (unsigned int i = 1; i<mle.size()-1; i++){
        std::cout << "|" << mle[i];
      }
      std::cout << ") = " << mle[mle.size()-1];
    }
  }
  MPI_Finalize();
  return 0;
}

// Inverse Problem Routine
double likelihoodRoutine(const std::vector<double> paramValues,
                         const void* functionDataPtr){
  const vector< vector< vector<double> > >& Full_Data = ((likelihoodRoutine_DataType *) functionDataPtr)->Full_Data;
  unsigned int ndim = 2;
  std::vector<double> Parameters(ndim, 0);
  for (unsigned int i = 0; i<ndim; i++)
    Parameters[i] = paramValues[i];
  unsigned int time_dim = Full_Data[0][0].size();
  std::vector<double> Small_Cells(time_dim,0);
  std::vector<double> Giant_Cells(time_dim,0);
  Small_Cells[0] = Full_Data[0][1][0];
  main_code(*((likelihoodRoutine_DataType *) functionDataPtr)->init,
            Small_Cells,
            Giant_Cells,
            Parameters,
            27,
            Full_Data[0][0]);
  double misfitValue = 0.;
  for(unsigned int i=1; i < time_dim; i++){
	double StdDevs  = paramValues[ndim-1];
	double ratio    = (Full_Data[0][1][i] - Small_Cells[i])/StdDevs;
	misfitValue += ratio*ratio
	             + log(StdDevs)
	             + 0.5*log(M_PI)
	             + 0.5*log(2.0);
  }
  return -1.0*misfitValue;
}

inline void read_full_data(std::vector< std::vector< std::vector<double> > >& v_data,
                           const std::string data_file){
  std::ifstream read;
  std::string line;
  read.open(data_file);
  if(!read.is_open())
    std::cout << "Error opening data file.\n";
  std::vector< std::vector<double> > data_vec;
  std::vector<double> aux_vec_t;
  std::vector<double> aux_vec_c;
  while(std::getline(read, line)){
    double t,c;
    std::istringstream(line) >> t >> c;
	aux_vec_t.push_back(t);
	aux_vec_c.push_back(c);
  }
  data_vec.push_back(aux_vec_t);
  data_vec.push_back(aux_vec_c);
  v_data.push_back(data_vec);
  read.close();
}
