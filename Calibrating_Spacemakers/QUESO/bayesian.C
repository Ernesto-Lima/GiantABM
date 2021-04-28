#include <queso/GenericVectorFunction.h>
#include <queso/GaussianVectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>
#include "general_libraries.h"
#include "main_abm_oxy.h"
#include "cell_abm.h"
#include "mpi.h"

using namespace std;

struct likelihoodRoutine_DataType{
  std::vector< std::vector< std::vector<double> > > Full_Data;
  LibMeshInit *init;
};

double likelihoodRoutine(const QUESO::GslVector& paramValues,
                         const QUESO::GslVector* paramDirection,
                         const void*             functionDataPtr,
                         QUESO::GslVector*       gradVector,
                         QUESO::GslMatrix*       hessianMatrix,
                         QUESO::GslVector*       hessianEffect);
			 
inline void read_full_data(std::vector< std::vector< std::vector<double> > >& v_data,
                           const std::string data_file);

int main(int argc, char* argv[]){
  MPI_Init(&argc,&argv);
  {
    LibMeshInit init (argc, argv, MPI_COMM_SELF);
    UQ_FATAL_TEST_MACRO(argc < 2,QUESO::UQ_UNAVAILABLE_RANK,"main()","input file must be argv[1]");
    QUESO::FullEnvironment env(MPI_COMM_WORLD,argv[1],"",NULL);
    if (env.fullRank() == 0) {
      std::cout << "\nBeginning run of 'Space Maker' code"
                << "\n my fullRank = "         << env.fullRank()
                << "\n my subEnvironmentId = " << env.subId()
                << "\n my subRank = "          << env.subRank()
                << "\n my interRank = "        << env.inter0Rank()
                 << std::endl << std::endl;
    }
    unsigned int N_par=3;
    //------------------------------------------------------
    // SIP Step 1 of 6: Instantiate the parameter space
    //------------------------------------------------------
    QUESO::VectorSpace<> paramSpace(env, "param_", N_par, NULL);
    //------------------------------------------------------
    // SIP Step 2 of 6: Instantiate the parameter domain
    //------------------------------------------------------
    QUESO::GslVector paramMins(paramSpace.zeroVector());
    QUESO::GslVector paramMaxs(paramSpace.zeroVector());
    paramMins[0] = 1.0;   // F_limit
    paramMaxs[0] = 10.0;
    paramMins[1] = 2.0;   // SM_timer
    paramMaxs[1] = 50.0;
    paramMins[2] = 0.0;   // Sigma_H
    paramMaxs[2] = 2.0;
    QUESO::BoxSubset<> paramDomain("param_",
                                   paramSpace,
                                   paramMins,
                                   paramMaxs);
    //------------------------------------------------------
    // SIP Step 3 of 6: Instantiate the likelihood function
    // object to be used by QUESO.
    //------------------------------------------------------
    likelihoodRoutine_DataType likelihoodRoutine_Data;
    likelihoodRoutine_Data.init = &init;
    read_full_data(likelihoodRoutine_Data.Full_Data,"../Data/full_small_cells.dat");
    read_full_data(likelihoodRoutine_Data.Full_Data,"../Data/full_giant_cells.dat");
    cout << "Number of datasets......" << likelihoodRoutine_Data.Full_Data.size() << endl;
    for (unsigned int i = 0; i<likelihoodRoutine_Data.Full_Data.size(); i++){
      cout << "Number of time points..." << likelihoodRoutine_Data.Full_Data[i].size() << endl;
      cout << "t = [" << likelihoodRoutine_Data.Full_Data[i][0][0];
      for (unsigned int j = 1; j<likelihoodRoutine_Data.Full_Data[i][0].size(); j++){
        cout << ", " << likelihoodRoutine_Data.Full_Data[i][0][j];
      }
      cout << "]\n";
    }
    QUESO::GenericScalarFunction<> likelihoodFunctionObj("like_",
                                                         paramDomain,
                                                         likelihoodRoutine,
                                                         (void *) &likelihoodRoutine_Data,
                                                         true);
    //------------------------------------------------------
    // SIP Step 4 of 6: Define the prior RV
    //------------------------------------------------------
    QUESO::UniformVectorRV<> priorRv("prior_", paramDomain);
    //------------------------------------------------------
    // SIP Step 5 of 6: Instantiate the inverse problem
    //------------------------------------------------------
    QUESO::GenericVectorRV<> postRv("post_", paramSpace);
    QUESO::StatisticalInverseProblem<> ip("",
                                          NULL,
                                          priorRv,
                                          likelihoodFunctionObj,
                                          postRv);
    //------------------------------------------------------
    // SIP Step 6 of 6: Solve the inverse problem, that is,
    // set the 'pdf' and the 'realizer' of the posterior RV
    //------------------------------------------------------
    
    //QUESO::GslVector paramInitials(paramSpace.zeroVector());
    //priorRv.realizer().realization(paramInitials);

    //QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());
    //proposalCovMatrix(0,0) = std::pow(std::abs(paramInitials[0]) / 20.0, 2.0);

    //ip.solveWithBayesMetropolisHastings(NULL, paramInitials);
    
    ip.solveWithBayesMLSampling();
  }
  MPI_Finalize();
  return 0;
}

// Inverse Problem Routine
double likelihoodRoutine(const QUESO::GslVector& paramValues,
                         const QUESO::GslVector* paramDirection,
                         const void*             functionDataPtr,
                         QUESO::GslVector*       gradVector,
                         QUESO::GslMatrix*       hessianMatrix,
                         QUESO::GslVector*       hessianEffect){
                         
  const vector< vector< vector<double> > >& Full_Data = ((likelihoodRoutine_DataType *) functionDataPtr)->Full_Data;
  const QUESO::BaseEnvironment& env = paramValues.env();
  MPI_Comm libMeshComm = env.subComm().Comm();
  unsigned int rank = env.subRank();
  int local_size;
  MPI_Comm_size(libMeshComm, &local_size);
  unsigned int ndim = 3;
  std::vector<double> Parameters(ndim, 0);
  for (unsigned int i = 0; i<ndim; i++)
    Parameters[i] = paramValues[i];
  unsigned int time_dim = Full_Data[0][0].size();
  std::vector<double> Small_Cells(time_dim,0);
  std::vector<double> Giant_Cells(time_dim,0);
  Small_Cells[0] = Full_Data[0][1][0];
  cout << paramValues[0] << " " << paramValues[1] << endl;
  main_code(*((likelihoodRoutine_DataType *) functionDataPtr)->init,
            Small_Cells,
            Giant_Cells,
            Parameters,
            27,
            Full_Data[0][0]);
  double misfitValue = 0.;
  for(unsigned int i=1; i < time_dim; i++){
	double StdDevs  = paramValues[ndim-1];
	//double ratio    = (Full_Data[0][1][i] - Small_Cells[i])/StdDevs;
	double ratio    = (Full_Data[1][1][i] - Giant_Cells[i])/StdDevs;
	//cout << Full_Data[0][1][i] << ", " << Small_Cells[i] << endl;
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
