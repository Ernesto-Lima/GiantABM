// Running the main function as a function to test for vegf_prod values.
/*

Number of intervals........11
Maximum number of levels...12
Number of parameters.......3
---------------------------------------------
	 Parameters being calibrated
spacemaker force threshold = (1|10)
spacemaker timer = (2|50)
standard deviation = (0.01|3)
---------------------------------------------
Number of datasets.........2
Number of time points......319
Number of time points......319
MLE: P(5.54545|42.9091|3) = -873.248
MLE: P(5.54545|42.9091|3) = -873.248
MLE: P(5.77273|44.0455|2.94182) = -819.886
MLE: P(5.77273|43.4773|2.97591) = -808.506
MLE: P(5.77273|43.1932|2.99295) = -808.397
MLE: P(5.77273|43.0511|3) = -808.359
MLE: P(5.77273|43.0511|3) = -808.359
MLE: P(5.77273|43.0511|3) = -808.359
MLE: P(5.77273|43.0511|3) = -808.359
MLE: P(5.77273|43.0511|3) = -808.359
MLE: P(5.77273|43.0511|3) = -808.359
MLE: P(5.77273|43.0511|3) = -808.359
MLE: P(5.77273|43.0511|3) = -808.359 
***************************************************************
* Done Running App  forward_model
***************************************************************

real	361m44.041s
user	1446m13.291s
sys	0m7.378s

*/
#include "general_libraries.h"
#include "main_abm_oxy.h"
#include "cell_abm.h"
#include "mpi.h" 

inline void read_full_data(std::vector< std::vector< std::vector<double> > >& v_data,
                           const std::string data_file);

int main(int argc, char** argv){
  
  MPI_Init(&argc, &argv);
  
  {
    LibMeshInit init (argc, argv, MPI_COMM_SELF);
    
    std::vector<double> Parameters(3, 0);
    
    Parameters[0]  = 50.77273;
    Parameters[1]  = 43.0511;
    Parameters[2]  = 3;
    
    std::vector< std::vector< std::vector<double> > > Full_Data;
    
    read_full_data(Full_Data,"../../Data/full_small_cells.dat");
    read_full_data(Full_Data,"../../Data/full_giant_cells.dat");
    cout << "Number of datasets........." << Full_Data.size() << endl;
    for (unsigned int i = 0; i<Full_Data.size(); i++)
      cout << "Number of time points......" << Full_Data[i][0].size() << endl;
    unsigned int time_dim = Full_Data[0][0].size();
    std::vector<double> Small_Cells(time_dim,0);
    std::vector<double> Giant_Cells(time_dim,0);
    Small_Cells[0] = Full_Data[0][1][0];
    main_code(init,
              Small_Cells,
              Giant_Cells,
              Parameters,
              27,
              Full_Data[0][0]);
    ofstream out_file;
    out_file.precision(16);
    out_file.open("forward_solution.txt");
    for (unsigned int i=0; i<time_dim; i++)
      out_file << scientific << Full_Data[0][0][i] << " " 
                             << Small_Cells[i] << " " 
                             << Giant_Cells[i] << endl;
    out_file.close();
  }
  MPI_Finalize();
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
