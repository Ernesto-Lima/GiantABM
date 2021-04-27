#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <functional>
#include <fstream>
#include <string>
#include "libapg.h"
#include <mpi.h>

std::vector< std::vector<double> > Parameter_Space::get_unique(std::vector< std::vector<double> > full_vec){
  std::vector< std::vector<double> > level_0;
  level_0.push_back(full_vec[0]);
  for(unsigned int i = 0; i<full_vec.size(); i++){
    bool add_it = true;
    for(unsigned int j = 0; j<level_0.size(); j++){
      bool equal_par = true;
      for(unsigned int p = 0; p<parameters.size(); p++){
        if(std::fabs(full_vec[i][p] - level_0[j][p])>1e-10)
          equal_par = false;
      }
      if(equal_par)
        add_it = false;
    }
    if(add_it)
      level_0.push_back(full_vec[i]);
  }
  return level_0;
}

std::vector< std::vector<double> > Parameter_Space::get_moore(std::vector<double> mle){
  int depth = parameters.size();
  const unsigned int new_NoI = 2;
  const int max_p = new_NoI+1;
  unsigned int grid_size = std::pow(max_p,depth)-1;
  std::vector<double> sample(depth,0);
  std::vector< std::vector<double> > level_0;
  for(unsigned int i = 0; i<grid_size; i++)
    level_0.push_back(sample);
  int* slots = (int*)alloca(sizeof(int) * depth);
  unsigned int pos = 0;
  for (int i = 0; i < depth; i++)
    slots[i] = 0;
  int index = 0;
  while (true){
    for (int i = 0; i < depth; i++){
      double range = (1.0/std::pow(2,GRID.size()))*(parameters[i].max/number_of_intervals);
      double p_low = mle[i]-range;
      if(p_low<parameters[i].low)
        p_low=parameters[i].low;
      double p_max = mle[i]+range;
      if(p_max>parameters[i].max)
        p_max=parameters[i].max;
      double value = p_low+slots[i]*range;
      if(value>p_max)
        level_0[pos][i] = p_max;
      else if (value<p_low)
        level_0[pos][i] = p_low;
      else
        level_0[pos][i] = value;
    }
    // Increment
    slots[0]++;
    bool add_value = false;
    for(unsigned int par = 0; par<depth; par++)
      if(level_0[pos][par] != mle[par])
        add_value = true;
    if(add_value)
      pos++;
    // Carry
    while (slots[index] == max_p){
      // Overflow, we're done
      if (index == depth - 1){
        return level_0;
      }
      slots[index++] = 0;
      slots[index]++;
    }
    index = 0;
  }
}

void Parameter_Space::grid_refine(){
  std::vector<double> mle;
  std::vector<double> mle2;
  bool first = true;
  bool second = true;
  const unsigned int g = GRID.size()-1;
  for(unsigned int l=0; l<GRID[g].size(); l++){
    if(parameters.size()<GRID[g][l].size()){
      if(first || second){
        if(first){
          mle = GRID[g][l];
          first = false;
        }
        else{
          if(mle[GRID[g][l].size()-1]<GRID[g][l][GRID[g][l].size()-1]){
            mle2 = mle;
            mle = GRID[g][l];
            second = false;
          }
          else{
            mle2 = GRID[g][l];
            second = false;
          }
        }
      }
      else{
        if(mle[GRID[g][l].size()-1]<GRID[g][l][GRID[g][l].size()-1]){
          mle2 = mle;
          mle = GRID[g][l];
        }
        else if(mle2[GRID[g][l].size()-1]<GRID[g][l][GRID[g][l].size()-1]){
          mle2 = GRID[g][l];
        }
      }
    }
  }
  std::vector< std::vector<double> > level_0 = get_moore(mle);
  level_0 = get_unique(level_0);
  std::vector< std::vector<double> > level_1 = get_moore(mle2);
  level_1 = get_unique(level_1);
  for(unsigned int i = 0; i<level_1.size(); i++){
    bool add_it = true;
    for(unsigned int j = 0; j<level_0.size(); j++){
      bool equal_par = true;
      for(unsigned int p = 0; p<parameters.size(); p++){
        if(std::fabs(level_1[i][p] - level_0[j][p])>1e-10)
          equal_par = false;
      }
      if(equal_par)
        add_it = false;
    }
    if(add_it)
      level_0.push_back(level_1[i]);
  }
  GRID.push_back(level_0);
}

void Parameter_Space::grid_refine_mle(){
  std::vector<double> mle;
  bool first = true;
  const unsigned int g = GRID.size()-1;
  for(unsigned int l=0; l<GRID[g].size(); l++){
    if(parameters.size()<GRID[g][l].size()){
      if(first){
        mle = GRID[g][l];
        first = false;
      }
      else{
        if(mle[GRID[g][l].size()-1]<GRID[g][l][GRID[g][l].size()-1])
          mle = GRID[g][l];
      }
    }
  }
  std::vector< std::vector<double> > level_0 = get_moore(mle);
  GRID.push_back(level_0);
}

std::vector<double> Parameter_Space::get_mle(){
  std::vector<double> mle;
  if(GRID.size()!=0){
    bool first = true;
    for(unsigned int g=0; g<GRID.size(); g++){
      for(unsigned int l=0; l<GRID[g].size(); l++){
        if(parameters.size()<GRID[g][l].size()){
          if(first){
            mle = GRID[g][l];
            first = false;
          }
          else{
            if(mle[GRID[g][l].size()-1]<GRID[g][l][GRID[g][l].size()-1])
              mle = GRID[g][l];
          }        
        }
      }
    }
    if(first){
      std::vector<double> tmp(1,0);
      std::cout << "ERROR: You must solve the model before.\n";
      return tmp;
    }
  }
  else{
    std::cout << "ERROR: You must initialize the grid before solving the problem.\n";
  }
  return mle;
}

void Parameter_Space::solve_mle(std::function<double(const std::vector<double> paramValues, const void *functionDataPtr)> func, const void* functionDataPtr){
  if(GRID.size()!=0){
    for(unsigned int l = 0; l < level_max-1; l++){
      const unsigned int g = GRID.size()-1;
      bool refine = false;
      for(unsigned int l=0; l<GRID[g].size(); l++){
        if(parameters.size()==GRID[g][l].size()){
          double value;
          value = func(GRID[g][l],functionDataPtr);
          GRID[g][l].push_back(value);
          refine = true;
        }
      }
      if(refine)
        grid_refine_mle();
    }
    const unsigned int g = GRID.size()-1;
    for(unsigned int l=0; l<GRID[g].size(); l++){
      if(parameters.size()==GRID[g][l].size()){
        double value;
        value = func(GRID[g][l],functionDataPtr);
        GRID[g][l].push_back(value);
      }
    }
  }
  else{
    std::cout << "ERROR: You must initialize the grid before solving the problem.\n";
  }
}

void Parameter_Space::solve(std::function<double(const std::vector<double> paramValues, const void *functionDataPtr)> func, const void* functionDataPtr){
  if(GRID.size()!=0){
    for(unsigned int l = 0; l < level_max-1; l++){
      const unsigned int g = GRID.size()-1;
      bool refine = false;
      int remainder = GRID[g].size() % world_size;
      int local_counts[world_size], offsets[world_size];
      int sum = 0;
      for (int i = 0; i < world_size; i++) {
        local_counts[i] = GRID[g].size() / world_size;
        if (remainder > 0) {
          local_counts[i] += 1;
          remainder--;
        }
        offsets[i] = sum;
        sum += local_counts[i];
      }
      int localArray[local_counts[world_rank]];
      int A[GRID[g].size()];
      double R[GRID[g].size()];
      double localLogLike[local_counts[world_rank]];
      if (!world_rank) {
        for (int i = 0; i < GRID[g].size(); i++) {
          A[i] = i;
        }
      }
      MPI_Scatterv(A, local_counts, offsets, MPI_INT, localArray, local_counts[world_rank], MPI_INT, 0, MPI_COMM_WORLD);
      for (unsigned int r = 0; r<local_counts[world_rank]; r++){
        unsigned int l = localArray[r];
        if(parameters.size()==GRID[g][l].size()){
          double value;
          value = func(GRID[g][l],functionDataPtr);
          refine = true;
          localLogLike[r] = value;
        }
        else{
          localLogLike[r] = GRID[g][l][GRID[g][l].size()-1];
        }
      }
      MPI_Allgatherv(localLogLike, local_counts[world_rank], MPI_DOUBLE, R, local_counts, offsets, MPI_DOUBLE, MPI_COMM_WORLD);
      for(unsigned int l=0; l<GRID[g].size(); l++)
        if(parameters.size()==GRID[g][l].size())
          GRID[g][l].push_back(R[l]);
      if(refine)
        grid_refine();
    }
    const unsigned int g = GRID.size()-1;
    int remainder = GRID[g].size() % world_size;
    int local_counts[world_size], offsets[world_size];
    int sum = 0;
    for (int i = 0; i < world_size; i++) {
      local_counts[i] = GRID[g].size() / world_size;
      if (remainder > 0) {
        local_counts[i] += 1;
        remainder--;
      }
      offsets[i] = sum;
      sum += local_counts[i];
    }
    int localArray[local_counts[world_rank]];
    int A[GRID[g].size()];
    double R[GRID[g].size()];
    double localLogLike[local_counts[world_rank]];
    if (!world_rank) {
      for (int i = 0; i < GRID[g].size(); i++) {
        A[i] = i;
      }
    }
    MPI_Scatterv(A, local_counts, offsets, MPI_INT, localArray, local_counts[world_rank], MPI_INT, 0, MPI_COMM_WORLD);
    for (unsigned int r = 0; r<local_counts[world_rank]; r++){
      unsigned int l = localArray[r];
      if(parameters.size()==GRID[g][l].size()){
        double value;
        value = func(GRID[g][l],functionDataPtr);
        localLogLike[r] = value;
      }
      else{
        localLogLike[r] = GRID[g][l][GRID[g][l].size()-1];
      }
    }
    MPI_Allgatherv(localLogLike, local_counts[world_rank], MPI_DOUBLE, R, local_counts, offsets, MPI_DOUBLE, MPI_COMM_WORLD);
    for(unsigned int l=0; l<GRID[g].size(); l++)
      if(parameters.size()==GRID[g][l].size())
        GRID[g][l].push_back(R[l]);
  }
  else{
    std::cout << "ERROR: You must initialize the grid before solving the problem.\n";
  }
}

void Parameter_Space::solve_serial(std::function<double(const std::vector<double> paramValues, const void *functionDataPtr)> func, const void* functionDataPtr){
  if(GRID.size()!=0){
    for(unsigned int l = 0; l < level_max-1; l++){
      const unsigned int g = GRID.size()-1;
      bool refine = false;
      for(unsigned int l=0; l<GRID[g].size(); l++){
        if(parameters.size()==GRID[g][l].size()){
          double value;
          value = func(GRID[g][l],functionDataPtr);
          GRID[g][l].push_back(value);
          refine = true;
        }
      }
      if(refine)
        grid_refine();
    }
    const unsigned int g = GRID.size()-1;
    for(unsigned int l=0; l<GRID[g].size(); l++){
      if(parameters.size()==GRID[g][l].size()){
        double value;
        value = func(GRID[g][l],functionDataPtr);
        GRID[g][l].push_back(value);
      }
    }
  }
  else{
    std::cout << "ERROR: You must initialize the grid before solving the problem.\n";
  }
}

void Parameter_Space::grid_init(){
  if(parameters.size()==0)
    std::cout << "WARNING: You must add a parameter before initializing the grid.\n";
  else{
    if(GRID.size()!=0)
      std::cout << "WARNING: Grid already initialized.\n";
    else{
      int depth = parameters.size();
      unsigned int grid_size = std::pow(number_of_intervals+1,depth);
      std::vector<double> sample(depth,0);
      std::vector< std::vector<double> > level_0;
      for(unsigned int i = 0; i<grid_size; i++)
        level_0.push_back(sample);
      int max = number_of_intervals+1;
      int* slots = (int*)alloca(sizeof(int) * depth);
      unsigned int pos = 0;
      for (int i = 0; i < depth; i++)
        slots[i] = 0;
      int index = 0;
      while (true){
        for (int i = 0; i < depth; i++){
          double value = parameters[i].low+slots[i]*parameters[i].max/number_of_intervals;
          if(value>parameters[i].max)
            level_0[pos][i] = parameters[i].max;
          else if (value<parameters[i].low)
            level_0[pos][i] = parameters[i].low;
          else
            level_0[pos][i] = value;
        }
        pos++;
        // Increment
        slots[0]++;
        // Carry
        while (slots[index] == max){
          // Overflow, we're done
          if (index == depth - 1){
            GRID.push_back(level_0);
            return;
          }
          slots[index++] = 0;
         slots[index]++;
        }
        index = 0;
      }
    }
  }
}

void Parameter_Space::grid_print(){
  if(GRID.size()!=0){
    for(unsigned int g=0; g<GRID.size(); g++)
      for(unsigned int l=0; l<GRID[g].size(); l++){
        std::cout << "[" << GRID[g][l][0];
        for(unsigned int p=1; p<GRID[g][l].size(); p++){
          std::cout << "|" << GRID[g][l][p];
        }
        std::cout << "]\n";
      }
  }
  else{
    int depth = parameters.size();
    if(depth==0){
      std::cout << "WARNING: You must add a parameter before printing the grid.\n";
    }
    else{
      std::cout << "WARNING: Initializing the grid inside the grid_print().\n";
      grid_init();
      grid_print();
    }
  }
}

void Parameter_Space::grid_save(std::string file_name){
  if(GRID.size()!=0){
    std::ofstream output_grid;
    std::string ext = ".txt";
    std::string full_name = file_name + ext;
    output_grid.open(full_name);
    for(unsigned int g=0; g<GRID.size(); g++){
      std::stringstream ss;
      ss << std::setw(10) << std::setfill('0') << g;
      std::string s = ss.str();
      std::string lvl = "_lvl";
      std::string level_name = file_name + lvl + s + ext;
      std::ofstream output_lvl;
      output_lvl.open(level_name);
      for(unsigned int l=0; l<GRID[g].size(); l++){
        output_grid << GRID[g][l][0];
        output_lvl << GRID[g][l][0];
        for(unsigned int p=1; p<GRID[g][l].size(); p++){
          output_grid << " " << GRID[g][l][p];
          output_lvl << " " << GRID[g][l][p];
        }
        output_grid << std::endl;
        output_lvl << std::endl;
      }
      output_lvl.close();
    }
    output_grid.close();
  }
  else{
    int depth = parameters.size();
    if(depth==0){
      std::cout << "WARNING: You must add a parameter before printing the grid.\n";
    }
    else{
      std::cout << "WARNING: Initializing the grid inside the grid_save().\n";
      grid_init();
      grid_save(file_name);
    }
  }
}

void Parameter_Space::print_parameters(){
  if(!world_rank){
    std::cout << "---------------------------------------------\n";
    std::cout << "\t Parameters being calibrated\n";
    for (unsigned i = 0; i<parameters.size(); i++)
      std::cout << parameters[i].name << " = (" << parameters[i].low << "|" << parameters[i].max << ")\n";
    std::cout << "---------------------------------------------\n";
  }
}

void Parameter_Space::add_parameter(std::string name,double min,double max){
  bool found = false;
  for (unsigned i = 0; i<parameters.size(); i++)
    if (parameters[i].name == name)
      found = true;
  if(found){
    std::cout << "ERROR: Parameter with the same name found.\n";
  }
  else{
    parameter new_parameter;
    new_parameter.name = name;
    if(min<max){    
      new_parameter.low = min;
      new_parameter.max = max;
    }
    else{
      std::cout << "WARNING: Lower bound was higher than upper bound.\n";
      new_parameter.low = max;
      new_parameter.max = min;
    }
    parameters.push_back(new_parameter);
  }
}

void Parameter_Space::add_parameter(std::string name){
  bool found = false;
  for (unsigned i = 0; i<parameters.size(); i++)
    if (parameters[i].name == name)
      found = true;
  if(found){
    std::cout << "ERROR: Parameter with the same name found\n";
  }
  else{
    parameter new_parameter;
    new_parameter.name = name;
    new_parameter.low = 0.0;
    new_parameter.max = 1.0;
    parameters.push_back(new_parameter);
  }
}

void Parameter_Space::set_maxlvl(int lvl){
  if(lvl>0)
    level_max = lvl;
  else{
    level_max = 1;
    std::cout << "WARNING: Max level should be higher than zero.\n";
  }
}

void Parameter_Space::set_inter(int inter){
  if(inter>0)
  number_of_intervals = inter;
  else{
    number_of_intervals = 1;
    std::cout << "WARNING: Number of intervals should be higher than zero.\n";
  }
}

unsigned int Parameter_Space::get_inter(){
  return number_of_intervals;
}

Parameter_Space::Parameter_Space(){
  level_max = 3;
  number_of_intervals = 4;
}

Parameter_Space::Parameter_Space(int argc, char* argv[]){
  level_max = 3;
  number_of_intervals = 4;
  int flag_i;
  MPI_Initialized(&flag_i);
  if ( ! flag_i)
    MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if(!world_rank)
    std::cout << std::endl;
}

unsigned int Parameter_Space::get_dim(){
  return parameters.size();
}

unsigned int Parameter_Space::get_maxlvl(){
  return level_max;
}

