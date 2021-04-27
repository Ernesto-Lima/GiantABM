#ifndef LIB_APG_H_
#define LIB_APG_H_

struct parameter {
    std::string name;
    double low;
    double max;
};

class Parameter_Space{
private:
  unsigned int level_max;
  unsigned int number_of_intervals;
  std::vector<parameter> parameters;
  std::vector< std::vector< std::vector<double> > > GRID;
  std::vector< std::vector<double> > get_moore(std::vector<double> mle);
  void grid_refine_mle();
  void grid_refine();
  std::vector< std::vector<double> > get_unique(std::vector< std::vector<double> > full_vec);
  Parameter_Space();
public:
  int world_rank;
  int world_size;
  Parameter_Space(int argc, char* argv[]);
  unsigned int get_dim();
  unsigned int get_maxlvl();
  void set_maxlvl(int lvl);
  void add_parameter(std::string name);
  void add_parameter(std::string name,double min,double max);
  void print_parameters();
  void set_inter(int inter);
  unsigned int get_inter();
  void grid_print();
  void grid_init();
  void grid_save(std::string file_name);
  void solve(std::function<double(const std::vector<double> paramValues, const void *functionDataPtr)> func, const void* functionDataPtr);
  void solve_serial(std::function<double(const std::vector<double> paramValues, const void *functionDataPtr)> func, const void* functionDataPtr);
  void solve_mle(std::function<double(const std::vector<double> paramValues, const void *functionDataPtr)> func, const void* functionDataPtr);
  std::vector<double> get_mle();
};

#endif
