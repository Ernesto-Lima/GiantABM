#!/bin/bash
rm *# *~ &> /dev/null
scale=16
#################### Create Options File ######################
input_file=options.in
cat << EOF > ${input_file}
mesh_name = circle_6384_lc400.msh
read_mesh = 1
ic_type = 27
initial_tum = 3
n_timesteps = 96
print_inter = 3
verbose = 1
print_sa = 1
domain_diameter = 6384
cluster_scale = ${scale}
con_b = 0.001
con_n = 0.001
time_step = 1.0
c_ccr = 10.0
c_tta = 0.488836
c_hha = 0.588836
max_outside = 2000
rand_seed = 5
prol_intens = 0.52
file_number = 1
k_var = 1
k_val = 0.1
nucleus_radius = 5.295
cell_radius = 9.953
action_prop = 1.214
ntri_ic = 0.5
EOF
#################### Create Python Mean Var Code ##############
python_file=mean_std.py
cat << EOF > ${python_file}
import sys
import pandas as pd
data = pd.read_csv('sbl.csv',header=None,delimiter=' ')
orig_stdout = sys.stdout
f = open('out_py.txt', 'w')
sys.stdout = f
for column in data:
    print (data[column].mean(),data[column].std())
sys.stdout = orig_stdout
f.close()
EOF
#################### Create Parameters File ###################
loop=1
echo "Realizations to add = ${loop}"
parameters_file=parameters.in
cat << EOF > ${parameters_file}
alpha_p = 4.9e-02
alpha_a = 4.1e-04
nu_diff = 50.0
n_lambd = 4.8e-02
live_ic = 0.5
dead_ic = 0.3
rate_pc = 0.0
t_death = 97.0
gamma_a = 2.4e-02
gamma_p = 50.0
sigma_h = 5.4e-02
size_loop = ${loop}
EOF
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#--------------- Run the model --------------------
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
FILE=sbl.csv
if [ -f "$FILE" ]; then
    cp $FILE tmp.csv
    make run
    cat tmp.csv >> sbl.csv
else 
    make run
fi
#################### Post-Process the Results ###################
number=$(more $FILE | wc -l | awk '{printf "%05d\n",$1+0}')
echo $number
pvpython mean_std.py
echo 
sed -n -e "1,33p" out_py.txt > live_${number}.dat
sed -n -e "34,66p" out_py.txt > dead_${number}.dat
rm *# *~ tmp.csv out_py.txt &> /dev/null
