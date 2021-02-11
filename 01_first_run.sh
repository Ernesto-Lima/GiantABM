#!/bin/bash
rm *# *~ &> /dev/null
scale=16
#################### Create Options File ######################
input_file=options.in
cat << EOF > ${input_file}
initial_tum = 3
n_timesteps = 1000
print_inter = 5
verbose = 1
print_sa = 1
domain_diameter = 6384
cluster_scale = ${scale}
time_step = 1.0
c_ccr = 10.0
c_tta = 0.488836
c_hha = 0.588836
rand_seed = 5
prol_intens = 0.52
k_var = 1
k_val = 0.1
nucleus_radius = 5.295
cell_radius = 9.953
action_prop = 1.214
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
alpha_p = 5.9e-02
alpha_a = 4.1e-04
live_ic = 0.1
dead_ic = 0.0
t_death = 10.0
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
lend=$(wc -l < out_py.txt | awk '{printf "%d\n",$1/2}')
echo 
sed -n -e "1,${lend}p" out_py.txt > live_${number}.dat
let lbgn=lend+1
lend=$(wc -l < out_py.txt | awk '{printf "%d\n",$1+0}')
sed -n -e "${lbgn},${lend}p" out_py.txt > dead_${number}.dat
#rm *# *~ tmp.csv out_py.txt &> /dev/null
