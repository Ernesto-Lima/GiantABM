#!/bin/bash
rm *# *~ &> /dev/null
cd Forward_Model
rm *# *~ saida1_*.txt job_output.txt cells_data*.txt sbl.txt live_*.txt dead_*.txt &> /dev/null
#make distclean
#################### Create Python Mean Var Code ##############
python_file=mean_std.py
cat << EOF > ${python_file}
import sys
import pandas as pd
data = pd.read_csv('sbl.txt',header=None,delimiter=' ')
orig_stdout = sys.stdout
f = open('out_py.txt', 'w')
sys.stdout = f
for column in data:
    print (data[column].mean(),data[column].std())
sys.stdout = orig_stdout
f.close()
EOF
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#--------------- Run the model --------------------
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
FILE=sbl.txt
if [ -f "$FILE" ]; then
    cp $FILE tmp.txt
    make
    time make run
    cat tmp.txt >> sbl.txt
else
    make
    time make run
fi
#################### Post-Process the Results ###################
number=$(more $FILE | wc -l | awk '{printf "%05d\n",$1+0}')
echo $number
pvpython mean_std.py
lend=$(wc -l < out_py.txt | awk '{printf "%d\n",$1/2}')
echo 
sed -n -e "1,${lend}p" out_py.txt > live_${number}.txt
let lbgn=lend+1
lend=$(wc -l < out_py.txt | awk '{printf "%d\n",$1+0}')
sed -n -e "${lbgn},${lend}p" out_py.txt > dead_${number}.txt
rm *# *~ tmp.txt out_py.txt &> /dev/null
cd ..
