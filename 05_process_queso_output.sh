#!/bin/bash
rm *# *~ &> /dev/null
if [ -d "./All_dataOutput/" ]; then
    if compgen -G "./Model_Calibration/outputData/" > /dev/null; then
        cd ./Model_Calibration/outputData/
            lines=$(grep zeros rawChain_mh.m | head -n1 | cut -d"(" -f2 | cut -d, -f1 | awk '{printf "%d\n",$1+1}')
            echo "Number of lines...${lines}"
            sed -n -e "2,${lines}p" rawChain_mh.m | cut -d"[" -f2 > ../../All_dataOutput/parameters_chain.dat
            sed -n -e "2,${lines}p" rawChain_mh_loglikelihood.m | cut -d"[" -f2 > ../../All_dataOutput/loglikelihood_chain.dat
        cd ..
    else
        echo "Error: Directory /path/to/dir does not exists."
    fi
else
    mkdir All_dataOutput
fi
