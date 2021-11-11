#!/bin/bash
rm *# *~ &> /dev/null
cd Forward_Model
    rm *# *~ saida1_*.txt job_output.txt cells_data*.txt sbl.txt live_*.txt dead_*.txt &> /dev/null
    #make distclean
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
cd ..
