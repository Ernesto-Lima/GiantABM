#!/bin/bash
cd Forward_Model
rm *# *~ merged*.pdf &> /dev/null
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#--------------- Plotting details -----------------
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
terminal="set terminal pdf font 'Times-New-Roman,7' lw 1"
ext=pdf
font=8
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#--------------- Plotting every scenario ----------
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
for file in cells_data*.txt; do
    number=$(echo ${file} | cut -da -f3 | cut -d. -f1)
    name=panel_${number}
    echo "Data ID......${number}"
    echo "Panel name...${name}"
    echo ${terminal} > figura.cmd
    echo "set output '${name}.${ext}'" >> figura.cmd
    echo "set multiplot layout 3,3" >> figura.cmd
    #--------------- Total confluence -----------------
    echo "set lmargin 10" >> figura.cmd
    echo "set bmargin 3" >> figura.cmd
    echo "set ylabel 'Total confluence' font ',${font}'" >> figura.cmd
    echo "set xlabel 'Time (hours)' font ',${font}'" >> figura.cmd
    echo "set xtics font ',${font}'" >> figura.cmd
    echo "set ytics 0.2 font ',${font}'" >> figura.cmd
    echo "plot [][0:1.1] '${file}' u 1:(\$2+\$3) t'' w l lc 'red' lw 2," >> figura.cmd
    #--------------- Total cell count -----------------
    echo "set lmargin 10" >> figura.cmd
    echo "set bmargin 3" >> figura.cmd
    echo "set ylabel 'Cell count' font ',${font}'" >> figura.cmd
    echo "set xlabel 'Time (hours)' font ',${font}'" >> figura.cmd
    echo "set xtics font ',${font}'" >> figura.cmd
    echo "set ytics 200 font ',${font}'" >> figura.cmd
    echo "plot [][0:800] '${file}' u 1:(\$4+\$5+\$6+\$7+\$8) t'' w l lc 'red' lw 2," >> figura.cmd
    #--------------- Proliferative count --------------
    echo "set lmargin 10" >> figura.cmd
    echo "set bmargin 3" >> figura.cmd
    echo "set ylabel 'Proliferative count' font ',${font}'" >> figura.cmd
    echo "set xlabel 'Time (hours)' font ',${font}'" >> figura.cmd
    echo "set xtics font ',${font}'" >> figura.cmd
    echo "set ytics 20 font ',${font}'" >> figura.cmd
    echo "plot [][0:] '${file}' u 1:5 t'' w l lc 'red' lw 2," >> figura.cmd
    #--------------- G1 count -------------------------
    echo "set lmargin 10" >> figura.cmd
    echo "set bmargin 3" >> figura.cmd
    echo "set ylabel 'G1 count' font ',${font}'" >> figura.cmd
    echo "set xlabel 'Time (hours)' font ',${font}'" >> figura.cmd
    echo "set xtics font ',${font}'" >> figura.cmd
    echo "set ytics 20 font ',${font}'" >> figura.cmd
    echo "plot [][0:] '${file}' u 1:6 t'' w l lc 'red' lw 2," >> figura.cmd
    #--------------- Quiescent count ------------------
    echo "set lmargin 10" >> figura.cmd
    echo "set bmargin 3" >> figura.cmd
    echo "set ylabel 'Quiescent count' font ',${font}'" >> figura.cmd
    echo "set xlabel 'Time (hours)' font ',${font}'" >> figura.cmd
    echo "set xtics font ',${font}'" >> figura.cmd
    echo "set ytics 200 font ',${font}'" >> figura.cmd
    echo "plot [][0:800] '${file}' u 1:4 t'' w l lc 'red' lw 2," >> figura.cmd
    #--------------- Dead cell count ------------------
    echo "set lmargin 10" >> figura.cmd
    echo "set bmargin 3" >> figura.cmd
    echo "set ylabel 'Dead count' font ',${font}'" >> figura.cmd
    echo "set xlabel 'Time (hours)' font ',${font}'" >> figura.cmd
    echo "set xtics font ',${font}'" >> figura.cmd
    echo "set ytics 2 font ',${font}'" >> figura.cmd
    echo "plot [][0:] '${file}' u 1:7 t'' w p pt 7 lc 'red' ps 0.1," >> figura.cmd
    #--------------- Spacemaker count -----------------
    echo "set lmargin 10" >> figura.cmd
    echo "set bmargin 3" >> figura.cmd
    echo "set ylabel 'Spacemaker count' font ',${font}'" >> figura.cmd
    echo "set xlabel 'Time (hours)' font ',${font}'" >> figura.cmd
    echo "set xtics font ',${font}'" >> figura.cmd
    echo "set ytics 5 font ',${font}'" >> figura.cmd
    echo "plot [][0:] '${file}' u 1:8 t'' w p pt 7 lc 'red' ps 0.1," >> figura.cmd
    #--------------- Force ----------------------------
    echo "set lmargin 10" >> figura.cmd
    echo "set bmargin 3" >> figura.cmd
    echo "set style fill transparent solid 0.4 noborder" | tee -a figura.cmd &> /dev/null
    echo "set ylabel '|Force|' font ',${font}'" >> figura.cmd
    echo "set xlabel 'Time (hours)' font ',${font}'" >> figura.cmd
    echo "set xtics font ',${font}'" >> figura.cmd
    #echo "unset ytics" >> figura.cmd
    echo "set ytics 2 font ',${font}'" >> figura.cmd
    echo -n "plot [][0:] '${file}' u 1:9 t'' w l lc 'red' lw 2," >> figura.cmd
    echo "'${file}' u 1:(\$9+\$10):(\$9-\$10) t'' with filledcu fill lc 'red'" | tee -a figura.cmd &> /dev/null
    gnuplot "figura.cmd" &> /dev/null
    pdfcrop ${name}.pdf ${name}t.pdf &> /dev/null
    mv ${name}t.pdf ${name}.pdf
done
cd ..
