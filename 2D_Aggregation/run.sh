#!/bin/bash

gnuplot -e "set xrange [2e-2:2];
                set yrange [0.5:1]; 
                set ytics 0.5; 
                set title 'Comparison with analytical zeroth moment'; 
                set logscale y 10; 
                set logscale x 10; 
                set xlabel 'Time';
                set ylabel 'Zeroth moment';
                plot 'Ntol.tmp' title 'Numerical', 'Ntol_ana.tmp' title 'Analytical' with lines; 
                pause 10
            "
gnuplot -e "set xrange [1e-4:1e1];
                set yrange [1e-1:1e3]; 
                set title 'Comparison with analytical number density'; 
                set logscale y 10; 
                set logscale x 10; 
                set xlabel 'Pivots';
                set ylabel 'Number density';
                plot 'numden.tmp' title 'Numerical', 'anaden_ini.tmp' title 'Inital' with lines, 'anaden_final.tmp' title 'Analytical' with lines;
                pause 10
            "