#!/bin/bash

# gnuplot -e "set xrange [1:11000];
#                 set yrange [1e-4:1e-2]; 
#                 set title 'Comparison with analytical zeroth moment'; 
#                 set xlabel 'Time';
#                 set ylabel 'Zeroth moment';
#                 plot 'Ntol.tmp' title 'Numerical', 'Ntol_ana.tmp' title 'Analytical' with lines; 
#                 pause 10
#             "
# gnuplot -e "set xrange [1e-4:1e1];
#                 set yrange [1e-1:1e3]; 
#                 set title 'Comparison with analytical number density'; 
#                 set logscale y 10; 
#                 set logscale x 10; 
#                 set xlabel 'Pivots';
#                 set ylabel 'Number density';
#                 plot 'numden.tmp' title 'Numerical', 'anaden_ini.tmp' title 'Inital' with lines, 'anaden_final.tmp' title 'Analytical' with lines;
#                 pause 10
#             "

gnuplot -e "set xrange [1e-5:1e4];
set size square 0.8,1;
                set xtics 1e-4,1e2,1e4;
                set yrange [1e-15:1e-1];
                set ytics 1e-15,1e5,35;
                set title '2d aggregation'; 
                set logscale y 10; 
                set logscale x 10; 
                set xlabel 'Pivots';
                set ylabel 'Number density';
                set key width -4;
                set key reverse Right;
                plot 'Nneww.tmp' title 'initial' with lines, 'numden2.tmp' title 'tau = 25', 'numden2_new.tmp' title 'tau = 5';
                pause 10
            "

# gnuplot -e "set xrange [1e-4:1e1];
#                 set yrange [1e-1:1e3]; 
#                 set title 'Comparison with matlab number density'; 
#                 set logscale y 10; 
#                 set logscale x 10; 
#                 set xlabel 'Pivots';
#                 set ylabel 'Number density';
#                 plot 'numden.tmp' title 'C++', 'matlab_result.txt' title 'Matlab' with lines;
#                 pause 10
#             "