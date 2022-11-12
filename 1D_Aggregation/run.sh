#!/bin/bash

gnuplot -e "set xrange [1:11000];
                set yrange [1e-4:1e-2]; 
                set title 'Comparison with analytical zeroth moment'; 
                set xlabel 'Time';
                set ylabel 'Zeroth moment';
                plot 'Ntol.tmp' title 'Numerical', 'Ntol_ana.tmp' title 'Analytical' with lines; 
                pause 10
            "
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

# gnuplot -e "set xrange [1e-5:1e4];
# set size square 0.8,1;
#                 set xtics 1e-4,1e2,1e4;
#                 set yrange [1e-15:35];
#                 set ytics 1e-15,1e5,35;
#                 set title 'Comparison with analytical number density'; 
#                 set logscale y 10; 
#                 set logscale x 10; 
#                 set xlabel 'Pivots';
#                 set ylabel 'Number density';
#                 set label 1 'N(t)/N(0) =0.74' at 1e-3,1e-2;
#                 set label 2 'N(t)/N(0) =8e-4' at 5e-3,5e-6;
#                 set key width -4;
#                 set key reverse Right;
#                 plot 'Nneww.tmp' title 'initial' with lines, 'numden2.tmp' title 's= 2 (0.074)', 'numden15.tmp' title 's=1.5 (0.074)', 'numden125.tmp' title 's=1.25 (0.074)','numden2_new.tmp' title 's=2 (8e-4)', 'numden15_new.tmp' title 's=1.5 (8e-4)', 'numden125_new.tmp' title 's=1.25 (8e-4)';
#                 pause 10
#             "

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