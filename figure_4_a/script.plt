set border 3
set xtics nomirror font ",16"
set ytics nomirror font ",16"
set grid 

set xlabel "x_1" font ",16"

max(x,y)=(x>y?x:y)

min(x,y)=(x<y?x:y)

f(x) = b*x*x+c*x+d

set colors classic

set key bottom center font ",14"
# set key top right font ",14"

set xrange [-1.5:1.5]
set yrange [0:1.2]

plot 'exact_solution' using ($1):($2) title 'Exact solution' with lines linestyle 3 lw 2, 'mc_estimates' using ($1):($2) title 'Numerical solution' with points pt 6 ps 0.7 lc -1

pause -1 

set terminal pdf enhanced 

set output "figure.pdf"

replot

set output 

set terminal x11
