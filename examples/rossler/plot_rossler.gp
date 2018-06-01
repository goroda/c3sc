set term png size 1000,1000 fontscale 2 linewidth 3

set output "rossler_trajs.png"

traj = "traj.dat"

load '/home/goroda/Software/gnuplot/gnuplot-colorbrewer/qualitative/Dark2.plt'
# load '/home/goroda/Software/gnuplot/gnuplot-colorbrewer/qualitative/Accent.plt'


set tmargin 0
set bmargin 0
set lmargin 3
set rmargin 3

set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set yrange [-1:1]
set xrange [0:10]

set multiplot layout 4,1 columnsfirst margins 0.2,0.97,.1,.99 spacing 0,0
set ylabel "x"
unset xtics
plot "traj.dat" u 1:2 w l ls 1 notitle,\
     "traj1.dat" u 1:2 w l ls 2 notitle,\
     "traj2.dat" u 1:2 w l ls 3 notitle,\
     "traj3.dat" u 1:2 w l ls 4 notitle,\
     "traj4.dat" u 1:2 w l ls 5 notitle,\
     "traj5.dat" u 1:2 w l ls 6 notitle
unset xlabel
set ylabel "y"
plot "traj.dat" u 1:3 w l ls 1 notitle,\
     "traj1.dat" u 1:3 w l ls 2 notitle,\
     "traj2.dat" u 1:3 w l ls 3 notitle,\
     "traj3.dat" u 1:3 w l ls 4 notitle,\
     "traj4.dat" u 1:3 w l ls 5 notitle,\
     "traj5.dat" u 1:3 w l ls 6 notitle

set ylabel "z"
plot "traj.dat" u 1:4 w l ls 1 notitle,\
     "traj1.dat" u 1:4 w l ls 2 notitle,\
     "traj2.dat" u 1:4 w l ls 3 notitle,\
     "traj3.dat" u 1:4 w l ls 4 notitle,\
     "traj4.dat" u 1:4 w l ls 5 notitle,\
     "traj5.dat" u 1:4 w l ls 6 notitle
set ylabel "u"

set xtics
set xlabel "t"
plot "traj.dat" u 1:5 w l ls 1 notitle,\
     "traj1.dat" u 1:5 w l ls 2 notitle,\
     "traj2.dat" u 1:5 w l ls 3 notitle,\
     "traj3.dat" u 1:5 w l ls 4 notitle,\
     "traj4.dat" u 1:5 w l ls 5 notitle,\
     "traj5.dat" u 1:5 w l ls 6 notitle

unset multiplot

set output "rossler_phase.png"
set xrange [-1:1]
set yrange [-1:1]
set zrange [-1:1]

set xlabel "x"
set ylabel "y"
set zlabel "z"
splot "traj.dat" u 2:3:4 w l ls 1 notitle,\
      "traj1.dat" u 2:3:4 w l ls 2 notitle,\
      "traj2.dat" u 2:3:4 w l ls 3 notitle,\
      "traj3.dat" u 2:3:4 w l ls 4 notitle,\
      "traj4.dat" u 2:3:4 w l ls 5 notitle,\
      "traj5.dat" u 2:3:4 w l ls 6 notitle,\


set output "rossler_cost_slice.png"
load '/home/goroda/Software/gnuplot/gnuplot-colorbrewer/sequential/YlOrRd.plt'

unset zrange
file = "costfunc_100.dat"
set contours
set cntrparam levels 20
set zlabel "Value function" rotate by 90
splot file u 1:2:3 w pm3d notitle
