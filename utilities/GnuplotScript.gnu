##gnuplot -e "a=1427191783" ../../utilities/GnuplotScript.gnu
reset
set key off

##set output
set term wxt

#set terminal pngcairo font ",10" enhanced
#set terminal pngcairo size 1024,1024
#set output sprintf("picture_%d.png",a)

#set terminal jpeg font ",10"
#set terminal jpeg size 1024,768
#set output sprintf("picture_%d.jpg",a)

#set terminal epscairo font ",5"
#set output sprintf("picture_%d.eps",a)

##set ranges and view
set xrange [0:*]
set zrange [0:*]
#set cbrange [0:2.5]

set xlabel  '{/Symbol r}'
set ylabel "z"

set view map
set size ratio 2

set pm3d interpolate 2,6

##plot
splot 'plot_data_'.a u 1:2:4 w pm3d, 'plot_horizon_common_'.a u 1:2:(0) w l lt rgb "green", 'plot_horizon_up_'.a u 1:2:(0) w l lt rgb "green", 'plot_horizon_down_'.a u 1:2:(0) w l lt rgb "green", 'plot_MITS_common_'.a u 1:2:(0) w l lt rgb "red", 'plot_MITS_up_'.a u 1:2:(0) w l lt rgb "red", 'plot_MITS_down_'.a u 1:2:(0) w l lt rgb "red"
pause -1
