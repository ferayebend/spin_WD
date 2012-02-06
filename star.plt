set terminal postscript portrait enhanced color solid "times" 20
set xlabel '{/=20 t (years)}'

set xrange [1:]
#set yrange[1e33:1e39]
set logscale x
set logscale y
set style data lines
set size 1,0.5
set mxtics 10
set mytics 10

set format x "10^{%T}"
set format y "10^{%T}"
set key top right
#------------------------------
#set label  1 "{/Symbol=15 m}_{/Text=12 30}{/Text=15 =2}" at 30000, 1.e34
#set label 10 "{/Symbol=15 m}_{/Text=12 30}{/Text=15 =10}" at 350, 1e34
#set label 12 "{/Symbol=15 w}_{cr}=1.05" at 200, 3e35
#set label 13 "{/Symbol=15 h}=0.01" at 200, 2e35
#-------------------------------
#luminosity
#set zero 1e-32
set style line 1 lt 1 lw 3 pt 3 ps 0.5
set style line 2 lt 2 lw 3 pt 3 ps 0.5
set style line 3 lt 3 lw 3 pt 3 ps 0.5
set style line 4 lt 4 lw 3 pt 3 ps 0.5


set output "J_star.ps"
set ylabel '{/=20 J_{star}/1E51'

#plot "star.out"  using 1:3 title "" ls 1
plot "star.out"  using 1:($9==1 ? $3 : 1/0) title "" ls 1,    \
     "star.out"  using 1:($9==0 ? $3 : 1/0) title "" ls 2


set output "period.ps"
set ylabel '{/=20 P (s)'

#plot "star.out"  using 1:5 title "" ls 1
plot "star.out"  using 1:($9==1 ? $5 : 1/0) title "" ls 1,    \
     "star.out"  using 1:($9==0 ? $5 : 1/0) title "" ls 2

set output "w_s.ps"
set ylabel '{/=20 w_s' 

#plot "star.out"  using 1:7 title "" ls 1
plot "star.out"  using 1:($9==1 ? $7 : 1/0) title "" ls 1,    \
     "star.out"  using 1:($9==0 ? $7 : 1/0) title "" ls 2

set yrange [1:1.5]
set mytics 1
set output "M_star.ps"
set ylabel '{/=20 M_{star}/M_{sun}'

plot "star.out"  using 1:($9==1 ? $2 : 1/0) title "" ls 1,    \
     "star.out"  using 1:($9==0 ? $2 : 1/0) title "" ls 2   
