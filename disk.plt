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

#set zero 1e-32
set style line 1 lt 1 lw 4 pt 3 ps 0.5
set style line 2 lt 2 lw 3 pt 3 ps 0.5
set style line 3 lt 3 lw 3 pt 3 ps 0.5
set style line 4 lt 4 lw 3 pt 3 ps 0.5



#***********************************************
set output "Mdot.ps"
set ylabel '{/=20 dM/dt g s^{-1}'

#plot "disk.out"  using 1:2 title "" ls 1

plot "disk.out"  using 1:($7==1 ? $2 : 1/0) title "" ls 1,    \
     "disk.out"  using 1:($7==0 ? $2 : 1/0) title "" ls 2   


#***********************************************

set output "R_in.ps"
set ylabel '{/=20 R_{in}/R_{star}'



plot "disk.out"  using 1:($7==1 ? $3 : 1/0) title "" ls 1,    \
     "disk.out"  using 1:($7==0 ? $3 : 1/0) title "" ls 2   

#***********************************************
set output "L.ps"
set ylabel '{/=20 L'

plot "disk.out"  using 1:4 title "L_{acc}" ls 1,   \
     "disk.out"  using 1:5 title "L_{disk}" ls 2

plot "disk.out"  using 1:($7==1 ? $4 : 1/0) title "" ls 1,    \
     "disk.out"  using 1:($7==0 ? $4 : 1/0) title "" ls 2,    \
     "disk.out"  using 1:($7==1 ? $5 : 1/0) title "" ls 1,    \
     "disk.out"  using 1:($7==0 ? $5 : 1/0) title "" ls 2      

#***********************************************

set output "f.ps"
set ylabel '{/=20 R_{in}/R_L'

plot "disk.out"  using 1:($7==1 ? $6 : 1/0) title "" ls 1,    \
     "disk.out"  using 1:($7==0 ? $6 : 1/0) title "" ls 2   
