set terminal png size 1200,1200
set output "Plot.png"

n=100 #number of intervals
max=0.81 #max value
min=-0.01 #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set boxwidth width*0.9
set style fill solid 0.5 # fill style

set multiplot
set size 1,0.33
set origin 0,0.66
set grid
set xlabel "Position"
set ylabel "Amplitude"
set title "Electrons"
set xrange[min:max]
#count and plot
plot "Output_e.txt" u (hist($1,width)):(1.0) smooth freq w boxes lc rgb"red" notitle

set origin 0,0.33
set title "Cations"
plot "Output_c.txt" u (hist($1,width)):(1.0) smooth freq w boxes lc rgb"blue" notitle

set origin 0,0
set title "Anions"
plot "Output_a.txt" u (hist($1,width)):(1.0) smooth freq w boxes lc rgb"green" notitle
