!./plotexample >plots.dat
set log xy
set style data l
plot 'plots.dat' us 1:(1/$2)
