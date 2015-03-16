unset key
set grid
set xlabel "Real part"
set ylabel "Imag part"
set title "Eigenvalues for cylinder in channel at Re = 150"
set term pdf
set out 'eig.pdf'
p 'eig.dat' w p pt 6
