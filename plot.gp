set term pdfcairo
set output "Interpolationspolynom.pdf"
set title "Interpolationspolynom"
 
set key left bottom

plot "polynom.dat" using 1: 2 title "Interpolation" w l, "int-pol.txt" using 1:2 title "Stützstellen", "int-pol.txt" smooth csplines title "cubic Spline"

set title "Trigonometrische Reihe"
plot "TrigReihe.dat" using 1:2 title "Reihe" w l , 3/(x**2 +1) 
unset term
unset outp