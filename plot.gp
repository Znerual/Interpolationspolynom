set term pdfcairo
set output "Interpolationspolynom.pdf"
set title "Interpolationspolynom"




plot "polynom.dat" using 1: 2 w l title "Interpolation", "int-pol.txt" using 1:2 title "Stützstellen", "int-pol.txt" smooth csplines title "cubic Spline"

unset term
unset outp