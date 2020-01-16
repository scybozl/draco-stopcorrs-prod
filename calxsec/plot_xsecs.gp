set term pdfcairo

set output 'xsecs.pdf'

set title 'Stop pair production runs'
set xlabel 'Stop mass m_{~t{.6\~}}'
set ylabel 'Neutralino mass m_{{/Symbol c}}'

set xrange [160:290]
set yrange [0:160]
set jitter

# define a function to map strings to palette indices
map_color(string) = ( \
  string eq 'green' ? 0x013220 : \
  string eq 'red' ? 0xff0000 : \
  0x000000)

plot 'badpoints_L.dat' using (column(1)-1):2:(map_color(stringcolumn(3))) notitle with points pt 7 ps 0.5 lc rgbcolor variable, \
     'badpoints_R.dat' using (column(1)+1):2:(map_color(stringcolumn(3))) notitle with points pt 7 ps 0.5 lc rgbcolor variable

