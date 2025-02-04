# plot individual frames
set term pngcairo size 600,200
set yrange [-1.0:1.0];
set xrange [0:30];
set size ratio 0.25

do for [i=0:13000] {
  fin = sprintf("output/1d_%d.dat", i)

  stats fin nooutput
  if (!GPVAL_ERRNO) {
    # extract time and set the title
    x = system("head -1 ".fin." | awk '{print $3}'")
    set title x

    fout = sprintf("output/plot_1d_%d.png",i)
    set output fout
    plot fin using 1:4 w l t "control" lt rgb "#FF0000", \
         fin using 1:($2 - 1) w l t "interface" lt rgb "#0000FF"
  }
}
set output

# plot animation
set term gif size 600,200 animate delay 12 loop 0 optimize
set output "output/plot_1d.gif"
do for [i=0:13000] {
  fin = sprintf("output/1d_%d.dat", i)

  stats fin nooutput
  if (!GPVAL_ERRNO) {
    # extract time and set the title
    x = system("head -1 ".fin." | awk '{print $3}'")
    set title x

    plot fin using 1:4 w l t "control" lt rgb "#FF0000", \
         fin using 1:($2 - 1) w l t "interface" lt rgb "#0000FF"
  }
}
set output
