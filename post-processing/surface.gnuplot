# plot individual frames
set term pngcairo size 600,600
set yrange [0:32];
set xrange [0:32];
set size ratio -1
set pm3d interpolate 0,0

do for [i=0:13000] {
  fin = sprintf("output/surface_%d.dat", i)

  stats fin nooutput
  if (!GPVAL_ERRNO) {
    # extract time and set the title
    x = system("head -1 ".fin." | awk '{print $3}'")
    set title x


    fout = sprintf("output/plot_surface_%d.png",i)
    set output fout
    # plot fin u 2:1:3 w image notitle
    splot fin u 2:1:($3 - 1) notitle with pm3d
  }
}
set output

set zrange [-1:1]

# plot animation
set term gif size 400,400 animate delay 5 loop 0 optimize
set output "output/plot_surface.gif"
do for [i=0:13000] {
  fin = sprintf("output/surface_%d.dat", i)

  stats fin nooutput
  if (!GPVAL_ERRNO) {
    # extract time and set the title
    x = system("head -1 ".fin." | awk '{print $3}'")
    set title x

    # plot fin u 2:1:3 w image notitle
    splot fin u 2:1:($3 - 1) notitle with pm3d
  }
}
set output
