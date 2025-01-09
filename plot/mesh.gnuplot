# get spine or elastic
stats "output/spine_output.dat" nooutput
if (GPVAL_ERRNO) {
  mode = "elastic"
} else {
  mode = "spine"
}

# plot individual frames
set term pngcairo size 600,600
set yrange [0:1.1];

do for [i=1:10] {
  fin = sprintf("output/%s_step_%d.dat", mode, i)
  stats fin nooutput
  if (!GPVAL_ERRNO) {
    fout = sprintf("plot/mesh_%d.png",i)
    set output fout
    plot fin using 1:2
  }
}
set output

# plot animation
set term gif size 600,600 animate delay 12 loop 0 optimize
set output "plot/mesh.gif"
do for [i=1:10] {
  fin = sprintf("output/%s_step_%d.dat", mode, i)
  stats fin nooutput
  if (!GPVAL_ERRNO) {
    plot fin using 1:2
  }
}
set output
