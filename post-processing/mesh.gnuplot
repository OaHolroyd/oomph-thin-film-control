# plot individual frames
set term pngcairo size 1200,400
set yrange [0.0:2.0];
set xrange [0:20];
set size ratio 0.25

do for [i=0:13000] {
  fin = sprintf("output/spine_step_%d.dat", i)

  stats fin nooutput
  if (!GPVAL_ERRNO) {
    fout = sprintf("output/mesh_%d.png",i)
    set output fout
    plot fin using 1:2 w l
  }
}
set output

# plot animation
set term gif size 600,200 animate delay 25 loop 0 optimize
set output "output/mesh.gif"
do for [i=0:13000] {
  fin = sprintf("output/spine_step_%d.dat", i)

  stats fin nooutput
  if (!GPVAL_ERRNO) {
    plot fin using 1:2 w l
  }
}
set output
