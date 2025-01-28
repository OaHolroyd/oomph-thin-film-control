# plot individual frames
set term pngcairo size 600,600
set yrange [0.9:1.1];
set xrange [0:20];

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
set term gif size 600,600 animate delay 100 loop 0 optimize
set output "output/mesh.gif"
do for [i=0:13000] {
  fin = sprintf("output/spine_step_%d.dat", i)

  stats fin nooutput
  if (!GPVAL_ERRNO) {
    plot fin using 1:2 w l
  }
}
set output
