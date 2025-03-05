# plot individual frames
set term pngcairo size 600,600
set yrange [0:16];
set xrange [0:32];
set size ratio -1
set pm3d interpolate 0,0

# do for [i=0:13000] {
#   fin = sprintf("output/surface_%d.dat", i)

#   stats fin nooutput
#   if (!GPVAL_ERRNO) {
#     # extract time and set the title
#     x = system("head -1 ".fin." | awk '{print $3}'")
#     set title x


#     fout = sprintf("output/plot_surface_%d.png",i)
#     set output fout
#     # plot fin u 2:1:3 w image notitle
#     splot fin u 2:1:($6 + 0.99):($6 + 1.0) notitle with image, \
#           fin u 2:1:3:3 notitle with pm3d
#     # splot fin u 2:1:6 notitle with pm3d
#   }
# }
# set output

MARGIN=0.1

set zrange [1 - 2 * MARGIN:1+MARGIN]
set cbrange [1-MARGIN:1+MARGIN]

# plot animation
set term gif size 500,400 animate delay 5 loop 0 optimize
set output "output/plot_surface.gif"
do for [i=0:13000] {
  fin = sprintf("output/surface_%d.dat", i)

  stats fin nooutput
  if (!GPVAL_ERRNO) {
    # extract time and set the title
    x = system("head -1 ".fin." | awk '{print $3}'")
    set title x

    # plot fin u 2:1:3 w image notitle

    splot fin u 2:1:($6 + 1 - 2 * MARGIN):($6 + 1.0) notitle with image, \
          fin u 2:1:3:3 notitle with pm3d
    # splot fin u 2:1:($6 + 0.99):($6 + 1.0) notitle with image
    # splot fin u 2:1:3:6 notitle with pm3d


    # splot fin u 2:1:3 notitle with pm3d

    # splot fin u 2:1:6 notitle with pm3d
  }
}
set output
