set term pngcairo size 600,600
set yrange [0.9:1.1];

do for [i=1:10] {
    fin = sprintf("output/spine_interface_%d.dat",i)
    fout = sprintf("plot/interface_%d.png",i)
    set output fout
    plot fin using 1:2 w l
}
set output


set term gif size 600,600 animate delay 12 loop 0 optimize
set output "plot/interface.gif"
do for [i=1:10] {
    fin = sprintf("output/spine_interface_%d.dat",i)
    plot fin using 1:2 w l
}
set output
