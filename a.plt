reset
set term gif animate
set output "animate.gif"
unset key
#do for [ii=1:100] { sp [-1:1][-1:1][-1:1 ]'paths.xyz' u 2:3:4 every :::1::ii w d, 'paths.xyz' u 2:3:4 every :::ii::ii w l lw 2}
do for [ii=1:1000] { sp [-1:1][-1:1][-1:1 ]'paths.xyz' u 3:4:5 every :::1::ii w d, 'paths.xyz' u 3:4:5 every :::ii::ii w l lw 2, 'paths.xyz' u 3:4:5:2 every :::ii::ii w labels}
