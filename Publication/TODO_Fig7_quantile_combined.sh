
# simply glue 4 simul vs observed slopes comarisons and name then Figure 5 ...
montage \
../ArchNew/BOOTSTRAP/R20_arch_qunatile_trend.original.trop.png \
../BOOTSTRAPS/R20_bact_qunatile_trend.original.trop.png \
-geometry +0-10 \
-tile 2x \
Figure7_quantiles.png
# montage `ls|grep R20_bact| sort -k3 -t'.'` -geometry +0-10 -tile 2x bact_R20.png
# montage `ls|grep Akashi_bact| sort -k3 -t'.'` -geometry +0-10 -tile 2x bact_Akashi.png
