
# simply glue 4 simul vs observed slopes comarisons and name then Figure 5 ...
montage `ls|grep IVYWREL_arch| sort -k3 -t'.'` -geometry +0-10 -tile 2x arch_IVYWREL.png
montage `ls|grep R20_arch| sort -k3 -t'.'` -geometry +0-10 -tile 2x arch_R20.png
montage `ls|grep Akashi_arch| sort -k3 -t'.'` -geometry +0-10 -tile 2x arch_Akashi.png
