
# simply glue 4 simul vs observed slopes comarisons and name then Figure 5 ...
montage `ls|grep IVYWREL_bact| sort -k3 -t'.'` -geometry +0-10 -tile 2x bact_IVYWREL.png
montage `ls|grep R20_bact| sort -k3 -t'.'` -geometry +0-10 -tile 2x bact_R20.png
montage `ls|grep Akashi_bact| sort -k3 -t'.'` -geometry +0-10 -tile 2x bact_Akashi.png
