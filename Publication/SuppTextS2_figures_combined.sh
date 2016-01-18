# # # generate data first (pictures to glue) ...
python DONE_brooms_MJ99_Argentina_bootstrap.py all all arch
python DONE_brooms_MJ99_Argentina_bootstrap.py all all bact


# # # simply glue 4 simul vs observed slopes comarisons and name then Figure 5 ...
montage \
exp_MTAD_all_all.arch.summary_SuppFig3.png \
exp_MTAD_all_all.bact.summary_SuppFig3.png \
-geometry +0+0 \
-tile 1x \
-border 0 \
-bordercolor white \
-background white \
SuppText2_Fig1.png


# # # simply glue 4 simul vs observed slopes comarisons and name then Figure 5 ...
montage \
exp_MTAD_all_all.arch.summary_MainFig6.png \
exp_MTAD_all_all.bact.summary_MainFig6.png \
-geometry +0+0 \
-tile 1x \
-border 0 \
-bordercolor white \
-background white \
SuppText2_Fig2.png

