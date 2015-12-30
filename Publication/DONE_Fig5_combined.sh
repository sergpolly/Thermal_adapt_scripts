
# simply glue 4 simul vs observed slopes comarisons and name then Figure 5 ...
montage \
exp_MTAD_all_all.arch.summary_Slopes_sim_vs_exp_protbact.png \
exp_MTAD_all_all.bact.summary_Slopes_sim_vs_exp_protbact.png \
exp_MTAD_cai_trop.arch.summary_Slopes_sim_vs_exp_protbact.png \
exp_MTAD_cai_trop.bact.summary_Slopes_sim_vs_exp_protbact.png \
-geometry +0-15 \
-tile 2x \
-border 0 \
-bordercolor white \
-background white \
Figure5.png

# glue 4 Figs together as is, leaving 0-natural hspace, and decreasing vspace by 15 (-15)
