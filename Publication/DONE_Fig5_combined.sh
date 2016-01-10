# generate data first (pictures to glue) ...
python DONE_brooms_MJ99_Argentina.py all all arch
python DONE_brooms_MJ99_Argentina.py all all bact
python DONE_brooms_MJ99_Argentina.py cai trop arch
python DONE_brooms_MJ99_Argentina.py cai trop bact
# simply glue 4 simul vs observed slopes comarisons and name then Figure 5 ...
montage \
exp_MTAD_all_all.arch.summary_Slopes_sim_vs_exp_protbact.png \
exp_MTAD_all_all.bact.summary_Slopes_sim_vs_exp_protbact.png \
exp_MTAD_cai_trop.arch.summary_Slopes_sim_vs_exp_protbact.png \
exp_MTAD_cai_trop.bact.summary_Slopes_sim_vs_exp_protbact.png \
-geometry +0+0 \
-tile 2x \
-border 0 \
-bordercolor white \
-background white \
Figure5.png

# -geometry +0-15 \
# glue 4 Figs together as is, leaving 0-natural hspace, and decreasing vspace by 15 (-15)
