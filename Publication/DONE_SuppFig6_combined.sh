# # generate data first (pictures to glue) ...
python DONE_brooms_MJ99_Argentina.py all all arch
python DONE_brooms_MJ99_Argentina.py all all bact
python DONE_brooms_MJ99_Argentina.py cai trop arch
python DONE_brooms_MJ99_Argentina.py cai trop bact
# # simply glue 4 simul vs observed slopes comarisons and name then Figure 5 ...
montage \
exp_MTAD_all_all.arch.summary_SuppFigure4.png \
exp_MTAD_all_all.bact.summary_SuppFigure4.png \
exp_MTAD_cai_trop.arch.summary_SuppFigure4.png \
exp_MTAD_cai_trop.bact.summary_SuppFigure4.png \
-geometry +0+0 \
-tile 2x \
-border 0 \
-bordercolor white \
-background white \
SuppFigure6.png

# -geometry +0-15 \
# glue 4 Figs together as is, leaving 0-natural hspace, and decreasing vspace by 15 (-15)


