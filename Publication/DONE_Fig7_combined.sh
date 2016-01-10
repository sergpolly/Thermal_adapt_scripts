python DONE_MainFig7_costs.py bact
python DONE_MainFig7_costs.py arch

convert \
-append \
-density 150 \
Fig7.arch.pdf \
Fig7.bact.pdf \
exp_MTAD_all_all.bact.summary_Fig7_costs.pdf \
Figure7.png








# # simply glue 4 simul vs observed slopes comarisons and name then Figure 5 ...
# montage \
# Fig7.arch.pdf \
# Fig7.bact.pdf \
# exp_MTAD_all_all.bact.summary_Fig7_costs.pdf \
# -geometry +0-1 \
# -tile 1x \
# -border 0 \
# -bordercolor white \
# -background white \
# Figure7.pdf

# # exp_MTAD_all_all.arch.summary_Fig7_costs.pdf \


# # exp_MTAD_all_all.arch.summary_Fig7_costs.pdf \
# # exp_MTAD_all_all.bact.summary_Fig7_costs.pdf \
