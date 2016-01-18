# # ~50 random proteins per organism per sample
# BOOTSTRAP_arch_org-all_cds-random.dat
# BOOTSTRAP_bact_org-all_cds-random.dat
python TODO_compare_bootstrap_slopes.py \
BOOTSTRAP_arch_org-all_cds-random.dat \
BOOTSTRAP_bact_org-all_cds-random.dat \
0.73 \
SuppS1_ribo_test.png

# # ~40% random organisms per sample, proteome-wide
# BOOTSTRAP_arch_org-random_cds-all.dat
# BOOTSTRAP_bact_org-random_cds-all.dat
python TODO_compare_bootstrap_slopes.py \
BOOTSTRAP_arch_org-random_cds-all.dat \
BOOTSTRAP_bact_org-random_cds-all.dat \
0.554 \
SuppS1_CUS_proteome_test.png



# # ~40% random organisms per sample, PHX proteins only (10th CAI decile)
# BOOTSTRAP_arch_org-random_cds-cai.dat
# BOOTSTRAP_bact_org-random_cds-cai.dat
python TODO_compare_bootstrap_slopes.py \
BOOTSTRAP_arch_org-random_cds-cai.dat \
BOOTSTRAP_bact_org-random_cds-cai.dat \
0.871 \
SuppS1_CUS_test_PHX.png


# # ~50 random proteins per organism (CUS-only) per sample
# BOOTSTRAP_arch_org-trop_cds-random.dat
# BOOTSTRAP_bact_org-trop_cds-random.dat
python TODO_compare_bootstrap_slopes.py \
BOOTSTRAP_arch_org-trop_cds-random.dat \
BOOTSTRAP_bact_org-trop_cds-random.dat \
0.871 \
SuppS1_CUS_PHX_test.png



# ##################################################
# # 	very strange CAI  as a organism criteria ....????.....
# #   most likely these are bad data ....
# # ~50 random proteins per organism (CUS only) per sample
# BOOTSTRAP_arch_org-cai_cds-random.dat
# BOOTSTRAP_bact_org-cai_cds-random.dat
python TODO_compare_bootstrap_slopes.py \
BOOTSTRAP_arch_org-cai_cds-random.dat \
BOOTSTRAP_bact_org-cai_cds-random.dat \
0.871 \
SuppS1_CUS_PHX_test_XXX.png


# glue PHX related pictures together ....
montage \
SuppS1_CUS_test_PHX.png \
SuppS1_CUS_PHX_test.png \
-geometry +0+0 \
-tile 2x \
-border 0 \
-bordercolor white \
-background white \
SuppS1_PHX_test.png










