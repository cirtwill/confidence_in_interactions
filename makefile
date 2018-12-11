# run by gmake -f makefile.mk
# output indicates if the rule was run

#$^ is for "all the prerequisites of this file"
#$......
# $@ is short for "the target of this rule"


# To calculate the alpha, beta, etc. presented in the main text, run code/Zillis_prior/calculating_prior_distributions.R
# Then enter calculate_distribution(calculate_parameters(gp_int_probs,0,0)) for the prior
# calculate_distribution(calculate_parameters(gp_int_probs,374,0)) for the most informed posterior


# Introduce a phony target to make everything run
all : manuscript/CI_interactions.pdf

# Zillis and Zillertal priors

# Make prior webs for both Zillis and Zillertal subwebs
data/Salix_example/Zillis/%_prior.csv : code/Zillis_prior/create_Zillis_priors.py data/Salix_example/Zillis/Zillis_web.csv data/Salix_example/Zillis/Zillertal_web.csv
  cd code/Zillis_prior && \
  python create_Zillis_priors.py && \
  cd ../../

# Calculate distributions that will be used to make posterior webs
data/Salix_example/Zillis/%/posterior_probabilities_%.tsv : code/Zillis_prior/calculating_prior_distributions.R data/Salix_example/cooccur_interact_galler_salix.csv data/Salix_example/cooccur_interact_galler_parasit.csv ../../data/Salix_example/Zillis/Zillis_SG_prior.csv ../../data/Salix_example/Zillis/oZillis_GP_prior.csv
  cd code/Zillis_prior && \
  Rscript calculating_prior_distributions.R && \
  cd ../../

# Create posterior webs
data/Zillis_webs/posterior/%.web : code/Zillis_prior/create_Zillis_posterior_webs.py data/Salix_example/Zillis/Salix_Galler/posterior_probabilities_*.tsv data/Salix_example/Zillis/Galler_Parasitoid/posterior_probabilities_*.tsv
  cd code/Zillis_prior && \
  python create_Zillis_posterior_webs.py && \
  cd ../../


# Calculate subweb properties
data/Zillis_webs/%_%_table_%.tsv : data/Zillis_webs/posterior/*.web code/calculate_Ziller_simweb_properties.py  
  cd code/Zillis_prior && \
  python calculate_Ziller_simweb_properties.py && \
  cd ../../
  

# Create subwebs, calculate properties of original webs
data/randomised_webs/%_NODF_table.tsv : code/calculate_NODF.R code/calculate_sim_network_properties.py code/create_posterior_webs.py data/Salix_example/Salix_Galler/posterior_probabilities.tsv
  cd code/ && \
  python create_posterior_webs.py && \
  python calculate_sim_network_properties.py && \
  Rscript calculate_NODF.R && \
  cd ../


# Calculate NODF of Zillis and Zillertal posterior webs
data/Zillis_webs/%_NODF_table_%.tsv : code/Zillis_prior/calculate_NODF.R data/Zillis_webs/*/*.web data/Zillis_webs/*/*/*.web
  cd code/Zillis_prior && \
  Rscript calculate_NODF.R && \
  cd ../


# FIGURES

# Histogram of interaction and cooccurence frequencies
manuscript/figures/Salix_Galler_histogram.eps : code/figure_generation/observation_interaction_histogram.py data/Salix_example/*/nonzero_cooccurances.tsv data/Salix_example/*/observed_interactions.tsv data/Salix_example/cooccur_interact_galler_*.csv
  # Datafiles generated with R code at top of fig file
  cd code/figure_generation && \
  python observation_interaction_histogram.py && \
  cd ../../

# Making figures
manuscript/figures/GPr_pdfs_increasing_N_Zillis.eps : code/figure_generation/pdfs_and_N.py data/Salix_example/Zillis/*/distfigure_*vals.tsv data/Salix_example/Zillis/*/distfigure_MLEs.tsv
  # Datafiles generated with R code at top of fig file
  cd code/figure_generation && \
  python pdfs_and_N.py && \
  cd ../../

# Posterior property figure, Barbour prior, shows SG and GP in one figure
manuscript/figures/Salix_Galler_posterior_properties.eps : code/figure_generation/posterior_property_dotplot.py data/randomised_webs/*_table.tsv
  cd code/figure_generation && \
  python posterior_property_dotplot.py &&\
  cd ../../

# Posterior property figure, Zillis prior
manuscript/figures/GP_posterior_properties_Zillis.eps : code/figure_generation/Zillis_posterior_plots.py data/Zillis_webs/*_table_*.tsv
  cd code/figure_generation && \
  python Zillis_posterior_plots.py && \
  cd ../../

# Samples to hit a threshold, in SI?
manuscript/figures/GP_samples_and_cdfs_Zillis.eps : code/figure_generation/samples_for_thresholds.py data/Salix_example/Zillis/*/samplefigure.tsv data/Salix_example/Zillis/*/samples_for_threshold.tsv
  # Datafiles generated with R code at top of fig file
  cd code/figure_generation && \
  python samples_for_thresholds.py && \
  cd ../../

# Plotting the Upper curve of CI, after Dom's R figure
manuscript/figures/upper_limit_DG.eps : code/figure_generation/upper_CI_after_DG.py data/R/CI_and_samples.tsv manuscript/figure_binomial.R
  cd code/figure_generation && \
  python upper_CI_after_DG && \
  cd ../../

# prior_web_para_only.csv - Barbour host-parasitoid prior
# binary_prior_web_SG.csv - binary Barbour salix-galler prior
# binary_prior_web_para_only.csv - binary version of Barbour host-parasitoid prior
# heatmap.py - not used
# posterior_property_plot.py - replaced with dotplot.


# Making the paper
manuscript/CI_interactions.pdf : manuscript/CI_interactions.tex manuscript/ecol_let.bst manuscript/manual_abbrev.bib manuscript/figures/Salix_Galler_histogram.eps manuscript/figures/upper_limit_DG.eps manuscript/figures/Salix_Galler_pdfs_increasing_N_Zillis.eps manuscript/figures/Salix_Galler_samples_and_cdfs_Zillis.eps manuscript/figures/Salix_Galler_posterior_properties_Zillis.eps
  cd manuscript && \
  latex CI_interactions.tex && \
  bibtex CI_interactions.aux && \
  latex CI_interactions.tex && \
  latex CI_interactions.tex && \
  dvips CI_interactions.dvi && \
  ps2pdf CI_interactions.ps && \
  cd ../

manuscript/supplemental.pdf : manuscript/supplemental.tex manuscript/ecol_let.bst manuscript/manual_abbrev.bib manuscript/figures/Salix_Galler_pdfs_increasing_N.eps manuscript/figures/Salix_Galler_samples_and_cdfs.eps manuscript/Figures/Salix_Galler_posterior_properties.eps
  cd manuscript && \
  latex supplemental.tex && \
  bibtex supplemental.aux && \
  latex supplemental.tex && \
  latex supplemental.tex && \
  dvips supplemental.dvi && \
  ps2pdf supplemental.ps && \
  cd ../

manuscript/coverletter.pdf : manuscript/coverletter.tex 
  cd manuscript && \
  latex coverletter.tex && \
  latex coverletter.tex && \
  dvips coverletter.dvi && \
  ps2pdf coverletter.ps && \
  cd ../