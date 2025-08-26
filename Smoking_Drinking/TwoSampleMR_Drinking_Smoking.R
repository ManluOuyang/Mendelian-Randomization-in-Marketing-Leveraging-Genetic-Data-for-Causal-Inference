install.packages('remotes')
library(remotes)
install_github("MRCIEU/TwoSampleMR")
remotes::install_github("mrcieu/ieugwasr")

library(TwoSampleMR)
library(dplyr)
library(ieugwasr)

#####======== Read Data ========#####
dw_exp_dat <- read_exposure_data(
  filename = "~/MR/DrinksPerWeek.csv",
  sep = ",",
  snp_col = "SNP_ID",
  beta_col = "Estimate.effect..direction.",
  se_col = "SE",
  effect_allele_col = "Effect.allele..EA.",
  other_allele_col = "other.allele",
  eaf_col = "Effect.allele.frequency..EAF.",
  pval_col = "P.value",
  samplesize_col = "Sample.size",
  chr_col = "Chromosome",
  pos_col = "Position"
)#11907805 SNPs

# Keep only significant SNPs
dw_exp_dat_sig <- dw_exp_dat %>% filter(pval.exposure <= 5e-8) #5204 SNPs

#####======== Check API KEY ========######
ieugwasr::get_opengwas_jwt()

ieugwasr::user()

#####======== Clumping ========#####
#options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

dw_exp_dat_clump <- clump_data(dw_exp_dat_sig) #37 SNPs

#####======== Import Outcome Data and match it with Exposure SNPs ========#####
si_out_dat <- read_outcome_data(
  snps = dw_exp_dat_clump$SNP,
  filename = "~/MR/SmokingInitiation.csv",
  sep = ",",
  snp_col = "SNP_ID",
  beta_col = "Estimate.effect..direction.",
  se_col = "SE",
  effect_allele_col = "Effect.allele..EA.",
  other_allele_col = "other.allele",
  eaf_col = "Effect.allele.frequency..EAF.",
  pval_col = "P.value",
  samplesize_col = "Sample.size",
  chr_col = "Chromosome",
  pos_col = "Position"
) #11634311 SNPs


#####======== Perform MR ========#####
# Harmonization
dat_ds <- harmonise_data(
  exposure_dat = dw_exp_dat_clump,
  outcome_dat = si_out_dat
)

dat_drinking_smoking <- subset(dat_ds, mr_keep)

# test outliers by using ivw_radial
outlier_ds <- ivw_radial(dat_drinking_smoking)

# remove outliers
outlier_ds <- outlier_ds$outliers
dat_ds_ready <- dat_drinking_smoking[which(!(dat_drinking_smoking$SNP %in% outlier_ds$SNP)),]#24 SNPs

# Run MR
res_ds <- mr(dat_ds_ready, method_list = c("mr_ivw_mre", "mr_raps", 
                                         "mr_weighted_median","mr_weighted_mode", "mr_egger_regression"))
res_ds

# Leave one out Analysis
res_loo_ds <- mr_leaveoneout(dat_ds_ready)

# Plot the result
setwd('~/MR/Plots')
Scatter_drinking_on_smoking <- mr_scatter_plot(res_ds, dat_ds_ready)

Scatter_drinking_on_smoking[[1]]

# Forest Plot
res_single <- mr_singlesnp(dat_ds_ready)
Forest_drinking_on_smoking <- mr_forest_plot(res_single)
Forest_drinking_on_smoking[[1]]

# Funnel Plot
res_single <- mr_singlesnp(dat_ds_ready)
Funnel_drinking_on_smoking <- mr_funnel_plot(res_single)
Funnel_drinking_on_smoking[[1]]

# Leave-one-out Plot
loo_ds <- mr_leaveoneout_plot(res_loo_ds)
loo_ds[[1]]