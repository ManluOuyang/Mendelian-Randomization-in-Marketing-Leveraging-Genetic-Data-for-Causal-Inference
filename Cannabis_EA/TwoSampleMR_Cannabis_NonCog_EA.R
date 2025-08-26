install.packages('remotes')
library(remotes)
install_github("MRCIEU/TwoSampleMR")
remotes::install_github("mrcieu/ieugwasr")

library(TwoSampleMR)
library(dplyr)
library(ieugwasr)

setwd('~/EA_Cannabis_NEW/GWAS')

#####======== Read Cannabis Data ========#####
CUD_exp_dat_new_raw <- read_exposure_data(
  filename = "CUD_new.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P.value",
  chr_col = "CHR",
  pos_col = "BP"
) #14780560 SNPs

# Keep only significant SNPs
CUD_exp_dat_strict_new <- CUD_exp_dat_new_raw %>% filter(pval.exposure <= 5e-8) #only 1274 SNPs left
CUD_exp_dat_linient_new <- CUD_exp_dat_new_raw %>% filter(pval.exposure <= 1e-6) #2762 SNPs left

#####======== Check API KEY ========######
ieugwasr::get_opengwas_jwt()

ieugwasr::user()

#####======== Clumping ========#####
#options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

CUD_exp_dat_clumped <- clump_data(CUD_exp_dat_strict_new) #21 SNPs left

#####======== Import Outcome Data and match it with Exposure SNPs ========#####
NonCog_EA_out_dat <- read_outcome_data(
  snps = CUD_exp_dat_clumped$SNP,
  filename = "NonCog_EA.csv",
  sep = ",",
  snp_col = "variant_id",
  beta_col = "est",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "MAF",
  pval_col = "p_value",
  chr_col = "chromosome",
  pos_col = "base_pair_location"
)#21 SNPs

#####======== Perform MR ========#####
# Harmonization
dat_NonCog <- harmonise_data(
  exposure_dat = CUD_exp_dat_clumped,
  outcome_dat = NonCog_EA_out_dat
)

# test outliers by using ivw_radial
outlier_NonCog <- ivw_radial(dat_NonCog)

# remove outliers
outlier_NonCog <- outlier_NonCog$outliers
dat_NonCog_ready <- dat_NonCog[which(!(dat_NonCog$SNP %in% outlier_NonCog$SNP)),]#13 SNPs

# Run MR
res_CUD_NonCog <- mr(dat_NonCog_ready, method_list = c("mr_ivw_mre", "mr_raps", 
                                                   "mr_weighted_median","mr_weighted_mode", "mr_egger_regression"))
res_CUD_NonCog

#####======== Visualization ========#####
setwd('~/EA_Cannabis_NEW/Plots')

# Scatter Plot
Scatter_CUD_on_NonCog <- mr_scatter_plot(res_CUD_NonCog, dat_NonCog_ready)
Scatter_CUD_on_NonCog[[1]]

# Forest Plot
res_single <- mr_singlesnp(dat_NonCog_ready)
Forest_CUD_on_EA <- mr_forest_plot(res_single)
Forest_CUD_on_EA[[1]]

# Funnel Plot
res_single <- mr_singlesnp(dat_NonCog_ready)
Funnel_CUD_on_EA <- mr_funnel_plot(res_single)
Funnel_CUD_on_EA[[1]]
