install.packages('remotes')
library(remotes)
install_github("MRCIEU/TwoSampleMR")
remotes::install_github("mrcieu/ieugwasr")

library(TwoSampleMR)
library(dplyr)
library(ieugwasr)
library(RadialMR)

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
EA_out_dat_new <- read_outcome_data(
  snps = CUD_exp_dat_clumped$SNP,
  filename = "EA_okbay.csv",
  sep = ",",
  snp_col = "rsID",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "Effect_allele",
  other_allele_col = "Other_allele",
  eaf_col = "EAF_HRC",
  pval_col = "P",
  chr_col = "Chr",
  pos_col = "BP"
)

#####======== Perform MR ========#####
# Harmonization
dat_EA_Cannabis <- harmonise_data(
  exposure_dat = CUD_exp_dat_clumped,
  outcome_dat = EA_out_dat_new
)

# test outliers by using ivw_radial
outlier <- ivw_radial(dat_EA_Cannabis)

# remove outliers
outliers <- outlier$outliers
dat_EA_Cannabis_ready <- dat_EA_Cannabis[which(!(dat_EA_Cannabis$SNP %in% outliers$SNP)),] #12 SNPs left

# Run MR
res_CUD_EA <- mr(dat_EA_Cannabis_ready, method_list = c("mr_ivw_mre", "mr_raps", 
                                        "mr_weighted_median","mr_weighted_mode", "mr_egger_regression"))
res_CUD_EA

#####======== Visualization ========#####
setwd('~/EA_Cannabis_NEW/Plots')

# Scatter Plot
Scatter_CUD_on_EA <- mr_scatter_plot(res_CUD_EA, dat_EA_Cannabis_ready)
Scatter_CUD_on_EA[[1]]

# Forest Plot
res_single <- mr_singlesnp(dat_EA_Cannabis_ready)
Forest_CUD_on_EA <- mr_forest_plot(res_single)
Forest_CUD_on_EA[[1]]

# Funnel Plot
res_single <- mr_singlesnp(dat_EA_Cannabis_ready)
Funnel_CUD_on_EA <- mr_funnel_plot(res_single)
Funnel_CUD_on_EA[[1]]
