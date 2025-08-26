install.packages('remotes')
library(remotes)
install_github("MRCIEU/TwoSampleMR")
remotes::install_github("mrcieu/ieugwasr")

library(TwoSampleMR)
library(dplyr)
library(ieugwasr)

setwd('~/EA_Cannabis/GWAS')

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
Cog_EA_out_dat <- read_outcome_data(
  snps = CUD_exp_dat_clumped$SNP,
  filename = "Cog_EA.csv",
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
) #21 SNPs


#####======== Perform MR ========#####
# Harmonization
dat_Cog <- harmonise_data(
  exposure_dat = CUD_exp_dat_clumped,
  outcome_dat = Cog_EA_out_dat
)

dat_Cog_keep <- subset(dat_Cog, mr_keep)

# test outliers by using ivw_radial
outlier_Cog <- ivw_radial(dat_Cog_keep)

# remove outliers
outlier_Cog <- outlier_Cog$outliers
dat_Cog_ready <- dat_Cog_keep[which(!(dat_Cog_keep$SNP %in% outlier_Cog$SNP)),]

# Run MR
res_CUD_Cog <- mr(dat_Cog_ready, method_list = c("mr_ivw_mre", "mr_raps", 
                                         "mr_weighted_median","mr_weighted_mode", "mr_egger_regression"))
res_CUD_Cog

#####======== Visualization ========#####
setwd('~/EA_Cannabis_NEW/Plots')

# Scatter Plot
Scatter_CUD_on_Cog <- mr_scatter_plot(res_CUD_Cog, dat_Cog_ready)
Scatter_CUD_on_Cog[[1]]

# Forest Plot
res_single <- mr_singlesnp(dat_Cog_ready)
Forest_CUD_on_Cog <- mr_forest_plot(res_single)
Forest_CUD_on_Cog[[1]]

# Funnel Plot
res_single <- mr_singlesnp(dat_Cog_ready)
Funnel_CUD_on_Cog <- mr_funnel_plot(res_single)
Funnel_CUD_on_Cog[[1]]