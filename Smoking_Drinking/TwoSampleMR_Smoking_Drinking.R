#file.edit("~/.Renviron")
usethis::edit_r_environ()

##Put API key into Renviron file (Key generated from 'https://api.opengwas.io/profile/')(key expired in 14 days)
#openGWAS_jwt=eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJvdXlhbmdtYW5sdUBvdXRsb29rLmNvbSIsImlhdCI6MTc0NDg2MDE4NCwiZXhwIjoxNzQ2MDY5Nzg0fQ.wroNEGYYpta16wiyqZsvjrXE71O6wumEIQWUluKKr4oeqFz8NbJygKVHjKkRJJxKaq1RKLLxFKBVj59Wy_KE_om6zuFt24nhd4Nc7RcDMeehtiMgf4Y2XYAIRtVkOTrKwHGQt2vpnE6SH0xyycHWtPBluTlU8_b7qXSIJw89g-xLAa8zd5Amzp4PdWC1IwbQ0VxOkUh9-meNmhU8VwipT-AECeebI7Ktkg752N748iKaYq8bi9amDrhsWggHm0kTy-zTlG8GYxK2vxMEF_83tMLcQh586mwU2xpx5chIPGYfFRBrTc4ffYfhihjDluPcIDU3O6MJpf5J5QtRzlgs4w

#####======== Restart R ========#####
install.packages('remotes')
library(remotes)
install_github("MRCIEU/TwoSampleMR")
remotes::install_github("mrcieu/ieugwasr")

library(TwoSampleMR)
library(dplyr)
library(ieugwasr)

#####======== Read Data ========#####
setwd('~')
si_exp_dat_raw <- read_exposure_data(
  filename = "MR/SmokingInitiation.csv",
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
) #11792288 SNPs

# Keep only significant SNPs
si_exp_dat_sig <- si_exp_dat_raw %>% filter(pval.exposure <= 5e-8) #7855 SNPs

#####======== Check API KEY ========######
ieugwasr::get_opengwas_jwt()

ieugwasr::user()

#####======== Clumping ========#####
#options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

si_exp_dat_clump <- clump_data(si_exp_dat_sig) #93 SNPs

#####======== Import Outcome Data and match it with Exposure SNPs ========#####
drinking_out_dat <- read_outcome_data(
  snps = si_exp_dat_clump$SNP,
  filename = "MR/DrinksPerWeek.csv",
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
)


#####======== Perform MR ========#####
# Harmonization
dat_sd <- harmonise_data(
  exposure_dat = si_exp_dat_clump,
  outcome_dat = drinking_out_dat
)

dat_smoking_drinking <- subset(dat_sd, mr_keep) #84 SNPs

# test outliers by using ivw_radial
outlier_sd <- ivw_radial(dat_smoking_drinking)

# remove outliers
outlier_sd <- outlier_sd$outliers
dat_sd_ready <- dat_smoking_drinking[which(!(dat_smoking_drinking$SNP %in% outlier_sd$SNP)),]#62 SNPs

# Run MR
res_sd <- mr(dat_sd_ready, method_list = c("mr_ivw_mre", "mr_raps", 
                                        "mr_weighted_median","mr_weighted_mode", "mr_egger_regression"))
res_sd

# Leave one out Analysis
res_loo_sd <- mr_leaveoneout(dat_sd_ready)

#####====== Plot the result =====#####
setwd('~/MR/Plots')
Scatter_smoking_on_drinking <- mr_scatter_plot(res_sd, dat_sd_ready)

Scatter_smoking_on_drinking[[1]]

# Forest Plot
res_single <- mr_singlesnp(dat_sd_ready)
p2 <- mr_forest_plot(res_single)
p2[[1]]

# Funnel Plot
res_single <- mr_singlesnp(dat_sd_ready)
p4 <- mr_funnel_plot(res_single)
p4[[1]]

# Leave-one-out Plot
loo_sd <- mr_leaveoneout_plot(res_loo_sd)
loo_sd[[1]]