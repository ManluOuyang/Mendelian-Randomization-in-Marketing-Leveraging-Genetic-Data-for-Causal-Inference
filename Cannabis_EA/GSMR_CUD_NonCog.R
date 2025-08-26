#####======== Estimate the LD correlation matrix ========#####
library("gsmr2")

#####======== EA -> Cannabis ========#####

#####======== Prepare Data ========#####
setwd('~/EA_Cannabis_NEW/GWAS')

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

#####======== Import Outcome Data and match it with Exposure SNPs ========#####
NonCog_out_dat_gsmr <- read_outcome_data(
  snps = CUD_exp_dat_strict_new$SNP,
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
)#1134 SNPs

# Harmonization
dat_gsmr_NonCog <- harmonise_data(
  exposure_dat = CUD_exp_dat_strict_new,
  outcome_dat = NonCog_out_dat_gsmr
)#135 SNPs dropped for being palindromic with intermediate allele frequencies

#####======== LD Matrix ========#####
setwd('~/EA_Cannabis_NEW/GSMR')

# Save the genetic variants and effect alleles in a text file using R 
write.table(dat_gsmr_NonCog[,c(1,2)], "gsmr_example_snps_NonCog.allele", col.names=F, row.names=F, quote=F) 

# Extract the genotype data from a GWAS dataset using GCTA 
# Prepare own reference data from 1000 genomes
system("~/tools/gcta/gcta-1.94.4-linux-kernel-3-x86_64/gcta64 --bfile EUR --extract gsmr_example_snps_NonCog.allele --update-ref-allele gsmr_example_snps_NonCog.allele --recode --out EUR")

# Estimate LD correlation matrix
snp_coeff_id = scan("EUR.xmat.gz", what="", nlines=1)
snp_coeff = read.table("EUR.xmat.gz", header=F, skip=2)

# Match the SNP genotype data with the summary data
snp_id = Reduce(intersect, list(dat_gsmr_NonCog$SNP, snp_coeff_id))
dat_new = dat_gsmr_NonCog[match(snp_id, dat_gsmr_NonCog$SNP),]
snp_order = match(snp_id, snp_coeff_id)
snp_coeff_id = snp_coeff_id[snp_order]
snp_coeff = snp_coeff[, snp_order]

# Calculate the LD correlation matrix
ldrho = cor(snp_coeff)

# Check the size of the correlation matrix and double-check if the order of the SNPs in the LD correlation matrix is consistent with that in the GWAS summary data
colnames(ldrho) = rownames(ldrho) = snp_coeff_id

dim(ldrho)

# Show the first 5 rows and columns of the matrix  
ldrho[1:5,1:5]


#####======== Perform GSMR ========#####
bzx = dat_new$beta.exposure   # SNP effects on the risk factor
bzx_se = dat_new$se.exposure    # standard errors of bzx
bzx_pval = dat_new$paval.exposure   # p-values for bzx
bzy = dat_new$beta.outcome    # SNP effects on the disease
bzy_se = dat_new$se.outcome    # standard errors of bzy
bzy_pval = dat_new$pval.outcome    # p-values for bzy
n_ref = 503    # Sample size of the reference sample used to calculate the LD matrix
gwas_thresh = 5e-8    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
multi_snps_heidi_thresh = 0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis
nsnps_thresh = 10   # the minimum number of instruments required for the GSMR analysis
heidi_outlier_flag = T    # flag for HEIDI-outlier analysis
ld_r2_thresh = 0.05    # LD r2 threshold to remove SNPs in high LD
ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
gsmr2_beta = 1     # 0 - the original HEIDI-outlier method; 1 - the new global HEIDI-outlier method 

# Run GSMR
gsmr_CUD_NonCog = gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis 
filtered_index=gsmr_CUD_NonCog$used_index
cat("The estimated effect of the exposure on outcome: ",gsmr_CUD_NonCog$bxy)
cat("Standard error of bxy: ",gsmr_CUD_NonCog$bxy_se)
cat("P-value for bxy: ", gsmr_CUD_NonCog$bxy_pval)
cat("Number of pleiotropic outliers: ", length(gsmr_CUD_NonCog$pleio_snps))


#####========= Bidirectional GSMR ========#####
gsmr_CUD_NonCog_bi = bi_gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta) 
cat("Effect of risk factor on disease: ",gsmr_EA_NonCog_bi$forward_bxy)

#####======== Visualization ========#####
effect_col = colors()[75]
vals = c(bzx[filtered_index]-bzx_se[filtered_index], bzx[filtered_index]+bzx_se[filtered_index])
xmin = min(vals); xmax = max(vals)
vals = c(bzy[filtered_index]-bzy_se[filtered_index], bzy[filtered_index]+bzy_se[filtered_index])
ymin = min(vals); ymax = max(vals)
par(mar=c(5,5,4,2))
plot(bzx[filtered_index], bzy[filtered_index], pch=20, cex=0.8, bty="n", cex.axis=1.1, cex.lab=1.2,
     col=effect_col, xlim=c(xmin, xmax), ylim=c(ymin, ymax),
     xlab=expression(Smoking~(italic(b[zx]))),
     ylab=expression(Drinking~(italic(b[zy]))))
abline(0, gsmr_CUD_NonCog$bxy, lwd=1.5, lty=2, col="dim grey")

nsnps = length(bzx[filtered_index])
for( i in 1:nsnps ) {
  # x axis
  xstart = bzx[filtered_index [i]] - bzx_se[filtered_index[i]]; xend = bzx[filtered_index[i]] + bzx_se[filtered_index[i]]
  ystart = bzy[filtered_index[i]]; yend = bzy[filtered_index[i]]
  segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
  # y axis
  xstart = bzx[filtered_index[i]]; xend = bzx[filtered_index[i]] 
  ystart = bzy[filtered_index[i]] - bzy_se[filtered_index[i]]; yend = bzy[filtered_index[i]] + bzy_se[filtered_index[i]]
  segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
}
