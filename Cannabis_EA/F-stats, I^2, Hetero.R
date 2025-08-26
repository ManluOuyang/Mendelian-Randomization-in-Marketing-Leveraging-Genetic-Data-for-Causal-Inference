df <- dat_sd_ready

# Per-SNP F and mean F
df$F_snp <- (df$beta.exposure / df$se.exposure)^2
mean_F <- mean(df$F_snp, na.rm = TRUE)

# I^2_GX (NOME check), orient exposure effects positive
bx <- abs(df$beta.exposure)
sx <- df$se.exposure
w  <- 1 / (sx^2)
bx_bar <- sum(w * bx, na.rm = TRUE) / sum(w, na.rm = TRUE)
Qx <- sum(w * (bx - bx_bar)^2, na.rm = TRUE)
L  <- sum(is.finite(bx))
I2_GX <- if (is.finite(Qx) && Qx > 0) max(0, (Qx - (L - 1)) / Qx) else 0

list(mean_F = mean_F, I2_GX = I2_GX)

# Heterogeneity test
heterogeneity <- mr_heterogeneity(df)
print(heterogeneity)
