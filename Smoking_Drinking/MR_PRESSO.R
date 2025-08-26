if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")

library(MRPRESSO)

# Run MR-PRESSO global method
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_sd_ready, NbDistribution = 1000,  SignifThreshold = 0.05)

mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_ds_ready, NbDistribution = 1000,  SignifThreshold = 0.05)
