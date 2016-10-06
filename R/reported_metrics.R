# Metrics cited in this manuscript derived from literature

# Specific growth rates
# =====================
# From Tanaka et al. 2007 Limnol Oceanogr: Table 1
C0 <- 343  # tissue carbon at time zero (µmolC/cm2)
C5 <- 495  # tissue carbon after 5 days of nutrient enrichment (µmolC/cm2)
Tanaka_gr <- log(C5/C0) / 5  # specific growth of tissue carbon (d^-1)
Tanaka_gr <- round(Tanaka_gr, 2)

# From Armoza-Zvuloni et al. 2014
Armoza_gr <- 0.2938  # Exponential growth rate reported on page 3 (wk^-1)
Armoza_gr <- 0.2938 / 7  # Convert to d^-1
Armoza_gr <- round(Armoza_gr, 2)

# From Shafir et al. 20016
Shafir_gr <- 0.0167  # d^-1; average of 5 species in a nursery setting, based on approximate volumes
