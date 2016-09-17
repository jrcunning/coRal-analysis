# Conversion values used in parameter value calculations
CmolR.m2 <- 3.8  # Edmunds et al. 2011
cm2.m2 <- 10000
mg.g <- 1000
h.d <- 24
s.h <- 3600
min.d <- 1440
cm2.m2 <- 10000
RuBP.CmolS <- 2.31e-7  # Brading et al. 2013
nmol.mol <- 1e9
mmol.mol <- 1e3
µmol.mol <- 1e6
CmolS.µm3 <- 2.29144e-14  # Brading et al. 2011
mm3.µm3 <- 1e-9
gchlA.CmolS <- 0.111  # Brading et al. 2013
gC.molC <- 12.011
mg.g <- 1000
µg.g <- 1e6
pg.g <- 1e12
gchlA.gC <- 0.01020408  # Brading et al. 2013
gchlA.molchlA <- 893.5
CmolR.m2_Spis <- 0.6  # Muscatine et al. 1989 in Edmunds et al. 2011 -- also 0.5 in Grover 2002
nmolchlA.cm2 <- 9.66  # Edmunds et al. 2011
µgchl.cm2_Spis <- 22  # Godinot et al. 2011 PLoS ONE
CmolR.m2_Gfas <- 2  # Schutter et al. 2010 [based on reported 6 mg afdw cm-2 (Fig. 9) and 41.4% Carbon by weight (Results section)]
gC.artemia <- 0.678e-6  # Wijgerde et al. 2011 JEB
pgC.cellS <- 234.5 # Average of two clade A phylotypes in Brading et al. 2013
#pgC.cellS <- 100 # Average of two types from Cassiopeia in Verde & McCloskey 1998
CmolS.cellS <- pgC.cellS / pg.g / gC.molC
