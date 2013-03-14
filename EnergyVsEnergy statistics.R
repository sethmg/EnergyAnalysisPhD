# a script to perform the energy difference statistics

gilchrist.energy <- read.table(
  '/home/seth/Documents/Research/Projects/Open Projects/12-035 Energy Analysis/EnergyComparisons.txt',
  header=T,row.names=1)

# set values of zero to Not Available and omit those entries for the stats
gilchrist.energy[gilchrist.energy == 0] <- NA
invisible( na.omit(gilchrist.energy) ) # output suppressed

# t test to see if the two are different
energy.different <- t.test(
  gilchrist.energy$DropTowerToMaxInstronForce.J.,
  gilchrist.energy$InstronToMaxForce.J.,
  paired=T) # paried, two tailed
print(energy.different)

# wilcoxon test to see if the two are different
energy.different_wilcox <- wilcox.test(
  gilchrist.energy$DropTowerToMaxInstronForce.J.,
  gilchrist.energy$InstronToMaxForce.J.,
  paired=T)
print(energy.different_wilcox)

# t test to see if the drop tower gives lower energy 
energy.dtLower <- t.test(
  gilchrist.energy$DropTowerToMaxInstronForce.J.,
  gilchrist.energy$InstronToMaxForce.J.,
  paired=T,alternative='less') # paired, one tailed
print(energy.dtLower)

# wilcoxon test to see if the dt gives lower energy
energy.dtLower_wilcox <- wilcox.test(
  gilchrist.energy$DropTowerToMaxInstronForce.J.,
  gilchrist.energy$InstronToMaxForce.J.,
  paired=T,alternative='less')
print(energy.dtLower_wilcox)

