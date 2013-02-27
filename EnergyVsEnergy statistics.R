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

# t test to see if the drop tower gives lower energy 
energy.dtLower <- t.test(
  gilchrist.energy$DropTowerToMaxInstronForce.J.,
  gilchrist.energy$InstronToMaxForce.J.,
  paired=T,alternative='less') # paired, one tailed
print(energy.dtLower)
