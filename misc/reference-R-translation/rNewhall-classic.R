
##
## use this as reference, re-factoring into /R/Newhall.classic.R
##


# Part 2: Newhall Simulation *** Only one month is running at a time in this test phase ***
#m = 3 # month for testing
# Processing steps

######################################################
###### Testing Function to run entire year ###########
######################################################
rNewhall.classic =  function(m){
  #get("soil.profile")
  #print(sum(soil.profile))

  # 2.1 Calculate precip and moisture conditions for both light Precip and heavy precip within the month

  # 2.1.1 Define monthly precipitation
  MP = precipitation[m]

  # 2.1.2 Calculate light Precipitation (LP) LP= MP/2
  LP = MP / 2

  # 2.1.3 potential evapotranspiration
  # Need to calculate potential evapotranspiration (PE) vis the Thornthwaite (1948) equation
  # For each month PE(no adjustment) = 16(10T/I)^a (in mm)
  # This value is then adjusted by a latitude constant (k) for each month  PE(Latitude adjusted) = k PE(no adjustment)
  # Define terms (from Thornthwaite 1948)
  I.m = (temperature/5)^1.514# need to account for negative temperatures
  I.m[is.nan(I.m)] <- 0
  I = sum(I.m) #Note that this value is calculated for the whole year, but then used in each monthly calculation
  a = (6.75*(10^-7)*(I^3)) - (7.71*(10^-5)*(I^2)) + (1.792*(10^-2)*I) + 0.49239
  if (temperature[m] > 0 &  temperature[m] < 26.5){
    raw.PE = 16 * (10 * temperature[m] / I)^a
  } else if (temperature[m] >= 38) { # If temp is over 38C, the PE value is fixed at 185.0mm
    raw.PE = 185.0
  } else if (temperature[m] < 0) { # If temp is over 38C, the PE value is fixed at 185.0mm
    raw.PE = 0
  } else {
    raw.PE =   pe.bins[findInterval(temperature[m], temp.bins)] #Use look up tables (temp.bins and pe.bins to assign a pe based on temperature)
  }
  # Asjust raw.PE using the K correction factor for daylight hours based on latitude
  if (nsHemisphere == 'N'){
    k = knorth[findInterval(latitude,as.numeric(dimnames(knorth)[[1]])),m] # look up k value by latitude and month in appropriate lookup table (n vs S hemis)
  } else {
    k = ksouth[findInterval(latitude,as.numeric(dimnames(ksouth)[[1]])),m]
  }

  PE = raw.PE * k

  # 2.1.4 Calculate effect of PE; called Net moisture activity (NMA) in Newhall and Berdanier (1996) and net potential
  # evapotranspiration (NPE) in Van Wambeke (2000). If NPE > 0, accretion will take place during this period; otherwise,
  # water will be extracted from the profile
  NMA = (LP - PE)


  # 2.1.5 Calculate heavy precipitation (HP). HP = MP/2 All HP enters the profile as accretion.
  HP = MP / 2

  # 2.1.6 Finally, calculate the montly water balance for the profile
  m.water.balance = MP - PE

  # 2.2 Model changes in water content during each month (both accretion and depletion)
  # 2.2.1 First half of the month (15 days) add or remove water from profile based on LP
  # If NPE > 0, apply NPE/2 to fill available slots; if NPE < 0, apply NPE/2 to deplete filled slots
  if (NMA > 0) {
    total.available.NMA = (NMA) / 2
    current.NMA = NMA / 2
    for (i in 1:length(accretion.order)){
      accretion.pos = accretion.order[i]
      water.in.slot = soil.profile[accretion.pos]
      water.to.fill.slot = water.per.slot - water.in.slot
      if (current.NMA > water.to.fill.slot) {
        soil.profile[accretion.pos] = water.in.slot + water.to.fill.slot
        current.NMA =  current.NMA - water.to.fill.slot
      } else {
        soil.profile[accretion.pos] = water.in.slot + current.NMA
        current.NMA =  current.NMA - current.NMA
      }
      if (current.NMA == 0) break # if there is no more water from preip, stop
      if (sum(soil.profile) == water.holding.capacity) break # if the soil profile is full, stop
    }
  }
  if (NMA < 0) {
    total.available.NMA = abs(NMA) / 2
    current.NMA = abs(NMA) / 2
    for (i in 1:length(depletion.order)){
      depletion.pos = depletion.order[i]
      water.in.slot = soil.profile[depletion.pos]
      energy.to.drain.slot = water.in.slot * depletion.req[i]
      if (soil.profile[depletion.pos] == 0) next
      if (current.NMA > energy.to.drain.slot) {
        soil.profile[depletion.pos] = 0
        current.NMA = current.NMA - energy.to.drain.slot
      } else {
        remaining.energy.proportion = current.NMA / energy.to.drain.slot
        soil.profile[depletion.pos] = remaining.energy.proportion * water.in.slot
        current.NMA = current.NMA - current.NMA
      }
      if (current.NMA == 0) break
      if (sum(soil.profile) == 0) break
    }
  }

  # # record mositure condition (add to moisture.states vector)
  if (sum(soil.profile[c(9,17,25)]) == 0){
    current.moisture.condition = 1
  }
  if (soil.profile[9] > 0 & soil.profile[17] > 0 & soil.profile[25] > 0){
    current.moisture.condition  = 3
  }
  if (soil.profile[9] == 0 | soil.profile[17] == 0 | soil.profile[25] == 0 & (sum(soil.profile[c(9,17,25)]) > 0)){
    current.moisture.condition  = 2
  }
  moisture.states <<- append(moisture.states,current.moisture.condition)

  #Duration of moisture state
  # m.diff = moisture.states[(length(moisture.states))] - moisture.states[(length(moisture.states)-1)]
  # if (m.diff == 0){
  #   moisture.calendar <<- append(moisture.calendar,rep(moisture.states[(length(moisture.states))], 15))
  # }
  # if (m.diff == -1){
  #   d.initial = round(15 * (NMA - total.available.NMA) / NMA)
  #   d.end = 15 - d.initial
  #   moisture.calendar <<- append(moisture.calendar,c(rep(moisture.states[(length(moisture.states)-1)], d.initial),
  #                                rep(moisture.states[(length(moisture.states))], d.end)))
  # }

  # 2.2.2 Create mid month graphs (for testing only)
  soil.profile.matrix = matrix(soil.profile, nrow = soil.profile.dims, ncol = soil.profile.dims, byrow = T)
  soil.profile.matrix.plot = melt(soil.profile.matrix)
  soil.profile.matrix.plot$value = soil.profile.matrix.plot$value/water.per.slot

  p1 = ggplot(soil.profile.matrix.plot, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    labs(x="", y="",
         title=paste("Month ", m, sep = "")) +
    scale_y_reverse(breaks=c(1:8),labels = c("Surface", "Moisture Control", "Moisture Control","","","","","Bottom"))+
    scale_x_continuous(breaks=c(1:8),labels = c("PWP", "", "","","","","","AWC")) +
    scale_fill_gradient(limits = c(0,1)) +#, colours=c("#132B43","#56B1F7" ))+
    theme_bw() + theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=10),
                       plot.title=element_text(size=14),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())+
    geom_hline(yintercept=.5) + geom_hline(yintercept=1.5) + geom_hline(yintercept=2.5) + geom_hline(yintercept=3.5) + geom_hline(yintercept=4.5) +
    geom_hline(yintercept=5.5) + geom_hline(yintercept=6.5) + geom_hline(yintercept=7.5) + geom_hline(yintercept=8.5) + geom_vline(xintercept=.5) +
    geom_vline(xintercept=1.5) + geom_vline(xintercept=2.5) + geom_vline(xintercept=3.5) + geom_vline(xintercept=4.5) + geom_vline(xintercept=5.5) + geom_vline(xintercept=6.5) +
    geom_vline(xintercept=7.5) + geom_vline(xintercept=8.5) +labs(fill = "% Filled")

  # 2.2.3 Heavy rain at the mid-point of the month all water from HP is added to the profile in available slots
  #print ("Big storm, filling up")
  available.HP = HP
  for (i in 1:length(accretion.order)){
    accretion.pos = accretion.order[i]
    water.in.slot = soil.profile[accretion.pos]
    water.to.fill.slot = water.per.slot - water.in.slot
    if (available.HP > water.to.fill.slot) {
      soil.profile[accretion.pos] = water.in.slot + water.to.fill.slot
      available.HP =  available.HP - water.to.fill.slot
    } else {
      soil.profile[accretion.pos] = water.in.slot + available.HP
      available.HP =  available.HP - available.HP
    }
    if (available.HP == 0) break # if there is no more water from preip, stop
    if (sum(soil.profile) == water.holding.capacity) break # if the soil profile is full, stop
  }
  # # record mositure condition 2
  if (sum(soil.profile[c(9,17,25)]) == 0){
    current.moisture.condition = 1
  }
  if (soil.profile[9] > 0 & soil.profile[17] > 0 & soil.profile[25] > 0){
    current.moisture.condition  = 3
  }
  if (soil.profile[9] == 0 | soil.profile[17] == 0 | soil.profile[25] == 0 & (sum(soil.profile[c(9,17,25)]) > 0)){
    current.moisture.condition  = 2
  }
  moisture.states <<- append(moisture.states,current.moisture.condition)

  soil.profile.matrix = matrix(soil.profile, nrow = soil.profile.dims, ncol = soil.profile.dims, byrow = T)
  soil.profile.matrix.plot = melt(soil.profile.matrix)
  soil.profile.matrix.plot$value = soil.profile.matrix.plot$value/water.per.slot


  p2 = ggplot(soil.profile.matrix.plot, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    labs(x="", y="",
         title=paste("Month ", m, sep = "")) +
    scale_y_reverse(breaks=c(1:8),labels = c("Surface", "Moisture Control", "Moisture Control","","","","","Bottom"))+
    scale_x_continuous(breaks=c(1:8),labels = c("PWP", "", "","","","","","AWC")) +
    scale_fill_gradient(limits = c(0,1))+#, colours=c("#132B43","#56B1F7" ))+
    theme_bw() + theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=10),
                       plot.title=element_text(size=14),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())+
    geom_hline(yintercept=.5) + geom_hline(yintercept=1.5) + geom_hline(yintercept=2.5) + geom_hline(yintercept=3.5) + geom_hline(yintercept=4.5) +
    geom_hline(yintercept=5.5) + geom_hline(yintercept=6.5) + geom_hline(yintercept=7.5) + geom_hline(yintercept=8.5) + geom_vline(xintercept=.5) +
    geom_vline(xintercept=1.5) + geom_vline(xintercept=2.5) + geom_vline(xintercept=3.5) + geom_vline(xintercept=4.5) + geom_vline(xintercept=5.5) + geom_vline(xintercept=6.5) +
    geom_vline(xintercept=7.5) + geom_vline(xintercept=8.5)+labs(fill = "% Filled")

  # 2.2.4 Second half of the month (15 days) add or remove water from profile based on LP
  # If NPE > 0, apply NPE/2 to fill available slots; if NPE < 0, apply NPE/2 to deplete filled slots
  if (NMA > 0) {
    total.available.NMA = (NMA) / 2
    current.NMA = NMA / 2
    for (i in 1:length(accretion.order)){
      accretion.pos = accretion.order[i]
      water.in.slot = soil.profile[accretion.pos]
      water.to.fill.slot = water.per.slot - water.in.slot
      if (current.NMA > water.to.fill.slot) {
        soil.profile[accretion.pos] = water.in.slot + water.to.fill.slot
        current.NMA =  current.NMA - water.to.fill.slot
      } else {
        soil.profile[accretion.pos] = water.in.slot + current.NMA
        current.NMA =  current.NMA - current.NMA
      }
      if (current.NMA == 0) break # if there is no more water from preip, stop
      if (sum(soil.profile) == water.holding.capacity) break # if the soil profile is full, stop
    }
  }
  if (NMA < 0) {
    total.available.NMA = abs(NMA) / 2
    current.NMA = abs(NMA) / 2
    for (i in 1:length(depletion.order)){
      depletion.pos = depletion.order[i]
      water.in.slot = soil.profile[depletion.pos]
      energy.to.drain.slot = water.in.slot * depletion.req[i]
      if (soil.profile[depletion.pos] == 0) next
      if (current.NMA > energy.to.drain.slot) {
        soil.profile[depletion.pos] = 0
        current.NMA = current.NMA - energy.to.drain.slot
      } else {
        remaining.energy.proportion = current.NMA / energy.to.drain.slot
        soil.profile[depletion.pos] = remaining.energy.proportion * water.in.slot
        current.NMA = current.NMA - current.NMA
      }
      if (current.NMA == 0) break
      if (sum(soil.profile) == 0) break
    }
  }
  # # record mositure condition 3
  if (sum(soil.profile[c(9,17,25)]) == 0){
    current.moisture.condition = 1
  }
  if (soil.profile[9] > 0 & soil.profile[17] > 0 & soil.profile[25] > 0){
    current.moisture.condition  = 3
  }
  if (soil.profile[9] == 0 | soil.profile[17] == 0 | soil.profile[25] == 0 & (sum(soil.profile[c(9,17,25)]) > 0)){
    current.moisture.condition  = 2
  }
  moisture.states <<- append(moisture.states,current.moisture.condition)


  #  m.diff = moisture.states[(length(moisture.states))] - moisture.states[(length(moisture.states)-1)]
  # if (m.diff == 0){
  #   moisture.calendar <<- append(moisture.calendar,rep(moisture.states[(length(moisture.states))], 15))
  # }
  # if (m.diff == -1){
  #   d.initial = round(15 * (NMA - total.available.NMA) / NMA)
  #   d.end = 15 - d.initial
  #   moisture.calendar <<- append(moisture.calendar,c(rep(moisture.states[(length(moisture.states)-1)], d.initial),
  #                                                    rep(moisture.states[(length(moisture.states))], d.end)))
  # }

  # 2.2.5 Calculate moisture condition duration
  # Steps needed:
  # 1) Determine how much of a difference there is  soil moisture between LP periods:
  #    ex. states 1 to 2, 2 to 3, or 1 to 3. Can be coded as 0, (-)1 or (-)2
  # 2) calcultae duration of initial state (15 * PE.trans / NPE )
  # 2a) calculate intermediate stage (if coded as (-)2)
  # 3) calculate duration of final moisture state (15 - duration of initial state )
  # 4) repeat for each 2 week period

  #moisture.conditions = c(moisture.condition.1, moisture.condition.2, moisture.condition.3)

  # 2.2.6 Create end month graphs (for testing only)
  soil.profile.matrix = matrix(soil.profile, nrow = soil.profile.dims, ncol = soil.profile.dims, byrow = T)
  soil.profile.matrix.plot = melt(soil.profile.matrix)
  soil.profile.matrix.plot$value = soil.profile.matrix.plot$value/water.per.slot
  p3 = ggplot(soil.profile.matrix.plot, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    labs(x="", y="",
         title=paste("Month ", m, sep = "")) +
    scale_y_reverse(breaks=c(1:8),labels = c("Surface", "Moisture Control", "Moisture Control","","","","","Bottom"))+
    scale_x_continuous(breaks=c(1:8),labels = c("PWP", "", "","","","","","AWC")) +
    scale_fill_gradient(limits = c(0,1))+#, colours=c("#132B43","#56B1F7" ))+
    theme_bw() + theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=10),
                       plot.title=element_text(size=14),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())+
    geom_hline(yintercept=.5) + geom_hline(yintercept=1.5) + geom_hline(yintercept=2.5) + geom_hline(yintercept=3.5) + geom_hline(yintercept=4.5) +
    geom_hline(yintercept=5.5) + geom_hline(yintercept=6.5) + geom_hline(yintercept=7.5) + geom_hline(yintercept=8.5) + geom_vline(xintercept=.5) +
    geom_vline(xintercept=1.5) + geom_vline(xintercept=2.5) + geom_vline(xintercept=3.5) + geom_vline(xintercept=4.5) + geom_vline(xintercept=5.5) + geom_vline(xintercept=6.5) +
    geom_vline(xintercept=7.5) + geom_vline(xintercept=8.5) +labs(fill = "% Filled")
  #
  # # 2.2.7 Create soil moisture condition plot
  moisture.conditions.dataframe = data.frame(Period=c(1, 2, 3), Condition = c(moisture.states[(length(moisture.states)-2)], moisture.states[(length(moisture.states)-1)], moisture.states[length(moisture.states)]))
  p4 = ggplot(moisture.conditions.dataframe, aes(Period,Condition)) +
    geom_col(fill = "#3182BD") +
    scale_y_discrete(limits = c(1,2,3),labels = c("Dry", "Moist/Dry", "Moist"))+
    scale_x_discrete(limits = c(1,2, 3),labels = c("Mid-Month 1","Mid-Month 2", "End of Month")) +
    labs(x="", y="",
         title=paste("Month ", m, sep = "")) +
    theme_bw() + theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=10),
                       plot.title=element_text(size=14))


  # 2.2.8 Create output
  soil.profile <<- soil.profile # pull out this value so it can be used in the next iteration
  #print(sum(soil.profile))
  return(new("newhalloutput", soilprofile = soil.profile, moiststates=moisture.states, midmonthplot = p1, midmonthplotHP = p2, endmonthplot = p3, moistplot = p4))
}



