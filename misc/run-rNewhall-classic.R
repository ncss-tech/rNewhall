
##
## first attempt at getting this working
##
library(soilDB)
library(aqp)
library(reshape2)
# library(rNewhall)

library(ggplot2)
library(grid)
library(gridExtra)


# ## DEB: not used in any of this code
# # but adjusted for above
# is.metric = FALSE

# name = 'Marena Mesonet Site, OK'
# latitude = 36.064166
# longitude = 97.2125

# in -> mm
precipitation <- c(0.21,2.27,4.36,4.36,2.94,1.95,5.64,6.29,3.63,6.29,0.40,0.61) * 25.4
# deg F -> deg C
temperature <- (c(35.4,29.6,36.2,42.7,54.9,64.8,62.8,65.2,59.9,45.5,41.6,23.6) - 32) * 5/9

# total water storage in mm
AWC <- 75


###

## near Sonora, CA
x <- fetchOSD('amador', extended = TRUE)
precipitation <- x$climate.monthly$q50[x$climate.monthly$variable == 'Precipitation (mm)']
temperature <- c(8, 10, 11, 14, 18, 22, 25, 25, 22, 18, 12, 8)
AWC <- 75


## must setup classic.env first
#
# # single month
# .Newhall.classic.month(
#   m = 1,
#   AWC = 71.3,
#   PPT = precipitation,
#   TAVG = temperature,
#   latitude = 36.064166,
#   longitude = 97.2125,
#   nsHemisphere = 'N',
#   ewHemisphere = 'W'
# )

# AWC units are mm

sim <- Newhall.classic(
  AWC = AWC,
  PPT = precipitation,
  TAVG = temperature,
  latitude = 36.064166,
  longitude = 97.2125,
  nsHemisphere = 'N',
  ewHemisphere = 'W'
)

str(sim[[1]], 1)

##
mcd <- do.call(
  'rbind',
  lapply(sim, '[[', 'moisture.conditions.dataframe')
)

mcd$Month <- factor(mcd$Month)
mcd$Period <- factor(mcd$Period)
mcd$Condition <- factor(mcd$Condition, levels = 1:3)
table(mcd$Month, mcd$Condition)


##
spm <- do.call(
  'rbind',
  lapply(sim, '[[', 'soil.profile.matrix.plot')
)


## ??
ms <- sapply(sim, '[[', 'moiststates')
ms

## TODO: check soil moisture
m <- sapply(sim, '[[', 'soilprofile')
colSums(m) / AWC


## p4
ggplot(mcd, aes(Period,Condition)) +
  facet_wrap(mcd$Month) +
  geom_col(fill = "#3182BD") +
  scale_y_discrete(limits = c(1,2,3), labels = c("Dry", "Moist/Dry", "Moist"))+
  scale_x_discrete(limits = c(1,2, 3), labels = c("MM 1","MM 2", "End")) +
  labs(x="", y="") +
  theme_bw() + theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=10),
                     plot.title=element_text(size=14))


## p3
ggplot(spm, aes(x = Var2, y = Var1)) +
  facet_wrap(spm$Month) +
  geom_raster(aes(fill=value)) +
  labs(x="", y="") +
  scale_y_reverse(breaks=c(1:8),labels = c("Surface", "Moisture Control", "Moisture Control","","","","","Bottom"))+
  scale_x_continuous(breaks=c(1:8),labels = c("PWP", "", "","","","","","AWC")) +
  scale_fill_gradient(limits = c(0,1))+#, colours=c("#132B43","#56B1F7" ))+
  theme_bw() + theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=10),
                     plot.title=element_text(size=14),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
  labs(fill = "% Filled")

  # geom_hline(yintercept=.5) + geom_hline(yintercept=1.5) + geom_hline(yintercept=2.5) + geom_hline(yintercept=3.5) + geom_hline(yintercept=4.5) +
  # geom_hline(yintercept=5.5) + geom_hline(yintercept=6.5) + geom_hline(yintercept=7.5) + geom_hline(yintercept=8.5) + geom_vline(xintercept=.5) +
  # geom_vline(xintercept=1.5) + geom_vline(xintercept=2.5) + geom_vline(xintercept=3.5) + geom_vline(xintercept=4.5) + geom_vline(xintercept=5.5) + geom_vline(xintercept=6.5) +
  # geom_vline(xintercept=7.5) + geom_vline(xintercept=8.5)






##
## old stuff
##




# 4.1 Run model for 1 year and display results graphically
#soil.profile = rep(1, times = 64) # this is 64 (water fills entire profile) * 12 (months); Maximum amount of water possible
# out = lapply(1:12, rNewhall.classic)

## DEB: currently commented-out in Newhall.classic()
# moisture.calendar.matrix = matrix(moisture.calendar, nrow = 12, ncol = 30, byrow = TRUE)
# moisture.calendar.matrix.plot = melt(moisture.calendar.matrix)
# ggplot(moisture.calendar.matrix.plot, aes(x = Var2, y = Var1)) +
#   geom_raster(aes(fill=value)) +
#   scale_y_reverse(breaks=c(1:12),labels = c(1,2,3,4,5,6,7,8,9,10,11,12))

#
# # Mid-month grid plot
# grid.arrange(out[[1]]@midmonthplot,out[[1]]@endmonthplot,
#              out[[2]]@midmonthplot,out[[2]]@endmonthplot,
#              out[[3]]@midmonthplot,out[[3]]@endmonthplot,
#              out[[4]]@midmonthplot,out[[4]]@endmonthplot,
#              out[[5]]@midmonthplot,out[[5]]@endmonthplot,
#              out[[6]]@midmonthplot,out[[6]]@endmonthplot,
#              out[[7]]@midmonthplot,out[[7]]@endmonthplot,
#              out[[8]]@midmonthplot,out[[8]]@endmonthplot,
#              out[[9]]@midmonthplot,out[[9]]@endmonthplot,
#              out[[10]]@midmonthplot,out[[10]]@endmonthplot,
#              out[[11]]@midmonthplot,out[[11]]@endmonthplot,
#              out[[12]]@midmonthplot,out[[12]]@endmonthplot,
#              top = textGrob(paste("rNewhall Soil Profile Moisture Content for ", name, sep = ""),
#                             gp=gpar(fontsize=16,font=2)))
#
# # Plot Moisture Conditions
# grid.arrange(out[[1]]@moistplot,
#              out[[2]]@moistplot,
#              out[[3]]@moistplot,
#              out[[4]]@moistplot,
#              out[[5]]@moistplot,
#              out[[6]]@moistplot,
#              out[[7]]@moistplot,
#              out[[8]]@moistplot,
#              out[[9]]@moistplot,
#              out[[10]]@moistplot,
#              out[[11]]@moistplot,
#              out[[12]]@moistplot,
#              top = textGrob(paste("rNewhall Soil Profile Moisture Conditions for ", name, sep = ""),
#                             gp=gpar(fontsize=16,font=2)))
#
#
# # Plot summer conditions
# grid.arrange( out[[6]]@midmonthplot,
#               out[[7]]@midmonthplot,
#               out[[8]]@midmonthplot,
#               out[[9]]@midmonthplot,
#               out[[6]]@midmonthplotHP,
#               out[[7]]@midmonthplotHP,
#               out[[8]]@midmonthplotHP,
#               out[[9]]@midmonthplotHP,
#               out[[6]]@endmonthplot,
#               out[[7]]@endmonthplot,
#               out[[8]]@endmonthplot,
#               out[[9]]@endmonthplot,
#               out[[6]]@moistplot,
#               out[[7]]@moistplot,
#               out[[8]]@moistplot,
#               out[[9]]@moistplot,
#               ncol = 4,
#               nrow = 4,
#               top = textGrob(paste("rNewhall Soil Profile Summer Conditions for ", name, sep = ""),
#                              gp=gpar(fontsize=16,font=2)))
#
# # Need to make sure water is filling up correctly with HP
#
# # end of month grid plot
# grid.arrange(out[[1]]@endmonthplot,
#              out[[2]]@endmonthplot,
#              out[[3]]@endmonthplot,
#              out[[4]]@endmonthplot,
#              out[[5]]@endmonthplot,
#              out[[6]]@endmonthplot,
#              out[[7]]@endmonthplot,
#              out[[8]]@endmonthplot,
#              out[[9]]@endmonthplot,
#              out[[10]]@endmonthplot,
#              out[[11]]@endmonthplot,
#              out[[12]]@endmonthplot,
#              top = textGrob(paste("rNewhall-Basic Soil Profile Moisture Content for ", name, sep = ""),
#                             gp=gpar(fontsize=16,font=2)))


