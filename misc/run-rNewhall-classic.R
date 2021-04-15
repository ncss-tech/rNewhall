
##
## first attempt at getting this working
##


library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)

# library(rNewhall)

# name = 'Marena Mesonet Site, OK'
# country = 'USA'
#
# latitude = 36.064166
# longitude = 97.2125
# nsHemisphere = 'N'
# ewHemisphere = 'W'

precipitation <- c(0.21,2.27,4.36,4.36,2.94,1.95,5.64,6.29,3.63,6.29,0.40,0.61) * 25.4
temperature <- (c(35.4,29.6,36.2,42.7,54.9,64.8,62.8,65.2,59.9,45.5,41.6,23.6) - 32) * 5/9

# start.year = 2009
# end.year = 2009

# ## DEB: not used in any of this code
# # but adjusted for above
# is.metric = FALSE
#

## TODO: these are modified as global variables in Newhall.classic
# set moisture calendar and moisture states to ()
moisture.states <- c(0)
moisture.calendar <- c()


# single month
Newhall.classic(
  m = 1,
  water.holding.capacity = 71.3,
  precipitation = precipitation,
  temperature = temperature,
  latitude = 36.064166,
  longitude = 97.2125,
  nsHemisphere = 'N',
  ewHemisphere = 'W'

)

# 12 months
sim <- lapply(
  1:12,
  Newhall.classic,
  water.holding.capacity = 71.3,
  precipitation = precipitation,
  temperature = temperature
)


# 3.1 Spin up model to create initial mositure conditions


diff = 1
iterations = 0
soil.profile.compare = rep(0, times = 64) # this is 64 (water fills entire profile) * 12 (months); Maximum amount of water possible
while (diff > 0.01) {
  out = lapply(1:12, rNewhall.classic)
  diff = abs(1 - (sum(soil.profile.compare) / sum(out[[12]]@soilprofile)))
  soil.profile.compare = out[[12]]@soilprofile
  iterations = iterations +1
}
print(paste("Spin up completed in", iterations, "iterations", sep = " "))


# 4.1 Run model for 1 year and display results graphically
#soil.profile = rep(1, times = 64) # this is 64 (water fills entire profile) * 12 (months); Maximum amount of water possible
out = lapply(1:12, rNewhall.classic)

# moisture.calendar.matrix = matrix(moisture.calendar, nrow = 12, ncol = 30, byrow = TRUE)
# moisture.calendar.matrix.plot = melt(moisture.calendar.matrix)
# ggplot(moisture.calendar.matrix.plot, aes(x = Var2, y = Var1)) +
#   geom_raster(aes(fill=value)) +
#   scale_y_reverse(breaks=c(1:12),labels = c(1,2,3,4,5,6,7,8,9,10,11,12))


# Mid-month grid plot
grid.arrange(out[[1]]@midmonthplot,out[[1]]@endmonthplot,
             out[[2]]@midmonthplot,out[[2]]@endmonthplot,
             out[[3]]@midmonthplot,out[[3]]@endmonthplot,
             out[[4]]@midmonthplot,out[[4]]@endmonthplot,
             out[[5]]@midmonthplot,out[[5]]@endmonthplot,
             out[[6]]@midmonthplot,out[[6]]@endmonthplot,
             out[[7]]@midmonthplot,out[[7]]@endmonthplot,
             out[[8]]@midmonthplot,out[[8]]@endmonthplot,
             out[[9]]@midmonthplot,out[[9]]@endmonthplot,
             out[[10]]@midmonthplot,out[[10]]@endmonthplot,
             out[[11]]@midmonthplot,out[[11]]@endmonthplot,
             out[[12]]@midmonthplot,out[[12]]@endmonthplot,
             top = textGrob(paste("rNewhall Soil Profile Moisture Content for ", name, sep = ""),
                            gp=gpar(fontsize=16,font=2)))

# Plot Moisture Conditions
grid.arrange(out[[1]]@moistplot,
             out[[2]]@moistplot,
             out[[3]]@moistplot,
             out[[4]]@moistplot,
             out[[5]]@moistplot,
             out[[6]]@moistplot,
             out[[7]]@moistplot,
             out[[8]]@moistplot,
             out[[9]]@moistplot,
             out[[10]]@moistplot,
             out[[11]]@moistplot,
             out[[12]]@moistplot,
             top = textGrob(paste("rNewhall Soil Profile Moisture Conditions for ", name, sep = ""),
                            gp=gpar(fontsize=16,font=2)))


# Plot summer conditions
grid.arrange( out[[6]]@midmonthplot,
              out[[7]]@midmonthplot,
              out[[8]]@midmonthplot,
              out[[9]]@midmonthplot,
              out[[6]]@midmonthplotHP,
              out[[7]]@midmonthplotHP,
              out[[8]]@midmonthplotHP,
              out[[9]]@midmonthplotHP,
              out[[6]]@endmonthplot,
              out[[7]]@endmonthplot,
              out[[8]]@endmonthplot,
              out[[9]]@endmonthplot,
              out[[6]]@moistplot,
              out[[7]]@moistplot,
              out[[8]]@moistplot,
              out[[9]]@moistplot,
              ncol = 4,
              nrow = 4,
              top = textGrob(paste("rNewhall Soil Profile Summer Conditions for ", name, sep = ""),
                             gp=gpar(fontsize=16,font=2)))

# Need to make sure water is filling up correctly with HP

# end of month grid plot
grid.arrange(out[[1]]@endmonthplot,
             out[[2]]@endmonthplot,
             out[[3]]@endmonthplot,
             out[[4]]@endmonthplot,
             out[[5]]@endmonthplot,
             out[[6]]@endmonthplot,
             out[[7]]@endmonthplot,
             out[[8]]@endmonthplot,
             out[[9]]@endmonthplot,
             out[[10]]@endmonthplot,
             out[[11]]@endmonthplot,
             out[[12]]@endmonthplot,
             top = textGrob(paste("rNewhall-Basic Soil Profile Moisture Content for ", name, sep = ""),
                            gp=gpar(fontsize=16,font=2)))


