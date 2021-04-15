library(raster)
library(sp)

r <- stack('e:/gis_data/prism/final_monthly_tavg_800m.tif')

p <- SpatialPoints(cbind(-120.69820, 38.02795), proj4string = CRS('+proj=longlat +datum=WGS84'))

e <- extract(r, p)

dput(as.vector(e))
