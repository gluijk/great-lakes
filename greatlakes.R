# Maps connecting the Great Lakes to the ocean
# www.overfitting.net
# https://www.overfitting.net/2024/01/conectando-los-grandes-lagos-al-mar-con.html


library(terra)  # read GeoTIFF, reprojection, crop and resample
library(tiff)  # save 16-bit TIFF's
library(png)  # save 8-bit PNG's


hillshademap=function(DEM, dx=25, dlight=c(0, 2, 3), gamma=1) {
    # hillshademap() inputs DEM data and outputs a hillshade matrix
    #
    # DEM: digital elevation map matrix
    # dx: DEM resolution (cell size in the same units as elevation values)
    # dlight: lighting direction (3D vector defined from observer to light source):
    #   dlight=c(0, 2, 3)  # sunrise
    #   dlight=c(0, 0, 1)  # midday
    #   dlight=c(0,-2, 3)  # sunset
    # gamma: optional output gamma lift

    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    # If array turn its first dimension into a matrix
    if (!is.matrix(DEM)) DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)

    dlightM=sum(dlight^2)^0.5
    
    # Vectorial product to calculate n (orthogonal vector)
    nx = 2*dx*(DEM[1:(DIMY-2), 2:(DIMX-1)] - DEM[3:DIMY,     2:(DIMX-1)])
    ny = 2*dx*(DEM[2:(DIMY-1), 1:(DIMX-2)] - DEM[2:(DIMY-1), 3:DIMX])
    nz = 4*dx^2
    nM = (nx^2 + ny^2 + nz^2)^0.5
    
    # Dot product to calculate cos(theta)
    dn = dlight[1]*nx + dlight[2]*ny + dlight[3]*nz  # (DIMY-2)x(DIMX-2) matrix
    
    # Reflectance (=cos(theta))
    hillshadepre=dn/(dlightM*nM)
    hillshadepre[hillshadepre<0]=0  # clip negative values
    
    # Add 1-pix 'lost' borders
    hillshademap=matrix(0, nrow=DIMY, ncol=DIMX)
    hillshademap[2:(DIMY-1), 2:(DIMX-1)]=hillshadepre
    rm(hillshadepre)
    hillshademap[c(1,DIMY),]=hillshademap[c(2,DIMY-1),]
    hillshademap[,c(1,DIMX)]=hillshademap[,c(2,DIMX-1)]
    
    return(hillshademap^(1/gamma))
}


###########################################################

# 1. PROCESS GEOTIFF DATA TO GET THE DEM AS A MATRIX

# https://download.gebco.net/
# The GEBCO_2023 Grid is a global terrain model for ocean and land,
# providing elevation data, in meters, on a 15 arc-second interval grid
# of 43200 rows x 86400 columns, giving 3,732,480,000 data points.
# The data values are pixel-centre registered i.e. they refer to elevations,
# in meters, at the centre of grid cells.
lakes=rast("geotiff_northamerica.tif")
lakes
plot(lakes)

# CROP Great Lakes (in long/lant degrees)
cropdef=ext(-100, -70, 38, 52)
lakes=crop(x=lakes, y=cropdef, threads=TRUE)
lakes
plot(lakes)

# REPROJECT raster from Longitude Latitude (+proj=longlat)/WGS84
# to Lambert Conic Conformal (+proj=lcc)/WGS84
# https://pygis.io/docs/d_understand_crs_codes.html
# https://stackoverflow.com/questions/36868506/how-to-change-a-lambert-conic-conformal-raster-projection-to-latlon-degree-r
CRS="+proj=lcc +ellps=WGS84 +lat_1=33 +lat_2=45 +lon_0=-84 +units=km"
lakes=project(x=lakes, y=CRS, method='bilinear', threads=TRUE)
lakes
plot(lakes)


# CROP raster to area of interest (in kms)
cropdef=ext(-675, 650, 4800, 5800)
lakes=crop(x=lakes, y=cropdef, threads=TRUE)
lakes
plot(lakes)
# resolution  : 0.3538445, 0.3538445  (x, y) in kms
# 0.3538445km grid resolution
resolution=res(lakes)[1]*1000  # 353.8445m grid resolution


# Convert to matrix and save as TIFF
DEM=matrix(as.array(lakes), nrow=nrow(lakes))
hist(DEM, breaks=1000)
abline(v=0, col='red')
DEM=DEM-min(DEM)
writeTIFF(DEM/max(DEM), "lakes.tif", compression='LZW', bits.per.sample=16)

DEM=matrix(as.array(lakes), nrow=nrow(lakes))
DIMY=nrow(DEM)
DIMX=ncol(DEM)


###########################################################

# 2. PROCESS MATRIX TO OBTAIN MAP CONTOURS AND HILLSHADE

# Calculate solid map contour
solid=DEM

# Sea level contours
solid[solid>=0]=1  # set >=0 areas to 1 (land)
solid[solid<0]=0  # set <0 areas to 0 (water)

# Superior, Michigan, Huron, Erie contours: ~179m above sea level
solid[solid<179]=0  # set <179 areas to 0 (water)
solid[solid>=179]=1  # set >=179 areas to 1 (land)

# Ontario contours: 74m above sea level
solid[solid<74]=0  # set <74 areas to 0 (water)
solid[solid>=74]=1  # set >=74 areas to 1 (land)

writePNG(solid, "mapsolid.png")


# Calculate outline map from solid map
outline=solid*0
# 1 pixel thickness outline
outline[2:(DIMY-1), 2:(DIMX-1)]=
    abs(solid[1:(DIMY-2), 2:(DIMX-1)] -
        solid[2:(DIMY-1), 2:(DIMX-1)]) +
    abs(solid[2:(DIMY-1), 1:(DIMX-2)] -
        solid[2:(DIMY-1), 2:(DIMX-1)])
# increase to 2 pixel thickness outline
outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
    outline[1:(DIMY-2), 2:(DIMX-1)]+outline[2:(DIMY-1), 3:(DIMX-0)]
# increase to 3 pixel thickness outline
outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
     outline[1:(DIMY-2), 2:(DIMX-1)]+outline[3:(DIMY-0), 2:(DIMX-1)]+
     outline[2:(DIMY-1), 1:(DIMX-2)]+outline[2:(DIMY-1), 3:(DIMX-0)]
outline[outline!=0]=1
writePNG(outline, "mapoutline.png")


# Calculate grayscale hillshade
MIX=0.85  # two light sources are mixed to fill dark areas a bit
hill=hillshademap(DEM, dx=resolution, dlight=c(0, 2, 3))
hillfill=hillshademap(DEM, dx=resolution, dlight=c(0, 3, 2))
hill=hill*MIX+hillfill*(1-MIX)
gamma=1/2.0
hill=(hill/max(hill))^(1/gamma)  # darken hillshade a bit

# Save hillshade
writeTIFF(hill, "hillshade.tif",
          bits.per.sample=16, compression="LZW")

# Display hillshade
image(t(hill[nrow(hill):1,]), useRaster=TRUE,
      col=c(gray.colors(256, start=0, end=1, gamma=0.5)),
      asp=nrow(hill)/ncol(hill), axes=FALSE)
