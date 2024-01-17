# Lake Michigan Horizons
# www.overfitting.net
# https://www.overfitting.net/2024/01/conectando-los-grandes-lagos-al-mar-con.html


library(terra)  # read GeoTIFF, reprojection, crop and resample
library(tiff)  # save 16-bit TIFF's
library(png)  # save 8-bit PNG's


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


# CROP raster to area of interest (in km)
cropdef=ext(-500, 0, 4800, 5500)
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

# 2. PROCESS MATRIX TO OBTAIN MAP CONTOURS

# Calculate solid map contour
solid=DEM

# Superior, Michigan, Huron, Erie contours: ~179m above sea level
solid[solid<179]=0  # set areas to 0 (water) 
solid[solid>=179]=1  # set areas to 1 (land)

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
# outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
#     outline[1:(DIMY-2), 2:(DIMX-1)]+outline[2:(DIMY-1), 3:(DIMX-0)]
# increase to 3 pixel thickness outline
# outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
#     outline[1:(DIMY-2), 2:(DIMX-1)]+outline[3:(DIMY-0), 2:(DIMX-1)]+
#     outline[2:(DIMY-1), 1:(DIMX-2)]+outline[2:(DIMY-1), 3:(DIMX-0)]
outline[outline!=0]=1

writePNG(outline, "mapoutline.png")


###########################################################

# 3. OBTAIN HORIZON VIEWS


indices=unname(which(outline==1, arr.ind=TRUE))  # contour pixels

Rearth=6371.23e3  # Earth average radius (m)
hmax=240  # max height (Madrid Towers)
dh=20  # 20m steps
h=seq(from=0, to=hmax, by=dh)

horizon=outline*0
NHORIZ=length(h)
for (j in 2:NHORIZ) {
    print(paste0(j, "/", NHORIZ, "..."))
    d=((Rearth+h[j])^2 - Rearth^2)^0.5
    R=round(d/resolution)  # radius in pixels
    horizontmp=outline*0
    for (i in 1:nrow(indices)) {
        x0=indices[i,2]
        y0=indices[i,1]
        if (x0>R & y0>R & x0<DIMX-R & y0<DIMY-R)
            for (x in round(x0-R):round(x0+R)) {
                for (y in round(y0-R):round(y0+R)) {
                    if ( ((x-x0)^2 + (y-y0)^2 )^0.5 < R) horizontmp[y,x]=horizontmp[y,x]+1
                }
            }
    }
    horizontmp[solid==1]=0  # remove calculated visible water on land
    horizontmp[horizontmp!=0]=1  # clip to 1 all visible water values
    horizon=horizon+horizontmp  # accumulate horizon ranges
}

writePNG(horizon/max(horizon), "horizon.png")

