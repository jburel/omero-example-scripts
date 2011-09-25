from omero.gateway import BlitzGateway
from omero.rtypes import *
user = 'will'
pw = 'ome'
host = 'localhost'
imageId = 101
conn = BlitzGateway(user, pw, host=host, port=4064)
conn.connect()


# Retrieve a given plane

# Use the pixelswrapper to retrieve the plane as 
# a 2D numpy array. See http://www.scipy.org/Tentative_NumPy_Tutorial

# Numpy array can be used for various analysis routines

image = conn.getObject("Image", imageId)
sizeZ = image.getSizeZ()
sizeC = image.getSizeC()
sizeT = image.getSizeT()
z,t,c = 0,0,0                       # first plane of the image
pixels = image.getPrimaryPixels()
plane = pixels.getPlane(z,c,t)      # get a numpy array.
print "\nPlane at zct: ", z, c, t
print plane
print "shape: ", plane.shape
print "min:", plane.min(), " max:", plane.max(), "pixel type:", plane.dtype.name


# Retrieve a given stack

# Get a Z-stack of tiles. 
# Using getTiles or getPlanes (see below) returns a generator of data (not all the data in hand)
# The RawPixelsStore is only opened once (not closed after each plane)
# Alternative is to use getPlane() or getTile() multiple times - slightly slower.

c,t = 0,0                   # First channel and timepoint
tile = (50, 50, 10, 10)     # x, y, width, height of tile
zctList = [(z, c, t, tile) for z in range(sizeZ)]     # list of [ (0,0,0,(x,y,w,h)), (1,0,0,(x,y,w,h)), (2,0,0,(x,y,w,h))....etc... ]
print "\nZ stack of tiles:"
planes = pixels.getTiles(zctList)
for i, p in enumerate(planes):
    print "Tile:", zctList[i], " min:", p.min(), " max:", p.max(), " sum:", p.sum()


#Retrieve a given hypercube

zctList = []
for z in range(sizeZ/2, sizeZ):     # get the top half of the Z-stack
    for c in range(sizeC):          # all channels
        for t in range(sizeT):      # all time-points
            zctList.append( (z,c,t) )
print "\nHyper stack of planes:"
planes = pixels.getPlanes(zctList)
for i, p in enumerate(planes):
    print "plane zct:", zctList[i], " min:", p.min(), " max:", p.max()