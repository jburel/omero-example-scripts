# -----------------------------------------------------------------------------
#  Copyright (C) 2018 University of Dundee. All rights reserved.
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# ------------------------------------------------------------------------------

# This Jython script uses ImageJ to analyse particles
# The purpose of the script is to be used in the Scripting Dialog
# of Fiji.
# Error handling is omitted to ease the reading of the script but
# this should be added
# if used in production to make sure the services are closed
# Information can be found at
# https://docs.openmicroscopy.org/omero/latest/developers/Java.html

from java.lang import Long
from java.lang import Float
from java.lang import System
from java.lang import Math
from java.util import ArrayList
from java.util import Arrays
from java.lang.reflect import Array
from jarray import zeros, array
import java

# Omero Dependencies
from omero.gateway import Gateway
from omero.gateway import LoginCredentials
from omero.gateway import SecurityContext
from omero.gateway.facility import BrowseFacility
from omero.gateway.facility import AdminFacility
from omero.gateway.facility import DataManagerFacility
from omero.gateway.facility import MetadataFacility
from omero.gateway.facility import ROIFacility
from omero.gateway.model import DatasetData
from omero.log import SimpleLogger

# this needs to be changed
from org.openmicroscopy.shoola.util.roi.io import ROIReader

from loci.formats import FormatTools, ImageTools
from loci.common import DataTools
from loci.plugins.util import ImageProcessorReader

from ij import IJ, ImagePlus, ImageStack
from ij.process import ByteProcessor, ShortProcessor
from ij.process import ImageProcessor
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable

from java.lang import Object
from omero.gateway.model import TagAnnotationData


# Setup
# =====

# OMERO Server details
HOST = "outreach.openmicroscopy.org"
PORT = 4064
#  parameters to edit
USERNAME = "username"
PASSWORD = "password"
# We want to process Images within this Dataset....
dataset_id = 974
# ...that are Tagged with this Tag
tag_text = "Mitosis"


# Connection method: returns a gateway object
def connect_to_omero():
    "Connect to OMERO"

    credentials = LoginCredentials()
    credentials.getServer().setHostname(HOST)
    credentials.getServer().setPort(PORT)
    credentials.getUser().setUsername(USERNAME.strip())
    credentials.getUser().setPassword(PASSWORD.strip())
    simpleLogger = SimpleLogger()
    gateway = Gateway(simpleLogger)
    user = gateway.connect(credentials)
    return gateway, user


# Convert omero Image object as ImageJ ImagePlus object
# (An alternative to OmeroReader)
def open_omero_image(ctx, image_id):
    browse = gateway.getFacility(BrowseFacility)
    print image_id
    image = browse.getImage(ctx, long(image_id))
    pixels = image.getDefaultPixels()
    size_z = pixels.getSizeZ()
    size_t = pixels.getSizeT()
    size_c = pixels.getSizeC()
    size_x = pixels.getSizeX()
    size_y = pixels.getSizeY()
    pixtype = pixels.getPixelType()
    pixels_type = FormatTools.pixelTypeFromString(pixtype)
    bpp = FormatTools.getBytesPerPixel(pixels_type)
    is_signed = FormatTools.isSigned(pixels_type)
    is_float = FormatTools.isFloatingPoint(pixels_type)
    is_little = False
    interleave = False
    store = gateway.getPixelsStore(ctx)
    pixels_id = pixels.getId()
    store.setPixelsId(pixels_id, False)
    stack = ImageStack(size_x, size_y)
    for t in range(0, size_t):
        for z in range(0, size_z):
            for c in range(0, size_c):
                plane = store.getPlane(z, c, t)

                channel = ImageTools.splitChannels(plane, 0, 1, bpp,
                                                   False, interleave)
                pixels = DataTools.makeDataArray(plane, bpp,
                                                 is_float, is_little)

                q = pixels
                if len(plane) != size_x*size_y:
                    tmp = q
                    q = zeros(size_x*size_y, 'h')
                    System.arraycopy(tmp, 0, q, 0, Math.min(len(q), len(tmp)))
                    if is_signed:
                        q = DataTools.makeSigned(q)
                if q.typecode == 'b':
                    ip = ByteProcessor(size_x, size_y, q, None)
                elif q.typecode == 'h':
                    ip = ShortProcessor(size_x, size_y, q, None)
                stack.addSlice('', ip)
    # Do something
    image_name = image.getName() + '--OMERO ID:' + str(image.getId())
    imp = ImagePlus(image_name, stack)
    imp.setDimensions(size_c, size_z, size_t)
    imp.setOpenAsHyperStack(True)
    imp.show()
    return imp


def list_images_in_dataset(ctx, datset_id):
    browse = gateway.getFacility(BrowseFacility)
    ids = ArrayList(1)
    ids.add(Long(dataset_id))
    return browse.getImagesForDatasets(ctx, ids)


def filter_images_by_tag(ctx, images, tag_value):
    metadata_facility = gateway.getFacility(MetadataFacility)
    tagged_image_ids = []
    for image in images:
        annotations = metadata_facility.getAnnotations(ctx, image)
        for ann in annotations:
            if isinstance(ann, TagAnnotationData):
                if ann.getTagValue() == tag_value:
                    tagged_image_ids.append(image.getId())
    return tagged_image_ids


def save_rois_to_omero(ctx, image_id, imp):
    # Save ROI's back to OMERO
    reader = ROIReader()
    roi_list = reader.readImageJROIFromSources(image_id, imp)
    roi_facility = gateway.getFacility(ROIFacility)
    return roi_facility.saveROIs(ctx, image_id, exp_id, roi_list)


# Prototype analysis example
gateway, user = connect_to_omero()
ctx = SecurityContext(user.getGroupId())
exp = gateway.getLoggedInUser()
exp_id = exp.getId()

images = list_images_in_dataset(ctx, dataset_id)
print "Number of images in Dataset", len(images)

ids = filter_images_by_tag(ctx, images, tag_text)

print "tagged_image_ids", ids

# Input ids in a comma seperated fashion
IJ.run("Set Measurements...", "area mean standard modal min centroid center \
        perimeter bounding fit shape feret's integrated median skewness \
        kurtosis area_fraction stack display redirect=None decimal=3")

for id1 in ids:
    # if target_user ~= None:
    # Switch context to target user and open omeroImage as ImagePlus object
    imp = open_omero_image(ctx, id1)
    # Some analysis which creates ROI's and Results Table
    IJ.setAutoThreshold(imp, "Default dark")
    IJ.run(imp, "Analyze Particles...",
           "size=50-Infinity display clear add stack")
    rm = RoiManager.getInstance()
    rm.runCommand(imp, "Measure")
    save_rois_to_omero(ctx, id1, imp)
    # Close the various components
    IJ.selectWindow("Results")
    IJ.run("Close")
    IJ.selectWindow("ROI Manager")
    IJ.run("Close")
    imp.close()

print "processing done"
