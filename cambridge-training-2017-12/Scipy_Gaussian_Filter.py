#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------------------
  Copyright (C) 2017 University of Dundee. All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

------------------------------------------------------------------------------

Run Gaussian Filter from Scipy on P/D/I/SPW.

@author Balaji Ramalingam
<a href="mailto:b.ramalingam@dundee.ac.uk">b.ramalingam@dundee.ac.uk</a>
"""

import omero
from omero.model import PlateI, ScreenI

import os
import sys
import subprocess
import re
import tempfile
import platform
import glob
import smtplib
# For guessing MIME type based on file name extension
import mimetypes
import datetime

import omero
import omero.scripts as scripts
from omero.rtypes import wrap, rlong
from omero.gateway import OriginalFileWrapper, MapAnnotationWrapper
from omero.gateway import BlitzGateway
from omero.rtypes import *  # noqa
from omeroweb.webgateway.marshal import imageMarshal

import json
import ConfigParser
from cStringIO import StringIO

import scipy
import scipy.ndimage as spi


def run(conn, params):
    """
    For each image defined in the script parameters run the correlation
    analyser and load the result into OMERO.
    Returns the number of images processed or -1 if there is a
    parameter error.
    @param conn   The BlitzGateway connection
    @param params The script parameters
    """

    print "Parameters = %s" % params

    if not params.get("Kernel_Window_Size"):
        print "Please enter a valid window size for the Gaussian kernel"
        return -1
    window_size = params.get("Kernel_Window_Size")

    if not params.get("Sigma"):
        print "Please enter a valid value for Sigma (>0)"
        return -1
    sigma = params.get("Sigma")

    if not params.get("IDs"):
        print "Please enter a valid omero_object (Project/Dataset/Image) Id"
        return -1

    images = []

    if params.get("Data_Type") == 'Dataset':
        for dsId in params["IDs"]:
            dataset = conn.getObject("Dataset", dsId)
            if dataset:
                for image in dataset.listChildren():
                    images.append(image)

    # Extract images
    image_ids = []
    for image in images:

        sizeZ = image.getSizeZ()
        sizeC = image.getSizeC()
        sizeT = image.getSizeT()
        zctList = []
        for z in range(sizeZ):
            for c in range(sizeC):
                for t in range(sizeT):
                    zctList.append((z, c, t))

        plane = planeGen(image, zctList, window_size, sigma)
        name = image.getName() + "_numpy_image"
        i = conn.createImageFromNumpySeq(plane, name, sizeZ, sizeC,
                                         sizeT, description="Gaussian Filter",
                                         dataset=dataset)

        image_ids.append(i.getId())

        add_map_annotation(conn, i, params)
    return image_ids


def planeGen(image, zctList, window, sigma):
    """generator will yield planes"""

    planes = image.getPrimaryPixels().getPlanes(zctList)
    for p in planes:
        truncated_param = (((window - 1)/2)-0.5)/sigma
        data = spi.filters.gaussian_filter(p, sigma=sigma,
                                           truncate=truncated_param)
        yield data


def add_map_annotation(conn, image, params):

    window_size = params.get("Kernel_Window_Size")
    sigma = params.get("Sigma")
    key_value_data = [["Kernel Window Size", str(window_size)],
                      ["Sigma", str(sigma)]]
    print "Adding MAP", key_value_data, image.getId()
    map_ann = MapAnnotationWrapper(conn)
    # Use 'client' namespace to allow editing in Insight & web
    map_ann.setNs(omero.constants.metadata.NSCLIENTMAPANNOTATION)
    map_ann.setValue(key_value_data)
    map_ann.save()
    # NB: only link a client map annotation to a single object
    image.linkAnnotation(map_ann)

# OMERO Figure Methods


def get_panel_json(image, x, y, width, height, channel=None):
    """
    Export the results as OMERO.figure
    """

    rv = imageMarshal(image)

    if channel is not None:
        for idx, ch in enumerate(rv['channels']):
            ch['active'] = idx == channel
            ch['color'] = 'ffffff'

    img_json = {
        "labels": [],
        "height": height,
        "channels": rv['channels'],
        "width": width,
        "sizeT": rv['size']['t'],
        "sizeZ": rv['size']['z'],
        "dx": 0,
        "dy": 0,
        "rotation": 0,
        "imageId": image.getId(),
        "name": image.getName(),
        "orig_width": rv['size']['width'],
        "zoom": 100,
        "shapes": [],
        "orig_height": rv['size']['height'],
        "theZ": rv['rdefs']['defaultZ'],
        "y": y,
        "x": x,
        "theT": rv['rdefs']['defaultT']
    }
    return img_json


def get_labels_json(panel_json, column, row):

    labels = []

    channels = panel_json['channels']
    imagename = panel_json['name']
    if row == 0:
        labels.append({"text": channels[column]['label'],
                       "size": 4,
                       "position": "leftvert",
                       "color": "000000"})
    if column == 0:
        labels.append({"text": imagename,
                       "size": 4,
                       "position": "top",
                       "color": "000000"})
    return labels


def create_figure_file(conn, image_ids):

    width = 512/10
    height = 512/10
    spacing_x = 512/50
    spacing_y = 512/50
    page_width = (width + spacing_x) * (5) * 1.25
    page_height = (height + spacing_y) * (len(image_ids)) * 1.25

    JSON_FILEANN_NS = "omero.web.figure.json"

    figure_json = {"version": 2,
                   "paper_width": page_width,
                   "paper_height": page_height,
                   "page_size": "mm",
                   "figureName": "from script",
                   }

    curr_x = 0
    curr_y = 0
    panels_json = []
    column_count = 2
    offset = 10

    gid = -1
    for z, image_id in enumerate(image_ids):
        image = conn.getObject('Image', image_id)
        curr_x = z * (width + spacing_x) + offset
        if z == 0:
            gid = image.getDetails().getGroup().getId()
        for c in range(image.getSizeC()):
            curr_y = c * (height + spacing_y) + offset
            j = get_panel_json(image, curr_x, curr_y, width, height, c)
            j['labels'] = get_labels_json(j, c, z)
            panels_json.append(j)

    figure_json['panels'] = panels_json

    figure_name = figure_json['figureName']
    if len(figure_json['panels']) == 0:
        raise Exception('No Panels')
    first_img_id = figure_json['panels'][0]['imageId']

    # we store json in description field...
    description = {}
    description['name'] = figure_name
    description['imageId'] = first_img_id

    # Try to set Group context to the same as first image
    conn.SERVICE_OPTS.setOmeroGroup(gid)

    json_string = json.dumps(figure_json)
    file_size = len(json_string)
    f = StringIO()
    json.dump(figure_json, f)

    update = conn.getUpdateService()
    orig_file = conn.createOriginalFileFromFileObj(f, '', figure_name,
                                                   file_size,
                                                   mimetype="application/json")
    fa = omero.model.FileAnnotationI()
    fa.setFile(omero.model.OriginalFileI(orig_file.getId(), False))
    fa.setNs(wrap(JSON_FILEANN_NS))
    desc = json.dumps(description)
    fa.setDescription(wrap(desc))
    fa = update.saveAndReturnObject(fa, conn.SERVICE_OPTS)
    return fa.getId().getValue()


if __name__ == "__main__":
    dataTypes = [rstring('Dataset')]
    client = scripts.client(
        'Scipy_Gaussian_Filter.py',
        """
    This script applies a gaussian filter to the selected images,
    uploads the generated images to OMERO and creates an OMERO.figure
        """,
        scripts.String(
            "Data_Type", optional=False, grouping="1",
            description="Choose source of images",
            values=dataTypes, default="Dataset"),

        scripts.List(
            "IDs", optional=False, grouping="2",
            description="Dataset IDs.").ofType(rlong(0)),

        scripts.Int(
            "Kernel_Window_Size", optional=False, grouping="3", default=20,
            description="Window size for the gaussian filter"),

        scripts.Int(
            "Sigma", optional=False, grouping="4", default=2,
            description="Sigma for the gaussian filter"),

        scripts.Bool(
            "Create_Omero_Figure", default=True, grouping="5",
            description="Create An OMERO.Figure from the resultant images"),

        authors=["Balaji Ramalingam", "OME Team"],
        institutions=["University of Dundee"],
        contact="ome-users@lists.openmicroscopy.org.uk",
    )

    try:
        # process the list of args above.
        scriptParams = {}
        for key in client.getInputKeys():
            if client.getInput(key):
                scriptParams[key] = client.getInput(key, unwrap=True)
        print scriptParams

        # wrap client to use the Blitz Gateway
        conn = BlitzGateway(client_obj=client)
        # # Call the main script - returns the number of images processed
        image_ids = run(conn, scriptParams)

        if scriptParams["Create_Omero_Figure"]:
            create_figure_file(conn, image_ids)

        message = "Done"
        client.setOutput("Message", rstring(message))

    finally:
        client.closeSession()
