"""
Makes kobuk_satellite.jpg, the satellite background kobuk_lake.jl plots under the floes: Esri
World Imagery of Kobuk Lake (Hotham Inlet, Alaska) resampled onto the simulation's
coordinates, covering the domain x in [-40, 40] km, y in [-35, 35] km at 50 m per pixel.

The coastline in kobuk_coast.csv has no geographic reference, so the map from simulation
coordinates (x, y) to UTM zone 4N (EPSG:32604) was found by fitting the coastline to the
lake's water in the image: a rotation by THETA counterclockwise, scaling by SCALE, then a
shift to (X_UTM, Y_UTM). Its error is a few hundred meters.

Requires numpy, scipy and Pillow. Usage: python3 make_kobuk_satellite.py
Imagery: Esri, Maxar, Earthstar Geographics, and the GIS User Community.
"""

import io
import os
import urllib.request

import numpy as np
from PIL import Image
from scipy import ndimage

# simulation coordinates -> UTM 4N
X_UTM, Y_UTM = 371079.0, 7404836.5
THETA = np.radians(4.140)
SCALE = 1.00137

# simulation domain and output resolution (m)
XMIN, XMAX, YMIN, YMAX = -40e3, 40e3, -35e3, 35e3
RES = 50.0

# source image in UTM 4N; finer Esri zoom levels mix in scenes with ice on the lake
SRC_X0, SRC_Y1, SRC_RES, SRC_NX, SRC_NY = 321477.0, 7455511.0, 100.0, 1100, 1000

url = (
    "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/export?"
    f"bbox={SRC_X0},{SRC_Y1 - SRC_NY * SRC_RES},{SRC_X0 + SRC_NX * SRC_RES},{SRC_Y1}"
    f"&bboxSR=32604&imageSR=32604&size={SRC_NX},{SRC_NY}&format=png&f=image"
)
with urllib.request.urlopen(url) as r:
    src = np.asarray(Image.open(io.BytesIO(r.read())).convert("RGB")).astype(float)

# pixel centers of the output image, top row first, in simulation coordinates
x = np.arange(XMIN + RES / 2, XMAX, RES)
y = np.arange(YMAX - RES / 2, YMIN, -RES)
X, Y = np.meshgrid(x, y)
U = SCALE * (np.cos(THETA) * X - np.sin(THETA) * Y) + X_UTM
V = SCALE * (np.sin(THETA) * X + np.cos(THETA) * Y) + Y_UTM
rows, cols = (SRC_Y1 - V) / SRC_RES - 0.5, (U - SRC_X0) / SRC_RES - 0.5
out = np.stack([ndimage.map_coordinates(src[..., k], [rows, cols], order=1) for k in range(3)], axis=-1)

fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "kobuk_satellite.jpg")
Image.fromarray(np.clip(out, 0, 255).astype(np.uint8)).save(fn, quality=90)
print(f"{fn}: {out.shape[1]} x {out.shape[0]} pixels")
