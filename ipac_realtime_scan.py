#!/usr/bin/env python
from __future__ import division
from flask import Flask, make_response, render_template
try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache
import functools
import glob
import logging
import os.path
import cStringIO as StringIO
import astropy.io.ascii
import astropy.io.fits
import astropy.stats
import numpy as np
import math

# PIL imports
import PIL.Image

app = Flask(__name__)
app.config['PROPAGATE_EXCEPTIONS'] = True
app.logger.setLevel(logging.INFO)

mirrordir = 'ptf.nersc.gov/project/deepsky/ptfvet/fermi'


@lru_cache(maxsize=20)
def psf_sub_catalog(exposure, chip):
    """Find PSF fit catalog for the subtraction image."""
    pattern = os.path.join(
        mirrordir,
        'e%d/c%d/p15/v1/OutDiffProds' % (exposure, chip),
        '*_flattened_pmtchscimrefpsffit.tbl')
    filename, = glob.glob(pattern)
    return astropy.io.ascii.read(filename)


@lru_cache(maxsize=20)
def sub_fits(exposure, chip):
    """Find PSF fit catalog for the subtraction image."""
    pattern = os.path.join(
        mirrordir,
        'e%d/c%d/p15/v1/OutDiffProds' % (exposure, chip),
        '*_flattened_pmtchscimref.fits')
    filename, = glob.glob(pattern)
    return astropy.io.fits.open(filename)


@lru_cache(maxsize=20)
def new_fits(exposure, chip):
    """Find PSF fit catalog for the subtraction image."""
    pattern = os.path.join(
        mirrordir,
        'e%d/c%d/p15/v1' % (exposure, chip),
        '*_flattened.fits')
    filename, = glob.glob(pattern)
    return astropy.io.fits.open(filename)


def _find_top_sources():
    for exposure in (407289, 407290, 407292, 407293, 407296, 407297, 407300, 407303):
        for chip in (0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11):
            try:
                catalog = psf_sub_catalog(exposure, chip)
            except ValueError:
                continue # Yuck! skipping fields w/o sub images
            keep = catalog['nneg'] < 30
            keep &= catalog['nneg'] != -999
            keep &= catalog['dnear'] > 4
            keep &= catalog['magpsf'] < 16
            keep_catalog = catalog[keep]
            if keep_catalog:
                keep_catalog.add_column(astropy.table.Column([exposure] * len(keep_catalog), 'exposure'))
                keep_catalog.add_column(astropy.table.Column([chip] * len(keep_catalog), 'chip'))
                yield keep_catalog
def find_top_sources():
    return astropy.table.vstack(list(_find_top_sources()))

@app.route('/')
def top_sources():
    catalog = find_top_sources()
    catalog.sort('magpsf')
    return render_template('index.html', catalog=catalog)


def plot_cutout(fits, i0, j0, width=101):
    """Plot a cutout image from a FITS file. Set the color scale
    automatically using sigma-clipping to extend one sigma below and three
    sigma above the median pixel value.

    The image is returned as an 8-bit grayscale PNG, with one pixel for each
    pixel in the cutout of the FITS image. We'll add the reticule as an SVG
    layer in the web page so that we don't need to encode color PNGs, which
    would use more bandwidth."""
    # Round image center to nearest pixel
    i0 = int(round(i0))
    j0 = int(round(j0))

    # Get dimensions of full FITS image
    imax, jmax = fits[0].data.shape

    # Compute sub-range of cutout, performing
    # bounds checking
    i1 = max(0, i0 - width // 2)
    i1_ = i1 - (i0 - width // 2)
    i2 = min(imax, i0 + (width + 1) // 2)
    i2_ = i2 - (i0 - width // 2)
    j1 = max(0, j0 - width // 2)
    j1_ = j1 - (j0 - width // 2)
    j2 = min(jmax, j0 + (width + 1) // 2)
    j2_ = j2 - (j0 - width // 2)

    # Get cutout from FITS image
    cutout = fits[0].data[i1:i2, j1:j2]

    # Get flattened view of cutout
    cutout_flattened = cutout.flatten()

    # Determine appropriate data range using sigma clipping
    loc = np.mean(astropy.stats.sigma_clip(cutout_flattened))
    scale = np.std(astropy.stats.sigma_clip(cutout_flattened))

    # Intialize image buffer
    data = np.empty((width, width), dtype=np.uint8)
    data[:, :] = 128

    data[i1_:i2_, j1_:j2_] = np.clip(np.round(255/4 * ((cutout - loc) / scale + 1)), 0, 255).astype(np.uint8)

    # Generate PNG output
    output = StringIO.StringIO()
    PIL.Image.fromarray(data[:, ::-1], 'L').save(output, 'png')
    response = make_response(output.getvalue())
    response.mimetype = 'image/png'
    return response


@app.route('/<int:exposure>/<int:chip>/<int:candidate>/sub.png')
def sub_png(exposure, chip, candidate):
    # Look up row in PSF fit catalog
    psf_row = psf_sub_catalog(exposure, chip)[candidate]

    i0 = psf_row['ypos']
    j0 = psf_row['xpos']

    # Open FITS file
    fits = sub_fits(exposure, chip)

    return plot_cutout(fits, i0, j0)


@app.route('/<int:exposure>/<int:chip>/<int:candidate>/new.png')
def new_png(exposure, chip, candidate):
    # Look up row in PSF fit catalog
    psf_row = psf_sub_catalog(exposure, chip)[candidate]

    i0 = psf_row['ypos']
    j0 = psf_row['xpos']

    # Open FITS file
    fits = new_fits(exposure, chip)

    return plot_cutout(fits, i0, j0)


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8888, threaded=True)

