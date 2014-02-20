#!/usr/bin/env python

from __future__ import division
from flask import Flask, make_response, render_template
try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache
import functools
import glob
import os.path
import cStringIO as StringIO
import astropy.io.ascii
import astropy.io.fits
import astropy.stats
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import math

# PIL imports
import PIL.Image

app = Flask(__name__)
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
            keep &= catalog['nneg'] != 999
            # keep &= catalog['dnear'] > 3
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
    #table_html = '<br>'.join('<a href="{exposure}/{chip}/{candidate}/sub.png">sub image</a>'.format(exposure=row['exposure'], chip=row['chip'], candidate=row['sourceid']) for row in catalog)
    #table_lines = catalog.pformat(html=True, max_width=np.inf, max_lines=np.inf)
    #for table_line in table_lines:
    #table_html = ''.join(catalog.pformat(html=True, max_width=np.inf, max_lines=np.inf))
    return render_template('index.html', catalog=catalog)

def plot_cutout(fits, i0, j0, width=101):
    # Get appropriate sub-ranges
    i1 = i0 - 0.5*width
    i2 = i0 + 0.5*width
    j1 = j0 - 0.5*width
    j2 = j0 + 0.5*width

    imax, jmax = fits[0].data.shape

    i1ind = max(0, int(np.floor(i1)))
    i2ind = min(imax, int(np.ceil(i2) + 1))
    j1ind = max(0, int(np.floor(j1)))
    j2ind = min(jmax, int(np.ceil(j2) + 1))

    imgdata = fits[0].data[i1ind:i2ind, j1ind:j2ind]

    # Determine appropriate data range using sigma clipping
    loc = np.mean(astropy.stats.sigma_clip(imgdata.flatten()))
    scale = np.std(astropy.stats.sigma_clip(imgdata.flatten()))

    # Set up figure
    fig = Figure(figsize=(2, 2), dpi=100)
    fig.subplots_adjust(top=1, left=0, bottom=0, right=1)
    ax = fig.add_subplot(111)

    # Plot reticule
    kwargs = dict(color=(0.5, 1, 0.5), linewidth=1)
    ax.plot([i0 - 0.125*width, i0 - 0.0625*width], [j0, j0], **kwargs)
    ax.plot([i0 + 0.125*width, i0 + 0.0625*width], [j0, j0], **kwargs)
    ax.plot([i0, i0], [j0 - 0.125*width, j0 - 0.0625*width], **kwargs)
    ax.plot([i0, i0], [j0 + 0.125*width, j0 + 0.0625*width], **kwargs)

    # Plot binary image (FIXME: make sure north up)
    ax.imshow(
        imgdata,
        extent=[i2ind, i1ind, j1ind, j2ind],
        interpolation='nearest',
        cmap='binary_r',
        vmin=loc - scale,
        vmax=loc + 3*scale)
    ax.set_xlim(i1, i2)
    ax.set_ylim(j1, j2)
    ax.axis('off')

    # Generate PNG output
    canvas = FigureCanvas(fig)
    output = StringIO.StringIO()
    canvas.print_png(output)
    response = make_response(output.getvalue())
    response.mimetype = 'image/png'
    return response

def plot_cutout_pil(fits, i0, j0, width=101):
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

    return plot_cutout_pil(fits, i0, j0)


@app.route('/<int:exposure>/<int:chip>/<int:candidate>/new.png')
def new_png(exposure, chip, candidate):
    # Look up row in PSF fit catalog
    psf_row = psf_sub_catalog(exposure, chip)[candidate]

    i0 = psf_row['ypos']
    j0 = psf_row['xpos']

    # Open FITS file
    fits = new_fits(exposure, chip)

    return plot_cutout_pil(fits, i0, j0)


if __name__ == '__main__':
    app.run(port=8888, debug=True, threaded=True)

