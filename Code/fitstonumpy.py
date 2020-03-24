
from pathlib import Path
from astropy.io import fits
import numpy as np
import glob
for filepath in glob.iglob('../Data/cutouts/**/*.fits', recursive=True):
    fp = Path(filepath)


    hdulist = fits.open(fp)
    scidata = hdulist[0].data

    parts = list(fp.parts)

    parts[2] = "tfcutouts"

    fp = Path(*parts)

    fp = fp.with_suffix('.npy')

    np.save(fp, scidata)

