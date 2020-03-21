import numpy as np
import matplotlib
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from astroquery.sdss import SDSS
from astropy import units as u
from astropy.nddata import Cutout2D
import glob
from astropy.io import fits
from astropy import wcs
from astropy.convolution import Gaussian2DKernel
from photutils import Background2D, MedianBackground
from photutils import detect_threshold, detect_sources
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import source_properties
from astropy.coordinates import SkyCoord
import time
galcounter = 0
starcounter = 0
for filepath in glob.iglob(r'../Data/newfits/*.fits', recursive=False):
    print(filepath)
    tic = time.perf_counter()
    fits_image_filename = filepath
    hdu_list = fits.open(fits_image_filename)
    hdu = fits.open(fits_image_filename)[1]
    w = wcs.WCS(hdu.header)
    image_data = fits.getdata(fits_image_filename,0)
    

    data = image_data
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (50, 50), filter_size=(3, 3), bkg_estimator=bkg_estimator)
    threshold = bkg.background + (10. * bkg.background_rms)
    
    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    npixels = 5
    segm = detect_sources(data, threshold, npixels=npixels,filter_kernel=kernel)
    #if ever we want to deblend the source
    #segm_deblend = deblend_sources(data, segm, npixels=npixels,filter_kernel=kernel, nlevels=32,contrast=0.001)
    

    cat = source_properties(data, segm, wcs=w)
    tbl = cat.to_table()
    df = tbl.to_pandas()
    
    #Setting parameters on data
    indexNames = df[df['xcentroid'] < 400].index
    dfsel = df.drop(indexNames)
    indexNames = dfsel[dfsel['xcentroid'] > 7776].index
    dfsel = dfsel.drop(indexNames)
    indexNames = dfsel[dfsel['ycentroid'] < 300].index
    dfsel = dfsel.drop(indexNames)
    indexNames = dfsel[dfsel['ycentroid'] > 5832].index
    dfsel = dfsel.drop(indexNames)
    indexNames = dfsel[dfsel['elongation'] > 10].index
    dfsel = dfsel.drop(indexNames)
    dfsel = dfsel[['xcentroid','ycentroid','sky_centroid.ra','sky_centroid.dec','elongation','equivalent_radius','area']]
    
    coords = list(zip(dfsel.xcentroid, dfsel.ycentroid))
    dfsel['coords'] = coords
    
    #Quereying database for galaxies
    lon1, lat1 = w.all_pix2world(0, 0, 0)
    lon2, lat2 = w.all_pix2world(6132, 8176, 0)
    maxlon = max(lon1,lon2)
    minlon = min(lon1,lon2)
    maxlat = max(lat1,lat2)
    minlat = min(lat1,lat2)
    
    query = f"""
    SELECT ra,dec
    FROM Galaxy
    WHERE ra between {minlon} and {maxlon}
    AND dec between {minlat} and {maxlat}
    AND g < 21
    """
    res = SDSS.query_sql(query)
    galaxies = res.to_pandas()
    
    

    c = SkyCoord(ra=dfsel['sky_centroid.ra']*u.degree, dec=dfsel['sky_centroid.dec']*u.degree)
    catalog = SkyCoord(ra=galaxies['ra']*u.degree, dec=galaxies['dec']*u.degree)
    
    #Cross referencing SDSS and sources
    max_sep = 3.0 * u.arcsec
    idx, d2d, d3d = c.match_to_catalog_3d(catalog)
    sep_constraint = d2d < max_sep
    c_matches = c[sep_constraint]
    catalog_matches = catalog[idx[sep_constraint]]
    res = SDSS.query_sql(query)
    galaxies = res.to_pandas()
    
    
    #Cutting out data
    size = u.Quantity([40, 40], u.pixel)
    for i in range(len(c_matches)):
            
        galcounter += 1
            
        # Make the cutout, including the WCS
        cutout = Cutout2D(image_data, position=c_matches[i], size=size, wcs=w)

        # Put the cutout image in the FITS HDU
        hdu = fits.PrimaryHDU(cutout.data)

        # Update the FITS header with the cutout WCS
        hdu.header.update(cutout.wcs.to_header(relax=False))


        # Write the cutout to a new FITS file
        cutout_filename = '../Data/galaxyfits/galaxyfits-'+str(galcounter)+'.fits'
        hdu.writeto(cutout_filename, overwrite=False)
    
    
    
    
    
    
    
    #Quereying database for stars
    lon1, lat1 = w.all_pix2world(0, 0, 0)
    lon2, lat2 = w.all_pix2world(6132, 8176, 0)
    maxlon = max(lon1,lon2)
    minlon = min(lon1,lon2)
    maxlat = max(lat1,lat2)
    minlat = min(lat1,lat2)
    
    query = f"""
    SELECT ra,dec
    FROM Star
    WHERE ra between {minlon} and {maxlon}
    AND dec between {minlat} and {maxlat}
    AND g < 21
    """
    res = SDSS.query_sql(query)
    stars = res.to_pandas()

    c = SkyCoord(ra=dfsel['sky_centroid.ra']*u.degree, dec=dfsel['sky_centroid.dec']*u.degree)
    catalog = SkyCoord(ra=stars['ra']*u.degree, dec=stars['dec']*u.degree)
    
    #Cross referencing SDSS and sources
    max_sep = 3.0 * u.arcsec
    idx, d2d, d3d = c.match_to_catalog_3d(catalog)
    sep_constraint = d2d < max_sep
    c_matches = c[sep_constraint]
    catalog_matches = catalog[idx[sep_constraint]]
    res = SDSS.query_sql(query)
    galaxies = res.to_pandas()
    
    
    #Cutting out data
    size = u.Quantity([40, 40], u.pixel)
    for i in range(len(c_matches)):
            
        starcounter += 1
            
        # Make the cutout, including the WCS
        cutout = Cutout2D(image_data, position=c_matches[i], size=size, wcs=w)

        # Put the cutout image in the FITS HDU
        hdu = fits.PrimaryHDU(cutout.data)

        # Update the FITS header with the cutout WCS
        hdu.header.update(cutout.wcs.to_header(relax=False))


        # Write the cutout to a new FITS file
        cutout_filename = '../Data/starfits/starfits-'+str(starcounter)+'.fits'
        hdu.writeto(cutout_filename, overwrite=False)
    
    toc = time.perf_counter()
    print(f"Runtime {toc - tic} seconds")