# author: Igor Andreoni - andreoni@caltech.edu

import subprocess
import os
import pdb

import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt

from astroquery.vizier import Vizier
Vizier.ROW_LIMIT = -1
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy import wcs

# For forced photometry
from astropy.io import fits
from astropy.visualization import simple_norm
from astropy.table import Table
from astropy.nddata import NDData
from photutils.psf import extract_stars
from photutils import EPSFBuilder
from photutils import find_peaks
from photutils.centroids import centroid_com

# Perform PSF photometry
from photutils.psf import BasicPSFPhotometry
from photutils.background import MMMBackground
from astropy.modeling.fitting import LevMarLSQFitter
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.stats import sigma_clip


class Target:
    """
    Class for the target source

    Parameters
    ----------
    ra float
        right ascension
    dec float
        declination
    """
    def __init__(self, ra, dec):
        self.ra = ra
        self.dec = dec


class Image:
    """
    Class for the image

    Parameters
    ----------
    filename str
        image file name
    ext int
        extension of interest in the fits file
    """
    def __init__(self, filename, ext=0, satlevel=None):
        self.filename = filename
        # zero point
        self.zp = 0.
        # header and data
        im = fits.open(filename)
        self.header = im[ext].header
        self.data = im[ext].data
        # filter
        self.filter = self.find_filter()
        # Saturation
        if satlevel is None:
            self.saturate = self.find_saturation()
        else:
            self.saturate = satlevel
        # Image size
        self.naxis1 = len(self.data[0])
        self.naxis2 =  len(self.data[:][:])
        # Pixel scale
        self.pixscale = self.header["PIXSCALE"]
        # right ascension and declination of the image center
        self.ra_center, self.dec_center = self.find_radec_center()
        self.box = np.max([self.naxis1, self.naxis2]) * self.pixscale

    def find_filter(self):
        # Find the filter info from the header
        try:
            return self.header["FILTER"]
        except KeyError:
            print(f"Cannot find FILTER card for image {self.filename}")
            return input("Enter the filter name manually:\n")

    def find_radec_center(self):
        # Find RA and Dec of the image center
        w = wcs.WCS(self.header, relax=True)
        coords = w.pixel_to_world(int(self.naxis1/2), int(self.naxis2/2))
        return coords.ra.deg, coords.dec.deg

    def find_saturation(self, sat=50000):
        # Find the saturation info from the header
        try:
            return self.header["SATURATE"]
        except KeyError:
            print(f"Cannot find FILTER card for image {self.filename}")
            print(f"Using default value of {sat} counts")
            return sat

    def find_size(self):
        # Find the image size
        try:
            return self.header["FILTERZ"]
        except KeyError:
            print(f"Cannot find FILTER card for image {self.filename}")
            return input("Enter the filter name manually:\n")


class Stars:
    #Class for the stars to be used for zeropoint calibration
    def __init__(self, im, target, args, sep_max=1):
        """
        Parameters:
        -----------
        im Image object
            input image
        target Target object
            Target of PSF photometry
        args arguments
        sep_max float
            maximum separation (arcsec)
        """
        self.gaia = self.find_stars_gaia(im, target,
                                         box_frac=args.bf)
        self.imsources = self.find_sources_img(im, n_max=600)

        # Crossmatch
        c_imsources = SkyCoord(ra=self.imsources['ra']*u.deg,
                               dec=self.imsources['dec']*u.deg)
        c_gaia = SkyCoord(ra=self.gaia['RA_ICRS'],
                               dec=self.gaia['DE_ICRS'])
        idx, d2d, _ = c_imsources.match_to_catalog_sky(c_gaia)
        sep_constraint = d2d < sep_max * u.arcsec
        self.stars = self.imsources[sep_constraint]

        # Write region files
        writeregionfile(im.filename.replace(".fits", "_gaia.reg"), c_gaia)
        writeregionfile(im.filename.replace(".fits", "_imsources.reg"),
                                            c_imsources)
        writeregionfile(im.filename.replace(".fits", "_stars.reg"), c_imsources[sep_constraint])


    def find_stars_gaia(self, im, target, box_frac=1.):
        """
        Find stars for calibration using Gaia

        Parameters:
        -----------
        im Image object
            input image
        box_frac float
            fraction of the box

        Returns
        -------
        t_star astropy table
            table with coordinates of stars
        """
        # Search box
        box = np.max([im.naxis1, im.naxis2]) * im.pixscale
        # Add 10% to the box size for the catalog search
        box += box * 0.1
        # FIXME
        box = box * box_frac

        # Update box size value for the image
        im.box = box
        # Query gaia centered at the target coordinates
        gaia_cat = prepare_gaia_catalog(target.ra, target.dec, box)
        # Select only stars. Following Luri+18, plx > 1.081
        #return gaia_cat[gaia_cat['Plx'] > 1.081]
        return gaia_cat[gaia_cat['Plx'] > 0.1]

    def find_sources_img(self, im, hsize=50, n_max=600):
        """
        Find sources in the image

        Parameters:
        -----------
        im Image object
            input image
        hsize int
            threshold dist from the edge
        n_max int
            maximum number of sources returned

        Returns
        -------
        t_im astropy table
            table with coordinates of sources in the image
        """
        # Find peaks

        peaks_tbl = find_peaks(im.data, centroid_func=centroid_com,
                               threshold=100., box_size=11)
        peaks_tbl['peak_value'].info.format = '%.8g'
        # Inverse sort by peak value
        peaks_tbl.sort('peak_value', reverse=True)
        # Remove saturating sources
        peaks_tbl = peaks_tbl[peaks_tbl['peak_value'] < im.saturate]

        # Make it short
        x = peaks_tbl['x_centroid']
        y = peaks_tbl['y_centroid']
        # Add 1 pixel to capture the center of the sources
        #x += 1
        #y += 1
        mask = ((x > hsize) & (x < (im.data.shape[1] -1 - hsize)) &
                (y > hsize) & (y < (im.data.shape[0] -1 - hsize)))
        # Create source table
        s_tbl = Table()
        # FIXME
        s_tbl['x'] = x[mask]
        s_tbl['y'] = y[mask]

        # WCS
        w = wcs.WCS(im.header, relax=True)
        s_tbl['ra'], s_tbl['dec'] = w.all_pix2world(s_tbl['x'], s_tbl['y'], 0)
        #s_tbl['dec'] = w.pixel_to_world(s_tbl['x'], s_tbl['y']).dec.deg
        s_tbl['peak_value'] = peaks_tbl[mask]['peak_value']
        s_tbl['name'] = [f"star_{n}" for n in np.arange(len(s_tbl))+1]
        # Select a subset of bright sources
        s_tbl = s_tbl[:n_max]
        return s_tbl


class GaiaAstrometry:
    """
    A class to provide a Gaia catalog for astrometric calibration.
    The query_gaia_astrom() function produces, for any requested field,
    a catalog of sources from Gaia DR2

    Standard usage:
     gaia = GaiaAstrometry( (ra, dec), field_size )  # with ra, dec 
                in degrees and field_size is the side of the field in arcsec
     t_gaia = gaia.query_gaia_astrom()
    """
    def __init__(self, field_center, field_size):
        self.field_ra = field_center[0]
        self.field_dec = field_center[1]
        self.boxsize = field_size/60. * u.arcmin
        # Catalog: Gaia DR2
        self.catID = "I/345/gaia2"

    def query_gaia_astrom(self):
        """
        Query Vizier server for Gaia sources found in a box with center
        self.field_ra, self.field_dec (in deg) and width
        self.field_width (in arcsec).
        Returns a table (if objects found) or None (if not)
        """
        # Create a SkyCoord object
        coords = SkyCoord(ra=self.field_ra*u.deg,
                          dec=self.field_dec*u.deg)
        # Select only those columns with relevant info
        columns_select = ['RA_ICRS', 'e_RA_ICRS', 'DE_ICRS', 'e_DE_ICRS',
                          'Source', 'Gmag', 'e_Gmag', 'Plx']
        # Vizier object
        v = Vizier(columns=columns_select, row_limit=-1)
        # Query Vizier
        t_result = v.query_region(coords, width=self.boxsize,
                                height=self.boxsize,
                                catalog=self.catID)
        # Check if the sources were found
        if len(t_result[self.catID]) == 0:
            return None
        else:
            return t_result[self.catID]

def prepare_gaia_catalog(ra, dec, field_size, write_region=True):
    """Query Gaia DR2 and save a table to do the astrometric calibration

    Parameters
    ----------
    ra float
        Right Ascension of the center of the field in degrees
    dec float
        Declination of the center of the field in degrees
    field_size float
        side of the field, in arcsec
    """
    # Create a GaiaAstrometry object
    gaia = GaiaAstrometry((ra, dec), field_size)
    # Gaia sources table for astrometry
    return gaia.query_gaia_astrom()


def writeregionfile(filename, objlist, color="green", system=''):
    """
    Creates region file for DS9 to box objiects from list

    Parameters
    ----------
    filename str
        output file name
    objlist
        list of ra and dec SkyCoord objects
    color str
        color of the region file
    """
    if system == '':
        system = 'wcs'
    out = open(filename, 'w')
    i = -1
    out.write(
        '# Region file format: DS9 version 4.0\nglobal color='+color +
        ' font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n'
    )
    if system == 'wcs':
        out.write('fk5\n')
        for ob in objlist:
            i += 1
            out.write("point(%.7f,%.7f) # point=box text={%i}\n" % (ob.ra.value, ob.dec.value, i))
    if system == 'img':
        out.write('image\n')
        for ob in objlist:
            i += 1
            out.write("point(%.3f,%.3f) # point=boxcircle text={%i}\n" % (ob.x, ob.y, i))
    out.close()


def query_ps1(ra, dec, box):
    """
    Query PS1 for photometric calibration

    Parameters:
    -----------
    ra float
        right ascension for the center of the query
    dec float
        declination for the center of the query
    box float
        side of the search box (arcsec)

    Returns
    -------
    result astropy table
        table 
    """
    # Query Pan-STARRS
    result = Vizier.query_region(SkyCoord(ra=ra, dec=dec,
                                          unit=(u.deg, u.deg)),
                                 radius=f"{box}s",
                                 catalog=["Pan-STARRS"])
    # Checks
    if len(result[0]) == 0:
        print(f"No match found in PS1!!")
        return None

    return result[0]


def get_psf(nddata, t_stars, size_box=21):

    # Extract stars
    ext_stars = extract_stars(nddata, t_stars, size=size_box)

    # Measure the PSF
    epsf_builder = EPSFBuilder(oversampling=4, maxiters=3,
                               progress_bar=False)
    epsf, fitted_stars = epsf_builder(ext_stars)
    norm = simple_norm(epsf.data, 'log', percent=99.)

    return epsf, norm


def get_psfphot(pos, image, psf_model, size_fit=21):
    # Fix the centroids
    psf_model.x_0.fixed = True
    psf_model.y_0.fixed = True

    sigma_psf = 2.
    daogroup = DAOGroup(2.0*sigma_psf*gaussian_sigma_to_fwhm)
    mmm_bkg = MMMBackground()
    photometry = BasicPSFPhotometry(group_maker=daogroup,
                                    bkg_estimator=mmm_bkg,
                                    psf_model=psf_model,
                                    fitter=LevMarLSQFitter(),
                                    fitshape=(size_fit,size_fit),
                                    finder=None,
                                    aperture_radius=5)
    result = photometry(image=image.data, init_guesses=pos)
    residual_image = photometry.get_residual_image()

    return result


def doCalibrate(result_t, stars_t, filt):
    # Mask out values which will result in problematic mag
    mask = (result_t['flux_fit'] > 0)
    result_t = result_t[mask]
    stars_t = stars_t[mask]
    mag_calib = np.array(stars_t[f'{filt}mag'])

    # Calibrate (remove last entry because it's the target)
    stars_t['flux'] = result_t['flux_fit']
    stars_t['flux_unc'] = result_t['flux_unc']
    stars_t['mag'] = -2.5 * np.log10(np.array(result_t['flux_fit']))
    
    # Diff calculated with sigma-clipping
    diff = np.ma.median(sigma_clip(np.abs(mag_calib - stars_t['mag']), sigma=2, maxiters=5))
    print(f"Zeropoint correction for {filt}-band: {diff}")
    filtered = sigma_clip(np.abs(mag_calib - stars_t['mag']), sigma=2, maxiters=5)
    print(f"The zeropoint was obtained using \
{np.size([x for x in filtered.mask if x==False])} stars")

    # Idicate in the table if the star was used for zpt
    idx = [1 if x==False else 0 for x in filtered.mask]
    stars_t['used_for_zp'] = idx

    stars_t['mag_calib'] = stars_t['mag'] + diff
    stars_t['mag_unc'] =  -2.5 * np.log10(np.array(result_t['flux_fit']-result_t['flux_unc'])) + 2.5 * np.log10(np.array(result_t['flux_fit']))
    stars_t[f'diff'] = stars_t[f'{filt}mag'] -  stars_t[f'mag_calib']

    # Correct the uncertainty by the error in the zp
    stars_t['mag_unc'] = np.abs(stars_t['mag_unc']) + np.std(stars_t[f'diff'])/np.sqrt(len(stars_t))
    stars_t['mag_unc_corr'] = np.abs(stars_t['mag_unc']) + np.std(stars_t[f'diff'])/np.sqrt(len(stars_t))
    print(f"Standard deviation from catalog: {np.std(np.abs(stars_t[f'diff']))}")

    return result_t, stars_t, diff

def getPos(stars_tbl, add1=False):
    # Create a table for the PSF photometry
    if 'used_for_zp' in stars_tbl.colnames:
        colnames = ["x", "y", "ra", "dec", "peak_value", "name", "rmag"]
        p = Table([stars_tbl[cn] for cn in colnames], names=colnames)
    else:
        p = Table(stars_tbl, copy=True)
    if add1 is True:
        p["x"] = [xx + 1 for xx in p["x"]]
        p["y"] = [xx + 1 for xx in p["y"]]
    p.rename_column("x", "x_0")
    p.rename_column("y", "y_0")
    # Add flux guess
    p["flux_0"] = stars_tbl['peak_value']

    return p


def do_forced(target, image, stars, args, use_calib_stars=False):
    # Convert ra, dec to xy
    w = wcs.WCS(image.header, relax=True)
    coords_target = SkyCoord(ra=target.ra*u.deg,
                             dec=target.dec*u.deg)
    # convert to pixel, add 1 to make it right
    x_target, y_target = w.all_world2pix(coords_target.ra, coords_target.dec, 0)
    #import pdb
    #pdb.set_trace()
    #x_target += 1
    #y_target += 1

    # Convert the image in NDData format
    nddata = NDData(data=image.data)

    # Get the PSF
    epsf, norm = get_psf(nddata, stars.stars, size_box=args.sbp)
    
    hdu = fits.PrimaryHDU(epsf.data)
    hdu.writeto('psf.fits', overwrite=True)

    # Stars for calibration
    stars_tbl = Table(stars.stars, copy=True)

    # FIXME generalize filter
    # Filter
    if image.filter in ["SDSS-R", "rp", "ZTF_r"]:
        filt = 'r'
    elif image.filter in ["SDSS-G", "gp", "ZTF_g"]:
        filt = 'g'
    elif image.filter in ["SDSS-I", "ip", "ZTF_i"]:
        filt = 'i'
 
    # Query the PS1 catalog
    index = []
    mag_calib = []

    # Query Pan-STARRS
    cat = query_ps1(target.ra, target.dec, image.box)

    # Crossmatch
    c_cat = SkyCoord(ra=cat['RAJ2000'],
                     dec=cat['DEJ2000'])
    c_stars = SkyCoord(ra=stars_tbl['ra']*u.deg,
                       dec=stars_tbl['dec']*u.deg)
    # Which stars have a crossmatch?
    idx, d2d, _ = c_stars.match_to_catalog_sky(c_cat)
    sep_constraint = d2d < 0.5 * u.arcsec
    stars_tbl = stars_tbl[sep_constraint]

    # Magnitudes for calibration
    stars_tbl[f'{filt}mag'] = [cat[i][f"{filt}mag"] for i in idx[sep_constraint]]

    # Create a table for the PSF photometry
    pos = getPos(stars_tbl, add1=True)

    # Get PSF photometry
    result_tbl = get_psfphot(pos, image, epsf, size_fit=args.sfp)

    # Find the best stars whose photometry best match the calibration catalog
    print("FIRST PHOTOMETRY ITER OF STARS")
    result_tbl, stars_tbl, diff = doCalibrate(result_tbl, stars_tbl, filt)

    # Obtain a new model for the PSF only for good stars
    epsf, norm = get_psf(nddata, stars_tbl[stars_tbl['used_for_zp'] == 1], size_box=args.sbp)

    # Select only those stars that were within the sigma clipping
    pos = getPos(stars_tbl[stars_tbl['used_for_zp'] == 1], add1=False)

    # Add the target at the bottom of the table
    pos.add_row([x_target, y_target,
                 coords_target.ra.deg,
                 coords_target.dec.deg, 1000,
                 'target', 99., 1000])
    # Perform again PSF phot
    result_tab = get_psfphot(pos, image, epsf, size_fit=args.sfp)
    # Separate the target
    result_tab_target = Table(result_tab[result_tab['name']=='target'], copy=True)

    print("SECOND PHOTOMETRY ITER OF STARS+TARGET")
    # Calibrate without the target
    result_tab, stars_tbl, diff = doCalibrate(result_tab[result_tab["name"] != 'target'],
                                    stars_tbl[stars_tbl['used_for_zp'] == 1],
                                    filt)

    # Fix the tables by separating the target again
    ##result_tab_target = Table(result_tab[-1], copy=True)
    #result_tab = Table(result_tab[result_tab['name']!='target'], copy=True)

    # Target table
    result_tab_target["mag_calib"] =  [-2.5 * np.log10(np.array(result_tab_target['flux_fit'])) + diff]
    result_tab_target["mag_unc"] = [-2.5 * np.log10(np.array(result_tab_target['flux_fit']-result_tab_target['flux_unc'])) + 2.5 * np.log10(np.array(result_tab_target['flux_fit']))]
    result_tab_target['mag_unc_corr'] = [np.abs(result_tab_target['mag_unc']) + np.std(stars_tbl[f'diff'])/np.sqrt(len(stars_tbl))]

    # Add possible upper limit
    if (result_tab_target['flux_fit']/result_tab_target['flux_unc']) < 5:
        result_tab_target['mag_lim'] = [-2.5 * np.log10(np.array(result_tab_target['flux_unc']*7)) + diff]
    print(result_tab_target)

    plt.plot(stars_tbl[stars_tbl['used_for_zp']==1][f'{filt}mag'], np.array(stars_tbl[stars_tbl['used_for_zp']==1]['mag_calib'])-np.array(stars_tbl[stars_tbl['used_for_zp']==1][f'{filt}mag']), 'bo', label=f"{len(stars_tbl[stars_tbl['used_for_zp']==1])} stars used for zp")
    plt.plot(stars_tbl[stars_tbl['used_for_zp']==0][f'{filt}mag'], np.array(stars_tbl[stars_tbl['used_for_zp']==0]['mag_calib'])-np.array(stars_tbl[stars_tbl['used_for_zp']==0][f'{filt}mag']), 'ro', label=f"{len(stars_tbl[stars_tbl['used_for_zp']==0])} stars sigma clipped")
    plt.plot([16, 24], [0,0], 'r--')
    plt.plot([16,24], [np.std(np.abs(stars_tbl['diff'])), np.std(np.abs(stars_tbl['diff']))], "k:")
    plt.plot([16,24], [-1*np.std(np.abs(stars_tbl['diff'])), -1*np.std(np.abs(stars_tbl['diff']))], "k:")
    plt.legend()
    plt.xlabel("magnitude")
    plt.ylabel("measured - cataloged mag")


#    plt.plot(np.arange(len(stars_tbl[f'diff'])), np.abs(stars_tbl[f'diff']), 'ro')
#    plt.plot([0,len(stars_tbl)], 2*np.array([np.std(np.abs(stars_tbl[f'diff'])), np.std(np.abs(stars_tbl[f'diff']))]), 'b--')
    plt.show()



def calib_sextractor(filename, doShow=True):
    if True:
        # Run SExtractor
        sexcommand = ["sex", filename, "-c config.sex",
                      f"-CATALOG_NAME {filename.replace('.fits','.cat')}",
                      "-DETECT_THRESH", f"{sex_thresh}",  "-ANALYSIS_THRESH",  f"{sex_thresh}",
                      "-CHECKIMAGE_TYPE SEGMENTATION", f"-CHECKIMAGE_NAME check_{f}.fits" ]

        print(" ".join(sexcommand))
        #subprocess.call(sexcommand) 
        os.system(" ".join(sexcommand))

        # Read in the catalog
        t = ascii.read(filename.replace('.fits','.cat'), format="sextractor")

        # Select bright stellar sources
        t_short = t[t["CLASS_STAR"] > sg_thresh]

        # Select against CRs
        t_short = t_short[t_short["FWHM_IMAGE"] > 3]

        t_short.sort("MAG_AUTO")
        #t_short = t_short[100:200]

        coords_list = SkyCoord(ra=t_short["X_WORLD"], dec=t_short["Y_WORLD"])
        selected_stars = []

        for coords, mag, x, y, cl in zip(coords_list, t_short["MAG_AUTO"],
                                     t_short["X_IMAGE"], t_short["Y_IMAGE"], t_short["CLASS_STAR"]):
            print(coords.ra.deg, coords.dec.deg, x, y, cl)
            result = Vizier.query_region(SkyCoord(ra=coords.ra, dec=coords.dec, unit=(u.deg, u.deg)),
                             radius="1s",
                             catalog=["Pan-STARRS"])
            if len(result) == 0:
                continue
            if result[0][0][f"{f}mag"] < 16:
                continue
            print("FOUND STAR")
            selected_stars.append({"RA": coords.ra, "Dec": coords.dec,
                                   "x": x, "y": y,
                                   "mag": mag, f: result[0][0][f"{f}mag"],
                                   "diff": mag - result[0][0][f"{f}mag"]})

        # Zero point
        zp = -1 * np.median(list(s["diff"] for s in selected_stars if np.ma.is_masked(s["diff"]) == False))

        # Find the limiting magnitude
        mag = -2.5 * np.log10(np.array(t["FLUX_AUTO"])) + zp
        snr = np.array(t["FLUX_AUTO"])/np.array(t["FLUXERR_AUTO"])

        mag_all = []
        snr_all = []
        for m, s in zip(mag, snr):
            if np.isnan(m) == False and np.isnan(s) == False: 
                if s < 50:
                    mag_all.append(m)
                    snr_all.append(s)

        fig, ax = plt.subplots(figsize=(8,6))
        ax.plot(snr_all, mag_all, 'bo', alpha=0.5)

        # Fit
        p = np.poly1d(np.polyfit(snr_all, mag_all, 3))
        x = np.arange(np.min(snr_all), np.max(snr_all))    
        ax.plot(x, p(x), "-r")

        values = [3, 5, 10]
        for v in values:
            ax.plot([v,v], [p(v)-2,p(v)+2], "--", label=f"{f}={p(v):.1f} mag for SNR={v}")
        ax.set_xlabel("SNR", fontsize=18)
        ax.set_ylabel(f"${f}$ (AB magnitude)", fontsize=18)
        ax.tick_params(axis='both',       # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                labelsize=16)
        ax.legend(fontsize=16)
        plt.savefig(f"plot_snr_mag_{f}.png")

        print(f"Zeropoint calculated using {len(selected_stars)} stars")
        print(f"Results obtained using {len(mag_all)} measurements from {len(t_short)} sources")
        for v in values:
            print(f"Limit of {p(v)} for SNR = {v}")
        if doShow is True:
            plt.show()

        return selected_stars


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Do PSF forced photometry')
    parser.add_argument('radec', metavar='RA, Dec', type=str, nargs='+',
                        help='RA and Dec (degrees)')
    parser.add_argument('--im', dest='images', type=str,
                        required=True,
                        help="Single image or file with a list of image \
names starting with @", default=None)
    parser.add_argument('--satlevel', '-sat', dest='satlevel', type=int,
                        required=False,
                        help="Manually provided saturation level in ADU",
                        default=None)
    parser.add_argument('--box-fraction', '-bf', dest='bf', type=float,
                        required=False,
                        help="Fraction of the box side used for the \
calibration star selection (centered on the target; the default is 10 percent\
more than the side of the detector)",
                        default=1.)
    parser.add_argument('--size-box-psf', '-sbp', dest='sbp', type=int,
                        required=False,
                        help="Box size (px) for PSF extraction",
                        default=21)
    parser.add_argument('--size-fit-psf', '-sfp', dest='sfp', type=int,
                        required=False,
                        help="Size of the PSF fitting (px)",
                        default=11)

    # Coordinates for photometry
    args = parser.parse_args()

    # RA and Dec
    ra, dec = float(args.radec[0]), float(args.radec[1])
    
    # Define the target object
    target = Target(ra, dec)

    # Images
    if args.images[0] == "@":
        with open(args.images[1:], 'r') as f:
            im_filenames = f.readlines()
            # Remove newlines
            im_filenames = [i.replace("\n", "") for i in im_filenames]
    else:
        im_filenames = [args.images]
    # Check that there is at least one image in the list
    if len(im_filenames) == 0:
        print("No images found as/in {args.images}")
        exit()

    # Iterate over the input images
    for filename in im_filenames:
        image = Image(filename, ext=0, satlevel=args.satlevel)
        stars = Stars(image, target, args)
        do_forced(target, image, stars, args, use_calib_stars=False)
