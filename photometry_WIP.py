# author: Igor Andreoni - andreoni@caltech.edu

import subprocess
import os
import pdb

import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt

from astroquery.vizier import Vizier
Vizier.ROW_LIMIT = 9999
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


def do_forced(ra, dec, filename, selected_stars=None, use_calib_stars=False):
    # Read in the image
    im = fits.open(filename)
    # Convert ra, dec to xy
    w = wcs.WCS(im[0].header)
    coords_target = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    x_target, y_target = w.world_to_pixel(coords_target)

    # Convert the image in NDData format
    nddata = NDData(data=im[0].data)

    # Size for the fit
    size_fit = 13
    size_box = 51

    if selected_stars is None:
        # Find peaks
        peaks_tbl = find_peaks(im[0].data, threshold=500.)
        peaks_tbl['peak_value'].info.format = '%.8g'
        # Only 100 brightest sources
        peaks_tbl.sort('peak_value', reverse=True)
        peaks_tbl = peaks_tbl[:100]
        hsize = (size_box - 1) / 2
        x = peaks_tbl['x_peak']  
        y = peaks_tbl['y_peak']  
        mask = ((x > hsize) & (x < (im[0].data.shape[1] -1 - hsize)) &
                (y > hsize) & (y < (im[0].data.shape[0] -1 - hsize)))  

        # Extract the stars
        stars_tbl = Table()
        stars_tbl['x'] = x[mask]  
        stars_tbl['y'] = y[mask]
        stars_tbl['ra'] = w.pixel_to_world(stars_tbl['x'], stars_tbl['y']).ra.deg
        stars_tbl['dec'] = w.pixel_to_world(stars_tbl['x'], stars_tbl['y']).dec.deg
        stars_tbl['name'] = [f"star_{n}" for n in np.arange(len(stars_tbl))+1]
        #stars_tbl['x'] = [s['x'] for s in selected_stars]
        #stars_tbl['y'] = [s['y'] for s in selected_stars]
        stars = extract_stars(nddata, stars_tbl, size=size_box)

    else:
        stars_tbl = Table()
        # Convert ra, dec to xy
        coords = SkyCoord(ra=np.array(selected_stars['ra'])*u.deg, dec=np.array(selected_stars['dec'])*u.deg)
        x_stars, y_stars = w.world_to_pixel(coords)
        stars_tbl['x'] = [round(s) for s in x_stars]
        stars_tbl['y'] = [round(s) for s in y_stars]
        stars_tbl['ra'] = selected_stars['ra']
        stars_tbl['dec'] = selected_stars['dec']
        stars_tbl['name'] = selected_stars['name']
        stars = extract_stars(nddata, stars_tbl, size=size_box)

    # REMOVE
    # Plot stars
    nrows = 2
    ncols = 2
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8, 12),
                           squeeze=True)
    ax = ax.ravel()
    for i in range(nrows*ncols):
        norm = simple_norm(stars[i], 'log', percent=99.)
        ##ax[i].imshow(stars[i], norm=norm, origin='lower', cmap='viridis')
    ##plt.show()

    # Measure the PSF
    epsf_builder = EPSFBuilder(oversampling=4, maxiters=3,
                               progress_bar=False)  
    epsf, fitted_stars = epsf_builder(stars)  
    norm = simple_norm(epsf.data, 'log', percent=99.)
    ##plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
    ##plt.colorbar()
    ##plt.show()
    
    hdu = fits.PrimaryHDU(epsf.data)
    hdu.writeto('psf.fits', overwrite=True)


    # Stars for calibration
        
    # Query the PS1 catalog
    with open(f"stars_{f}.reg", "w") as tf:
        tf.write("# Region file format: DS9 version 4.1 \n")
        tf.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n")
        tf.write("fk5 \n")
        tf.write("\n")
    index = []
    mag_calib = []
    for i, x_star, y_star, ra_star, dec_star, name_star in zip(np.arange(len(stars_tbl)),
                                                               stars_tbl['x'], stars_tbl['y'],
                                                               stars_tbl['ra'],
                                                               stars_tbl['dec'],
                                                               stars_tbl['name']):
        result = Vizier.query_region(SkyCoord(ra=ra_star, dec=dec_star, unit=(u.deg, u.deg)),
                             radius="1s",
                             catalog=["Pan-STARRS"])
        if len(result) == 0:
            print(f"No match found for {name_star}")
            print(ra_star, dec_star)
            continue
        print(f"FOUND match for {name_star}")
        with open(f"stars_{f}.reg", "a") as tf:
            for l in result[0]:
                tf.write(f"circle({l['RAJ2000']}, {l['DEJ2000']}, 1.5" + '") \n')
        index.append(i)
        mag_calib.append(result[0][0][f"{f}mag"])

    print(f"Found {len(index)} stars for calibration")
    if len(index) < 2:
        print("Too few!")
        exit()

    stars_tbl = stars_tbl[index]
    stars_tbl[f'{f}mag'] = mag_calib

    # Perform PSF photometry
    from photutils.psf import BasicPSFPhotometry
    from photutils.background import MMMBackground
    from astropy.modeling.fitting import LevMarLSQFitter
    from photutils.psf import IntegratedGaussianPRF, DAOGroup
    from astropy.stats import gaussian_sigma_to_fwhm
    from astropy.stats import sigma_clip

    # Photometry of stars (note: with un-fixed centroids it does not work)
    psf_model = epsf
    psf_model.x_0.fixed = True
    psf_model.y_0.fixed = True
    pos = stars_tbl
    pos.rename_column("x", "x_0")
    pos.rename_column("y", "y_0")
    # Add flux guess
    pos["flux_0"] = 100

    sigma_psf = 2.
    daogroup = DAOGroup(2.0*sigma_psf*gaussian_sigma_to_fwhm)
    mmm_bkg = MMMBackground()
    photometry = BasicPSFPhotometry(group_maker=daogroup,
                                    bkg_estimator=mmm_bkg,
                                    psf_model=psf_model,
                                    fitter=LevMarLSQFitter(),
                                    fitshape=(size_fit,size_fit))
    result_tab = photometry(image=im[0].data, init_guesses=pos)
    residual_image = photometry.get_residual_image()
    hdu = fits.PrimaryHDU(residual_image)
    hdu.writeto(f'residual_stars_{f}.fits', overwrite=True)

    # Photometry of the target
    psf_model.x_0.fixed = True
    psf_model.y_0.fixed = True
    pos = Table([[],[],[],[],[]],
                names=("x_0", "y_0", "ra", "dec", "flux_0"))
    pos.add_row([x_target, y_target, coords_target.ra.deg, coords_target.dec.deg, 99.])
    result_tab_target = photometry(image=im[0].data, init_guesses=pos)
    residual_image = photometry.get_residual_image()
    hdu = fits.PrimaryHDU(residual_image)
    hdu.writeto(f'residual_target_{f}.fits', overwrite=True)

    # Calibrate
    stars_tbl['flux'] = result_tab['flux_fit']
    stars_tbl['flux_unc'] = result_tab['flux_unc']
    stars_tbl['mag'] = -2.5 * np.log10(np.array(result_tab['flux_fit']))
    
    # Diff calculated with sigma-clipping
    diff = np.median(sigma_clip(np.abs(mag_calib - stars_tbl['mag']), sigma=2, maxiters=5))
    print(f"Zeropoint correction for {f}-band: {diff}")
    stars_tbl['mag_calib'] = stars_tbl['mag'] + diff
    stars_tbl['mag_unc'] =  -2.5 * np.log10(np.array(result_tab['flux_fit']-result_tab['flux_unc'])) + 2.5 * np.log10(np.array(result_tab['flux_fit']))
    stars_tbl[f'diff'] = stars_tbl[f'{f}mag'] -  stars_tbl[f'mag_calib']

    # Correct the uncertainty by the error in the zp
    stars_tbl['mag_unc'] = np.abs(stars_tbl['mag_unc']) + np.std(stars_tbl[f'diff'])/np.sqrt(len(stars_tbl))
    stars_tbl['mag_unc_corr'] = np.abs(stars_tbl['mag_unc']) + np.std(stars_tbl[f'diff'])/np.sqrt(len(stars_tbl))
    print(stars_tbl)
    print(f"Standard deviation from catalog: {np.std(np.abs(stars_tbl[f'diff']))}")

    # Target table
    result_tab_target["mag_calib"] =  -2.5 * np.log10(np.array(result_tab_target['flux_fit'])) + diff
    result_tab_target["mag_unc"] = -2.5 * np.log10(np.array(result_tab_target['flux_fit']-result_tab_target['flux_unc'])) + 2.5 * np.log10(np.array(result_tab_target['flux_fit']))
    result_tab_target['mag_unc_corr'] = np.abs(result_tab_target['mag_unc']) + np.std(stars_tbl[f'diff'])/np.sqrt(len(stars_tbl))

    # Add possible upper limit
    if (result_tab_target['flux_fit']/result_tab_target['flux_unc']) < 5:
        result_tab_target['mag_lim'] = -2.5 * np.log10(np.array(result_tab_target['flux_unc']*7)) + diff
    print(result_tab_target)

    plt.close()
    plt.plot(np.arange(len(stars_tbl[f'diff'])), np.abs(stars_tbl[f'diff']), 'ro')
    plt.plot([0,len(stars_tbl)], 2*np.array([np.std(np.abs(stars_tbl[f'diff'])), np.std(np.abs(stars_tbl[f'diff']))]), 'b--')
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

        # Query the PS1 catalog
        with open(f"stars_{f}.reg", "w") as tf:
            tf.write("# Region file format: DS9 version 4.1 \n")
            tf.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n")
            tf.write("fk5 \n")
            tf.write("\n")

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
            with open(f"stars_{f}.reg", "a") as tf:
                for l in result[0]:
                    tf.write(f"circle({l['RAJ2000']}, {l['DEJ2000']}, 1.5" + '") \n')

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

    # Coordinates for photometry
    ra, dec = 316.2947045, 12.8235382
    for f in ["i"]:
        if f == "g":
            filename = "ztf21aaquyjp_g.fits"
            sg_thresh = 0.8
            sex_thresh = 10
        if f == "r":
            filename = "ztf21aaquyjp_r.fits"
            sg_thresh = 0.5
            sex_thresh = 80
        if f == "i":
            filename = "ztf21aaquyjp_i.fits"
            sg_thresh = 0.5
            sex_thresh = 80

        #selected_stars = calib_sextractor(filename, doShow=False)
        stars_for_psf = ascii.read('stars.csv', format='csv')
        #stars_for_psf = [
        #                 (716.57218, 1133.473),
        #                 (1542.5084, 371.23749),
        #                 (828.18237, 342.6856),
        #                 (1640.2606, 1595.5735)
        #                 #(42.164519, 12.158084)
        #                 ]
 

        do_forced(ra, dec, filename, selected_stars=stars_for_psf, use_calib_stars=False)
