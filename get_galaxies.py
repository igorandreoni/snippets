# Author: Igor Andreoni
# email: andreoni@caltech.edu


import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astroquery.vizier import Vizier
from astropy.table import Table, vstack
from astropy.cosmology import Planck15 as cosmo


def query_glade(ra, dec, rad=1, dist_min=0, dist_max=np.inf,
                sep_max_kpc=np.inf, catalog='VII/281/glade2'):
    """
    Fetch galaxies from the GLADE catalog

    ---
    Parameters

    ra float
        Right ascension (deg)
    dec float
        Declination (deg)
    rad float
        Search radius (deg)
    dist_min float
        Min distance of the galaxies (Mpc)
    dist_max float
        Max distance of the galaxies (Mpc)
    sep_max_kpc float
        Min projected distance of the galaxies (kpc)
    catalog str
        Exact version of the catalog

    ---
    Return

    glade_select list
        Selected galaxies, all info
    sep_select list
        Separations for the selected galaxies (arcsec)
    dist_kpc_select list
        Separations for the selected galaxies (kpc)

    """

    # Create a SkyCoord object
    coord = SkyCoord(ra=ra, dec=dec,unit=(u.deg, u.deg))
    # Query GLADE via Vizier
    glade = Vizier.query_region(coord, radius=rad*u.deg,
                                   catalog="glade")

    # Empty result
    if len(glade) == 0:
        return None, None, None

    # Check if any of the galaxies found are within an
    # acceptable distance
    glade_select = []
    dist_kpc_select = []
    sep_select = []

    for g in glade[0]:
        # Ignore galaxies too nearby
        if g['Dist'] < dist_min or g['Dist'] > dist_max:
            continue
        sep = SkyCoord(ra=g['RAJ2000']*u.deg,
                       dec=g['DEJ2000']*u.deg).separation(coord)
        # Projected distance (kpc)
        dist_kpc = g['Dist']*(10**3)*np.sin(sep)/np.cos(sep)
        if dist_kpc >= 0 and dist_kpc < sep_max_kpc:
            glade_select.append(g)
            dist_kpc_select.append(dist_kpc)
            sep_select.append(sep)
    # No selected galaxies
    if len(glade_select) == 0:
        return None, None, None
    else:
        return glade_select, sep_select, dist_kpc_select


def query_6df(ra, dec, rad=1, dist_min=0, dist_max=np.inf,
                sep_max_kpc=np.inf, catalog='VII/259/spectra'):
    """
    Fetch galaxies from the 6df catalog

    ---
    Parameters

    ra float
        Right ascension (deg)
    dec float
        Declination (deg)
    rad float
        Search radius (deg)
    dist_min float
        Min distance of the galaxies (Mpc)
    dist_max float
        Max distance of the galaxies (Mpc)
    sep_max_kpc float
        Min projected distance of the galaxies (kpc)
    catalog str
        Exact version of the catalog

    ---
    Return

    gal_select list
        Selected galaxies, all info
    sep_select list
        Separations for the selected galaxies (arcsec)
    dist_kpc_select list
        Separations for the selected galaxies (kpc)

    """

    # Create a SkyCoord object
    coord = SkyCoord(ra=ra, dec=dec,unit=(u.deg, u.deg))
    # Query GLADE via Vizier
    gal = Vizier.query_region(coord, radius=rad*u.deg,
                                   catalog=catalog)

    # Empty result
    if len(gal) == 0:
        return None, None, None

    # Check if any of the galaxies found are within an
    # acceptable distance
    gal_select = []
    dist_kpc_select = []
    sep_select = []

    # Get distances from redshifts
    print("Note: distances are luminosity distances (Mpc) \
using Planck+15 cosmology")

    for g in gal[0]:
        # Ignore galaxies too nearby
        dist_g = cosmo.luminosity_distance(g['z']).value
        if dist_g < dist_min or dist_g > dist_max:
            continue
        ra_deg = Angle(g['RAJ2000'] + " hours").deg
        dec_deg = Angle(g['DEJ2000'] + " degrees").deg
        sep = SkyCoord(ra=ra_deg*u.deg,
                       dec=dec_deg*u.deg).separation(coord)
        # Projected distance (kpc)
        dist_kpc = dist_g*(10**3)*np.sin(sep)/np.cos(sep)
        if dist_kpc >= 0 and dist_kpc < sep_max_kpc:
            gal_select.append(g)
            dist_kpc_select.append(dist_kpc)
            sep_select.append(sep)
    # No selected galaxies
    if len(gal_select) == 0:
        return None, None, None
    else:
        return gal_select, sep_select, dist_kpc_select


def query_generic(ra, dec, rad=1, catalog=None):
    """
    Fetch galaxies from a generic input catalog

    ---
    Parameters

    ra float
        Right ascension (deg)
    dec float
        Declination (deg)
    rad float
        Search radius (deg)
    catalog str
        Exact version of the catalog

    ---
    Return

    gal_select list
        Selected galaxies, all info

    """

    # Create a SkyCoord object
    coord = SkyCoord(ra=ra, dec=dec,unit=(u.deg, u.deg))
    # Query GLADE via Vizier
    gal = Vizier.query_region(coord, radius=rad*u.deg,
                                   catalog=catalog)

    # Empty result
    if len(gal) == 0:
        return None
    else:
        return gal


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Query galaxy catalogs')
    parser.add_argument('--ra', dest='ra', type=float,
                        required=True,
                        help='Right Ascension of the center of \
                        the query (deg)')
    parser.add_argument('--dec', dest='dec', type=float,
                        required=True,
                        help='Declination of the center of the \
                        query (deg)')    
    parser.add_argument('--r', dest='rad', type=float,
                        required=False,
                        help='Radius of the query (deg), \
                        default=1deg', default=1.)
    parser.add_argument('--c', dest='catalog', type=str,
                        required=False, default="glade",
                        help='Catalog name. \
                        The default is GLADE2.3 (VII/281/glade2). \
                        Available catalogs with distances & separations: \n \
                         GLADE v2.3 (VII/281/glade2); \
                         6dF DR3 spec (VII/259/spectra). \
                        Other catalogs will not have calculated separations. \n \
                        Some examples are: \
                          2MASS extended sources (VII/233/xsc); \
                          HYPERLEDA (VII/237/pgc); \
                          SDSS DR12 (V/147/sdss12); \
                          GAIA S/G class (VII/285/gdr2ext)')
    parser.add_argument('--dist-min', dest='dist_min', type=float,
                        required=False, default=0,
                        help='Minimum distance (Mpc)')
    parser.add_argument('--dist-max', dest='dist_max', type=float,
                        required=False, default=np.inf,
                        help='Maximum distance (Mpc)')
    parser.add_argument('--sep-max', dest='sep_max_kpc', type=float,
                        required=False, default=np.inf,
                        help='Maximum projected separation (kpc)')
    parser.add_argument('--out', dest='out', type=str,
                        required=False, default=None,
                        help='Output file name (CSV)')
    args = parser.parse_args()

    # Query the GLADE catalog
    Vizier.ROW_LIMIT = -1
    print("Querying the GLADE catalog")
    print(f"Central coordinates: RA, Dec = {args.ra}, {args.dec} (deg)")
    print(f"Search radius: {args.rad} deg")
    print(f"Minimum distance: {args.dist_min} Mpc")
    print(f"Maximum distance: {args.dist_max} Mpc")
    print(f"Maximum projected separation: {args.sep_max_kpc} kpc")
    print("")

    if args.catalog == 'VII/281/glade2':
        galaxies, sep, dist_kpc = query_glade(args.ra, args.dec,
                                              rad=args.rad,
                                              dist_min=args.dist_min,
                                              dist_max=args.dist_max,
                                              sep_max_kpc=args.sep_max_kpc,
                                              catalog=args.catalog)

    elif args.catalog == 'VII/259/spectra':
        galaxies, sep, dist_kpc = query_6df(args.ra, args.dec,
                                              rad=args.rad,
                                              dist_min=args.dist_min,
                                              dist_max=args.dist_max,
                                              sep_max_kpc=args.sep_max_kpc,
                                              catalog=args.catalog)
    else:
        galaxies = query_generic(args.ra, args.dec,
                                  rad=args.rad, catalog=args.catalog)
        if galaxies is not None:
            galaxies = galaxies[0]
#        print("Un-recognized catalog.")
#        print("Available catalogs:")
#        print("VII/281/glade2")
#        print("VII/259/6dfgs")

    if galaxies is None:
        print("No galaxies found with the given search parameters.")
        print("(double-check: does the input catalog exist in Vizier?)")
        exit()
    else:
        print(f"Found {len(galaxies)} galaxies")
        print("")

    # Create a table with the results
    tbl = vstack(galaxies)
    if args.catalog == 'VII/281/glade2' or args.catalog == 'VII/259/spectra':
        tbl["sep_arcsec"] = sep
        tbl["dist_kpc"] = dist_kpc
    print(tbl)

    if args.out is None:
        out = f"galaxies_{str(args.ra)}_{str(args.dec)}_{args.catalog.split('/')[-1]}.csv"
    else:
        out = args.out
    print(f"The results were written in CSV format in {out}")
    tbl.write(out, format='csv', overwrite=True)
