# Author: Igor Andreoni
# To install dl: pip install --ignore-installed --no-cache-dir astro-datalab

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

from dl import queryClient as qc
import dl


def query_coords_ls(ra, dec, radius_deg,
                    catalog='ls_dr10',
                    outfile="datalab_query.csv"):
    '''Query datalab'''

    #Crossmatch with tractor table
    query = qc.query(sql=f"SELECT objid, ra, dec, type, \
mag_g, flux_g, flux_ivar_g, \
mag_r, flux_r, flux_ivar_r, \
mag_i, flux_i, flux_ivar_i, \
mag_z, flux_z, flux_ivar_z \
from {catalog}.tractor \
where 't' = Q3C_RADIAL_QUERY(ra, dec, {ra}, {dec}, {radius_deg})")

    # Write the result in CSV format
    with open(outfile, "w") as f:
        f.write(query)


if __name__ == "__main__":
    # WLM coordinates
    ra, dec = 0.49208333, -15.46083333

    # Search radius
    radius_deg = 2.0

    # Output filename (CSV)
    outfile = f"ls_query_ra{ra}_dec{dec}_rad{radius_deg}.csv"
    query_coords_ls(ra, dec, radius_deg, catalog='ls_dr10',
                    outfile=outfile)
