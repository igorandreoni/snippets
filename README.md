Collection of useful astronomy snippets. All work with Python 3.

## Query Legacy Survey DR8 photoz

Query photometric redshifts from Legacy Survey DR8<br>

```
usage: query_photoz_datalab.py [-h] [-r RADIUS] RA, Dec [RA, Dec ...]

positional arguments:
  RA, Dec     RA and Dec (degrees)

optional arguments:
  -h, --help  show this help message and exit
  -r RADIUS   Search radius (arcsc)
```

**Notes:**
* Make sure that the latest version of the datalab package is installed
`pip install --upgrade noaodatalab`
 

**Example:**

```
python query_photoz_datalab.py 199.651874 30.3243492 -r 2
```
```
z_phot_median, z_phot_std, z_phot_l95, ra, dec, type, flux_z, sep_arcsec
0.362533, 0.086208, 0.258176, 199.65187593677, 30.324329766903, EXP, 7.06425, 0.07021754722965218
```

## Find sources in ZTF data

Query Kowalski and obtain a summary table of all sources in a cone <br>
```
usage: cone_search_ztf.py [-h] [-r RADIUS] [--date-start DATE_START]
                          [--date-end DATE_END] [--ndethist NDETHIST_MIN]
                          [--pid PROGRAMID_LIST [PROGRAMID_LIST ...]]
                          [--drb DRB_MIN] [--out OUT]
                          RA, Dec [RA, Dec ...]
```

ZTF alert cone search using Kowalski
```
positional arguments:
  RA, Dec               RA and Dec (degrees)

optional arguments:
  -h, --help            show this help message and exit
  -r RADIUS             Search radius (arcmin)
  --date-start DATE_START
                        Start date of the query, in ISO format. Example:
                        '2017-08-17 12:41:04.4'
  --date-end DATE_END   End date of the query, in ISO format. Example:
                        '2017-08-18 12:00:00.0'
  --ndethist NDETHIST_MIN
                        Minimum number of detections
  --pid PROGRAMID_LIST [PROGRAMID_LIST ...]
                        ZTF program IDs
  --drb DRB_MIN         Minimum drb score
  --out OUT             Output filename: if given, a CSV file will be created
```

**Example:**

```
python cone_search_ztf.py 272.312512 -9.62908 -r 10 --date-start '2021-06-01' --date-end '2021-06-12' --pid 1 --out my_sources.csv
```

## Upload spectra to Fritz

usage: upload_spectra_fritz.py [-h] [-d] [--date DATE] [--inst INST_ID]
                               ID [ID ...] <br>
```
Upload spectra to the Fritz marshal, individually or in bulk.

positional arguments:
  ID              Spectra filenames

optional arguments:
  -h, --help      show this help message and exit
  --date DATE     Date of the observations,
                  for example 2020-11-10T00:00:00
                  or 2021-01-20
  --inst INST_ID  Instrument ID, e.g. inst_id = 3 for DBSP. Instrument IDs can
                  be found here: https://fritz.science/api/instrument
  -d              Use default reducer ID and observer - please customize the
                  code if you want to use this (unrequired) option
```

**Notes:**
* If the name of the ZTF source is in the filename, it will be recognized automatically; filenames without the ZTF source name are fine, too
* Make sure that 'ZTF' appears only once at most in the filename
* The spectra are assumed to have wavelength in the first column (in angstrom), flux in the second column, possibly flux error in the third column
* The code will automatically let you know under which group IDs the source was saved, but they can be changed interactively when uploading 
* Please copy your API token before using the code. You can find the API code in your Fritz account page
* Headers will be uploaded along with the spectra; DBSP_DRP output headers will be merged automatically before the upload

**Useful links:**
* Instrument IDs: https://fritz.science/api/instrument
* User IDs: https://fritz.science/api/user
* Group IDs: https://fritz.science/api/groups
* Skyportal API guide: https://skyportal.io/docs/api.html
* Fritz marshal user's guide https://fritz-marshal.org/doc/

**Examples:**
```
python upload_spectra_fritz.py ZTF*fits --inst 3 --date 2020-10-28T07:00:00
python upload_spectra_fritz.py lris20201017_ZTF20thebest.spec --inst 7 --date 2020-10-17
```

## Crop and display DBSP spectra

```
usage: dbsp_crop_spec.py [-h] [--doPlot] [-n NAMES [NAMES ...]] <br>
                         [--suffix OUT_SUFFIX]
```
Crop DBSP spectra and plot them up. <b>This is NOT a necessary step <br>
for the upload on Fritz</b>, if FITS files from DBSP_DRP are available. <br>
<br>
If no option is given, all spectra in the format ./ZTF\*.fits <br>will be processed.<br>
<br>
```
optional arguments:
  -h, --help            show this help message and exit
  --doPlot              Plot up the spectra
  -n NAMES [NAMES ...]
                        Names of the spectra; if not provided, all spectra in
                        the format ./ZTF\*.fits will be processed
  --suffix OUT_SUFFIX   suffix for the output; default = '_crop.txt'
```

**Example:**

```
python dbsp_crop_spec.py --n ZTF19cool.fits ZTF20lesscool.fits --doPlot
```

**Limitations:** Headers are not yet transferred



## Get galaxies
usage: get_galaxies.py [-h] --ra RA --dec DEC [--r RAD] [--c CATALOG]<br>
                       [--dist-min DIST_MIN] [--dist-max DIST_MAX]<br>
                       [--sep-max SEP_MAX_KPC] [--out OUT]<br>
```
Query galaxy catalogs

optional arguments:
  -h, --help            show this help message and exit
  --ra RA               Right Ascension of the center of the query (deg)
  --dec DEC             Declination of the center of the query (deg)
  --r RAD               Radius of the query (deg), default=1deg
  --c CATALOG           Catalog name. The default is GLADE2.3
                        (VII/281/glade2). Available catalogs with distances &
                        separations: GLADE v2.3 (VII/281/glade2); 6dF DR3 spec
                        (VII/259/spectra). Other catalogs will not have
                        calculated separations. Some examples are: 2MASS
                        extended sources (VII/233/xsc); HYPERLEDA
                        (VII/237/pgc); SDSS DR12 (V/147/sdss12); GAIA S/G
                        class (VII/285/gdr2ext)
  --dist-min DIST_MIN   Minimum distance (Mpc)
  --dist-max DIST_MAX   Maximum distance (Mpc)
  --sep-max SEP_MAX_KPC
                        Maximum projected separation (kpc)
  --out OUT             Output file name (CSV)
```


**Example:**
```
python get_galaxies.py --ra 219.51950000 --dec -60.34800000 --r 3. --out my_galaxies.csv --dist-max 100
```

**Example** for a generic catalog (say HYPERLEDA, see the help for suggestions):

```
python get_galaxies.py --ra 329.41950000 --dec -80.3580 --r 3.0 --c VII/237/pgc --out galaxies_FRB190711_hyperleda.csv
```


## Create JSON

Create JSON files for DECam given an observing sequence table
```
positional arguments:
  filename              Sequence file name (CSV)

optional arguments:
  -h, --help            show this help message and exit
  -oh OVERHEAD, --overhead OVERHEAD
                        Overhead between exposures (s)
  -max MAX_TIME, --max-time MAX_TIME
                        Maximim time per sequence (hr)
  -pi PI, --principal-investigator PI
                        Program Principal Investigator
  -prog PROGRAM, --program PROGRAM
                        Program name
  -id PROPID, --proposal-id PROPID
                        Proposal ID
  -d OUTDIR, --directory OUTDIR
                        Path to the directory where the JSON files will be
                        written
```

The input CSV file must have these columns:
```
fieldname, ra, dec, filter, exptime
CDFS, 52.5, -28.1, g, 140
CDFS, 52.5,-28.1, i, 170
...
```

Example:
```
python create_json.py observing_sequence.csv -max 2 -pi "Albert Einstein" -prog "Fun With Relativity" -id NOAO-1908A --directory JSON_files
```

Example output:
```
Observing sequence found with 148 exposures
Assuming 30s overhead between exposures
(this can be changed by passing the --overhead argument)
Sequence 1 has 43 exposures, total time 1.96hr
Sequence 2 has 39 exposures, total time 2.00hr
Sequence 3 has 39 exposures, total time 2.00hr
Sequence 4 has 27 exposures, total time 1.34hr
Expected total run time for the night: 7.36hr
```

## Get a new transient name

Given a list of transients, a base for the name, and the MJD of a discovery, output the new ordered label for the transient

```
options:
  -h, --help            show this help message and exit
  -n NAMES [NAMES ...], --name NAMES [NAMES ...]
                        <Required> List of transient names e.g. AT2018gfo AT2022cmc
  --mjd MJD             MJD of the new transient detection (default: todays MJD)
  -b BASE, --base BASE  Base name for the transients, e.g. AT
```
**Example:**
```
python nameObject.py -n ZTF23ab ZTF22a ZTF23azz --base ZTF --mjd 60093.55
```
Example output:
```
ZTF23baa
```
