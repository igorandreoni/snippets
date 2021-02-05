Collection of useful astronomy snippets. All work with Python 3.

## Upload spectra to Fritz

usage: upload_spectra_fritz.py [-h] [-d] [--date DATE] [--inst INST_ID]
                               ID [ID ...] <br>

Upload spectra to the Fritz marshal, individually or in bulk. <br>

positional arguments:<br>
  ID              Spectra filenames <br>

optional arguments:<br>
  -h, --help      show this help message and exit<br>
  --date DATE     Date of the observations, <br>
                  for example 2020-11-10T00:00:00
                  or 2021-01-20 <br>
  --inst INST_ID  Instrument ID, e.g. inst_id = 3 for DBSP. Instrument IDs can<br>
                  be found here: https://fritz.science/api/instrument<br>
  -d              Use default reducer ID and observer - please customize the
                  code if you want to use this (unrequired) option<br>

**Notes:**
* If the name of the ZTF source is in the filename, it will be recognized automatically; filenames without the ZTF source name are fine, too
* Make sure that 'ZTF' appears only once at most in the filename
* The spectra are assumed to have wavelength in the first column (in angstrom), flux in the second column, possibly flux error in the third column
* The code will automatically let you know under which group IDs the source was saved, but they can be changed interactively when uploading 
* Please update your API before using the code. You can find the API code in your Fritz account page
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
usage: dbsp_crop_spec.py [-h] [--doPlot] [-n NAMES [NAMES ...]] <br>
                         [--suffix OUT_SUFFIX]<br>
<br>
Crop DBSP spectra and plot them up. <b>This is NOT a necessary step <br>
for the upload on Fritz</b>, if FITS files from DBSP_DRP are available. <br>
<br>
If no option is given, all spectra in the format ./ZTF\*.fits <br>will be processed.<br>
<br>
optional arguments:<br>
  -h, --help            show this help message and exit<br>
  --doPlot              Plot up the spectra<br>
  -n NAMES [NAMES ...]<br>
                        Names of the spectra; if not provided, all spectra in <br>
                        the format ./ZTF\*.fits will be processed <br>
  --suffix OUT_SUFFIX   suffix for the output; default = '_crop.txt' <br>

**Example:**

```
python dbsp_crop_spec.py --n ZTF19cool.fits ZTF20lesscool.fits --doPlot
```

**Limitations:** Headers are not yet transferred



## Get galaxies
usage: get_galaxies.py [-h] --ra RA --dec DEC [--r RAD] [--c CATALOG]<br>
                       [--dist-min DIST_MIN] [--dist-max DIST_MAX]<br>
                       [--sep-max SEP_MAX_KPC] [--out OUT]<br>
<br>
Query galaxy catalogs<br>
<br>
optional arguments:<br>
  -h, --help            show this help message and exit<br>
  --ra RA               Right Ascension of the center of the query (deg)<br>
  --dec DEC             Declination of the center of the query (deg)<br>
  --r RAD               Radius of the query (deg), default=1deg<br>
  --c CATALOG           Catalog name. The default is GLADE2.3
                        (VII/281/glade2). Available catalogs with distances &
                        separations: GLADE v2.3 (VII/281/glade2); 6dF DR3 spec
                        (VII/259/spectra). Other catalogs will not have
                        calculated separations. Some examples are: 2MASS
                        extended sources (VII/233/xsc); HYPERLEDA
                        (VII/237/pgc); SDSS DR12 (V/147/sdss12); GAIA S/G
                        class (VII/285/gdr2ext)<br>
  --dist-min DIST_MIN   Minimum distance (Mpc)<br>
  --dist-max DIST_MAX   Maximum distance (Mpc)<br>
  --sep-max SEP_MAX_KPC<br>
                        Maximum projected separation (kpc)<br>
  --out OUT             Output file name (CSV)<br>



**Example:**

```
python get_galaxies.py --ra 219.51950000 --dec -60.34800000 --r 3. --out my_galaxies.csv --dist-max 100
```

**Example** for a generic catalog (say HYPERLEDA, see the help for suggestions):

```
python get_galaxies.py --ra 329.41950000 --dec -80.3580 --r 3.0 --c VII/237/pgc --out galaxies_FRB190711_hyperleda.csv
```
