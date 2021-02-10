# Author: Igor Andreoni
# email: andreoni@caltech.edu

import glob

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1', 'Yes', 'True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', 'No', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Crop DBSP spectra, \
and plot them up. If no option is given, all spectra in the \
format ./ZTF*fits will be processed')

    parser.add_argument("--doPlot",  action="store_true",
                        default=False, help='Plot up the spectra')
    parser.add_argument('-n', dest='names', nargs='+', required=False,
                        help='Names of the spectra; if not provided, \
all spectra in the format ./ZTF*fits will be processed', default=None)
    parser.add_argument('--suffix', dest='out_suffix', type=str,
                        required=False, help="suffix for the output; \
default = '_crop.txt'", default="_crop.txt")

    args = parser.parse_args()

    # Input files
    if args.names is None:
        list_files = glob.glob("ZTF*fits")
    else:
        list_files = args.names

    print("The spectra currently range from 2,788A to 10,805A.")
    docrop = str2bool(input("Would you like to crop the spectra? \n"))

    # Define the ranges to crop the spectra
    if docrop is True:
        print("The default crop is from 3,200A to 10,000A and \
the interval [5190, 5572] is also removed.")
        default_crop = str2bool(input("Are you happy with this? \n"))
        if default_crop is True:
            ranges = [(3200, 5190), (5572, 10000)]
        else:
            happy = False
            while happy is False:
                nranges = int(input("Desired number of ranges: \n"))
                ranges = []
                for nr in np.arange(nranges):
                    raw_range = input(f"Range {nr+1} (two numbers, \
comma-separated): \n").split(",")
                    new_range = tuple([int(s) for s in raw_range])
                    ranges.append(new_range)
                print(f"Your rages for cropping are {ranges}")
                happy = str2bool(input(f"Are you happy with these ranges?\n"))
    else:
        ranges = None

    # Iterate over the spectra
    for filename in list_files:
        hdulist = fits.open(filename)
        # No cropping
        if docrop is False:
            wav = [p[0] for p in hdulist[1].data]
            flux = [p[1] for p in hdulist[1].data]
            fluxerr = [p[2] for p in hdulist[1].data]
        else:
            wav, flux, fluxerr = [], [], []
            for r in ranges:
                wav_crop = [p[0] for p in hdulist[1].data 
                            if (p[0] > r[0] and p[0] < r[1])]
                flux_crop = [p[1] for p in hdulist[1].data 
                             if (p[0] > r[0] and p[0] < r[1])]
                fluxerr_crop = [p[2] for p in hdulist[1].data 
                                if (p[0] > r[0] and p[0] < r[1])]
                wav += wav_crop
                flux += flux_crop
                fluxerr += fluxerr_crop

        # Plot
        if args.doPlot is True:
            plt.errorbar(wav, flux, yerr=fluxerr, color='grey', alpha=0.5)
            plt.plot(wav, flux, 'b')
            plt.xlabel("Wavelength")
            plt.ylabel("Flux")
            plt.show()
            plt.close()

        # Write the chopped spectra
        with open(filename.replace(".fits", args.out_suffix), 'w') as out:
            for w, f, fe in zip(wav, flux, fluxerr):
                out.write(f"{w}, {f}, {fe}\n")
