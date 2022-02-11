# Read an observing sequence and write JSON files for DECam
# Author: Igor Andreoni

import json

import numpy as np
from astropy.io import ascii


def get_exp_data(fieldname, ra, dec, exptime, filt,
                 program, propID, seqID, seqnum, seqtot=None,
                 expType='object', comment='', count=1, doWait=False):
    """Create a dictionary for each exposure

    Parameters
    ----------
    fieldname str
        Field name
    ra float
        Right Ascension of the pointing
    dec float
        Declination of the pointing
    exptime float
        Exposure time
    filt str
        Filter
    program str
        Name of the program
    propID str
        Proposal ID
    seqID str
        Name of the sequence
    seqnum int
        Number of this exposure in the sequence
    seqtot int
        Total number of exposures in the sequence
    expType str
        Exposure type (default 'object')
    comment str
        Comments
    count int
        Number of time that the same setup sould be
        repeated (default 1)
    doWait bool
        Any wait time needed? (default False)

    Returns
    -------
    data_exp dict
        Dictionary with relevant data for each exposure
    """
    data_exp = {'object': fieldname,
                'RA': str(ra),
                'dec': str(dec),
                'expType': expType,
                'exptime': str(exptime),
                'filter': filt,
                'comment': comment,
                'count': str(count), # Number of times an exposure is repeated
                'program': program,
                'propid': propID,
                'seqid': seqID,
                'seqnum': str(seqnum), # index of an exposure within a sequence
                'seqtot': str(seqtot), # total number of exposures in some sequence
                'wait': str(doWait)
                }
    return data_exp


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Create JSON for observing \
sequence')

    parser.add_argument(dest='filename',
                        nargs=1, type=str,
                        help='Sequence file name (CSV)')
    parser.add_argument('-oh', '--overhead', dest='overhead',
                        type=float, default=30,
                        help='Overhead between exposures (s)')
    parser.add_argument('-max', '--max-time', dest='max_time',
                        type=float, default=1.,
                        help='Maximim time per sequence (hr)')
    parser.add_argument('-pi', '--principal-investigator',
                        dest='pi', type=str, default='Zhang',
                        help='Program Principal Investigator')
    parser.add_argument('-prog', '--program',
                        dest='program', type=str, default='KNTraP',
                        help='Program name')
    parser.add_argument('-id', '--proposal-id',
                        dest='propID', type=str, default='2022A-679480',
                        help='Proposal ID')
    parser.add_argument('-d', '--directory',
                        dest='outdir', type=str, default='./',
                        help='Path to the directory where the JSON files\
will be written')
    args = parser.parse_args()

    # Read the observing sequence file as a table
    t = ascii.read(args.filename[0], format='csv')
    print(f"Observing sequence found with {len(t)} exposures")

    # Create the sequences, keep track of the time per sequence
    sequence = 1
    tot_time = 0
    tot_time_all = 0
    sequences = {str(sequence): []}
    seqID = f"{args.program}_seq{sequence:02d}"
    seqnum = 0

    print(f"Assuming {args.overhead}s overhead between exposures")
    print(f"(this can be changed by passing the --overhead argument)")

    for l in t:
        fieldname = l['fieldname']
        ra = l['ra']
        dec = l['dec']
        filt = l['filter']
        exptime = l['exptime']
        # Keep track of the total time for all sequences
        tot_time_all += exptime + args.overhead
        # Check that the exposure time is small enough
        if exptime > 3600*args.max_time:
            print("WARNING: exposure time greater than the maximum time \
allowed in the sequence!!")
            print(l)
        tot_time += exptime + args.overhead
        # Same sequence or new sequence?
        if tot_time > 3600*args.max_time:
            print(f"Sequence {sequence} has {seqnum} exposures, \
total time {(tot_time - exptime - args.overhead)/3600:.2f}hr")
            # Initialize a new sequence
            sequence += 1
            sequences[str(sequence)] = []
            tot_time = l['exptime'] + args.overhead
            seqnum = 0
        # Increment the number in the sequence
        seqnum += 1
        data = get_exp_data(fieldname, ra, dec, exptime, filt,
                 args.program, args.propID, seqID, seqnum, seqtot=None,
                 expType='object', comment=f"PI {args.pi}",
                 count=1, doWait=False)
        sequences[str(sequence)].append(data)
    # Verbose for the last sequence
    print(f"Sequence {sequence} has {seqnum} exposures, \
total time {(tot_time - exptime - args.overhead)/3600:.2f}hr")

    # Total time for all sequences
    print(f"Expected total run time for the night: {tot_time_all/3600:.2f}hr")

    # For each sequence..
    for k in sequences.keys():
        # Add seqtot info
        seqtot = np.max([int(d["seqnum"]) for d in sequences[k]])
        for d in sequences[k]:
            d["seqtot"] = str(seqtot)
        # create json files
        j_seq = json.dumps(sequences[k], indent=4)
        out_filename = f"{args.outdir}/{args.program.replace(' ', '_')}_seq{int(k):02d}.json"
        with open(out_filename, "w") as f:
            f.write(j_seq)
