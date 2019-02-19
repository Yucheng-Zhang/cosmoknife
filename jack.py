import argparse
import knife
import pandas as pd
import numpy as np


def load_data_pd(fn):
    '''Load random data file.'''
    print('>> Loading data: {}'.format(fn))
    tb = pd.read_table(fn, delim_whitespace=True, comment='#', header=None)
    rand = tb.values
    # RA, DEC, weight
    return np.column_stack((rand[:, 0], rand[:, 1], rand[:, 3]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make jackknife regions.')
    parser.add_argument('-njr', type=int, default=50,
                        help='Number of jackknife regions.')
    parser.add_argument('-nra', type=int, default=10,
                        help='Number of slices in RA.')
    parser.add_argument('-rand', type=str, default='',
                        help='Random file, with columns: RA, DEC, Z, weight')
    parser.add_argument('-nside', type=int, default=256,
                        help='Jackknife regions will be a Healpix map.')
    parser.add_argument('-save', type=int, default=1,
                        help='Save to file or not.')
    parser.add_argument('-fmap', type=str, default='out_jk_map.fits',
                        help='Output jackknife regions in Healpix map.')
    parser.add_argument('-fbounds', type=str, default='out_jk_bounds.dat',
                        help='Output jackknife regions in bounds: [ra_min, ra_max, dec_min, dec_max]')
    args = parser.parse_args()

if __name__ == "__main__":

    rand = load_data_pd(args.rand)
    jk_map, jk_bounds = knife.knife(rand, args.njr, args.nra, args.nside,
                                    fmap=args.fmap, fbounds=args.fbounds, plot=True)
