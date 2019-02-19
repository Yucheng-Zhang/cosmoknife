import argparse
import knife
import pandas as pd
import numpy as np
import label


def load_data_pd(fn, tp=''):
    '''Load data file.'''
    print('>> Loading data: {}'.format(fn))
    tb = pd.read_table(fn, delim_whitespace=True, comment='#', header=None)
    tb = tb.values
    if tp == 'knife':
        # RA, DEC, weight
        return np.column_stack((tb[:, 0], tb[:, 1], tb[:, 3]))
    else:
        return tb


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
    parser.add_argument('-plotmap', type=bool, default=True,
                        help='Plot jackknife healpix map.')

    parser.add_argument('-fmap', type=str, default='out_jk_map.fits',
                        help='Output jackknife regions in Healpix map.')
    parser.add_argument('-fbounds', type=str, default='out_jk_bounds.dat',
                        help='Output jackknife regions in bounds: \
                            [ra_min, ra_max, dec_min, dec_max]')

    parser.add_argument('-data', type=str, default='',
                        help='Data file, with columns: RA, DEC, Z, weight')
    parser.add_argument('-fodata', type=str, default='out_data_jk.dat',
                        help='Output data file.')
    parser.add_argument('-forand', type=str, default='out_rand_jk.dat',
                        help='Output random file.')

    parser.add_argument('-tp', type=str, default='bounds',
                        help='map or bounds')

    args = parser.parse_args()

if __name__ == "__main__":

    print('====== Making jackknife regions ======')
    rand = load_data_pd(args.rand, tp='knife')
    jk_map, jk_bounds = knife.knife(rand, args.njr, args.nra, args.nside,
                                    fmap=args.fmap, fbounds=args.fbounds,
                                    plot=args.plotmap)

    print('====== Labeling data points ======')
    data = load_data_pd(args.data)
    label.label(data, jk_bounds, tp=args.tp, f_data=args.fodata)

    print('====== Labeling random points ======')
    data = load_data_pd(args.rand)
    label.label(data, jk_bounds, tp=args.tp, f_data=args.forand)
