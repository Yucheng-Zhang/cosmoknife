import argparse
import sys
from kniferand import main_knife
from knifemask import main_mask

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Make jackknife regions.')

    #--- General options ---#
    parser.add_argument('-me', type=str, default='',
                        help='[rand, mask]')
    parser.add_argument('-njr', type=int, default=20,
                        help='Number of jackknife regions.')
    parser.add_argument('-nra', type=int, default=5,
                        help='Number of slices in RA.')
    parser.add_argument('-rra', type=float, default=0.,
                        help='Rotation in RA if cross 0.')

    #--- Common options ---#
    parser.add_argument('-fmap', type=str, default='',
                        help='Output jackknife regions in Healpix map.')
    parser.add_argument('-plotmap', type=int, default=1,
                        help='Plot jackknife healpix map.')
    parser.add_argument('-sf', type=int, default=0,
                        help='Shuffle(1) or not(0) when plotting.')

    #--- For jackknife with randoms ---#
    parser.add_argument('-rand', type=str, default='',
                        help='Random file, with columns: RA, DEC, Z, weight')

    parser.add_argument('-nside', type=int, default=1024,
                        help='Jackknife regions will be a Healpix map.')

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

    parser.add_argument('-jk0', type=int, default=0,
                        help='Initial label of jackknife regions. \
                            All the labels will be [jk0, jk0+1, ..., jk0+njr-1]')
    parser.add_argument('-lb', type=int, default=1,
                        help='0: just get jk bounds with random; \
                              1: 0 + label data and random;')

    parser.add_argument('-bdf', type=str, default='',
                        help='Input jackknife bounds, if provied, \
                            will do the label directly')

    parser.add_argument('-rand_lbed', type=str, default='', help='')

    #--- For jackknife with mask ---#
    parser.add_argument('-mask', type=str, default='',
                        help='Input mask file for jackknife.')

    args = parser.parse_args()

if __name__ == "__main__":

    if args.me == 'rand':
        main_knife(args)

    elif args.me == 'mask':
        main_mask(args)

    else:
        sys.exit('Wrong -me option.')
