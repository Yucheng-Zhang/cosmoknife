import argparse
import knife
import label
import sys
import miscfuncs
import numpy as np


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Make jackknife regions.')

    parser.add_argument('-njr', type=int, default=20,
                        help='Number of jackknife regions.')
    parser.add_argument('-nra', type=int, default=5,
                        help='Number of slices in RA.')
    parser.add_argument('-rra', type=float, default=0.,
                        help='Rotation in RA if cross 0.')

    parser.add_argument('-rand', type=str, default='',
                        help='Random file, with columns: RA, DEC, Z, weight')

    parser.add_argument('-nside', type=int, default=256,
                        help='Jackknife regions will be a Healpix map.')
    parser.add_argument('-plotmap', type=int, default=1,
                        help='Plot jackknife healpix map.')

    parser.add_argument('-fmap', type=str, default='',
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

    args = parser.parse_args()


if __name__ == "__main__":

    # make jackknife regions
    if args.bdf == '' and args.rand_lbed == '':
        print('====== Making jackknife regions ======')
        rand = miscfuncs.load_data_pd(args.rand, tp='knife')
        jk_map, jk_bounds = knife.knife(
            rand, args.njr, args.nra, args.nside, args.rra)

        if args.fmap != '':
            miscfuncs.save_jk_map(jk_map, args.fmap)

        if args.plotmap:
            print('-- note: not labeled yet, just demo of regions')
            miscfuncs.plot_jk_map(jk_map, shuffle=True, njr=args.njr)

        if args.fbounds != '':
            miscfuncs.save_jk_bounds(jk_bounds, args.fbounds)

        if args.tp == 'bounds':
            jkr = jk_bounds
        elif args.tp == 'map':
            jkr = jk_map
        else:
            print('>> Error: wrong tp option!')
            sys.exit()

    # load bounds file if provided
    if args.bdf != '' and args.rand_lbed == '':
        print('>> Loading bounds file: {}'.format(args.bdf))
        jkr = np.loadtxt(args.bdf)

    # label data and random points
    if (args.lb == 1 or args.bdf != '') and args.rand_lbed == '':
        print('====== Labeling data points ======')
        data = miscfuncs.load_data_pd(args.data)
        label.label(data, jkr, tp=args.tp,
                    f_data=args.fodata, jk0=args.jk0)

        print('====== Labeling random points ======')
        data = miscfuncs.load_data_pd(args.rand)
        label.label(data, jkr, tp=args.tp,
                    f_data=args.forand, jk0=args.jk0)

    # analyze labeled random points
    if args.rand_lbed != '':
        rand = miscfuncs.load_data_pd(args.rand_lbed)
        miscfuncs.analyze_rand(rand)
