'''
Label the jackknife region for each point.
'''
import numpy as np
import healpy as hp
import argparse
import miscfuncs as mf
import sys


def label_map(data, jkr, ll):
    '''Get the label for data with jackknife regions marked in Healpix map.'''
    jkl = np.zeros(len(data[:, 0])) - 1.

    rot = hp.Rotator(coord=['C', 'G'])
    nside = hp.npix2nside(len(jkr))
    # convert to theta, phi used by default in Healpy
    theta_equ, phi_equ = np.deg2rad(90.-data[:, 1]), np.deg2rad(data[:, 0])
    theta_gal, phi_gal = rot(theta_equ, phi_equ)
    ipix = hp.ang2pix(nside, theta_gal, phi_gal)  # pixel number for each point

    jkl = ll[jkr[ipix]]

    return jkl


def cat_jkl(jkl, me):
    '''Have a look at the jackknife labels.'''
    n_tot = len(jkl)
    n_lost = np.count_nonzero(jkl == -1)
    print('-- Label info, method: {}'.format(me))
    print('-- # total points: {0:d}'.format(n_tot))
    print('-- # points not covered in jk regions: {0:d}'.format(n_lost))
    print('-- percent: {0:f} %'.format(100. * n_lost / n_tot))

    return n_lost


def label_bounds(data, jkr, ll):
    '''Get the label for data with jackknife regions marked in bounds.'''
    jkl = np.zeros(len(data[:, 0])) - 1.
    for i, p in enumerate(data):
        for j, r in enumerate(jkr):
            if p[1] >= r[2] and p[1] <= r[3]:  # DEC
                if r[0] < r[1]:
                    if p[0] >= r[0] and p[0] <= r[1]:
                        jkl[i] = ll[j]
                        break
                else:  # regions crossing zero
                    if p[0] >= r[0] or p[0] <= r[1]:
                        jkl[i] = ll[j]
                        break

    return jkl


def save_data(data, fn):
    '''Save the output data with jk label.'''
    header = 'RA   DEC   redshift   weight   jackknife'
    np.savetxt(fn, data, header=header)
    print(':: Data written to file: {}'.format(fn))


def rm_lost_points(data):
    '''Remove the points not covered in jackknife regions.'''
    return data[data[:, 4] != -1]


def label(data, jkr, tp, f_data='', jk0=0):
    '''Main function.'''
    # make label list
    label_list = [jk0+i for i in range(len(jkr))]

    if tp == 'map':
        print('>> Labeling with map {}'.format(jkr))
        jkr = hp.read_map(jkr)
        jkl = label_map(data, jkr, label_list)
        n_lost = cat_jkl(jkl, tp)

    elif tp == 'bounds':
        print('>> Labeling with bounds {}'.format(jkr))
        jkr = np.loadtxt(jkr)
        jkl = label_bounds(data, jkr, label_list)
        n_lost = cat_jkl(jkl, tp)

    else:
        sys.exit('Wrong -tp.')

    # Jackknife label will be added to the last column of data.
    data = np.column_stack((data, jkl))

    # get rid of lost points
    if n_lost > 0:
        print('>> Removing points not covered in jackknife regions')
        data = rm_lost_points(data)

    if f_data != '':
        save_data(data, f_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Label the jackknife region for each point.')

    parser.add_argument('-data', type=str, default='',
                        help='Data file with points to be labeled.')
    parser.add_argument('-tp', type=str, default='',
                        help='map or bounds')
    parser.add_argument('-jkr', type=str, default='',
                        help='file with jackknife regions.')
    parser.add_argument('-fo', type=str, default='', help='output file.')

    args = parser.parse_args()

    print('>> Loading file: {}'.format(args.data))
    data = mf.load_data_pd(args.data)

    if args.tp == 'map':
        print('>> Labeling with map {}'.format(jkr))
        jkr = hp.read_map(jkr)

    elif args.tp == 'bounds':
        print('>> Labeling with bounds {}'.format(jkr))
        jkr = np.loadtxt(jkr)

    else:
        sys.exit('Wrong -tp.')

    label(data, jkr, args.tp, f_data=args.fo)
