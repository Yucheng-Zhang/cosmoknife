'''
Label points with jackknife regions.
'''
import utils
import numpy as np
import healpy as hp


def cat_jkl(jkl, me=None):
    '''Have a look at the jackknife labels.'''
    n_tot = len(jkl)
    n_lost = np.count_nonzero(jkl == -1)
    if me is not None:
        print('-- Label info, method: {}'.format(me))
    print('-- # total points: {0:d}'.format(n_tot))
    print('-- # points not covered in jk regions: {0:d}'.format(n_lost))
    print('-- percent: {0:f} %'.format(100. * n_lost / n_tot))

    return n_lost


def rm_lost_points(data):
    '''Remove the points not covered in jackknife regions.'''
    return data[data[:, -1] != -1]


def save_labeled_data(data, fn):
    '''Save the output data with jk label.'''
    header = 'RA   DEC   redshift   weight   jackknife'
    fmt = '% .15e   % .15e   % .15e   % .15e   %8d'
    np.savetxt(fn, data, fmt=fmt, header=header)
    print(':: Data written to file: {}'.format(fn))


def label_w_map(data, jkr):
    '''Label data points with jackknife regions given in Healpix map.'''
    theta_gal, phi_gal = utils.get_theta_phi(data[:, 0], data[:, 1])

    nside = hp.npix2nside(len(jkr))
    ipix = hp.ang2pix(nside, theta_gal, phi_gal)  # pixel number for each point

    jkl = jkr[ipix]
    np.where(jkl != hp.UNSEEN, jkl, -1)  # label the points not covered with -1

    n_lost = cat_jkl(jkl, me='map')

    # Jackknife label will be added to the last column of data.
    data = np.column_stack((data, jkl))

    # get rid of lost points
    if n_lost > 0:
        print('>> Removing points not covered in jackknife regions')
        data = rm_lost_points(data)

    return data


def label_w_bounds(data, jkr):
    '''Label data points with jackknife regions given in bounds.'''
    jkl = np.zeros(len(data[:, 0])) - 1
    for i, p in enumerate(data):
        for j, r in enumerate(jkr):
            if p[1] >= r[2] and p[1] <= r[3]:  # DEC
                if r[0] < r[1]:
                    if p[0] >= r[0] and p[0] <= r[1]:
                        jkl[i] = j
                        break
                else:  # regions crossing zero
                    if p[0] >= r[0] or p[0] <= r[1]:
                        jkl[i] = j
                        break

    n_lost = cat_jkl(jkl, me='bounds')

    # Jackknife label will be added to the last column of data.
    data = np.column_stack((data, jkl))

    # get rid of lost points
    if n_lost > 0:
        print('>> Removing points not covered in jackknife regions')
        data = rm_lost_points(data)

    return data
