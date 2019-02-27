'''
Label the jackknife region for each point.
'''
import numpy as np
import healpy as hp


def label_map(data, jkr, ll):
    '''Get the label for data with jackknife regions marked in Healpix map.'''
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
            if p[0] >= r[0] and p[0] <= r[1] and p[1] >= r[2] and p[1] <= r[3]:
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


def label(data, jkr, tp='bounds', f_data='', jk0=0):
    '''Main function.'''
    # make label list
    label_list = [jk0+i for i in range(len(jkr))]

    if tp == 'map':
        print('>> Labeling with map')
        jkl = label_map(data, jkr, label_list)
        n_lost = cat_jkl(jkl, tp)

    elif tp == 'bounds':
        print('>> Labeling with bounds')
        jkl = label_bounds(data, jkr, label_list)
        n_lost = cat_jkl(jkl, tp)

    # Jackknife label will be added to the last column of data.
    data = np.column_stack((data, jkl))

    # get rid of lost points
    if n_lost > 0:
        print('>> Removing points not covered in jackknife regions')
        data = rm_lost_points(data)

    if f_data != '':
        save_data(data, f_data)
