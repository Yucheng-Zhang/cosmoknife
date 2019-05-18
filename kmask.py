'''
Make jackknife Healpix masks.
'''
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from . import utils


def jk_masks_w_bounds(mask, bound, froot, test=True):
    '''Make jackknife masks with bounds.'''
    nside = hp.get_nside(mask)
    npix = hp.nside2npix(nside)

    ipix = np.arange(npix)
    # get RA, DEC for all pixels
    theta_gal, phi_gal = hp.pix2ang(nside, ipix)
    ra, dec = utils.get_ra_dec(theta_gal, phi_gal)

    idx = np.where(mask > 0.)
    dec, ra, ipix = dec[idx], ra[idx], ipix[idx]
    npix_c = len(dec)

    if test:
        test_mask = np.ones(npix)

    for j, r in enumerate(bound):  # loop over bounds
        jk_mask = np.copy(mask)

        for i in range(npix_c):  # loop over pixels covered
            if dec[i] >= r[2] and dec[i] <= r[3]:  # DEC
                if r[0] < r[1]:
                    if ra[i] >= r[0] and ra[i] <= r[1]:
                        jk_mask[ipix[i]] = 0.  # set pixel in the jk region to zero
                else:  # regions corssing zero
                    if ra[i] >= r[0] or ra[i] <= r[1]:
                        jk_mask[ipix[i]] = 0.

        fn = froot + '_jk_{0:d}.fits'.format(j)
        hp.write_map(fn, jk_mask)
        print('>> jk mask {0:d} : {1:s}'.format(j, fn))
        if test:
            test_mask = test_mask * jk_mask

    if test:
        hp.mollview(test_mask, coord='GC')
        plt.show()
