'''
Make jackknife regions based on the mask.
'''
import numpy as np
import healpy as hp
from kniferand import cut_in_dec, cut_in_ra, make_jk_map
import miscfuncs as mf


def knife(data, njr, nra, nside, rra):
    '''Knife function. data includes 3 columns [ra, dec, mask].'''
    w_total = np.sum(data[:, 2])
    w_dec = w_total / njr  # weight for final jackknife regions

    # make weights for RA pieces
    if njr % nra == 0:
        w_ra = np.full(nra, w_total/nra)
    else:
        res = njr % (nra-1)
        ndec = (njr - res) / (nra - 1)
        w_ra = np.full(nra - 1, ndec * w_dec)
        w_ra = np.append(w_ra, res * w_dec)

    d_ra = cut_in_ra(data, w_ra, rra)

    d_dec = cut_in_dec(d_ra, w_dec)

    jk_map = make_jk_map(d_dec, nside)

    return jk_map


def main_mask(args):
    '''Main function for jackknife with mask.
    Value on each pixel should be in range [0, 1] or hp.UNSEEN+(0,1].'''
    print('>> Loading mask: {}'.format(args.mask))
    mask = hp.read_map(args.mask)
    nside = hp.get_nside(mask)
    npix = hp.nside2npix(nside)
    ipix = np.array([i for i in range(npix)])

    # cut off the pixels with value 0 or UNSEEN
    idx = np.where((mask != 0.) & (mask != hp.UNSEEN))
    mask, ipix = mask[idx], ipix[idx]

    theta, phi = hp.pix2ang(nside, ipix)
    ra, dec = mf.get_ra_dec(theta, phi)
    data = np.column_stack((ra, dec, mask))

    jk_map = knife(data, args.njr, args.nra, nside, args.rra)

    if args.plotmap:
        mf.plot_jk_map(jk_map, shuffle=args.sf, njr=args.njr)

    if args.fmap != '':
        mf.save_jk_map(jk_map, args.fmap)
