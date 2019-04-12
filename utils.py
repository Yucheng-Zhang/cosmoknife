'''
Script to use some useful functions.
'''
import argparse
import miscfuncs as mf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Useful funcs for jackknife analysis.')

    parser.add_argument('-func', type=str, default='',
                        help='''Function to use.''')
    #--- merge_masks ---#
    parser.add_argument('-fmasks', type=str, nargs='+',
                        help='Masks to merge.')
    parser.add_argument('-nside', type=int, help='nside of the masks.')
    parser.add_argument('-fo0', type=str, default='merged-mask.fits')
    #--- sep_mask_ra ---#
    parser.add_argument('-fmask', type=str, help='Mask to seperate.')
    parser.add_argument('-RAs', type=float, nargs='+',
                        help='Two seperations in RA.')
    parser.add_argument('-fo1', type=str, nargs='+',
                        default=['sep-mask1.fits', 'sep-mask2.fits'])

    args = parser.parse_args()


if __name__ == "__main__":

    if args.func == 'merge_masks':
        mf.merge_masks(args.fmasks, args.nside, args.fo0)
    elif args.func == 'sep_mask_ra':
        mf.sep_mask_ra(args.fmask, args.RAs, args.fo1)
    else:
        pass
