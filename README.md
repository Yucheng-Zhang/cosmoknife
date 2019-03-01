# Jackknife

(Need update)

- Make jackknife regions based on the weights of the random points.
  - First sort and cut into pieces in the RA direction.
  - Then for each RA piece, sort and cut into sub-pieces in the DEC direction.

- Mark the jackknife regions.
  - Bounds in RA and DEC for each region: [RA_min, RA_max, DEC_min, DEC_max].
  - Healpix map may also be used to mark the angular jackknife regions.

- For now, input random and data file should have four columns: `RA, DEC, redshift, weight`.
- Output file have one more columns with jackknife regions label, starting from `0`. (New: You can specify with jk0 option.)

## Usage
- `python jack.py --help` to see all options.
- Example
```bash
#!/bin/bash

frand=rand_DR14_QSO_S_radec_z.dat # random file
fdata=data_DR14_QSO_S_radec_z.dat # data file

n_jk=20 # number of jackknife regions
n_ra=5 # number of slices in the RA direction

f_bounds=out_jk_bounds.dat # bounds for jk regions

fo_rand=rand_DR14_QSO_S_rdzwj.dat # output random file
fo_data=data_DR14_QSO_S_rdzwj.dat # output data file

python jack.py -rand $frand -data $fdata -forand $fo_rand -fodata $fo_data -njr $n_jk -nra $n_ra -fbounds $f_bounds

```
