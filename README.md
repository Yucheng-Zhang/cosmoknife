# Jackknife

- Make jackknife regions based on the weights of the random points.
  - First sort and cut into pieces in the RA direction.
  - Then for each RA piece, sort and cut into sub-pieces in the DEC direction.
- Mark the jackknife regions.
  - Healpix map is used to mark the angular jackknife regions.
  - Bounds in RA and DEC for each region: [RA_min, RA_max, DEC_min, DEC_max].