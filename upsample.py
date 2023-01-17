# This is a function designed to take in a square flux map and upscale it to the size of the sub-sampled model velocity grid.
# The flux in each subpixel element is simply the total flux divided by the amount of subpixels within an original ALMA pixel.
# Written by Kyle K. M. Kabasares
# Test commit

def upsample(A,ssf):
    B = np.repeat(A,ssf,axis=0)
    C = np.repeat(B,ssf,axis=1)
    C = C/(ssf**2)
    return C