def freq_axis(data_cube):
      # Open and extract the data from the input data cube 
    hdu = fits.open(data_cube)
    hdu.info()
    data = hdu[0].data
    
    # Speed of light in km/s
    c = 299792.458
    # Pick out the 3 dimensions of the data cube and swap axes to have the axes ordered 
    # as x,y,z
    data = data[0,:,:,:]
    #data[data < 0] = 0
    # At this point, the order of the axes is: (z,x,y). 
    # We have to interchange z with y to get (y,x,z), and then switch x and y to get it in regular (x,y,z) -> (RA,dec,freq)
    data = np.swapaxes(np.swapaxes(data,0,2),0,1)
    

    # Dimensions of each axis in the data cube 
    dx = int(np.size(data,0))
    dy = int(np.size(data,1))
    dz = int(np.size(data,2))
    
    restfreq = hdu[0].header['RESTFRQ']/1e9
    f_init = hdu[0].header['CRVAL3']
    f_spacing = hdu[0].header['CDELT3'] # frequency step in the data cube
    f_final  = (f_init+(dz*f_spacing))

    f_range = np.arange(f_init,f_final,f_spacing)/1e9
    
#     # Generate velocity grid using optical definition of radial velocity 
#     v_range_optical = c*(1-(restfreq/f_range))
#     # Generate velocity grid with the radio definition of radial velocity
#     v_range_radio = c*((restfreq - f_range))/(restfreq)

    return f_range