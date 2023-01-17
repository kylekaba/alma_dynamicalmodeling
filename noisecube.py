def noisecube(pbcordata,filename,paramfile,primbeam=None):
    
    """ 
    Purpose
    -------
    This program will create two noise cubes, 1 that takes into account the primary beam correction and 1 that
    doesn't, using an ALMA primary-beam corrected data cube and the primary beam itself. 
    
    Input parameters
    ----------------
    
    pbcordata: String of the file name of the original ALMA data cube that the noise will be calculated from.
    If the cube's name is 'cube.fits', then type that string (with quotes) into the first argument.
    
    paramfile: .txt parameter file that contains the information regarding the spatial regions the noise will 
    be calculated on and determines its output size. The required parameters that need to be listed are as follows:
        
            xlo_data: Leftmost x position of the spatial fit region in the data cube 
                      as determined by the user.
            
            xhi_data: Rightmost x position of the spatial fit region in the data cube 
                      as determined by the user.
                      
            ylo_data: Bottommost y position of the spatial fit region in the data cube 
                      as determined by the user.
            
            
            yhi_data: Topmost y position of the spatial fit region in the data cube 
                      as determined by the user.
                      
            xlo_noise: Leftmost x position of the noise region in the data cube 
                       as determined by the user.
            
            xhi_noise: Rightmost x position of the noise region in the data cube 
                       as determined by the user.
                       
            ylo_noise: Bottommost y position of the noise region in the data cube 
                       as determined by the user.
                       
            yhi_noise: Topmost y position of the noise region in the data cube 
                       as determined by the user.
                       
            rebin: Rebinning factor that you want to downsample the data to.
            
    OPTIONAL INPUTS:
            prim_beam: String of the primary beam data cube 

            
    Written by Kyle K. M. Kabasares
    """
    
    # Import the ALMA data and correctly orient it.
    hdul = fits.open(pbcordata)
    data = hdul[0].data
    data = data[0,:,:,:]
    data = swap_cube_axes(data)
    
    # Create the non-primary-beam corrected data cube by multiplying each slice of the .pbcor cube
    # with each slice of the primary beam cube
    
    dx = np.size(data,1)
    dy = np.size(data,0)
    dz = np.size(data,2)
    flatcube = np.zeros((data.shape))
    
    # Identify the square region that the dynamical models will be built on and the region the noise will be 
    # extracted on
    
    FILE = open(paramfile,"r")
    fixedparameters = defaultdict(str)
    for line in FILE:
        paramval = line.strip().split('=')
        fixedparameters[paramval[0].strip()] = paramval[1].strip()

    xlo_data = int(fixedparameters['xlo_data'])
    xhi_data = int(fixedparameters['xhi_data'])
    ylo_data = int(fixedparameters['ylo_data'])
    yhi_data = int(fixedparameters['yhi_data'])
    
    xlo_noise = int(fixedparameters['xlo_noise'])
    xhi_noise = int(fixedparameters['xhi_noise'])
    ylo_noise = int(fixedparameters['ylo_noise']) 
    yhi_noise = int(fixedparameters['yhi_noise']) 
    
    rebin = int(fixedparameters['rebin'])  
    # Identify the spatial region the data is extracted on by creating a spatial array 
    # Make the map 0 everywhere except in the region where the data is taken 
    # This will be rebinned to produce a noise cube on the same scale 
    # As the final ALMA data model in the dynamical program 
    
    # Rebin the data_map array by a factor of rebin
    # Noise region 
    noise_map = np.zeros((dy,dx))
    noise_map[ylo_noise:yhi_noise,xlo_noise:xhi_noise] = 1
    
    # Rebin this noise map by the factor of the rebin scale
    noise_maprebin = block_reduce(noise_map,rebin,np.mean)
    
    # This identifies the location of where we want to extract the noise on the rebinned scale
    noise_mapindicator = np.where(noise_maprebin == 1)
    noise_x = noise_mapindicator[1]
    noise_y = noise_mapindicator[0]
    
    # for loop that will iterate the process of creating 
    
    data_rebin = np.zeros((int(dx/rebin),int(dy/rebin),dz))
    noise_vector_pbcor = np.zeros(dz)
    
    # Create the noise cubes using a for loop by multiplying the measurement of the standard deviation of the 
    # noise on both the flat cube scale and the .pbcor scale by the ones arrays I've generated.
    
    # Pre-allocate the final noise cube arrays
    data_subset = data[ylo_data:yhi_data,xlo_data:xhi_data,:]
    
    # Create a ones cube that will be used to store the noise information per slice on the rebinned scale 
    subdx = np.size(data_subset,1)
    subdy = np.size(data_subset,0)
    
    # This noise cube is on the scale of the final rebinned model cube
    uniformnoise = np.ones((int(subdy/rebin),int(subdx/rebin),dz))
    
    for k in range(dz):
        
        # Rebin the primary beam corrected data cube and the non-corrected primary beam cube
        data_rebin[:,:,k] = block_reduce(data[:,:,k],rebin,np.mean)
   
        # Take a noise measurement from the primary beam corrected cube
        noise_vector_pbcor[k] = np.std(data_rebin[noise_y,noise_x,k])
        
        # Store the uniform noise in a cube of the correct shape
        uniformnoise[:,:,k] = noise_vector_pbcor[k]*uniformnoise[:,:,k]
    
    # Save the cubes in a .fits file
    hdul_staticnoise = fits.PrimaryHDU(np.swapaxes(uniformnoise,0,2))
    filename1 = filename + '.fits'
    hdul_staticnoise.writeto(filename1,overwrite=True)
    
    if primbeam is not None:
        
    # Import the primary beam cube and correctly orient it.
        hdul2 = fits.open(primbeam)
        pbeam = hdul2[0].data
        pbeam = pbeam[0,:,:,:]
        pbeam = swap_cube_axes(pbeam)

        
        pbeam_rebin = np.zeros((int(subdx/rebin),int(subdy/rebin),dz))
        flatcube_rebin = np.zeros((int(dx/rebin),int(dy/rebin),dz))
        noise_vector_uniform = np.zeros(dz)
        uniformnoise = np.ones((int(subdy/rebin),int(subdx/rebin),dz)) 
        
    # Create the flat-cube from the .pbcor data cube and the primary beam data cube.
        for k in range(dz):
            flatcube[:,:,k] = data[:,:,k]*pbeam[:,:,k]
        
    # This noise cube is on the scale of the final rebinned model cube
        varyingnoise = np.ones((int(subdy/rebin),int(subdx/rebin),dz))
        
        flatcube_final = np.ones((int(subdy/rebin),int(subdx/rebin),dz))
        for k in range(dz):
        # Rebin the primary beam corrected data cube, primary beam itself, and the non-corrected primary beam cube
            data_rebin[:,:,k] = block_reduce(data[:,:,k],rebin,np.mean)
            pbeam_rebin[:,:,k] = block_reduce(pbeam[ylo_data:yhi_data,xlo_data:xhi_data,k],rebin,np.mean)
            flatcube_rebin[:,:,k] = block_reduce(flatcube[:,:,k],rebin,np.mean)
        
        # Take a noise measurement from the flat (non primary beam corrected) data cube and the primary beam 
        # corrected cube
            noise_vector_uniform[k] = np.std(flatcube_rebin[noise_y,noise_x,k])
            noise_vector_pbcor[k] = np.std(data_rebin[noise_y,noise_x,k])
        
        # Store the flat noise and the varying noise in cubes of the correct shape
            uniformnoise[:,:,k] = noise_vector_pbcor[k]*uniformnoise[:,:,k]
            varyingnoise[:,:,k] = noise_vector_uniform[k]*(flatcube_final[:,:,k]/pbeam_rebin[:,:,k])
          
        # Save the cubes in a .fits file
        hdul_staticnoise = fits.PrimaryHDU(np.swapaxes(uniformnoise,0,2))
        filename1 = filename + 'uniformcube' + '.fits'
        hdul_staticnoise.writeto(filename1,overwrite=True)
    
        hdul_varyingnoise = fits.PrimaryHDU(np.swapaxes(varyingnoise,0,2))
        filename2 = filename + 'varyingcube' + '.fits'
        hdul_varyingnoise.writeto(filename2,overwrite=True)
    
        return uniformnoise,varyingnoise
    return uniformnoise