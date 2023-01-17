def dynamical(freeparams,
              textparamfile,
              galaxyname,
              f_range,
              FluxMap,
              datacube,
              eps,
              method=None,
              counter=None,
              r_int=None,
              vc_st=None,
              r_gas_int=None,
              vc_gas=None,
              modelcube=None,
              fit_region=None,
              SaveFITS=False,
              SaveAny=False,
              verbose=False,
              loud=False):

    """ 
    Purpose
    -------

    This program creates model integrated Gaussian line profiles and optimizes them to ALMA data cubes for the 
    purpose of measuring the masses of supermassicve black holes. The optimization makes use of the open-source
    LMFIT package: https://lmfit.github.io/lmfit-py/. 
    
    Input parameters
    ----------------
    
    freeparams: LMFIT parameters that will be used as initial guesses to the model fit and will be optimized with
    each successive iteration.
    
    textparamfile: .txt parameter file that contains the fixed parameter information for the particular galaxy the
    fit will be optimized 
    
         Required parameters
         
         --------------------
        xi: Leftmost x position of the spatial region that will be used in creating a dynamical model. This region
        is typically a subset of the entire spatial region of the original data cube.
        
        xf: Rightmost x position of the spatial region that will be used in creating a dynamical model. This region
        is typically a subset of the entire spatial region of the original data cube.
        
        yi: Bottommost y position of the spatial region that will be used in creating a dynamical model. This region
        is typically a subset of the entire spatial region of the original data cube.
        
        yf: Topmost y position of the spatial region that will be used in creating a dynamical model. This region
        is typically a subset of the entire spatial region of the original data cube.
        
        *** IMPORTANT NOTE *** 
        This region needs to be square, so xf - xi = yf - yi, and the difference between the two regions needs to be
        evenly divisible by the value of the rebin parameter.
        
        nu_i: First frequency channel with visible emission in the data cube.
        
        nu_f: Last frequency channel with visible emission in the data cube.
        
        N: The amount of radial rings used in the v_extended model. 
        
        r_init: The starting radius (in parsecs) where the v_extended model will place its first ring. This number
        is unimportant if the v_extended model is not being used.
        
        G: Netwon's gravitational constant in units where masses are in solar masses, velocities are in km/s, and radii
        are in parsecs. (0.00430237)
        
        dpc: Distance to the galaxy in parsecs.
        
        cs: Speed of light in km/s
        
        sub: Sub-sampling factor. Recommended minimum value is 2.
        
        rebin: The rebinning factor that will set how much the final model is binned down. i.e. If the original
        data cube size is (dx,dy), the rebinned model has a size of (dx/rebin, dy/rebin)
        
        gridsize: Integer (that needs to be odd) that represents how big the primary beam/PSF image will be on
        each side. One good rule of thumb is to make the length larger than 5x the primary beam's average size 
        (take the major and minor axis average) in units of pixels.
        
        a: Semimajor axis of the ellipse that will be used to define the fitting region where the data and model 
        will be compared on.
        
        q: Axial ratio (b/a) of the minor to major axis ratio of the fitting ellipse.
        
        xc_ellipse: x coordinate of the center of the fitting ellipse
        
        yc_ellipse: y coordinate of the center of the fitting ellipse
        
        GammaEllipse: Position angle of the fitting ellipse. This typically needs to be as close to the gas disk
        position angle as possible.
    
    galaxyname: String extension that will be added to all the output files of the program.
    
    f_range: The frequency axis that the model cube will be built on. This is the same frequency axis as the
    original data cube, and can be computed from using the user-defined function "freq_axis".

    FluxMap: A map of integrated flux that will be used to weight the Gaussian line profiles (set the area under
    their curves). This should have the same spatial scale as the original ALMA data cube.
    
    datacube: String of the file name of the original ALMA data cube that models will be built to be compared with.
    
    eps: The noise cube that is used in the chi-squared calculation of the model optimization. The noise cube needs
    to be the same size as the final rebinned model size. 
    
    OPTIONAL INPUTS:
    ----------------
    
    method: Used to choose how we parameterize the stellar mass/stellar circular velocity profile if no
    Multi-Gaussian expansion is available. Currently, the 2 supported options are "vext" and "power" which 
    indicate using either the v_extended model or a simple 2 parameter power law describing the mass profile: 
    M(r) = c1*r^beta
    
    r_int: Radius (in parsecs) of the stellar circular velocity profile. Needs to be the same
    length as "vc_st"
    
    vc_st: Stellar circular velocity (in km/s) profile typically derived from a Multi-Gaussian (MGE) expansion, which is used
    as a proxy for enclosed stellar mass at a given radius. Needs to be the same length as r_int.
    
    r_gas_int: Radius (in parsecs) of the gas disk circular velocity profile. Needs to be the same length as 
    "vc_gas".
    
    vc_gas: Molecular gas disk circular velocity (in km/s) profile. This is typically derived from observations
    of the molecular gas distribution in the ALMA data.
    
    modelcube: An ALMA model cube created from a previous dynamical modeling run. This typically is only used when
    one wants to perform a Monte Carlo simulation to better understand the statistics of the fitting parameters.
    
    fit_region: A 3D array the same dimensions as the final model in the dynamical modeling process that indicates
    what pixels will be included in the dynamical modeling fit. A value of "1" represents pixels to be included
    whereas a value of "0" represents pixels to be excluded.
    
    SaveFITS: Set to either to True if you wish to save additional .fits files in your working directory (list
    is at the end of the dynamical modeling program), and False if you don't.
    
    SaveAny: Set to either to True if you wish to save all of the key .fits files while running the optimization or
    False if you simply want to optimize. Key fits files include: 
    
    verbose: Setting this to "True" will show all the plots, while "False" hides them.
    
    quiet: Setting this to "True" will hide all the print text displayed in the optimization.

    Written by Kyle K. M. Kabasares 
     """
    time_start = time.time()
    clear_output()    
## Clear the previous output
    ### INPUT GRID (FIXED) PARAMETERS FROM TEXT PARAMETER FILE
    FILE = open(textparamfile,"r")
    parameters = defaultdict(str)
    
    for line in FILE:
        paramval = line.strip().split('=')
        parameters[paramval[0].strip()] = paramval[1].strip()
    
  
    xi = int(parameters['xi'])
    xf = int(parameters['xf'])
    yi = int(parameters['yi'])
    yf = int(parameters['yf'])
    z_i = int(parameters['nu_i'])
    z_f = int(parameters['nu_f'])
    N = int(parameters['N'])
    r_init = float(parameters['r_init'])
    G = float(parameters['G']) 
    D = float(parameters['dpc']) 
    cs = float(parameters['cs'])
    vsys = float(parameters['vsys'])
    ssf = int(parameters['sub'])
    rebin = int(parameters['rebin'])
    gridsize = int(parameters['gridsize'])
    a = float(parameters['a'])
    q = float(parameters['q'])
    xc_ellipse = float(parameters['xc_ellipse'])
    yc_ellipse = float(parameters['yc_ellipse'])
    Gamma = float(parameters['GammaEllipse'])
    RL = int(parameters['RL'])

    # Import the ALMA Data Cube and extract key information from the header
    hdul = fits.open(datacube)
    data = hdul[0].data
    data = data[0,:,:,:]
    data = swap_cube_axes(data)
    if verbose == True:
        print('The shape of the data cube is',data.shape)
    
    # Key information from the header
    # Beam Major and Minor Axis
    x_std = hdul[0].header['BMAJ']*3600
    y_std = hdul[0].header['BMIN']*3600
    if verbose == True:
        print('The major and minor axis of the primary beam in arcseconds is',x_std,y_std)
    
    # Primary Beam Rotation
    resbeamPA = hdul[0].header['BPA']
    PSF_PA = (90.+resbeamPA)*(np.pi/180)
    
    # Pixel Scale 
    res = abs(hdul[0].header['CDELT1']*3600)
    if verbose == True:
        print('The pixel scale is', res)
    
    # Frequency spacing
    f_spacing = hdul[0].header['CDELT3']
    
    # Rest frequency of the molecular gas
    restfreq = hdul[0].header['RESTFRQ']
    if verbose == True: 
        print('The rest frequency of the molecular gas line in Hz is',restfreq)

    ### FREE PARAMETERS OF THE MODEL
    mbh = freeparams['mbh']
    MtoL = freeparams['MtoL']
    xc = freeparams['xc']
    yc = freeparams['yc']
    z = freeparams['z']
    theta = freeparams['theta']
    incl = freeparams['incl']
    F_0 = freeparams['F_0']
    sigma_0 = freeparams['sigma_0']
    sigma_1 = freeparams['sigma_1']
    sigma_2 = freeparams['sigma_2']
    mu = freeparams['mu']
    r_0 = freeparams['r_0']
    alpha = freeparams['alpha']
    PSF_scale = freeparams['PSF_scale']
    
    #### PRINT THE CURRENT PARAMETER VALUE 
    if counter is not None:
        print('The iteration count is at',int(counter))
    
    #### PRINT THE CURRENT PARAMETER VALUE 
    if loud == True:
        print('The black hole mass in solar masses is', float(mbh))
        print('The central x pixel is',float(xc))
        print('The central y pixel is',float(yc))
        print('The position angle is at',float(theta))
        print('The inclination angle is at',float(incl))
        print('The redshift is',float(z))
        print('The Mass to Light ratio is',float(MtoL))
        print('The flux multiplier constant is',float(F_0))
        print('The constant velocity dispersion term is',float(sigma_0))
        print('The amplitude of the velocity dispersion Gaussian is',float(sigma_1))
        print('The amplitude of the velocity dispersion exponential is',float(sigma_2))
        print('The radius offset of the velocity dispersion Gaussian is',float(r_0))
        print('The standard deviation of the velocity dispersion Gaussian is',float(mu))
        print('The value of alpha is',float(alpha))
        print('The value of PSF scale is',float(PSF_scale))
    
    ### Import the FIXED Parameters from the ALMA Parameter File     
    # Start a timer 
    start = time.time()

    # Disk inclination angle and rotation angle (PA)
    # in degrees, where i=0 is face-on and i=90 deg.
    # is edge-on and th=0 is north (up) and th=90 deg. is
    # east (left) of the receding disk major axis.
    # Afterwards, transform these angles to radians and
    # from the inclination angle compute the
    # disk minor/major axis ratio, qv
          
    incl=incl/180.*np.pi
    theta=(270.-theta)/180.*np.pi
    qv=np.cos(incl)
    
    # Construct matrices that define the native x/y positions
    # of an array with (ndx,ndy) dimensions. The native x/y positions
    # are shifted by the (xc,yc) disk centers at this stage.

    if modelcube is None:
        if verbose == True:
            print('No model cube is present')
        
        # Save the original dimensons of the full-sized ALMA grid
        X,Y,Z = data.shape
            
        # Truncate the data to work on a smaller grid 
        data = data[yi:yf,xi:xf,:]
        ndx = np.size(data,1)
        ndy = np.size(data,0)
        ndz = np.size(data,2)
        
        # Pre-Allocate Data Rebin Array
        data_rebin = np.zeros((int(ndy/rebin),int(ndx/rebin),ndz))
        for i in range(ndz):
            data_rebin[:,:,i] = block_reduce(data[:,:,i],rebin,np.mean)
            
        if verbose == True:
            print('ndx,ndy,and ndz are',(ndx,ndy,ndz))
            
        # Model cube is for fitting to an already made model cube, typically in the context of a Monte Carlo 
        # simulation

    elif modelcube is not None:
        X,Y,Z = data.shape
        
        if verbose == True:
            print('Model cube is present')
        
        data = data[yi:yf,xi:xf,:]
        ndx = np.size(data,1)
        ndy = np.size(data,0)
        ndz = np.size(data,2)
        
        # Pre-Allocate Data Rebin Array
        data_rebin = np.zeros((int(ndy/rebin),int(ndx/rebin),ndz))
        for i in range(ndz):
            data_rebin[:,:,i] = block_reduce(data[:,:,i],rebin,np.mean)
        
        if verbose == True:
            print('ndx,ndy,and ndz are',(ndx,ndy,ndz))
            

    # SAVE THE RE-BINNED DATA FOR FUTURE USE
    if verbose == True:
        print('The filename is',galaxyname)

    figtitle = galaxyname + 'RebinnedData' + '.fits'
    
    if SaveFITS == True:
        hdu=fits.PrimaryHDU(np.swapaxes(np.swapaxes(data_rebin,0,2),1,2))
        hdu.writeto(figtitle,overwrite=True)
    
    # Now, we must shift the xc and yc coordinates by the xi,xf, yi, and yf shift.
    xc = xc - xi
    yc = yc - yi
    
    xva=np.linspace(1,ssf*ndx,ssf*ndx)/ssf-xc
    yva=(np.linspace(1,1,ssf*ndy)/ssf)*np.linspace(1,ssf*ndy,ssf*ndy)-yc
    xva, yva = np.meshgrid(xva,yva)    
    
    # Transform the native x/y positions into
    # physical (disk) x/y positions using the
    # disk inclination and rotation angles.
    # Construct array that maps observed pixel
    # to a physical disk radius
    xv=(xva*np.cos(theta)-yva*np.sin(theta))
    yv=(yva*np.cos(theta)+xva*np.sin(theta))/qv
    
    # Radius at a projected location (x',y') on the sky plane
    rv=np.sqrt((xv)**2+(yv)**2)
    
    # Radius map
    if verbose == True:
        plt.figure(1)
        plt.imshow(rv,origin='lower')
        plt.title('Radius (pixels)')
        plt.xlabel('x pixels')
        plt.ylabel('y pixels')
        cb = plt.colorbar()
        cb.set_label('Radius ',fontsize = 16)
        plt.show()
    
    # Compute the fraction of the maximum line-of-sight (LOS)
    # velocity (along the disk major axis) at each observed
    # pixel position
    los_frac=(xva*np.cos(theta)-yva*np.sin(theta))/rv
    
    ## CREATE a los_frac free fall
    los_frac_ff = yv/rv

    # Convert the radius array from pixels to parsecs (pc)
    # Determine the parsec to arcsecond scale
    
    # Calculate angular diameter distance
    D_A = D/(1+(vsys/cs))**2
    
    pc2arcsec = D_A*np.tan(1/206265)
    
    # Determine the pc to pixel scale
    scale = pc2arcsec*res
    
    if loud == True:
        print('The pc to arcsec scale is',pc2arcsec)
        print('The pc to pixel scale is',scale)
    
    rv_pc=rv*scale
    
    ### Stellar and Gas contributions to the circular velocity 
    ### These will be added in quadrature with the black hole mass contribution
    
    if r_int is None and vc_st is None:
        if method == 'vext':
         ### New Free Parameters for the Extended Velocity Distribution with logarithmic bins
        # We start with 10 bins to start
            delta0 = freeparams['delta0']
            delta1 = freeparams['delta1']
            delta2 = freeparams['delta2']
            delta3 = freeparams['delta3']
            delta4 = freeparams['delta4']
            delta5 = freeparams['delta5']
            delta6 = freeparams['delta6']
            delta7 = freeparams['delta7']
            delta8 = freeparams['delta8']
            
            if verbose == True:
                print('The value of delta0 is',float(delta0))
                print('The value of delta1 is',float(delta1))
                print('The value of delta2 is',float(delta2))
                print('The value of delta3 is',float(delta3))
                print('The value of delta4 is',float(delta4))
                print('The value of delta5 is',float(delta5))
                print('The value of delta6 is',float(delta6))
                print('The value of delta7 is',float(delta7))
                print('The value of delta8 is',float(delta8))
    
            ### Create the logarithmic bins for the extended stellar circular velocity profile 
            logxbins = np.logspace(np.log10(r_init),np.log10(a*scale),N)
        
            # Set the first entry to correspond to r = 0.
            logybins = np.zeros((len(logxbins)+1))
            logybins[0] = 0
            logybins[1::] = logxbins[:]
            
            if verbose == True:
                print(logybins)
            # Create a mass array with N entries, one for each radii
            mguess = np.zeros(N+1)
            mguess[0] = 0
        
            # Create free parameter placeholders
            for i in range(0,N):
                varname = 'm' + str(i)
                mguess[i+1] = freeparams[varname]
            
            if verbose == True:
                print(mguess)
            
        # Accumulate the mass profile in order for it to be monotonically increasing
        # Interpolate the data with a cubic spline on the coarse scale
            interpfunc = interpolate.PchipInterpolator(logybins,mguess)

            # Create a more finely fampled radius vector that we can linearly interpolate over the entire rv_pc grid
            rnew = np.linspace(0,np.rint(a*scale),np.rint(np.round(a*scale)))
        
            # Re-interpolate with the monotonically increasing spline
            newinterpfunc = interpolate.PchipInterpolator(rnew,interpfunc(rnew))
        
            massinterp = np.array(newinterpfunc(rnew))
        
            # Plot the mass as a function of radius with the interpfunc
            plt.plot(rnew,massinterp)
            plt.title('Mass Profile')
            plt.xlabel('r (pc)')
            plt.ylabel('$M_{\odot}$')
            plt.show()
        
            # Convert the enclosed mass into a circular velocity at each radii 
            # Set the circular velocity to be 0 at r = 0
            vcirc = np.zeros(np.size(massinterp))
            vcirc[0] = 0
            vcirc[1::] = np.sqrt((G*massinterp[1::])/rnew[1::])
            
            if verbose == True:
                print('The velocity of stars at r = 0 is',vcirc[np.where(rnew == 0)])
        
            # Plot the circular velocity curve as a function of radius
            plt.plot(rnew,vcirc)
            plt.title('Interpolated Velocity')
            plt.xlabel('r (pc)')
            plt.ylabel('V (km/s)')
            plt.show()
    
            vcirctextfile = galaxyname + 'vext' + '.txt'
            masstextfile = galaxyname + 'mass' + '.txt'
            np.savetxt(vcirctextfile, np.c_[rnew, vcirc], fmt='%1.3f')
            np.savetxt(masstextfile, np.c_[rnew, newinterpfunc(rnew)], fmt='%1.3f')
        
            vc2ml = np.interp(rv_pc,rnew,vcirc)
            
            if verbose == True:
                plt.imshow(vc2ml,origin='lower')
                plt.title('Stellar Circular Velocity (km/s)')
                cb = plt.colorbar()
                plt.show()
        
            hdu=fits.PrimaryHDU(vc2ml)
            vcname = galaxyname + 'vext' + '.fits'
            
            if SaveFITS == True:
                hdu.writeto(vcname,overwrite=True)
            
        # Assuming that the mass goes as M(r) = c_1 r^beta
        elif method == 'power':
            c1 = freeparams['c1']
            beta = freeparams['beta']
            vc2ml = np.sqrt(G*c1*(rv_pc**(beta-1)))
            
            print('The value of c1 is',float(c1))
            print('The value of beta is',float(beta))
            
            if verbose == True:
                plt.imshow(vc2ml,origin='lower',cmap='jet')
                plt.title('Stellar Mass Profile - Power Law')
                cb = plt.colorbar()
                plt.show()
                
        else: 
            return print('None of the possible options were selected.')
        
    else:
        vc2ml = np.interp(rv_pc,r_int,vc_st)
        
        if verbose == True:
            plt.imshow(vc2ml,origin='lower')
            plt.title('Stellar Circular Velocity Contribution (km/s)')
            cb = plt.colorbar()
            plt.show()
        
    
    if r_gas_int is None and vc_gas is None:
        vc_gas = 0
    else:
        vc_gas = np.interp(rv_pc,r_gas_int,vc_gas)
        
        if verbose == True:
            plt.imshow(vc_gas,origin='lower')
            plt.title('Gas Circular Velocity Contribution (km/s)')
            cb = plt.colorbar()
            plt.show()
        
    # Determine Total and LOS circular velocity using the stellar mass-to-light
    # ratio, ml, and Newton's constant G
    vctotal=np.sqrt(MtoL*vc2ml**2 +(G*mbh/rv_pc)+vc_gas**2)
    
    ## Calculate a free-fall velocity, that is square root of 2 times the total circular speed
    vff = vctotal*np.sqrt(2)

    vlosrot = alpha*vctotal
    vlosff = np.sqrt(2*(1-alpha**2))*vctotal
    
    if freeparams['alpha'].vary is True:
        # Plot the total free-fall velocity and the component that we see as an inflow velocity
        
        # Total Inflow Velocity 
        if verbose == True:
            plt.imshow(vlosff,origin='lower')
            plt.title('Total Free-Fall Velocity')
            cb = plt.colorbar()
            plt.show()
                 
        # Save to a fits file
        if SaveFITS == True:
            hdu=fits.PrimaryHDU(vlosff)
            vlosffname = galaxyname + 'vlosff' + '.fits'
            hdu.writeto(vlosffname,overwrite=True)
        
        # LOS Inflow Velocity
        vlosff_frac = -np.sqrt(2*(1-alpha**2))*vctotal*los_frac_ff*np.sin(incl)
        
        if verbose == True:
            plt.imshow(vlosff_frac,origin='lower',cmap='jet')
            plt.title('LOS Inflow Velocity ')
            cb = plt.colorbar()
            plt.show()
          
        # Save to a fits file
        if SaveFITS == True:
            hdu = fits.PrimaryHDU(vlosff_frac)
            vlosff_fracname = galaxyname + 'vlosff_frac' + '.fits'
            hdu.writeto('vlosff',overwrite=True)

    vlostotal = alpha*vctotal*los_frac*np.sin(incl) - np.sqrt(2*(1-alpha**2))*vctotal*los_frac_ff*np.sin(incl)
    
    
    end_grid = time.time()
    print('The model grid construction time in seconds is',end_grid-start)

    sigmaturb = sigma_0 + sigma_1*np.exp(-(rv_pc-r_0)**2/(2*mu**2)) + sigma_2*np.exp(-rv_pc/r_0)

    # Plotting total circular velocity at each spatial position on the grid
    if verbose == True:
        plt.figure(2) 
        plt.imshow(vctotal,cmap='viridis',origin='lower')
        plt.xlabel('X Disk')
        plt.ylabel('Y Disk')
        cb = plt.colorbar()
        cb.set_label('$V_{Circular}$ (km/s)',fontsize = 16) 
        plt.show()
    
    # Plotting the LOS velocity 
        plt.figure(3)
        plt.imshow(vlostotal,cmap='jet',origin='lower')
        plt.xlabel('X Disk')
        plt.ylabel('Y Disk')
        plt.title('$V_{LOS}$ Total')
        cb = plt.colorbar()
        cb.set_label('$V_{LOS}$ (km/s)',fontsize = 16)  
        plt.show()
    
        plt.imshow(sigmaturb,origin='lower')
        plt.title('Turbulent Velocity Dispersion')
        plt.xlabel('X Disk')
        plt.ylabel('Y Disk')
        cb = plt.colorbar()
        plt.show()
    
    # Determine the corresponding frequency centroid and frequency width at each spatial position
    # f1 and f2 define the starting and ending frequency channels that define the line profile.
    f_0 = restfreq/1e9
    f_obs = (f_0/(1+z))*(1-(vlostotal/cs))
    df_obs = (f_0/(1+z))*(sigmaturb/cs)
    
    # Create a velocity range based on the optical definition of radial velocity 
    # c(f-f_0)/f
    v_range = cs*(f_range-f_0)/(f_range)
    
    # Model the PSF as an elliptical gaussian with a mean = 0 and standard deviation proportional to the 
    # FWHM of the synthesized beam 
    # Set the center of the Gaussian at (ndx/2,ndy/2)
 
    # In degrees
    PSF_PA = (90.+resbeamPA)*(np.pi/180)
    
    x_std *= PSF_scale
    y_std *= PSF_scale
    
    # If the the PSF size is made to vary, the grid size will change to the nearest odd integer to 8*sigma_x
    if PSF_scale.vary is False:
        gridsize = gridsize
    
    if PSF_scale.vary is True:
        # Ensure the grid size is still an odd number
        gridsize = 2*math.floor(np.rint(5*(x_std/res)*PSF_scale)/2) + 1
        print('The grid size is now',gridsize)
        
    # Create the ALMA Primary Beam
    PSF_sub = np.array(Gaussian2DKernel((x_std/res)/2.3548,(y_std/res)/2.3548,PSF_PA,x_size=gridsize,y_size=gridsize))
    PSF_sub = PSF_sub/np.max(PSF_sub)
    if SaveAny == True:
        hdu=fits.PrimaryHDU(PSF_sub)
        PSFtitle = galaxyname + 'PSF' + '.fits'
        hdu.writeto(PSFtitle,overwrite=True)
    
    # Create line profiles from the f_centroid and f_width arrays
    # Sigma is the dispersion, which is assumed to be flat for now. 
    # Centroid velocity is the systemic velocity of the galaxy

    delta_f = f_spacing/1e9
    # if statement to carry the correct minus sign in the event the frequency axis is increasing or decreasing
    # from the start
    if f_range[0] < f_range[1]:
        delta_f *= -1
        print('Delta frequency is', delta_f)
    
    elif f_range[0] > f_range[0]:
        delta_f *= 1
        print('Delta frequency is', delta_f)
    
    # Create an integrated gaussian line profile following the methodology of B. Boizelle
    # The integrated line profile will be weighted by the deconvolved flux map in the following step
    glineflux = np.swapaxes(np.swapaxes(-0.5*np.array([-scipy.special.erf((i-(delta_f/2)-f_obs)/(np.sqrt(2)*df_obs))+scipy.special.erf(((i+(delta_f/2)-f_obs)/(np.sqrt(2)*df_obs))) for i in f_range]),2,0),0,1)
    #glineflux[glineflux < 1e-6*np.max(glineflux)] = 0
    
    if verbose == True:
        print('The line profile array shape is currently',glineflux.shape)
    end_lineprofile = time.time()
    
    if verbose == True:
        print('The time to construct a model line profile in seconds is',end_lineprofile-start)
    
    unweightedlineprof_title = galaxyname + 'unweightedlineprofile' + '.fits'
    if SaveAny == True:
        hdu=fits.PrimaryHDU(np.swapaxes(np.swapaxes(glineflux,0,2),1,2))
        hdu.writeto(unweightedlineprof_title,overwrite=True)
    
    # Upscale the flux map by (ssf x ssf) 
    # Deconvolve the flux map with the Richardson-Lucy algorithm 
    # It is important that the flux map has been made to be strictly-non negative before the deconvolution
    
    FluxMap = restoration.richardson_lucy(FluxMap,PSF_sub,num_iter=RL)
    FluxMap = upsample(FluxMap[yi:yf,xi:xf],ssf)
 
    
    if SaveFITS == True:
        hdu=fits.PrimaryHDU(FluxMap)
        fluxmaptitle = galaxyname + 'fluxmapbeforeRL' + '.fits'
        hdu.writeto(fluxmaptitle,overwrite=True)
    
    if verbose == True:
        plt.imshow(FluxMap,origin='lower')
        plt.title('Flux Map (Jy/Beam) (Post RL Deconvolution)')
        cb = plt.colorbar()
        plt.xlabel('x pixels')
        plt.ylabel('y pixels')
        plt.show()
    
    # Weight the line profile by multiply each slice by the flux map
    for i in range(ndz):
        glineflux[:,:,i] = glineflux[:,:,i]*FluxMap[:,:]*F_0 
       
    if verbose == True:
        plt.imshow(np.sum(glineflux,2),origin='lower')
        cb=plt.colorbar()
        plt.title('Collapsed Normalized Line Profile')
        plt.xlabel('x pixels')
        plt.ylabel('y pixels')
        plt.show()
        
    ZZZ = np.sum(glineflux,2)
    
    if SaveFITS == True:
        hdu=fits.PrimaryHDU(ZZZ)
        modelMoment0 = galaxyname +'ModelMoment0' + '.fits'
        hdu.writeto(modelMoment0,overwrite=True)
        
        if SaveAny == True:
            modelpreconvolution = galaxyname + 'modelpreconvolution' + '.fits'
            hdu=fits.PrimaryHDU(np.swapaxes(np.swapaxes(glineflux,0,2),1,2))
            hdu.writeto(modelpreconvolution,overwrite=True)
    
    #glineflux[glineflux < np.max(glineflux)*1e-6] = 0
    
    if verbose == True:
        print('The shape of the model line profile array is',glineflux.shape)
      
    #Display the PSF 
    if verbose == True:
        plt.figure(4)
        plt.imshow(PSF_sub,interpolation='none',origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        cb = plt.colorbar()
        cb.set_label('PSF',fontsize = 16) 
        plt.show()
    
    # Display the flux map
        plt.figure(5)
        plt.imshow(FluxMap,origin='lower')
        plt.xlabel('x ')
        plt.ylabel('y ')
        cb = plt.colorbar()
        cb.set_label('Flux Map',fontsize = 16)
        plt.show()
        
    # Re-bin the integrated gaussian line profile to the scale of the original ALMA data for convolution efficiency
    # Pre-allocate the array first
    rebinned_glineflux = np.zeros((ndx,ndy,ndz))
    for i in range(ndz):
        rebinned_glineflux[:,:,i] = block_reduce(glineflux[:,:,i],ssf,np.sum)
        rebinned_glineflux[np.isnan(rebinned_glineflux)] = 0

    rebinned_glineflux[rebinned_glineflux < 1e-6*np.max(rebinned_glineflux)] = 0    
    
    # Convolve the PSF with the integrated Gaussian line profile 
    # First pre-allocate arrays to be filled at both the original scale and at the block-averaged scale.
    convolvetest = np.zeros((ndx,ndy,ndz))
    convolvetest_sub = np.zeros((int(ndx/(rebin)),int(ndy/(rebin)),ndz))
    
    # Define the fitting ellipse that determines where the model optimizations occur

    semimaj = a
    semimin = q*a
    Gamma=(90.+Gamma)/180.*np.pi
    e = Ellipse2D(amplitude=1., x_0=(xc_ellipse-xi), y_0=(yc_ellipse-yi), a=semimaj, b=semimin,theta=Gamma)
    y, x = np.mgrid[0:ndx,0:ndy]

    # Select the regions of the ellipse we want to fit
    # Create a fitting-cube that will contain the regions where the fits will occur
    fitting_ellipse = np.array(e(x,y))
    
    # Plot the elliptical region on the scale of the ALMA data 
    if verbose == True:
        plt.figure(6)              
        plt.imshow(fitting_ellipse,origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.title('Elliptical Fitting Region')
        cb = plt.colorbar()
        plt.show()
    
    # Save a 2D fits image on the full scale of the fitting ellipse
    if SaveFITS == True:
        fitellipsename = galaxyname + 'FittingEllipse_OriginalALMA_Scale' + '.fits'
        hdu=fits.PrimaryHDU(fitting_ellipse)
        hdu.writeto(fitellipsename,overwrite=True)
    
    ### FIND THE HIGHEST AND LOWEST X AND Y VALUES TO USE FOR THE CONVOLUTION BOX
    ### STORE THE X AND Y POINTS OF THE ELLIPSE IN ARRAYS
    good_ellipse_pixels = np.where(fitting_ellipse == 1)
    y_ellipse = good_ellipse_pixels[0]
    x_ellipse = good_ellipse_pixels[1]
    
    ### FIND THE MAX AND MIN VALUES OF BOTH THE X AND Y ARRAYS 
    x_ellipse_max = np.max(x_ellipse)
    x_ellipse_min = np.min(x_ellipse)
    y_ellipse_max = np.max(y_ellipse)
    y_ellipse_min = np.min(y_ellipse)
    
    if verbose == True:
        print('The minimum and maximum of the ellipse in the x direction is',x_ellipse_min,x_ellipse_max)
        print('The minimum and maximum of the ellipse in the y direction is', y_ellipse_min,y_ellipse_max)
    
    if verbose == True:
        print('The x length of the un-rebinned ellipse is',x_ellipse_max-x_ellipse_min)
        print('The y length of the un-rebinned ellipse is',y_ellipse_max-y_ellipse_min)
    
    ### ADD A BUFFER OF ABOUT 2*FWHM ON BOTH ENDS
    ### USE THESE VALUES TO PROPERLY CHOOSE THE CORRECT VALUES OF THE CUBE TO CONVOLVE
    
    # Take the average of x_std and y_std
    beamavg = np.mean((x_std/res,y_std/res))
    
    box_xlo = int(np.rint(x_ellipse_min-(beamavg)))
    box_xhi = int(np.rint(x_ellipse_max+(beamavg)))
    box_ylo = int(np.rint(y_ellipse_min-(beamavg)))
    box_yhi = int(np.rint(y_ellipse_max+(beamavg)))
    
    # If these limits go past the edges of the spatial dimension of the ALMA, set strict limits
    if box_xlo <= 0:
        box_xlo = 0
    if box_xhi >= ndx:
        box_xhi = ndx
    if box_ylo <= 0:
        box_ylo = 0
    if box_yhi >= ndy:
        box_yhi = ndy
    
    if verbose == True:
        print('The low and high and for the convolution box x dimension is',box_xlo,box_xhi)
        print('The low and high for the convolution box y dimension is',box_ylo,box_yhi)
    
    ### GENERATE THE CONVOLUTION BOX'S X WIDTH AND Y WIDTH VALUES
    box_x_width = box_xhi - box_xlo
    box_y_width = box_yhi - box_ylo
    
    if box_x_width >= ndx:
        box_x_width = ndx
    
    if box_y_width >= ndy:
        box_y_width = ndy
    
    if verbose == True:
        print('The convolution box width in the x-direction in pixels is',box_x_width)
        print('The convolution box width in the y-direction in pixels is',box_y_width)
    
    ### THE REPLACEMENT CONVOLUTION METHOD IS TO CREATE A RECTANGULAR REGION THAT ENCAPUSLATES THE FITTING ELLIPSE
    box_region = Box2D(amplitude=1.,x_0=(xc_ellipse-xi),y_0=(yc_ellipse-yi),x_width=box_x_width,y_width=box_y_width)
    convolve_box = np.array(box_region(x,y))
        
    ### PLOT THE CONVOLUTION BOX ON THE SCALE OF THE ALMA DATA
    if verbose == True:
        plt.figure(7)              
        plt.imshow(convolve_box,origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.title('Convolution Box Region')
        cb = plt.colorbar()
        plt.show()
    
    preconvolution_small = galaxyname + 'smallerregionmodelbeforeconvolution' + '.fits'
    if SaveFITS == True:
        hdu=fits.PrimaryHDU(np.swapaxes(np.swapaxes(rebinned_glineflux[box_ylo:box_yhi,box_xlo:box_xhi,:],0,2),1,2))
        hdu.writeto(preconvolution_small,overwrite=True)
    
    ### MULTIPROCESSING CONVOLUTION 
    mp1 = time.time()
    
    # Now select the portion of convolvetest to be convolved properly
    cubepsflist = []
    
    # No boundary='extend' as this crops the cube somehow 
    for i in list(range(z_i-1,z_f-1)):
        cubepsflist.append((rebinned_glineflux[box_ylo:box_yhi,box_xlo:box_xhi,i],PSF_sub,'extend','normalize_kernel=True'))
    
    pool = multiprocessing.Pool(processes=2,maxtasksperchild=1)
    if __name__ == '__main__':
        mappedarray = np.array(pool.starmap(convolve,cubepsflist))
        pool.close()
        pool.join()
        
    mappedarray = swap_cube_axes(mappedarray)
    
    if loud == True:
        print('The shape of the mapped array is',mappedarray.shape)
    
    mp2 = time.time()
    
    if loud == True:
        print('The time multiprocessing convolution takes is',mp2-mp1)
    
    # Convolved Model
    for i in list(range(z_i-1,z_f-1)):
        convolvetest[box_ylo:box_yhi,box_xlo:box_xhi,i] = mappedarray[:,:,i-(z_i-1)]
        #convolvetest[convolvetest < 1e-6*np.max(convolvetest)] = 0 
    
    if loud == True:
        print('The shape of the convolved array is',convolvetest.shape)
    
    # Save the model on the ALMA scale (cropped spatial axis)
    if SaveFITS is True:
        hdu=fits.PrimaryHDU(np.swapaxes(np.swapaxes(convolvetest,0,2),1,2))
        ALMAsizemodeltitle = galaxyname + 'ModelCubeALMAscale' + '.fits'
        hdu.writeto(ALMAsizemodeltitle,overwrite=True)
        
    end_convolve = time.time()
    
    if loud == True:
        print('The time to perform the convolution in seconds has taken',end_convolve-start)
    
    # Save the model on the ALMA scale (full spatial axis)
    modelFS = np.zeros((Y,X,Z))

    # Save the model cube 
    modelFS[yi:yf,xi:xf,:] = convolvetest[:,:,:]
    if SaveFITS is True:
        hdu=fits.PrimaryHDU(np.swapaxes(np.swapaxes(modelFS,0,2),1,2))
        ALMAsizemodeltitle_full = galaxyname + 'ModelCubeALMAscale_Full' + '.fits'
        hdu.writeto(ALMAsizemodeltitle_full,overwrite=True)
    
    
    # Re-Bin the Model by averaging over m x m blocks in the spatial domain.
    for i in range(ndz):
        convolvetest_sub[:,:,i] = block_reduce(convolvetest[:,:,i],rebin,func=np.mean)
        
    if SaveFITS is True:
        rebinned_finalmodeltitle = galaxyname + 'FinalModelRebinned' + '.fits'
        hdu=fits.PrimaryHDU(np.swapaxes(np.swapaxes(convolvetest_sub,0,2),1,2))
        hdu.writeto(rebinned_finalmodeltitle,overwrite=True)
        
    model = convolvetest_sub
    
    if verbose == True:
        print('The new noise cube''s shape is',eps.shape)
        print('The shape of the model ALMA data cube is',model.shape)
    
    if SaveFITS is True:
        noise_title = galaxyname + 'noisecubesmallregion' + '.fits'
        hdu=fits.PrimaryHDU(np.swapaxes(np.swapaxes(eps,0,2),1,2))
        hdu.writeto(noise_title,overwrite=True)
    
    # Sub-sampling factor
    # Create a smaller fitting ellipse that will be used to perform fit on the sub-sampled scale
    fitting_ellipse_small = np.array(block_reduce(fitting_ellipse,rebin,np.mean))
    
    if verbose == True:
        plt.imshow(fitting_ellipse_small,origin='lower')
        plt.title('Smaller Elliptical Region for Chi-Squared Fits')
        plt.show()
    
    # Create a smaller fit cube
    fit_cube_small = np.ones((int(ndx/rebin),int(ndy/rebin),ndz))
    for i in range(ndz):
        if i < (z_i-1):
            fit_cube_small[:,:,i] = 0
        elif i > (z_f-1):
            fit_cube_small[:,:,i] = 0
        else: 
            fit_cube_small[:,:,i] = fitting_ellipse_small*fit_cube_small[:,:,i]
            
    ### SAVE THE FIT CUBE REGION
    if SaveFITS is True:
        fitcubetitle = galaxyname + 'FITCUBE' + '.fits'
        hdu=fits.PrimaryHDU(np.swapaxes(np.swapaxes(fit_cube_small,0,2),1,2))
        hdu.writeto(fitcubetitle,overwrite=True)

    # Save the Fitting Region Cube (Down-sampled version)
    
    # Select the Fit region to be the regions that only contain a value greater than 0.5
    if fit_region is None:
        fittingregion = np.where(fit_cube_small >= 0.5)
        if loud == True:
            print('The amount of points in the fitting region is',np.size(fittingregion))
    
        x_fit = fittingregion[1]
        y_fit = fittingregion[0]
    
    elif fit_region is not None:
        if verbose == True:
            plt.imshow(fit_region[:,:,z_i+5],origin='lower')
        fittingregion = np.where(fit_region >= 0.5)
        
        if SaveFITS == True:
            fithdu = fits.PrimaryHDU(swap_cube_axes(fit_region))
            fithdutitle = galaxyname + 'FITCUBE_IMPORTED' + '.fits'
            fithdu.writeto(fithdutitle,overwrite=True)
            
        if verbose == True:
            print('Imported fit cube is used')
            print('The amount of points in imported fitting region is',np.size(fittingregion))
        
    
    # Create a residual cube on the scale of the rebinned ALMA data
    residual_cube = model-data_rebin
    
    ### PRINT OUT SLICES OF BOTH THE MODEL AND DATA CUBES FOR A QUICK, QUALITATIVE INSPECTION
    if verbose == True:
        plt.figure(8)              
        plt.imshow(model[:,:,(z_i+10)],origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        cb = plt.colorbar()
        plt.title('Model Frequency Slice z_i + 10')
        plt.show()
    
        plt.figure(9)              
        plt.imshow(data_rebin[:,:,(z_i+10)],origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        cb = plt.colorbar()
        plt.title('Data Frequency Slice z_i + 10')
        plt.show()
    
        plt.figure(10)              
        plt.imshow(model[:,:,(z_i+20)],origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        cb = plt.colorbar()
        plt.title('Model Frequency Slice z_i + 20')
        plt.show()
    
        plt.figure(11)              
        plt.imshow(data_rebin[:,:,(z_i+20)],origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        cb = plt.colorbar()
        plt.title('Data Frequency Slice z_i + 20')
        plt.show()
    
    
    if modelcube is not None:
        if verbose == True:    
            plt.figure(12)              
            plt.imshow(modelcube[:,:,(z_i+10)],origin='lower')
            plt.xlabel('x [pixels]')
            plt.ylabel('y [pixels]')
            cb = plt.colorbar()
            plt.title('Imported Model Cube Slice z_i + 10')
            plt.show()
        
            plt.figure(13)              
            plt.imshow(modelcube[:,:,(z_i+5)],origin='lower')
            plt.xlabel('x [pixels]')
            plt.ylabel('y [pixels]')
            cb = plt.colorbar()
            plt.title('Imported Model Cube Slice z_i + 20')
            plt.show()
        
    
    # Save the Synthetic Data Cube Before Flattening 
    if SaveAny is True:
        final_save = galaxyname + 'FinalModelBeforeFlattening' + '.fits'
        hdu=fits.PrimaryHDU(np.swapaxes(np.swapaxes(convolvetest_sub,0,2),1,2))
        hdu.writeto(final_save,overwrite=True)
    
    # Identify the fitting regions in the data, model, and noise cube and flatten them to a 1D array
    # This must be done because LMFIT only accepts 1D arrays to perform chi-squared minimization.  
    model_vector =  model[fittingregion].flatten('C')
    if verbose == True:
        print('The number of data points in the chi-squared fit is',np.size(model_vector))
    end_chi = time.time()
    if verbose == True:
        print('The time it has taken to calculate chi in seconds is',end_chi-start)
    
  
    # Create a Chi-Square Map
    chisquarecube = np.zeros((np.size(model,0),np.size(model,1),np.size(model,2)))
#     for i in range(z_i-1,z_f-1):
#         chisquarecube[:,:,i] = (data_rebin[:,:,i] - model[:,:,i])**2/(eps[:,:,i]**2)

    chisquarecube[np.where(fit_cube_small == 1)] = (data_rebin[np.where(fit_cube_small == 1)] - model[np.where(fit_cube_small == 1)])**2/(eps[np.where(fit_cube_small == 1)]**2)
        
    chisquarecubetitle = galaxyname + 'chisquarecube' + '.fits'
    
    if SaveFITS == True:
        hdu = fits.PrimaryHDU(np.swapaxes(chisquarecube,0,2))
        hdu.writeto(chisquarecubetitle,overwrite=True)
        
    chisquaremap = np.nansum(chisquarecube,2)  
    
    if verbose == True:
        plt.imshow(chisquaremap,origin='lower')
        plt.title('Chi-Square Map')
        plt.show()
    
    if SaveFITS is True:
        hdu=fits.PrimaryHDU(chisquaremap)
        chi_square_map_name = galaxyname + 'chi_square_map' + '.fits'
        hdu.writeto(chi_square_map_name,overwrite=True)
    
    data_vector = data_rebin[fittingregion].flatten('C')
    
    eps_vector = eps[fittingregion].flatten('C')

    if modelcube is None:
        if loud == True:
            print('No model cube is present')
        chi = (data_vector-model_vector)/(eps_vector)
    elif modelcube is not None:
        print('Model cube is present')
        modelcubevector = modelcube[fittingregion].flatten('C')
        chi = (modelcubevector - model_vector)/(eps_vector)
    
    if verbose == True:
        print('The reduced chi-squared value is',np.sum(chi**2)/(np.size(chi)-np.size(freeparams)-2))

    # Save all the FITS files with the usage of boolean argument
    if SaveFITS is True:
    # Create FITS (Flexible Image Transport System) files 
    # from scratch to save the centroid frequencies and frequency widths
        hdu=fits.PrimaryHDU(f_obs)
        fobstitle = galaxyname +'fobsarray' + '.fits'
        hdu.writeto(fobstitle,overwrite=True)
    
        hdu=fits.PrimaryHDU(df_obs)
        dfobstitle = galaxyname + 'dfobsarray' + '.fits'
        hdu.writeto(dfobstitle,overwrite=True)
    
    # Save the deconvovled flux map as a FITS file for future use
        hdu=fits.PrimaryHDU(FluxMap)
        deconvolvedfluxmap = galaxyname + 'deconvolvedfluxmap' + '.fits'
        hdu.writeto(deconvolvedfluxmap,overwrite=True)
    
    # Write the Gaussian Line Profile to a FITS file
        if SaveAny is True:
            hdu=fits.PrimaryHDU(np.swapaxes(np.swapaxes(glineflux,0,2),0,1))
            glinefluxtesttitle = galaxyname + 'gaussianlineprofiletest' + '.fits'
            hdu.writeto(glinefluxtesttitle,overwrite=True)
    
    # Create FITS (Flexible Image Transport System) files 
    # from scratch to save the radial position and LOS velocity arrays.
        hdu=fits.PrimaryHDU(rv_pc)
        rv_pctitle = galaxyname + 'rv_pc' + '.fits'
        hdu.writeto(rv_pctitle,overwrite=True)
        
        hdu=fits.PrimaryHDU(los_frac)
        VLOSfractitle = galaxyname + 'vlosfrac' + '.fits'
        hdu.writeto(VLOSfractitle,overwrite=True)
        
        hdu=fits.PrimaryHDU(vlostotal)
        VLOStotaltitle = galaxyname + 'vlostotal' + '.fits'

        hdu.writeto(VLOStotaltitle,overwrite=True)
        
        hdu=fits.PrimaryHDU(vctotal)
        vctotaltitle = galaxyname + 'vctotaltitle' + '.fits'
        hdu.writeto(vctotaltitle,overwrite=True)
        
    # Save the residual cube to a FITS file.
        residualcubetitle = galaxyname + 'residualcube' + '.fits'
        hdu = fits.PrimaryHDU(np.swapaxes(residual_cube,0,2))
        hdu.writeto(residualcubetitle,overwrite=True)
    time_end = time.time()
    print('The time to optimize this program in seconds is',time_end-time_start)
    return chi