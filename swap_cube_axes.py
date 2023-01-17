# Import both the primary beam corrected data cube and the primary beam cube itself to create a noise cube.
# Also import the swap cube function to properly swap the data cube

def swap_cube_axes(data_cube):
    # This assumes the cube has 3 dimensions and currently has the orientation (z,x,y)
    data_cube = np.swapaxes(np.swapaxes(data_cube,0,2),0,1)
    
    return data_cube