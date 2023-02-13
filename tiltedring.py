# Define Tilted Ring Class

# Import standard numpy libraries
import numpy as np
import astropy
from astropy.io import fits
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.interpolate import SmoothBivariateSpline
from scipy.interpolate import griddata
from jampy.mge_vcirc import mge_vcirc



class TiltedRing:
	"""
	This is an implementation of a tilted-ring model for a rotating
	and inclined disk. 

	Attributes
	----------
	n_rings: Int 
			The number of rings to model the velocity field
			typically given as the number of resolution elements
			across the major axis of the disk

	x_c: float 
		The disk centroid x-position in pixels

	y_c: float
		The disk centroid y-position in pixels

	incl: float
		The ring's inclination angle in degrees

	PA: float
		The ring's position angle in degrees

	v_circ: float
		

	"""

	def __init__(self,x_c,y_c,incl,PA,v_circ,v_sys):
		""" Initialize attributes to describe a single tilted ring which
		is described by the 6 parameters described above """

		self.x_c = x_c
		self.y_c = y_c
		self.incl = incl
		self.PA = PA
		self.v_circ = v_circ
		self.v_sys = v_sys


