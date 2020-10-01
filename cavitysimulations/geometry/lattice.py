#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# ---------------------------------------------------------------------#
#  Lattice.py
#
#  Copyright 2018 Andreas Grisch <andreas.gritsch@mpq.mpg.de>
#
# This class represents a 1d waveguide.
# Its most important attribute is (m = Lx, n = Ly) matrix.
# The matrix columns represent (x, y, z, hole_r_x, hole_r_y)

# ---------------------------------------------------------------------#
# LIBRARIES
# ---------------------------------------------------------------------#
import numpy as np
import math
from copy import copy
import h5py
import meep as mp# CHANGE CHANGE CHANGE
# ---------------------------------------------------------------------#
# del_w, del_a, del_hy, del_hx = 0.05, 0.001, 0.025, 0.025   
# w_max ,a_max = 0.7, 0.45
# w_min, a_min, hy_min, hx_min = 0.65, 0.25, 0.1, 0.05

del_w, del_a, del_hy, del_hx = 0.05, 0.001, 0.025, 0.025   
w_max ,a_max = 0.7, 0.45
w_min, a_min, hy_min, hx_min = 0.65, 0.25, 0.1, 0.05
class OneDLattice():
    
    #p_data = "cavitysimulations/geometry/bandstructure_data/sweep_data.hdf5"
    p_data = "bandstructure_data/sweep_data.hdf5"


    def __init__(self, Lx, filename = p_data):
        """
        Creates an numpy array (5,Lx) representing the cavity.
        The first 3 entries are the x,y,z coordinate.
        The last 2 represent hx and hy.
        """
    # -------------------CHANGE CHANGE CHANGE----------------------------------
        
        self.coordinates = np.zeros((Lx, 5))
        self.filename = filename
        self.Lx = Lx
        for i in range(Lx):
            self.coordinates[i, 0] = i
            self.coordinates[i, 3] = 1
            self.coordinates[i, 4] = 1
        self.data = self.load_data(self.filename)

    def set_y(self, y):
        """
        Set the y value for all holes. Analog to set_z()
        """
        for i in range(self.Lx):
            self.coordinates[i, 1] = y

    def set_z(self, z):
        """
        Set the z value for all holes. Analog to set_y()
        """
        for i in range(self.Lx):
            self.coordinates[i, 2] = z

    def set_hx(self, hx):
        """
        Set hx value for all holes. Analog to set_hy()
        """
        for i in range(self.Lx):
            self.coordinates[i, 3] = hx

    def set_hy(self, hy):
        """
        Set hy value for all holes. Analog to set_hx()
        """
        for i in range(self.Lx):
            self.coordinates[i, 4] = hy
    # -------------------------------

    def modify_xposition(self, i, deltax):
        """
        Modify x coordinate by adding deltax at position i.
        """
        self.coordinates[i, 0] = self.coordinates[i, 0] + deltax
    # -------------------------------

    def modify_hx(self, i, hx_new):
        """
        Set hx to hx_new at position i. Analog to modify_hy()
        """
        self.coordinates[i, 3] = hx_new

    def modify_hy(self,i,hy_new):
        """
        Set hx to hy_new at position i. Analog to modify_hx()
        """
        self.coordinates[i, 4] = hy_new

    def modify_diameter(self, i, diameter_new):
        """
        Set hx and hy to diameter at position i.
        """
        self.coordinates[i, 3] = diameter_new
        self.coordinates[i, 4] = diameter_new
    # -------------------------------

    def remove_hole(self, i):
        """
        Remove a hole by setting its x_coordinate to -1
        """
        self.coordinates[i, 0] = -1
    # -------------------------------
    
    def normal_spacing(self, a, hx, hy):
        self.spacing = []
        for i in range(self.Lx):
            self.spacing.append(a)
        self.set_hx(hx)
        self.set_hy(hy)
        self.coordinates[0,0] = self.spacing[0] / 2
        
        for i in range(1, self.Lx):
            self.coordinates[i,0] = sum(self.spacing[:i]) + self.spacing[i]/2 
        
        return self.spacing
    def output_full(self, lattice_constant=1, hx=1, hy=1):
        """
        Return 'full' x,y,z,hx,hy, meaning the whole np.array plus its mirrored values at x=0
        """
        idxs = np.any(self.coordinates >= 0, axis=1)
        output = self.coordinates[idxs, :]*np.array([lattice_constant, lattice_constant, lattice_constant, hx, hy]) + np.array([lattice_constant/2, 0, 0, 0, 0])
        output = np.append((output*np.array([-1, 1, 1, 1, 1]))[::-1, :], output, axis=0)
        return output

    def output_full_negative_axis(self, lattice_constant=1, hx=1, hy=1):
        """
        Return 'full' x,y,z,hx,hy, only on positive x_axis
        """
        idxs = np.any(self.coordinates >= 0, axis=1)
        output = self.coordinates[idxs, :]*np.array([lattice_constant, lattice_constant, lattice_constant, hx, hy]) + np.array([lattice_constant/2, 0, 0, 0, 0])
        output = output*np.array([-1, 1, 1, 1, 1])
        return output

    def output_full_positive_axis(self, lattice_constant=1, hx=1, hy=1):
        """
        Return 'full' x,y,z,hx,hy, only on negative x_axis
        """
        idxs = np.any(self.coordinates >= 0, axis=1)
        output = self.coordinates[idxs, :]*np.array([lattice_constant, lattice_constant, lattice_constant, hx, hy]) + np.array([lattice_constant/2, 0, 0, 0, 0])
        return output

    def output_coordinates_positive_axis(self, lattice_constant=1):
        """
        Return x,y,z, only on positive x_axis, analog to output_coordinates_negative_axis()
        """
        idxs = np.any(self.coordinates >= 0, axis=1)
        output = self.coordinates[idxs, 0:3]*lattice_constant + np.array([lattice_constant/2, 0, 0])
        return output

    def output_coordinates_negative_axis(self, lattice_constant=1):
        """
        Return x,y,z, only on negative x_axis, analog to output_coordinates_positive_axis()
        """
        idxs = np.any(self.coordinates >= 0, axis=1)
        output = self.coordinates[idxs, 0:3]*lattice_constant + np.array([lattice_constant/2, 0, 0])
        output = output*np.array([-1, 1, 1])
        return output

    def output_coordinates_full(self, lattice_constant=1):
        """
        Return x,y,z, on both axis.
        """
        idxs = np.any(self.coordinates >= 0, axis=1)
        self.output = self.coordinates[idxs, 0:3]*lattice_constant + np.array([lattice_constant/2, 0, 0])
        self.output = np.append((self.output*[-1, 1, 1])[::-1], self.output, axis=0)
        return self.output

    def output_holes_diameters_positive_axis(self, hx=1, hy=1):
        """
        Return hx,hy, on positive axis.
        """
        idxs = np.any(self.coordinates >= 0, axis=1)
        self.holes = self.coordinates[idxs,3:4]*[hx,hy]
        return self.holes

    def output_holes_diameters_negative_axis(self, hx=1, hy=1):
        """
        Return hx,hy, on negative axis.
        """
        idxs = np.any(self.coordinates >= 0, axis=1)
        self.holes = self.coordinates[idxs,3:4]*[hx,hy]
        return self.holes

    def output_holes_diameters_full(self, hx=1, hy=1):
        """
        Return hx,hy, on both sides.
        """
        idxs = np.any(self.coordinates >= 0, axis=1)
        self.holes = self.coordinates[idxs,3:4]*[hx,hy]
        self.holes = np.append(self.holes,self.holes*[1,1], axis=0)
        return self.holes

    def tapered_spacings(self,a_center,a_lin_center, a_lin_mirror, n_lin ):
        """
        Taper spacing accroding to Loncar from 2008/simpetus tutorial 'experimental'
        """
        self.spacing = np.zeros(self.Lx)
        self.spacing[0] = a_center/2.
        delta_a=(a_lin_mirror-a_lin_center)/(n_lin)
        for i in range(n_lin):
            self.spacing[i+1]=a_lin_center+delta_a*i
        for i in range(self.Lx-n_lin-1):
            self.spacing[i+n_lin+1]=a_lin_mirror
        return self.spacing
    
    def pow_degree_a_taper(self, 
                           N_taper,
                           hx , 
                           hy ,
                           w  ,
                           a_center, 
                           a_mirror, 
                           pow = 2):
        '''
        N_taper = Total tapered segments ( including centre and mirror segments)
        pow = degree of tapering polynomial 
        '''
        N_taper = N_taper
        self.pow_spacing = []
        
        for i in range(N_taper):
            self.pow_spacing.append(a_center - ((i/(N_taper-1)) ** pow) * (a_center - a_mirror)) 
            
        for i in range(self.Lx - N_taper):
            self.pow_spacing.append(a_mirror)
        
        if (hx != 0 and hy != 0):
            self.set_hx(hx)
            self.set_hy(hy)

        return self.pow_spacing
    
    def apply_pow_spacing(self):        
        
        self.coordinates[0,0] = self.pow_spacing[0] / 2
        for i in range(1, len(self.pow_spacing)):
            self.coordinates[i,0] = sum(self.pow_spacing[:i]) + self.pow_spacing[i]/2  
            
            '''
            |<--------------- line of symmetry
            |______________________
            | o | o | o | o | o |o |
            |______________________
            |   
            ---------> x
                    3rd Hole
                    
            The ith hole has x = a_center/2(=self.polyspacing[0]/2) + self.polyspacing[1] +.... self.poly_spacing[i-1]
            
            '''

    def apply_tapered_spacings(self):
        """
        Apply these :Taper spacing accroding to Loncar from 2008/simpetus tutorial 'experimental'
        """
        self.coordinates[0,0] = self.spacing[0]
        for i in range(1,len(self.spacing-1)):
            self.coordinates[i,0] = self.spacing[i]+self.coordinates[i-1,0]

    def create_background_pattern(self):
        """
        Create background pattern. I dont know why this is here.
        """
        for i in range(self.Lx):
            self.remove_hole(i)

    def lin_hole_taper(self, relative_diameter_at_end):
        """
        Make a linear hole taper. you can specify the relative_diameter_at_end
        """
        delta_radius = (1-relative_diameter_at_end)/self.Lx
        for i in range(self.Lx):
            self.modify_diameter(i,1-delta_radius*i)

    def quadratic_hole_taper(self, relative_diameter_at_end, number_of_tapered_holes=100):
        """
        Make a quadratic hole taper. you can specify the relative_diameter_at_end.
        Specify number_of_tapered_holes if you want to have 'mirror'-holes at the end.
        """
        N = number_of_tapered_holes # use this abbrevation against eye cancer
        if N>self.Lx:
            N=self.Lx
        y_off = (N-1)**2
        for i in range(N):
            self.modify_diameter(i, (y_off-i**2*(1-relative_diameter_at_end))/y_off)
        for i in range(N, self.Lx):
            self.modify_diameter(i, (y_off-N**2*(1-relative_diameter_at_end))/y_off)

    def quartic_elliptical_hole_taper(self, relative_diameter_at_end, number_of_tapered_holes=100):
        """
        Make a quadratic hole taper. you can specify the relative_diameter_at_end.
        Specify number_of_tapered_holes if you want to have 'mirror'-holes at the end.
        """
        N = number_of_tapered_holes  # use this abbrevation against eye cancer
        if N > self.Lx:
            N = self.Lx
        y_off = (N - 1) ** 4
        for i in range(N):
            self.modify_hx(i, 1)
            self.modify_hy(i, (y_off - i ** 4 * (1 - relative_diameter_at_end)) / y_off)
        for i in range(N, self.Lx):
            self.modify_hx(i, 1)
            self.modify_hy(i, (y_off - N ** 4 * (1 - relative_diameter_at_end)) / y_off)
            
    # -------------------CHANGE CHANGE CHANGE----------------------------------#
    
    def polynomial_elliptical_hole_taper(self, number_of_tapered_holes, hx , hy , 
                                         w = 0.7 ,
                                         a_center = 0.390 ,  
                                         a_mirror = 0.449):
                                                              
        '''
        
        a_center(gamma_center) : lattice constant ( gamma ) at the cavity segment
        a_mirror(gamma_mirror) : lattice constant ( gamma ) for non-tapered mirror segment
        data                   : sweep data
        polynomial             : curve thata fits gamma vs a data
        poly_coeff             : degrees of freedom for the polynomial, array = size (degree+1,)
        gamma_arr              : array of N_taper equispaced gammas between (and including) gamma_center and gamma_mirror
        a_arr                  : array of N_taper values of a for each mirror segment
        
        
        '''

#         del_w, del_a, del_hy, del_hx = 0.05, 0.001, 0.025, 0.025
       
#         w_max ,a_max = 0.7, 0.45
        
#         w_min, a_min, hy_min, hx_min = 0.65, 0.25, 0.1, 0.05
        
        
        index_a_center = int((a_center - a_min)/ del_a + 0.1)
        index_a_mirror = int((a_mirror - a_min)/del_a + 0.1)
        index_w = int((w - w_min) / del_w + 0.1)
        index_hy = int((hy - hy_min) / del_hy + 0.1)
        index_hx = int((hx - hx_min) / del_hx + 0.1)
        
        gamma_center = self.data[index_w, index_a_center, index_hy, index_hx]
        gamma_mirror = self.data[index_w, index_a_mirror, index_hy, index_hx]
        
        for i in range(self.Lx):
            self.modify_hx(i, hx)
            self.modify_hy(i, hy)
            
        N_taper = number_of_tapered_holes
        N_mirror = self.Lx - N_taper
        
        data = self.data
        
        gamma_arr = np.linspace(gamma_center, gamma_mirror, N_taper) 
    
        poly_coeff = polynomial_fit(self.data, w , hy, hx , gamma_mirror, degree = 4 )
        polynomial = np.poly1d(poly_coeff)
        
        a_arr = polynomial(gamma_arr)
        
        self.poly_spacing = np.append(a_arr, np.zeros(N_mirror))
            
        for i in range(N_mirror):
            self.poly_spacing[N_taper+i] = a_mirror
        
        self.poly_spacing[0] = a_center
        
        return self.poly_spacing
    
    # -------------------CHANGE CHANGE CHANGE----------------------------------#
    
    def apply_poly_spacing(self):        
        
        self.coordinates[0,0] = self.poly_spacing[0] / 2
        for i in range(1, len(self.poly_spacing)):
            self.coordinates[i,0] = sum(self.poly_spacing[:i]) + self.poly_spacing[i]/2  
            
            '''
            |<--------------- line of symmetry
            |______________________
            | o | o | o | o | o |o |
            |______________________
            |   
            ---------> x
                    3rd Hole
                    
            The ith hole has x = a_center/2(=self.polyspacing[0]/2) + self.polyspacing[1] +.... self.poly_spacing[i-1]
            
            '''
        
    def load_data(self, filename = "sweep_data(1).hdf5"):    # CHANGE CHANGE CHANGE
        '''
        Loads data from hdf5 file
        '''
        hf = h5py.File(filename, 'r')
        data = np.array( hf.get("gamma"))
        self.data = data
        hf.close()
        
        x = np.where(data ==-1)
        for i in range(len(x[0])):
            data[x[0][i], x[1][i], x[2][i], x[3][i] ] = 0    # data cleaning 
        
        return data
    
    def get_index(self, w, a, hy, hx):
#         del_w, del_a, del_hy, del_hx = 0.05, 0.001, 0.025, 0.025
#         w_max ,a_max = 0.7, 0.45
#         w_min, a_min, hy_min, hx_min = 0.65, 0.25, 0.1, 0.05
        
        index_a = int((a - a_min)/del_a + 0.1)
        index_w = int((w - w_min) / del_w + 0.1)
        index_hy = int((hy - hy_min) / del_hy + 0.1)
        index_hx = int((hx - hx_min) / del_hx + 0.1)
        
        return index_w, index_a, index_hy, index_hx
    
 # -------------------CHANGE CHANGE CHANGE---------------------------------- #   
    
def polynomial_fit(data, w, hy, hx, gamma_mirror, degree = 10, a_mirror = 0.347 ):            
        
        
    '''
    NOT a class function 
    
    Fits a polynomial curve of mentioned degree to gamma vs lattice constant(a) data ( for fixed w, hy and hx)
    '''
    
    
        
#     del_w, del_a, del_hy, del_hx = 0.05, 0.001, 0.025, 0.025
#     w_max ,a_max = 0.7, 0.45
#     w_min, a_min, hy_min, hx_min = 0.65, 0.25, 0.1, 0.05
       
#     w_min, a_min, hy_min, hx_min = 0.65, 0.25, 0.1, 0.05


    index_a_mirror = int((a_mirror - a_min)/del_a + 0.1)
    index_w = int((w - w_min) / del_w + 0.1)
    index_hy = int((hy - hy_min) / del_hy + 0.1)
    index_hx = int((hx - hx_min) / del_hx + 0.1)
    
    
    a_list = np.arange(a_min, a_max, del_a)                        # list of a
        
    gamma_arr = data[index_w, :, index_hy, index_hx]               # array of gamma for fixed w, hy and hx
    gamma_NZ = gamma_arr[ gamma_arr > 0]
    
    gamma_interest = gamma_NZ[ : np.where( (gamma_NZ == gamma_mirror))[0][0]]
    a_interest = a_list[ np.where(gamma_arr > 0)][ : len(gamma_interest) ]
        
    p_coeff = np.polyfit(gamma_interest, a_interest, degree)       # degress of freedom for the polynomial
    
    return p_coeff
    
