#!/usr/bin/env python
# coding: utf-8

## This script is developed for the purpose of visualization 
## of the FARGO output
## This script combines the FARGO recommended script and the multifluid script developed by me

# Author: sayantan
# Date : 9 September 2019
# modified : 20 November 2019 to modify the polar plots with cartesian coordinates
# moddified: 26 December for generating the training set data


####### Importing the required modules

from pylab import*
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rc('text', usetex=True)
import os

## optional import for multiprocessing (under developement)
# import time
# import glob
# from matplotlib import ticker, cm
# import multiprocessing as mp
# import concurrent.futures
# import cv2



## Making the Parameter class 
class Parameters(object):
    """
    Class for reading the simulation parameters.
    input: string -> name of the parfile, normally variables.par
    """
    def __init__(self, output_folder):

        if len(output_folder) > 1:
            if output_folder[-1] != '/':
                output_folder += '/'
        try:
            params = open(output_folder +"variables.par",'r') #Opening the parfile
        except IOError:                  # Error checker.
            print  (paramfile + " not found.")
            return
        lines = params.readlines()     # Reading the parfile
        params.close()                 # Closing the parfile
        par = {}                       # Allocating a dictionary
        for line in lines:             #Iterating over the parfile
            name, value = line.split() #Spliting the name and the value (first blank)
            try:
                float(value)           # First trying with float
            except ValueError:         # If it is not float
                try:
                    int(value)         # we try with integer
                except ValueError:     # If it is not integer, we know it is string
                    value = '"' + value + '"'
            par[name] = value          # Filling the dictory
        self._params = par             # A control atribute, actually not used, good for debbuging
        for name in par:               # Iterating over the dictionary
            exec("self."+name.lower()+"="+par[name]) #Making the atributes at runtime




## Making the Mesh Class
class Mesh():
    """
    Mesh class, for keeping all the mesh data.
    Input: directory [string] -> this is where the domain files are.
    """
    def __init__(self, output_folder=""):
        if len(output_folder) > 1:
            if output_folder[-1] != '/':
                output_folder += '/'
        try:
            domain_x = np.loadtxt(output_folder+"domain_x.dat")
        except IOError:
            print("IOError with domain_x.dat")
        try:
            #We avoid ghost cells
            domain_y = np.loadtxt(output_folder+"domain_y.dat")[3:-3]
        except IOError:
            print ("IOError with domain_y.dat")
        self.xm = domain_x #X-Edge
        self.ym = domain_y #Y-Edge
        
        self.xc = 0.5*(domain_x[:-1] + domain_x[1:]) #X-Center
        self.yc = 0.5*(domain_y[:-1] + domain_y[1:]) #Y-Center
        
#         #(Surfaces computed from edge to edge)
#         #We first build 2D arrays for x & y (theta,r)
#         T,R = meshgrid(self.xm, self.ym)
#         R2  = R*R
#         self.surf = 0.5*(T[:-1,1:]-T[:-1,:-1])*(R2[1:,:-1]-R2[:-1,:-1])


## Making a Visualization Class 

class fargo_visualization(object):
    
    """
    Class for reading the simulation parameters.
    input: string -> name of the output folder, normally fargo_multifluid
          : string -> path the output folder
    """
    ## Global definitation
    Earth_mass  = 3e-6
    
    def __init__(self,output_folder):
        Mesh.__init__(self, output_folder)       #All the Mesh attributes inside Field (## Setting up the coordinates )
        Parameters.__init__(self, output_folder) #All the Parameters attributes inside Field    
        
        
        self.output_folder = output_folder
        self.number_of_outputs = int(self.ntot/self.ninterm) ## Total number of outputs
        self.sample_number = int(output_folder.split('_')[-1]) ## this picks the file number

        # self.number_of_orbit = int(self.ntot/20) ## Total number of orbits (ninterm=20 is 1 orbit)
#         print(self.output_folder[-3:-1])
        print (self.nx, self.ny, self.omegaframe,self.alpha, self.epsilon1, self.sigmaslope, self.sigma0,self.planetmass)
                

    def _initialize(self,fluid,fluid_property,output=None):
        
        '''  
        Returns: initial_density,mean_intial_density,density,mean_density
        
        1. Method to read the intial density and calculate the azimuathally averaged density
        2. Also to read the density and it's azimuthal average at each output
        
        '''
        
        intial_file = self.output_folder +'/' + fluid + fluid_property +str(0)+".dat"
        
        initial_density = np.fromfile(intial_file).reshape(self.ny,self.nx) ## initial gas density P.ny * P.nx
        mean_intial_density = initial_density.mean(axis=1)
        
        if output == None:
            output = self.number_of_outputs ## this automatically picks the last output    
            
        orbit = int(output*self.ninterm/20) #(ninterm=20 is 1 orbit)
        
        filename = self.output_folder +'/' + fluid + fluid_property + str(output)+".dat"
            
        density = np.fromfile(filename).reshape(self.ny,self.nx) ## evolving density 
        mean_density=density.mean(axis=1)
        
        return(initial_density,mean_intial_density,density,mean_density)
        
    def _density_plot(self,fluid,fluid_property,output= None,plot_type=None,fraction=None,ax=None):
        
        '''
         Input: 1. fluid type(dust or gas) 
                2. field type(dens, vy, vx) 
                Optional input
                3. Output at which orbit, by default it takes the last
                4. plot type (1D, 2D or 2D polar), by default it takes 1D
                5. fraction= to display the gapwidth if needed
                
         Output: density plot 
        '''
        
        ############## Creating output folders to save the output images##########
        current_directory = os.getcwd() 
        # print("the current directory in", current_directory)
        path= current_directory + '/movie_outputs_'+ fluid +'_'+ fluid_property
        path2 = current_directory + '/movie_folder'
                
        try:      
            os.makedirs(path)
        except OSError:  
            pass
#              print ("Creation of the directory %s failed/ not needed as it already exit" % path)
        else:  
            print ("Successfully created the directory %s" % path)

        #######################################################
               
        if output == None:
            output = self.number_of_outputs ## this automatically picks the last output    
            
        orbit = int(output*self.ninterm/20) #(ninterm=20 is 1 orbit)
        
        print("currently plotting output= ",output)
        
        ## Reading the intial value of the density and the evolving density
        initial_density,mean_intial_density,density,mean_density = self._initialize(fluid,fluid_property,output=output)
        
        
        file_number = self.output_folder[-3:-1]

        ## this is used for outputs from diffent models
#         outputfilename_figfile = os.path.join(path, fluid +'_dens_' +'%s'%file_number+ ".jpg")
        ## this is used for outputs from differnt orbits but same model
        outputfilename_figfile = os.path.join(path, fluid +'_dens_' +'%s'%output+ ".jpg")
        
        Planet_mass = int(self.planetmass/self.Earth_mass)
        
        if plot_type == None or plot_type == "1D":
            fig = plt.figure(figsize=(5, 5))
        
            plt.plot(self.yc,mean_intial_density,label= "Initial %s"%fluid + " Density" + " res: %s"%self.ny +r"$\times$ %s"%self.nx )  
            plt.plot(self.yc,mean_density,linestyle='--',label=  r" $M_{P}$= %s"%Planet_mass+r" $M_{\rm {Earth}}$" + " Orbit = %s"%orbit)  
#             plt.plot(yc,mean_intial_density,linestyle='-',label= "Evolving %s"%fluid + " Density"+" Orbit = %s"%number) ## plotting the radial profile
            plt.title ("Azimuthally averaged %s"%fluid+ " density profile")       
            plt.legend(loc=4)
            plt.yscale('log')
            
            ## Setting the limits of the y-axis 
            y_upper = 10*self.sigmagas
            y_lower = self.sigmagas/10
            
            if fluid == "gas":
                plt.ylim(y_lower,y_upper)
            else:
                plt.yscale('log')
#                 plt.ylim(y_lower*self.epsilon1,y_upper*self.epsilon1)    
            plt.ylabel(r" $\Sigma(r)$", fontsize=15)
        #         plt.ylabel(r"Log $\Delta \Sigma/\Sigma_{0}$", fontsize=15)
            plt.xlabel(r"Radius", fontsize=15)
            plt.tight_layout()
            plt.savefig(outputfilename_figfile,format='jpg',dpi=300)
            plt.close()
            
        if plot_type == "2D":
            fig = plt.figure(figsize=(5, 5))
            imshow(np.log10(density),origin='lower', aspect='auto',extent=[self.xmin,self.xmax,self.ymin,self.ymax])
#             plt.title("Fluid = %s"%fluid + "  Orbit = %s"%number+ " res: %s"%P.ny +r"$\times$ %s"%P.nx, fontsize=12,color='Black')
            plt.ylabel(r"Radius", fontsize=15)
            plt.xlabel(r"$\varphi$")
            plt.yscale("log")
            plt.colorbar()
#             cbar.set_ylabel(r"Log $\Sigma$")
            plt.savefig(outputfilename_figfile,format='jpg',dpi=300)
            plt.close()
            
            
        if plot_type == "2D_polar":            
            r, theta = np.meshgrid(self.xc, self.yc)
            self.data = np.log10(density)
            norm = matplotlib.colors.Normalize(vmin = np.min(self.data), vmax = np.max(self.data), clip = False)
            x = self.yc[:, np.newaxis] * np.cos(self.xc)
            y = self.yc[:, np.newaxis] * np.sin(self.xc)

            # add one additional column to make the contour plot complete
            x = np.hstack([x, x[:, :1]])
            y = np.hstack([y, y[:, :1]])
            data = np.hstack([self.data, self.data[:, :1]])

            xlabel, ylabel = "x", "y"
            if ax is None: ## for single plot
                plt.figure(figsize=(5, 5))
                ax = plt.gca()  
            
#             ax = plt.subplot(111)
#             cax = ax.contourf(x, y, data,100,rstride=1, cstride=1, norm=norm, alpha=1, cmap=plt.cm.Spectral)
            fig = ax.contourf(x, y, data,100, norm=norm, alpha=1, cmap='YlOrBr')  
            cbar = plt.colorbar(fig,orientation = "vertical", fraction=0.046, pad=0.04,ax=ax)
            ax.axis('equal')
            ax.axis('square')
            
            ax.set_xlabel(xlabel,fontsize=16)
            ax.set_ylabel(ylabel,fontsize=16)
            
            ## additonal info if needed for plotting optional
            Planet_mass = int(self.planetmass/self.Earth_mass)
            gap_list_, number ,gap_depth  = self._get_the_disk_gap(fluid,fluid_property,fraction=fraction,output=output)
            # print(gap_list_)
            ax.set_title(str(fluid[0:4]).upper() + " Gap-widths =%s"% np.round(gap_list_,3)+ " $M_{P}$= %s"%Planet_mass+ r" $M_{\rm E}$"+ " Orbit = %s"%orbit, fontsize=10,color='Black')
            plt.tight_layout()

            current_directory = os.getcwd() 
            if fluid == 'gas':
                path = current_directory + '/analysis_output/' + 'Disk_gas_plots' # introduced for the traning sets
            else:                
                path = current_directory + '/analysis_output/' + 'Disk_dust_plots'
            fig_D = os.path.join(path, fluid + '_gap_' + '%s' % self.sample_number + '.jpg')
            # ax.legend(loc=4, fontsize=6)
            # plt.tight_layout()            
            plt.savefig(fig_D,format='jpg',dpi=300)
            plt.close()
            # plt.savefig(outputfilename_figfile,format='jpg',dpi=300)
#             plt.close()

        return fig
    
   
        
    def _get_the_disk_gap(self,fluid,fluid_property,fraction=None,output= None,plot=None, ax=None):
        
        initial_density, mean_intial_density, density, mean_density = self._initialize(fluid, fluid_property, output=output)

        if output == None:
            output = self.number_of_outputs  # this automatically picks the last output  
        orbit = int(output * self.ninterm / 20)  # (ninterm=20 is 1 orbit)

        if fraction == None:
            fraction = 0.5  # by default the gap is calcualted at 0.5 * the intial density
            
        fraction_with_respect_to_initial_density = mean_intial_density * fraction 
#         print(fraction_with_respect_to_initial_density)
        index_above_radius1 = np.asarray(np.nonzero(self.yc >= 1))
        dens = mean_density 
        
        ## Checking for multiplicity, if the density at any radius (here the position of the planet)
        ## is greater than the comparison density then multiple gaps:
#         radius_in_list =[]
#         radius_out_list = []

        
#         ## if there multiple gaps gap the following is used     
#         if dens[min(min(index_above_radius1))] >= fraction_with_respect_to_initial_density[min(min(index_above_radius1))]:  
# #             print(dens[min(min(index_above_radius1))])
#             radius_list = []# defining empty lists  
#             for index in range(min(min(index_above_radius1)), 0, -1):
#                 if dens[index] <= fraction_with_respect_to_initial_density[index]:                    
#                     radius_list.append(self.yc[index])
#             if len(radius_list) != 0: ## if there are no values below the given density fraction               
#                 radius_in_list.append(min(radius_list)) 
#                 radius_out_list.append(max(radius_list))
#             radius_list = []                     
#             for index in range(min(min(index_above_radius1)), len(self.yc)):                
#                 if dens[index] <= fraction_with_respect_to_initial_density[index]:                    
#                     radius_list.append(self.yc[index])
#             if len(radius_list) != 0: ## if there are no values below the given density fraction               
#                 radius_in_list.append(min(radius_list)) 
#                 radius_out_list.append(max(radius_list))
            
        
#         else: ## if there is only one gap the following is used     
#             for index in range(min(min(index_above_radius1)), len(self.yc)): # index above r = 1 for density above cutoff
#                 if dens[index] >= fraction_with_respect_to_initial_density[index]:  # for density greater than the given fraction the intial density                    
#                     radius_out = self.yc[index]
# #                     print("radius_out",radius_out)
#                     break
#             radius_out_list.append(radius_out)
#             for index in range(min(min(index_above_radius1)), 0, -1):        # index below r = 1 for density above cutoff            
#                 if dens[index] >= fraction_with_respect_to_initial_density[index]:  # for density greater than the given fraction the intial densit                   
#                     radius_in = self.yc[index]
# #                     print("radius_in",radius_in)
#                     break       
#             radius_in_list.append(radius_in) 
        radius_in_list =[] # radius inner wall of gap
        radius_out_list = [] # radius outer wall of gap
        gap_depth = [] # to store the values of gap depth
        gap_depth_fraction = [] # to store the depth fraction wrt to the initial density
        radius_at_gap_depth = [] # to store the values radius of gap depth
        
        '''Algorithm One
           Used to calculate all the gap widths and depth
        '''
        
        index_start = 0 # Default starting index
        index_intersections = [] # empty list to store the radius index of gaps outer and inner edge
        while index_start <= len(self.yc):   # a while loop to find all the gaps 
        ## the condition to identify a gap is when the mean density goes below a defined value
            for index in range(index_start, len(self.yc), 1):
                if dens[index] <= fraction_with_respect_to_initial_density[index]: # then the density below the fraction                
                    radius_in_list.append(self.yc[index]) # we identify the radius as radius_in
                    index_intersections.append(index) # we also store the index value for later use
                    index_start = index
#                     print(index)
                    break 
                else :  # if density is above the fraction the i.e no radius value the while loop is broken
                    index_start = len(self.yc) + 1 # to break the loop as there are no gaps (the staring index in modified)
            for index in range(index_start, len(self.yc), 1):                 
                if dens[index] >= fraction_with_respect_to_initial_density[index]:                    
                    radius_out_list.append(self.yc[index])
                    index_intersections.append(index)
                    index_start = index
                    break                  
                else :  # if index_start is beyong range we continue to break the loop
                    index_start = len(self.yc) + 1 # to break the loop as there are no gaps
        
        ## Calculation the gap depth and the radius at gap depth
        for gap_numbers in range(0, len(index_intersections),2): ## identifing the radius index at gap edges
            den_list = [] # empty list o store the density values in the gap
            local_depth = [] # empty list to store the local gap depth
            index_at_gap_depth = [] # empty list to store the index at the gap depth
            ## finding the density between the gaps and findng its min
            for index in range(index_intersections[gap_numbers], index_intersections[gap_numbers+1],1): 
                den_list.append(dens[index]) ## all the density values in the gap
            gap_depth.append(min(den_list)) #  the gap depth minimum value
            local_depth.append(min(den_list)) # storing the local value for comparison
            index_at_gap_depth = np.asarray(np.nonzero(mean_density == local_depth)) # index at gap_depth
            radius_at_gap_depth.append(self.yc[min(min(index_at_gap_depth))]) # radius at depth
            gap_depth_fraction.append(local_depth/mean_intial_density[min(min(index_at_gap_depth))])

        Planet_mass = int(self.planetmass/self.Earth_mass)        
#         print(radius_in_list,radius_out_list)
        gap_list = np.asarray(radius_out_list) - np.asarray(radius_in_list)
        multiplicity = len(radius_in_list)
        if multiplicity == 0: ## of there are no gaps
            gap_list =np.zeros(2) ## fill zero
            gap_depth = np.zeros(1)
            # radius_at_gap_depth = np.zeros(1)
        
        ## for plotting the density profile along with the gaps
        if plot!=None:
            if ax is None: ## for single plot
                plt.figure(figsize=(5, 5))
                ax = plt.gca()                
            fig = ax.plot(self.yc,mean_intial_density,label= "Initial %s"%fluid + " Density" + " res: %s"%self.ny +r"$\times$ %s"%self.nx )
            ax.plot(self.yc,mean_density,linestyle='--',label=  r" $M_{P}$= %s"%Planet_mass+r" $M_{\rm {E}}$" + " Orbit = %s"%orbit)
            for x in radius_in_list+radius_out_list:
                ax.axvline(x,linestyle='--',linewidth=0.5, color='k')
            if multiplicity >=1: 
                ax.scatter(radius_at_gap_depth, gap_depth,marker='x')

            ax.plot(self.yc,fraction_with_respect_to_initial_density,linestyle='-.',linewidth=0.5,label= str(fraction*100)+r" percent of Initial %s"%fluid + " Density")
            ax.set_yscale('log')
            ax.set_ylabel(r" $\Sigma(r)$", fontsize=15)
            ax.set_xlabel(r"$R/R_{\rm 0}$", fontsize=15)
            ax.set_title(str(fluid[0:4]).upper() + " Gap-widths =%s"% np.round(gap_list,3)+r"$R_{\rm 0}$", fontsize=15,color='Black')
            current_directory = os.getcwd() 
            if fluid == 'gas':
                path = current_directory + '/analysis_output/' + 'gas_gap' # introduced for the traning sets
            else:                
                path = current_directory + '/analysis_output/' + 'dust_gap'
            fig1 = os.path.join(path, fluid + '_gap_' + '%s' % self.sample_number + '.jpg')
            ax.legend(loc=4, fontsize=6)
            plt.tight_layout()            
            plt.savefig(fig1,format='jpg',dpi=300)
            plt.close()
        return gap_list, multiplicity, gap_depth 
        
    def _movie_density_plot(self,fluid,fluid_property,output=None,plot_type=None,frequency=None):
        
        if frequency == None:
            interval = 10 ## the default frequency
        else: 
            interval = frequency
        start = time.time()
        print("Plotting now...")
        
        
        for number in range(0,self.number_of_outputs+1,interval):
            self._density_plot(fluid,fluid_property,output=number,plot_type=plot_type)
        print("End of simulation")    
        end = time.time()
        T = end-start
        print("Total time for parallel run %.2f"%T, 'sec') 
        
# #         start = time.perf_counter()    
#         numbers = [number for number in range(0,self.number_of_outputs,interval)]      
#         def do_something(number):
#              return number +1            
# #             self._density_plot(fluid,fluid_property,output=number,plot_type=plot_type)
#             with concurrent.futures.ProcessPoolExecutor() as executor:                
#                 results = executor.map(do_something, numbers)        
# #             for result in results:
# #                 print(result)
#         finish =time.perf_counter()

#         print(f'Finished in {round(finish-start,2)}seconds(s)')
        

    def _time_evol_disk_gap(self,fluid,fluid_property,fraction=None,frequency=None,ax=None):
        
        Planet_mass = int(self.planetmass/self.Earth_mass)
        if frequency == None:
            interval = 10 ## the default frequency
        else: 
            interval = frequency
        gap_list_1 = []
        gap_list_2 = []
        orbit_number = []
#         orbit = int(output*self.ninterm/20) #(ninterm=20 is 1 orbit)
        for number in range(0,self.number_of_outputs+1,interval):
            gap, num,gap_depth  = self._get_the_disk_gap(fluid,fluid_property,fraction=fraction,output=number)
            gap_list_1.append(gap[0])
            if num >1:
                gap_list_2.append(gap[1])
            else:
                gap_list_2.append(0)
            orbit_number.append(int(number*self.ninterm/20))
        if ax is None:
            ax = plt.gca()
        fig=ax.plot(orbit_number,gap_list_1,label= r" $M_{P}$= %s"%Planet_mass+r" $M_{\rm {E}}$") 
        ax.plot(orbit_number,gap_list_2,label= "2nd Gap width "+str(fraction)+r"$\times$ Initial %s"%fluid +" Density") 
        ax.legend(loc=4)
        ax.set_ylabel(r"$\Delta_{gap}/R_{\rm P}$", fontsize=12)
        ax.set_xlabel(r"No of orbits", fontsize=12)
        current_directory = os.getcwd() 
        if fluid == 'gas':
            path = current_directory + '/analysis_output/' + 'gas_gap_evol' # introduced for the traning sets
        else:                
            path = current_directory + '/analysis_output/' + 'dust_gap_evol'
        figname = os.path.join(path, fluid + '_gap_' + '%s' % self.sample_number + '.jpg')
        plt.tight_layout() 
        plt.savefig(figname,format='jpg',dpi=300)
        plt.close()
        return fig







