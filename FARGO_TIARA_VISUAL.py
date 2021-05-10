#!/usr/bin/env python
# coding: utf-8

# In[1]:


###### This script is used to vizualise the data directly in the TIARA ##### 

# Author: sayantan
# Date : 22 December 2020

## importing the class module developed by Sayantan
if __name__ ==  '__main__': ## This was implemented to fix the error while running parallel
	import Fargo_visual_class as FV

	################################### User defined inputs ##################

	#folder="Training_25April_21/TS_001_018/outputs/fargo_multifluid_new_1/"
	folder = "Training_25April_21/TS_001_018/outputs/fargo_multifluid_new_10/"
	#folder ="Training_3May_21/TS_000_020/outputs/fargo_multifluid_new_1/"
	fluid_list = ['gas','dust1']
	fluid_property_list = ['dens','vy','vx']

	#######please select the particular fluid and the field from the list##########

	f=1;p=0
	fluid = fluid_list[f]
	fluid_property= fluid_property_list[p]
	fraction = 0.5
	V = FV.fargo_visualization(output_folder=folder)



	# V._density_plot(fluid,fluid_property,plot_type="1D",output=None,fraction=fraction,ax=None)

	V._time_lapse(fluid,fluid_property,plot_type="2D_polar",frequency=1,gap_type=None,Parallel=True)

#	V._time_lapse(fluid,fluid_property,plot_type=None,frequency=10,gap_type=None,Parallel=True) 
#	V._make_movie(fluid,fluid_property,plot_type="2D_polar",frequency=None,gap_type=None,Parallel=True)







