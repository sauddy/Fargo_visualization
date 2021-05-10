#!/usr/bin/env python
# coding: utf-8

# ## Modified code to generate simple images for CNN analysis

# In[1]:


## This script is dev for the purpose of generating images for the CNN analysis. 
## The important addtions are that that 2D polar images are now very simple.
## The other alternative is to generate the images simple 2D
## Also more images a simple alternative is to generate time lapse images when the system is stable

## modification : 23 September 2020 : We are generating the time lapse images and the corresponding 
## CSV files


####### Importing the required modules

import os
import glob
import csv
import numpy as np

## importing the class module developed by Sayantan
import Fargo_visual_class as FV


# In[2]:


## Please Specify then output file ## 
output_no = 96
output_orbit = 2400
current_directory = os.getcwd() 
# print("the current directory in", current_directory)
path= current_directory + '/analysis_output_' + str(output_orbit)
gas_gap = path + '/gas_gap'
dust_gap = path + '/dust_gap'
gas_gap_evol = path + '/gas_gap_evol'
dust_gap_evol = path + '/dust_gap_evol'
Disk_gas_plots = path + '/Disk_gas_plots'
Disk_dust_plots = path + '/Disk_dust_plots'

output_folder_list = [path,gas_gap,dust_gap,gas_gap_evol,dust_gap_evol,Disk_gas_plots,Disk_dust_plots]
for file in output_folder_list:
    try:      
        os.makedirs(file)
    except OSError:  
        print ("Creation of the directory %s failed/ not needed as it already exit" % file)
    else:  
        print ("Successfully created the directory %s" % file)
        


# In[3]:


########### User input is the folder address ############### the rest is automatic #########
# folder_address= '/Users/sayantan/Desktop/Programming/Dusty_Disk_Gap/Training_set/'
#folder_address = '/Volumes/Auddy_Akiti/ML_data/DISK_gap_correction_21June_20/Training_set_18June20/'
# folder_address='/Volumes/My_Seagate/Training_set/'
#list_of_training_set = glob.glob(folder_address+'TS_*')
#list_TS_sorted_folder =sorted(list_of_training_set, key=lambda x: int(x.replace(folder_address+'TS_', '')))
# print(list_TS_sorted_folder)
#list_TS_sorted_folder = ['High_mass_runs'] ## added 31 dec 2020 for test data fo M_Earth > 120


folder_address= 'Training_22Jan-21/'
list_of_training_set = glob.glob(folder_address+'TS_*')
list_TS_sorted_folder =sorted(list_of_training_set, key=lambda x: int(x.replace(folder_address+'TS_', '')))
print(list_TS_sorted_folder)

print(list_TS_sorted_folder)

# In[4]:


def data_analysis(sample):
    '''
    Input : the folder name of the traning set
    output: csv with the gap properties
            folder with dust 1D, 2D images at the final orbit
            folder with gas 1D, 2D images
            folder for dust gap profile evolution
            folder with gas gap profile evolution
    
    
    '''
    output_directory = sample+'/outputs/'
    list_of_folder = glob.glob(output_directory+'fargo_multifluid_new*')
    list_sorted_folder =sorted(list_of_folder, key=lambda x: int(x.replace(output_directory+'fargo_multifluid_new_', '')))
#     print(list_sorted_folder)
    ## for naming the .csv files according to the Training folder
    file_start = sample.split('_')[-2]
#
    file_end = sample.split('_')[-1]

    output_filename = os.path.join(path,'Disk_gap_param_'+ file_start+'_'+file_end +'.csv')
    ## for the high mass data the new name : 31 December 2020
#    output_filename = os.path.join(path,'Disk_gap_param_' +'.csv')
    csv_file = open(output_filename,'w')
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(['Sample#','Planet_Mass','Epsilon', 'Alpha', 'Stokes','Aspect_Ratio','SigmaSlope',
                         'Dust_gap_1', 'Gas_gap_1','Dust_depth_1', 'Gas_depth_1','Dust_gap_2','Dust_depth_2',
                         '#_DG','#_GG'])


    for folder in list_sorted_folder:
        print(folder)
        dust_fraction = 0.5 # for dong this only helps  guide what contrast is a gap
        gas_fraction = 0.5  # for dong this only helps  guide what contrast is a gap
        gap_type = None
        sample_number = int(folder.split('_')[-1])

        ## calling the fargo visualiazation class
        V = FV.fargo_visualization(output_folder=folder,output_no = output_orbit)


        ## estimating the dust and the gas gap and multiplicity(if any)
        dust_gap, dust_num, dust_depth  =  V._get_the_disk_gap('dust1','dens',fraction=dust_fraction,output=output_no,gap_type=gap_type, plot=True)
        gas_gap, gas_num,gas_depth  =  V._get_the_disk_gap('gas','dens',fraction=gas_fraction,output=output_no,gap_type=gap_type, plot=True)
        
        print(dust_gap, dust_num, dust_depth)
        print(gas_gap, gas_num,gas_depth)
        ## for creating the list with the gap widths
        dust_gap_list = np.zeros(2) ## filling zeros for filling the c.sv
        gas_gap_list =  np.zeros(2)
        dust_gap_list[0] = dust_gap[0]        
        gas_gap_list[0] = gas_gap[0]

        if len(dust_gap)!=1:
            dust_gap_list[1] = dust_gap[1]
        if len(gas_gap)!=1:
            gas_gap_list[1] = gas_gap[1]

#         print(dust_depth)
        ## for the moment we are interested in saving depth of one gap
        dust_depth_list = np.zeros(2)
        gas_depth_list = np.zeros(2)
        dust_depth_list[0] = dust_depth[0]
#         gas_depth_list[0] = gas_depth[0]

        if len(dust_depth)!=1:
            dust_depth_list[1] = dust_depth[1] 
        if len(gas_depth)!=0:
            gas_depth_list[0] = gas_depth[0]



#         ## plotting the gap evolution for user defined fraction
#         V._time_evol_disk_gap('dust1' ,'dens',fraction=dust_fraction,gap_type=gap_type,frequency=1)
#         V._time_evol_disk_gap('gas','dens',fraction=gas_fraction,gap_type=gap_type,frequency=1)

#         ## plotting the gas and dust disk images 
        V._density_plot('dust1','dens',output=output_no,plot_type="2D_polar", fraction=dust_fraction,gap_type=gap_type)
#         V._density_plot('gas','dens',output=output_no,plot_type="2D_polar", fraction=gas_fraction,gap_type=gap_type)


        print("dust_gap=",dust_gap)
        print("gas_gap=",gas_gap)
        print("dust_depth=",dust_depth)
        print("gas_depth=",gas_depth)

        csv_writer.writerow([sample_number,'%.2e'%V.planetmass,'%.2e'%V.epsilon1,'%.2e' %V.alpha,'%.2e'%V.ts1,'%.2e'%V.aspectratio
                              ,'%.2e'%V.sigmaslope,'%.3f' %dust_gap_list[0],'%.3f' %gas_gap_list[0],'%.2e' %dust_depth_list[0],
                             '%.2e' %gas_depth_list[0],'%.3f' %dust_gap_list[1],'%.2e'%dust_depth_list[1],'%.1f' %dust_num,'%.1f' %gas_num])

    csv_file.close()


# In[5]:


## this 1st try to run parallel. If not possible then it does in series
try: 
    import time
    import concurrent.futures ## intentionally changed futures to future
    run ="Parallel"
    print("running", run)
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
        sample = [sample for sample in list_TS_sorted_folder]
#         sample = [sample for sample in list_TS_sorted_folder[:12]]
#         print(sample)
        start = time.time()
        results = executor.map(data_analysis,sample)
except ImportError :
    run ="Series"
    print("running", run)
    start = time.time()
    for sample in list_TS_sorted_folder:
        data_analysis(sample)   

print("End of simulation")    
end = time.time()
T = end-start
print("Total time for",run,"run %.2f"%T, 'sec') 


# In[ ]:




