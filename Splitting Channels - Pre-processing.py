
# coding: utf-8

# ## Gastruloid Pipeline

# ### Importing modules

# In[ ]:


get_ipython().magic(u'matplotlib inline')

import os
import tifffile as tiff
from scipy import ndimage
import glob
import numpy as np
import pandas as pd
import sys


# In[ ]:


## Put the location of the 16-bit or 8-bit(following Fiji processing) images below; it will create new subfolders
#with each channel name as the name, and split each channel into its equivalent folder 

data_locn = "/Users/alexandrabaranowski/Desktop/Cells for pipeline analysis/5 Background Normalisation/Signalling/"
# it doesn't accept file names that have dots in them so changed the 3.25 & 0.5 to dashes 
# Make sure all file names with single digits are named 00, 01, etc otherwise they will be ordered wrong 

folders = False
blur = 0
channels = [0, 1,2,3]
files_list = []
tiff_list = []
############################# generating file list ##############################

if folders == True:  
    files_list.append(glob.glob(data_locn+'*/*.tif'))
else:
    files_list.append(glob.glob(data_locn+'*tif'))
    

for i, v in enumerate(files_list[0]):
    tiff_list.append(tiff.imread(v))


# In[ ]:


filenames = []
# no blur here #
for channel in channels:
    #print channel
    blurred_list = []
    for i in tiff_list:
        blurred_list.append(ndimage.gaussian_filter(i[channel], sigma=blur)) 

############################### saving the image ################################
    
    for i in files_list[0]:
        base = os.path.basename(i)
        base2 = os.path.splitext(base)[0]
        filenames.append(os.path.splitext(base2)[0])
    
    if not os.path.exists(data_locn+'Ch' + str(channel) + '/'):
        os.makedirs(data_locn+'Ch' + str(channel) + '/')
    
    for i,v in enumerate(blurred_list):
        tiff.imsave(data_locn+'Ch' + str(channel) + '/' +str(filenames[i]) + 'ch_' + str(channel) + ' .tif', v) 

