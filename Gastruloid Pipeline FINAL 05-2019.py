
# coding: utf-8

# In[ ]:


get_ipython().magic(u'matplotlib inline')

import os
import tifffile as tiff
from scipy import ndimage
from scipy import stats
from scipy.integrate import trapz
import glob
import cv2
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
import sys
from PIL import Image
import seaborn as sns
import mahotas


# ### Retrieving the files (16 and 8bit, unblurred, separated channels)
# ##### 16-bit are used to calculate fluorescence, 8-bit for all other measurements 

# In[ ]:


## location of the 8-bit split channels, per channel

data_locn_mask = ["/Users/alexandrabaranowski/Desktop/Cells for pipeline analysis/3 hGlds Signalling modulation/2019-03-04_Wnt3aPT_05and3Chi_400c_20x_72h_1_8bit/Ch0/", 
                  "/Users/alexandrabaranowski/Desktop/Cells for pipeline analysis/3 hGlds Signalling modulation/2019-03-04_Wnt3aPT_05and3Chi_400c_20x_72h_1_8bit/Ch1/", 
                 "/Users/alexandrabaranowski/Desktop/Cells for pipeline analysis/3 hGlds Signalling modulation/2019-03-04_Wnt3aPT_05and3Chi_400c_20x_72h_1_8bit/Ch2/", 
                  "/Users/alexandrabaranowski/Desktop/Cells for pipeline analysis/3 hGlds Signalling modulation/2019-03-04_Wnt3aPT_05and3Chi_400c_20x_72h_1_8bit/Ch3/"]

## location of the 16-bit split channels, per channel
data_locn_full = ["/Users/alexandrabaranowski/Desktop/Cells for pipeline analysis/3 hGlds Signalling modulation/2019-03-04_Wnt3aPT_05and3Chi_400c_20x_72h_1/Ch0/",
                 "/Users/alexandrabaranowski/Desktop/Cells for pipeline analysis/3 hGlds Signalling modulation/2019-03-04_Wnt3aPT_05and3Chi_400c_20x_72h_1/Ch1/",
                  "/Users/alexandrabaranowski/Desktop/Cells for pipeline analysis/3 hGlds Signalling modulation/2019-03-04_Wnt3aPT_05and3Chi_400c_20x_72h_1/Ch2/",
                 "/Users/alexandrabaranowski/Desktop/Cells for pipeline analysis/3 hGlds Signalling modulation/2019-03-04_Wnt3aPT_05and3Chi_400c_20x_72h_1/Ch3/"]

## location of the 16-bit empty image used for background subtraction, per channel
data_locn_bg = ["/Users/alexandrabaranowski/Desktop/Cells for pipeline analysis/5 Background Normalisation/Signalling/Ch0/",
                "/Users/alexandrabaranowski/Desktop/Cells for pipeline analysis/5 Background Normalisation/Signalling/Ch1/",
               "/Users/alexandrabaranowski/Desktop/Cells for pipeline analysis/5 Background Normalisation/Signalling/Ch2/",
               "/Users/alexandrabaranowski/Desktop/Cells for pipeline analysis/5 Background Normalisation/Signalling/Ch3/"]


# In[ ]:


channels = [0,1,2,3]


# In[ ]:


#### Retrieving the images from the location and ordering them in ascending order (N.B. Make sure single digit 
# numbers are labeled 00, 01, 02 etc otherwise they will not be ordered correctly)

files_list0 = []; files_list1 = []; files_list2 = []; files_list3 = []
files_list_full0 = []; files_list_full1 = []; files_list_full2 = []; files_list_full3 = []
bg_list0 = []; bg_list1 = []; bg_list2 = []; bg_list3 = []; 

files_list = [files_list0, files_list1, files_list2, files_list3]
files_list_full = [files_list_full0, files_list_full1, files_list_full2, files_list_full3]
bg_list = [bg_list0, bg_list1, bg_list2, bg_list3]

for i in range(0, len(channels)):
    files_list[i].append(glob.glob(data_locn_mask[i]+'*tif'))
    files_list[i][0].sort()
    
    files_list_full[i].append(glob.glob(data_locn_full[i]+'*tif'))
    files_list_full[i][0].sort()
        
    bg_list[i].append(glob.glob(data_locn_bg[i]+'*tif'))
    bg_list[i][0].sort()


# In[ ]:


filenames0 = []; filenames1 = []; filenames2 = []; filenames3 = []
filenames_full0 = []; filenames_full1 = []; filenames_full2 = []; filenames_full3 = []

filenames = [filenames0, filenames1, filenames2, filenames3]
filenames_full = [filenames_full0, filenames_full1, filenames_full2, filenames_full3]

for j in range(0,len(channels)):
    for i in files_list[j][0]:
        base = os.path.basename(i)
        base2 = os.path.splitext(base)[0]
        filenames[j].append(os.path.splitext(base2)[0])

for j in range(0,len(channels)):
    for i in files_list_full[j][0]:
        base = os.path.basename(i)
        base2 = os.path.splitext(base)[0]
        filenames_full[j].append(os.path.splitext(base2)[0])


# In[ ]:


## Check that the list names are int the correct order - if the files were not named 00, 01 etc. they will not be

print files_list[0]
files_list_full[1]


# #### Every channel is processed separately, starting with the empty channel; the pipeline follows the format of: 
# 1) Removal of noise <br>
# 2) Thresholding <br>
# 3) Feature Extraction <br>
# 4) QC to ensure all feature lists are of the correct length 
# #### Fluorescence measurements are taken together for all channels 
# 5) Fluorescence Noise Removal, Normalisation and Measurements <br>
# 6) Proportion of gastruloid covered by each channel
# #### Then, general QC is performed to remove images without gastruloids, with debris etc.
# 7) Different methods are used and evaluated <br>
# 8) Based on the one with best results, those indices are removed from all other lists <br>
# 9) Extract proportions and proportion of overlap measurements <br>
# 10) Features are compiled into a dataframe, ready for export to R <br>

# ## Empty Channel - Removing noise & Thresholding

# 1) De-noise

# In[ ]:


get_ipython().magic(u'matplotlib inline')
matplotlib.rcParams.update({'figure.max_open_warning': 0})

blurred_0 = []
image_copies_0 = []
masks_0 = []
blurred_masks_0 = []
threshold_0 = []

for i, v in enumerate(files_list[0][0]):
    # read the original image
    image_0 = np.array(Image.open(files_list[0][0][i]))

    
    #Make a copy to prevent distortion of the original, bc all the transformations change the image
    image_copy = image_0.copy()
    image_copies_0.append(image_copy)
    image_copy2 = image_0.copy()

    
    #Gaussian Blur 

    blur = cv2.GaussianBlur(image_copy,(25,25),3) # change this depending on the image 
    #blur2 = cv2.GaussianBlur(image_copy,(35,35),0) 
    #diff = blur - blur2
    #imgf = image_copy + diff 
    blurred_0.append(blur)
    
    #Set initial threshold 
    ret, th0 = cv2.threshold(blur, 100, 255, cv2.THRESH_BINARY_INV+cv2.THRESH_OTSU)
    threshold_0.append(th0)
    mask_inv= cv2.bitwise_not(th0)
    
    
    #blurred mask of the gastruloid - for the fluorescence
    blurred_masked_img = cv2.bitwise_and(blur,blur, mask = th0)

    # bluured mask of the background - to remove bg noise 
    #blurred_bg_mask = cv2.bitwise_and(blur,blur, mask = mask_inv)      
    blurred_masks_0.append(blurred_masked_img)

    #fig, ax = plt.subplots(1, 4, figsize = (10,3)) 
    #ax[0].imshow(image_0)
    #ax[1].imshow(th0)
    #ax[2].imshow(th1)
    #ax[3].imshow(imgf)
    #ax[3].imshow(diff)


# In[ ]:


### Histogram-based determination of the correct threshold for segmentation ###
### If the image after Gaussian blur looks bimodal, then use Otsu's ###
# Otsu's will set threshold automatically; if not, use historgrams to determine threshold to use #

# https://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_imgproc/py_thresholding/py_thresholding.html
for i, v in enumerate(blurred_0[0:6]): # choose a subset of the images or all of them
    
    ## finding the range of the grayscale values so we can do the thresholding properly 
    # https://mmeysenburg.github.io/image-processing/05-creating-histograms/ 
    hist = cv2.calcHist([v], [0], None, [256], [0, 256])
    histogram2=cv2.calcHist(image_copies_0[i], [0], None, [256], [0, 256])
    fig, ax = plt.subplots(1,2,figsize=(10,3))
    
    ax[0].plot(histogram2)
    ax[0].set_title("Grayscale Histogram") # Before de-noise
    ax[0].set_xlabel("Grayscale Value")
    ax[0].set_ylabel("Pixels")
    ax[1].plot(hist)
    ax[1].set_title("Grayscale Blurred Histogram") # After Gaussian Filter
    ax[1].set_xlabel("Grayscale Value")
    ax[1].set_ylabel("Pixels")

plt.tight_layout()
plt.show()

save_locn = "/Users/alexandrabaranowski/Desktop/Images for report /"
#fig.savefig(save_locn+'1. 2015BraGFP 96h Historgram Gaussian3.png')


# 2) Thresholding

# when you have identified the threshold, identify the gastruloid as the entity within the frame with the largest area

# In[ ]:


from scipy.interpolate import splprep, splev
im_out_ = [] # stores the initial floodfilled mask - this will suffice for some gastruloids
contours_0rough = [] # only needed if contour smoothing is applied
contours_0 = [] # stores the contours for downstream use
erosions = [] # stores the mask after erosions to remove debris surrounding gastruloids 
smoothed_0 = [] # stores the smoothed contours

for i, v in enumerate(image_copies_0):
    blurred_copy = v.copy()
    blurred_copy1 = v.copy()
    image_copy = v.copy()
    v2 = v.copy()
    
    #blur = cv2.GaussianBlur(image_copy, (35,35), 0)
    
    # Try DoG (imgf) instead of simple Gaussian when Glds have debris around them
    blur = cv2.GaussianBlur(image_copy,(25,25),3) 
    #blur2 = cv2.GaussianBlur(image_copy,(35,35),0) 
    #diff = blur - blur2
    #imgf = image_copy + diff
    
    
    # The threshold can be kept at 255 bc the Otsu's binarisation algorithm finds the optimal threshold for each image
    # this is returned as the ret value 
    ret, th01 = cv2.threshold(blur, 0, 255, cv2.THRESH_BINARY_INV+cv2.THRESH_OTSU)
    mask_inv = cv2.bitwise_not(th0)
    th0 = cv2.erode(th01, (5,5) ,iterations =1)
    
    im_floodfill_0 = th0.copy()
    
    # Mask used for flood filling; Notice the size needs to be 2 pixels larger than the image.
    h, w, = th0.shape[:2]
    mask = np.zeros((h+2, w+2), np.uint8)

    # Floodfill from point (0, 0)
    cv2.floodFill(im_floodfill_0, mask, (0,0), 255)

    # Invert floodfilled image
    im_floodfill_inv_0 = cv2.bitwise_not(im_floodfill_0)

    # Combine the two images to get the foreground
    im_out_0 = th0 | im_floodfill_inv_0
    im_out_.append(im_out_0)
    
    #### Cleaning the mask - some gastruloids will NOT require this (e.g., if they are clearly distinguishable from
    # background, have smooth edges etc)
    kernel = np.ones((25,25),np.uint8)
    erosion = cv2.erode(im_out_0, kernel ,iterations =1)
    #erosion = cv2.morphologyEx(im_out, cv2.MORPH_OPEN, kernel, iterations = 1)
    #erosion = cv2.morphologyEx(erosion, cv2.MORPH_OPEN, kernel, iterations = 1)
    erosions.append(erosion)
    
    #alternative method to erosion - erosion followed by dilation 
    #opening = cv2.morphologyEx(th0, cv2.MORPH_OPEN, kernel, iterations = 2)
    #opening = cv2.erode(opening, kernel, iterations = 1)
    
    # Finding the contours for each method - this is for comparison; choose the one that best delineates the gastruloids
    (_, contours, _) = cv2.findContours(erosion, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    (_, contours2, _) = cv2.findContours(erosion, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    (_, contours3, _) = cv2.findContours(im_out_0, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    (_, contours4, _) = cv2.findContours(th0, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    contours2 = sorted(contours2, key = cv2.contourArea, reverse = True)[:1]
    contours3 = sorted(contours3, key = cv2.contourArea, reverse = True)[:1]
    contours = sorted(contours, key = cv2.contourArea, reverse = True)[:1]
    
    if not contours2: 
        contours2 = np.array([0])
    contours_0rough.append(contours2)  
    
    if not contours3: 
        contours3 = np.array([0])
    smoothed_0.append(contours3)
    
    # ---- add smoothing code excerpt here ------
    # -------------------------------------------
    
    v3 = v.copy()
    cv2.drawContours(v3, contours2, -1, (255, 0, 0), 10)     # draw contours over original image
    cv2.drawContours(v2, contours3, -1, (255, 0, 0), 10)     # draw contours over original image
    #cv2.drawContours(v, contours4, -1, (255, 0, 0), 10)     # draw contours over original image

    #plot it 
    fig, ax = plt.subplots(1, 6, figsize = (16,5))
    ax[0].imshow(th0)
    ax[0].set_xlabel("Distance / Pixels")
    ax[0].set_ylabel("Distance / Pixels")
    ax[0].set_title("Before Floodfilling")
    
    ax[1].imshow(im_out_0)
    ax[1].set_xlabel("Distance / Pixels")
    #ax[1].set_ylabel("Distance / Pixels")
    ax[1].set_title("After Floodfilling")
    
    ax[2].imshow(erosion) 
    ax[2].set_xlabel("Distance / Pixels")
    #ax[2].set_ylabel("Distance / Pixels")
    ax[2].set_title("After Floodfilling & Erosion ")
    
    ax[3].imshow(v)
    ax[3].set_xlabel("Distance / Pixels")
    #ax[3].set_ylabel("Distance / Pixels")
    ax[3].set_title("No Floodfill contour")
    
    ax[4].imshow(v2)
    ax[4].set_xlabel("Distance / Pixels")
    #ax[4].set_ylabel("Distance / Pixels")
    ax[4].set_title("Floodfill contour")

    ax[5].imshow(v3)
    ax[5].set_xlabel("Distance / Pixels")
    #ax[4].set_ylabel("Distance / Pixels")
    ax[5].set_title("Floodfill & Erosion contour")
    
    
    plt.tight_layout()
    plt.show()
    
    
#save_locn = "/Users/alexandrabaranowski/Desktop/Images for report /"
#fig.savefig(save_locn+'1. 2015BraGFP 72h Empty.png')
 


# In[ ]:



################# this is for gastruloids that have unclear borders (e.g., dead cells) and so normal contouring 
#doesn't capture the edges well ###########
# adapted from: https://agniva.me/scipy/2016/10/25/contour-smoothing.html
# if used, ADD IT TO THE ABOVE CODE at the point of ------- add smoothing code excerpt here ---------
smoothed_0 = []
for i, v in enumerate(erosions):
v2 = v.copy()
(_, contours, _) = cv2.findContours(v, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
contours = sorted(contours, key = cv2.contourArea, reverse = True)[:1]
cv2.drawContours(v3, contours2, -1, (255, 0, 0), 10)     # draw contours over original image

for contour in contours:
    x,y = contour.T
    # Convert from numpy arrays to normal arrays
    x = x.tolist()[0]
    y = y.tolist()[0]
    
    # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.splprep.html
    tck, u = splprep([x,y], u=None, s=1.0, per=1)
    
    # https://docs.scipy.org/doc/numpy-1.10.1/reference/generated/numpy.linspace.html
    u_new = np.linspace(u.min(), u.max(), 25)
   
    # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.splev.html
    x_new, y_new = splev(u_new, tck, der=0)
    
    # Convert it back to numpy format for opencv to be able to display it
    res_array = [[[int(i[0]), int(i[1])]] for i in zip(x_new,y_new)]
    
    smoothed.append(np.asarray(res_array, dtype=np.int32))
        
    # Draw the smoothed contours on the original image
    cv2.drawContours(v2,smoothed, -1, (255, 0, 0), 10)
    if not smoothed: 
        smoothed = np.array([0])
    smoothed_0.append(smoothed)
            
    fig, ax = plt.subplots(1, 2, figsize = (10,3))
    ax[0].imshow(v)
    ax[0].imshow(v2)
    
    #msk = cv2.bitwise_and(blurred_copy,blurred_copy, mask = contours)


# 3) Feature Extraction

# In[1]:


areas_0 = []
areas0 = [] # this will only be used for QC
for i, v in enumerate(erosions):
    #find areas
    filledareas = cv2.countNonZero(v)
    areas0.append(filledareas)


for i, c in enumerate(smoothed_0):
    if np.size(c)>1:
        cnt = c[0]
        area = cv2.contourArea(cnt)
        areas_0.append(area)
    else:
        areas_0.append(0)

print areas_0
print areas0


# ### Find different shape attributes of the empty channel 

# In[ ]:


get_ipython().magic(u'matplotlib inline')
matplotlib.rcParams.update({'figure.max_open_warning': 0})

perimeters_0 = []
circularity_0 = []
minc_centres_0 = []
minc_radii_0 = []
min_ellipse_0 = []
boxes_0 = []


# In[ ]:


for (i,c) in enumerate(contours_0rough): # input is the list of contours, depending on which was used above 
    
    if np.size(c)>1:
        cnnt = c[0]
        # perimeters
        perimeters = cv2.arcLength(cnnt,True)
        perimeters_0.append(perimeters)

        #circularity
        circularity = ((4*np.pi*areas_0[i])/perimeters**2)
        circularity_0.append(circularity)

        #minimum enclosing circle, radius and centre of the circle are used as parameters 
        (x,y),radius = cv2.minEnclosingCircle(cnnt)
        centre = (int(x),int(y))
        radius = int(radius)
        #cv2.circle(image_copies_0[i],centre,radius,(0,255,0),2)
        minc_centres_0.append(centre)
        minc_radii_0.append(radius)

        # fitting minimum rotated rectangle
        rect = cv2.minAreaRect(cnnt)
        box = cv2.boxPoints(rect)
        box = np.int0(box)
        boxes_0.append(box)
        #box_contour = cv2.drawContours(v,[box],0,(0,0,255),2)
    else:
        #min_ellipse_0.append(0)
        perimeters_0.append(0)
        circularity_0.append(0)
        minc_centres_0.append(0)
        minc_radii_0.append(0)
        boxes_0.append(0)
        
        


# #### Calculating moments and aspect ratio

# In[ ]:


moments_0 = []
hu_moments_0 = []

cx_0 = [] # the x co-ordinate of the centroid
cy_0 = [] # the y co-ordinate of the cetnroid
for i, v in enumerate(contours_0rough):
    cnt = v[0]
    M = cv2.moments(cnt)
    
    #Centroid
    cx = int(M['m10']/M['m00'])
    cy = int(M['m01']/M['m00'])
    
    #Calculate the Hu moments - these are scale, rotation and translation invariant so these will be used 
    hu = cv2.HuMoments(M).flatten()
    
    hu_moments_0.append(hu)
    moments_0.append(M)   
    cx_0.append(cx)
    cy_0.append(cy)


# #### Aspect ratios normalised for orientation of the gastruloids

# In[ ]:


aspect_ratios_0 = []

width = []
height = []

for i, v in enumerate(contours_0rough):
    cnt = v[0]    
    #Find the aspect ratio, which is the ratio of the width to height of bounding rectangle of the contour
    x,y,w,h = cv2.boundingRect(cnt)
    width.append(w)
    height.append(h)
    
for i, v in enumerate(width): # normalise the AR so that all values are > 1 so the actual distance can be measured
    k = height[i]
    #print k
    if v > k:
        aspect_ratios_0.append(np.float(v) / k)
    else:
        aspect_ratios_0.append(np.float(k) / v)
    

print aspect_ratios_0


# #### Separate the 7 Hu moments 

# In[ ]:


# the Hu moments are a way to look at the shape - if an image is rotated but is the same image, then they are 
#expected to be the same; it is used a lot for object DETECTION 

# These values are proved to be invariants to the image scale, rotation, 
#and reflection except the seventh one, whose sign is changed by reflection.
# https://docs.opencv.org/3.0-beta/modules/imgproc/doc/structural_analysis_and_shape_descriptors.html

I1 = []; I2 = []; I3 = []; I4 = []; I5 = []; I6 = []; I7 = []

for i in range(0,len(hu_moments_0)):
    moment = hu_moments_0[i]
    for j in range(0,2): 
        i1 = moment[0]
        i2 = moment[1]
        i3 = moment[2]
        i4 = moment[3]
        i5 = moment[4]
        i6 = moment[5]
        i7 = moment[6]
        #print moment
    
    I1.append(i1); I2.append(i2); I3.append(i3); I4.append(i4); I5.append(i5); I6.append(i6); I7.append(i7)


# #### Log-transform the Hu moments - this might not always be necessary but in some cases will show the features better

# In[ ]:


# first the array with all 7 together 
log_hu = []
for i, v in enumerate(hu_moments_0):
    log_hu_i = -np.sign(v) * np.log10(np.abs(v))
    log_hu.append(log_hu_i)


# In[ ]:


# Also for each separate moment 
log_i1 = []; log_i2=[]; log_i3=[]; log_i4=[]; log_i5=[]; log_i6=[]; log_i7=[]
I=[I1, I2, I3, I4, I5, I6, I7]

for j in range(0,7):
    for i,v in enumerate(I[j]):
        log_hu_i = -np.sign(v) * np.log10(np.abs(v))
        if j == 0:
            log_i1.append(log_hu_i)
        if j == 1:
            log_i2.append(log_hu_i)
        if j == 2:
            log_i3.append(log_hu_i)
        if j == 3: 
            log_i4.append(log_hu_i)
        if j == 4:
            log_i5.append(log_hu_i)
        if j == 5:
            log_i6.append(log_hu_i)
        elif j ==6 :
            log_i7.append(log_hu_i)

log_i = [log_i1, log_i2, log_i3, log_i4, log_i5, log_i6, log_i7]


# #### Make a mask to get the Haralick features from inside the Gastruloid only

# In[ ]:


blurred_masks_0 = []
for i, v in enumerate(image_copies_0):  
    img = v.copy()
    
    #Gaussian Blur - same as that used to de-noise when making the original Ch0 masks 
    blur = cv2.GaussianBlur(image_copy,(5,5),0)
    
    # Create a mask using the area of the gastruloid found in the empty channel 
    blurred_masked_img = cv2.bitwise_and(blur,blur, mask = erosions[i])
    blurred_masks_0.append(blurred_masked_img)


# In[ ]:


#https://gogul09.github.io/software/image-classification-python
### Get Haralick features - texture / more global than intensity values ###
from mahotas import features

haralick_features = []

for i, v in enumerate(blurred_masks_0): 
    #img = cv2.GaussianBlur(v, (3,3), 0)
    img = v.copy()
    haralick = features.haralick(img).mean(axis=0)
    haralick_features.append(haralick)


# In[ ]:


# split each into its own variable
H1ang2ndmoment = []; H2contrast = []; H3correlation = []; H4sumsqvar = []; H5invdiffmoment = []
H6sumavg = []; H7sumvar = []; H8sumentropy = []; H9entropy = []; H10diffvar = []; H11diffentropy = []
H12imcorr1 = []; H13imcorr2 = []


for i in range(0,len(haralick_features)):
    f = haralick_features[i]
    for j in range(0,2):
        H1 = f[0]
        H2 = f[1]
        H3 = f[2]
        H4 = f[3]
        H5 = f[4]
        H6 = f[5]
        H7 = f[6]
        H8 = f[7]
        H9 = f[8]
        H10 = f[9]
        H11 = f[10]
        H12 = f[11]
        H13 = f[12]
        
    H1ang2ndmoment.append(H1); H2contrast.append(H2); H3correlation.append(H3); H4sumsqvar.append(H4); 
    H5invdiffmoment.append(H5); H6sumavg.append(H6); H7sumvar.append(H7)
    H8sumentropy.append(H8); H9entropy.append(H9); H10diffvar.append(H10); H11diffentropy.append(H11); H12imcorr1.append(H12); 
    H13imcorr2.append(H13)


# #### Post-processing - convexity hull

# In[ ]:


hull_areas_0 = []
hull_0 = []
convexity_defects_0 = []

for (i,c) in enumerate(contours_0rough): 
    cnt = c[0]

    #find the hull for for the main contour 
    hull = cv2.convexHull(cnt,returnPoints = False) # this is used for the convexity defects, but is not kept
    hull2 = cv2.convexHull(cnt,returnPoints = True)
    hull_0.append(hull2)
    
    #Calculating the hull area 
    hullarea = cv2.contourArea(hull2)
    hull_areas_0.append(hullarea)   

    #convexity defects 
    defects = cv2.convexityDefects(cnt,hull)
    convexity_defects_0.append(defects)


# #### Calculating solidity / convexity - i.e., gld area / area of the convex hull

# In[ ]:


convexity_0 = []
for i, v in enumerate(areas_0):
    convexity = v / hull_areas_0[i]
    #print convexity
    convexity_0.append(convexity)


# #### From the ellipse: keep major and minor axis lengths as these don't change by rotation, but the angle does (so not saving that)

# In[ ]:


major_ellipse_axis = [] # length of the major axis of the ellipse 
minor_ellipse_axis = []
for i, c in enumerate(contours_0rough):
    if np.size(c) > 5:
        cnt = c[0]
        (x,y),(MA,ma),angle = cv2.fitEllipse(cnt)
        #print x,y, MA, ma, angle
        major_ellipse_axis.append(ma)
        minor_ellipse_axis.append(MA)
        
    else:
        major_ellipse_axis.append(0)
        minor_ellipse_axis.append(0)


# 1) Ensuring all lists are of equal length 

# In[ ]:


print "1) Area list:", len(areas_0)
print "2) Perimeter list:", len(perimeters_0)
print "3) Circularity list:", len(circularity_0)
print "4) Hull list:", len(hull_0)
print "5) Convexity defects list:", len(convexity_defects_0)
print "6) Hull area list:", len(hull_areas_0)
print "7) Minimum enclosing circle radii list:", len(minc_radii_0)
print "8) Minimum enclosing rotated rectangle list:", len(boxes_0)

print "\n9) Cy list:", len(cy_0)
print "10) Cx list:", len(cx_0)
print "11) Moments list:", len(moments_0)
print "12) Hu Moments list:", len(hu_moments_0)
print "13) Hu-I1:", len(I1),',', "14) Hu-I2:", len(I2),',', "15) Hu-I3:", len(I3),',', "16) Hu-I4:", len(I4),',', "17) Hu-I5:", len(I5),',',"18) Hu-I6:", len(I6),',', "19) Hu-I7:", len(I7)
print "20) Log Hu Moments list:", len(log_hu)
print "21) LogHu-I1:", len(log_i1),',', "22) LogHu-I2:", len(log_i2),',', "23) LogHu-I3:", len(log_i3),',', "24) LogHu-I4:", len(log_i4),',', "25) LogHu-I5:", len(log_i5),',',"26) LogHu-I6:", len(log_i6),',', "27) LogHu-I7:", len(log_i7)

print "\n28) H1ang2ndmoment:", len(H1ang2ndmoment), ',', "29) H2contrast:", len(H2contrast), ',', "30) H3correlation:" len(H3correlation),',', "31) H4sumsqvar:", len(H4sumsqvar), ',', "32) H5invdiffmoment:", len(H5invdiffmoment), ',', "33) H6sumavg:", len(H6sumavg), ',',"34) H7sumvar:", len(H7sumvar), ',', "35) H8sumentropy:", len(H8sumentropy), ',', "36) H9entropy:", len(H9entropy), ',',"37) H10diffvar:", len(H10diffvar), ',', "38) H11diffentropy:", len(H11diffentropy), ',', "39) H12imcorr1:", len(H12imcorr1), ',',"40) H13imcorr2:", len(H13imcorr2)


print "\n41) AR list:", len(aspect_ratios_0)
print "42) Convexity / Solidity list:", len(convexity_0)
print "43) Major Ellipse list:", len(major_ellipse_axis)
print "44) Minor Ellipse list:", len(minor_ellipse_axis)

if len(areas_0) == len(perimeters_0) == len(circularity_0) == len(hull_0) == len(convexity_defects_0) ==len(hull_areas_0) == len(minc_radii_0) == len(boxes_0) == len(moments_0) == len(cy_0) == len(cx_0) == len(hu_moments_0)== len(I1) ==len(I2) ==len(I3) ==len(I4) ==len(I5)==len(I6)==len(I7)==len(log_hu) ==len(log_i1)==len(log_i2) ==len(log_i3) ==len(log_i4) ==len(log_i5) ==len(log_i6) ==len(log_i7) == len(H1ang2ndmoment) == len(H2contrast) len(H3correlation) ==len(H4sumsqvar)==len(H5invdiffmoment) ==len(H6sumavg)==len(H7sumvar) ==len(H8sumentropy) ==len(H9entropy) ==len(H10diffvar)== len(H11diffentropy) ==len(H12imcorr1)==len(H13imcorr2) == len(aspect_ratios_0) == len(convexity_0) == len(major_ellipse_axis) ==len(minor_ellipse_axis) :
    print "\nCheck complete - All lists are of same length"
else: 
    print "\nError: all list lengths are not equal"


# ## Fluorescent Channels

# ##### Ch1 (here Sox2)

# 1) De-noise

# In[ ]:


get_ipython().magic(u'matplotlib inline')
matplotlib.rcParams.update({'figure.max_open_warning': 0})

blurred_1 = []
image_copies_1 = []
blurred_masks_1 = []

for i, v in enumerate(files_list[1][0]):  
    # readthe original image
    image_1 = np.array(Image.open(files_list[1][0][i]))
    #image_colour0 = np.array(Image.open(files_list_full[0][0][i]))
    
    #print image_1.shape
    
    #Make a copy bc all the transformations change the image, and we don't want to change the original
    image_copy = image_1.copy()
    img_copy = image_1.copy()
    image_copy1 = image_1.copy()
    image_copies_1.append(image_copy)
    
    #Gaussian Blur 
    blur = cv2.GaussianBlur(image_copy,(5,5),0) 
    blurred_1.append(blur)
    
    # Create a mask using the area of the gastruloid found in the empty channel 
    blurred_masked_img = cv2.bitwise_and(blur,blur, mask = erosions[i])
    blurred_masks_1.append(blurred_masked_img)

    #fig, ax = plt.subplots(1, 3, figsize = (10,3))
    #ax[0].imshow(image_1)
    #ax[1].imshow(blur) 
    #ax[2].imshow(blurred_masked_img)
    
    #ax[2].imshow(masked_img)
    #ax[3].imshow(th1)
    #ax[1].imshow(blurred_masked_img)


# In[ ]:


for i, v in enumerate(blurred_1):

    ## finding the range of the grayscale values so we can do the thresholding properly 
    # https://mmeysenburg.github.io/image-processing/05-creating-histograms/ 
    histogram = cv2.calcHist([v], [0], None, [256], [0, 256])
    plt.plot(histogram)

    plt.title("Grayscale Histogram")
    plt.xlabel("Grayscale Value")
    plt.ylabel("Pixels")
plt.show()


# 2) Thresholding

# In[ ]:


matplotlib.rcParams.update({'figure.max_open_warning': 0})

im_out_1 = []
contours_1 = []
masks_1 = []
th_in_use1 = []

for i, v in enumerate(blurred_masks_1):
    ret, th1 = cv2.threshold(v, 160, 255, cv2.THRESH_BINARY)
    mask_inv = cv2.bitwise_not(v)
    th_in_use1.append(th1)

    im_floodfill_1 = th1.copy()
    
    # Mask used for flood filling; Notice the size needs to be 2 pixels larger than the image.
    h, w, = th1.shape[:2]
    mask = np.zeros((h+2, w+2), np.uint8)

    # Floodfill from point (0, 0)
    cv2.floodFill(im_floodfill_1, mask, (0,0), 255)

    # Invert floodfilled image
    im_floodfill_inv_1 = cv2.bitwise_not(im_floodfill_1)

    # Combine the two images to get the foreground.
    im_out_1s = th1 | im_floodfill_inv_1

    im_out_1.append(im_out_1s)
    
    (_, contours, _) = cv2.findContours(th1, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    # If only the largest contour is required, add [:1] at the end of the expression; when expression is speckled, 
    # measure all contours to count the number of speckles 
    contours = sorted(contours, key = cv2.contourArea, reverse = True)[:1]
   
    cv2.drawContours(v, contours, -1, (255, 0, 0), 10)     # draw contours over original image
    if not contours:
        contours = np.array([0])
    contours_1.append(contours)
    
    #create mask 
    masked_img = cv2.bitwise_and(v,v, mask = im_out_1s)
    masks_1.append(masked_img)
    
    #plot it 
    fig, ax = plt.subplots(1, 3, figsize = (13,3))
    ax[0].imshow(th1)
    ax[0].set_xlabel("Distance / Pixels")
    ax[0].set_ylabel("Distance / Pixels")
    ax[0].set_title("Mask before floodilling")
    ax[1].imshow(im_out_1s)
    ax[1].set_xlabel("Distance / Pixels")
    ax[1].set_ylabel("Distance / Pixels")
    ax[1].set_title("Mask After Floodfilling")
    ax[2].imshow(v) 
    ax[2].set_xlabel("Distance / Pixels")
    ax[2].set_ylabel("Distance / Pixels")
    ax[2].set_title("Contour drawn around Fluorescence \nArea In Gastruloid Mask")
    
    plt.tight_layout()
    plt.show()
    
    #save_locn = "/Users/alexandrabaranowski/Desktop/Images for report /"
    #fig.savefig(save_locn+'1. 2015BraGFP 72h Ch1.png')


# 3) Feature Extraction

# In[ ]:


areas1 = [] # this will be used to see the total ch1 fluo
for i, v in enumerate(im_out_1):
    #find areas
    filledareas = cv2.countNonZero(v)
    areas1.append(filledareas)

areas_1=[] # this will be used to find the area of the biggest contour, so to calculate the perimeter, circ etc 
for i, c in enumerate(contours_1):
    if np.size(c)>1:
        cnt = c[0]
        area = cv2.contourArea(cnt)
        areas_1.append(area)
    else:
        areas_1.append(0)


# In[ ]:


matplotlib.rcParams.update({'figure.max_open_warning': 0})

perimeters_1 = []
circularity_1 = []
minc_centres_1 = []
minc_radii_1 = []
min_ellipse_1 = []
boxes_1 = []


# In[ ]:


for (i,c) in enumerate(contours_1): 
    
    if np.size(c)>3:
        cnt = c[0]

        # perimeters
        perimeters = cv2.arcLength(cnt,True)
        perimeters_1.append(perimeters)

        #circularity
        circularity = ((4*np.pi*areas_1[i])/perimeters**2)
        circularity_1.append(circularity)

        #minimum enclosing circle, radius and centre of the circle are used as parameters 
        (x,y),radius = cv2.minEnclosingCircle(cnt)
        centre = (int(x),int(y))
        radius = int(radius)
        #cv2.circle(image_copies_0[i],centre,radius,(0,255,0),2)
        minc_centres_1.append(centre)
        minc_radii_1.append(radius)

        # fitting minimum rotated rectangle
        rect = cv2.minAreaRect(cnt)
        box = cv2.boxPoints(rect)
        box = np.int0(box)
        boxes_1.append(box)
    else:
        perimeters_1.append(0)
        circularity_1.append(0)
        minc_centres_1.append(0)
        minc_radii_1.append(0)
        #min_ellipse_1.append(0)
        boxes_1.append(0)


# #### Calculating moments and aspect ratio

# In[ ]:


moments_1 = []
hu_moments_1 = []

cx_1 = [] # the x co-ordinate of the centroid
cy_1 = [] # the y co-ordinate of the cetnroid
for i, v in enumerate(contours_1):
    cnt = v[0]
    M = cv2.moments(cnt)
    
    #Calculate the Hu moments 
    hu = cv2.HuMoments(M).flatten() 
    hu_moments_1.append(hu)
    
    try:
        #Centroid
        cx = int(M['m10']/M['m00'])
        cy = int(M['m01']/M['m00'])
     
        moments_1.append(M)   
        cx_1.append(cx)
        cy_1.append(cy)
        
    except ZeroDivisionError:
        #hu_moments_1.append(0)
        moments_1.append(0)   
        cx_1.append(0)
        cy_1.append(0)


# #### Separate the 7 Hu moments 

# In[ ]:


I1_1 = []; I2_1 = []; I3_1 = []; I4_1 = []; I5_1 = []; I6_1 = []; I7_1 = []

for i in range(0,len(hu_moments_1)):
    moment = hu_moments_1[i]
    for j in range(0,2): 
        i1 = moment[0]
        i2 = moment[1]
        i3 = moment[2]
        i4 = moment[3]
        i5 = moment[4]
        i6 = moment[5]
        i7 = moment[6]
    
    I1_1.append(i1); I2_1.append(i2); I3_1.append(i3); I4_1.append(i4); I5_1.append(i5); I6_1.append(i6); I7_1.append(i7)


# #### Log-transform the Hu moments - this might not always be necessary but in some cases will show the features better

# In[ ]:


# first the array with all 7 together 
log_hu1 = []
np.seterr(divide='ignore', invalid='ignore')

for i, v in enumerate(hu_moments_1):
    try: 
        log_hu_i = -np.sign(v) * np.log10(np.abs(v))
        log_hu1.append(log_hu_i)
    except v==0:
        log_hu1.append(0)



# In[ ]:


# Also for each separate moment 
log_i11 = []; log_i21=[]; log_i31=[]; log_i41=[]; log_i51=[]; log_i61=[]; log_i71=[]
I_1=[I1_1, I2_1, I3_1, I4_1, I5_1, I6_1, I7_1]

for j in range(0,7):
    for i,v in enumerate(I_1[j]):
        log_hu_i = -np.sign(v) * np.log10(np.abs(v))
        if j == 0:
            log_i11.append(log_hu_i)
        if j == 1:
            log_i21.append(log_hu_i)
        if j == 2:
            log_i31.append(log_hu_i)
        if j == 3: 
            log_i41.append(log_hu_i)
        if j == 4:
            log_i51.append(log_hu_i)
        if j == 5:
            log_i61.append(log_hu_i)
        elif j ==6 :
            log_i71.append(log_hu_i)

log_ich1 = [log_i11, log_i21, log_i31, log_i41, log_i51, log_i61, log_i71]


# #### Hull & Convexity

# In[ ]:


hull_areas_1 = []
hull_1 = []
convexity_defects_1 = []

for (i,c) in enumerate(contours_1): 
    if np.size(c)>1:
        cnt = c[0]

    #find the hull for for the main contour 
        hull = cv2.convexHull(cnt,returnPoints = False) # this is used for the convexity defects, but is not kept
        hull2 = cv2.convexHull(cnt,returnPoints = True)
        hull_1.append(hull2)
    
        #Calculating the hull area 
        hullarea = cv2.contourArea(hull2)
        hull_areas_1.append(hullarea)   

        #convexity defects 
        defects = cv2.convexityDefects(cnt,hull)
        convexity_defects_1.append(defects)
    
    else:
        hull_1.append(0)
        hull_areas_1.append(0)
        convexity_defects_1.append(0)


# #### Convexity / Solidity 

# In[ ]:


convexity_1 = []
for i, v in enumerate(areas_1):
    try:
        convexity = v / hull_areas_1[i]
    #print convexity
        convexity_1.append(convexity)
    except ZeroDivisionError:
        convexity_1.append(0)
print len(convexity_1)


# #### Major and Minor ellipses of the fluo area - don't need this if expression isn't homogeneous

# In[ ]:


aspect_ratios_1 = []
width1 = []
height1 = []

for i, c in enumerate(contours_1):
    if np.size(c) >1:
        cnt = c[0]
        #Find the aspect ratio, which is the ratio of the width to height of bounding rectangle of the contour
        x,y,w,h = cv2.boundingRect(cnt)
        width1.append(w)
        height1.append(h)
    else:
        width1.append(0)
        height1.append(0)
    
for i, v in enumerate(width1):
    k = height1[i]
    #print k
    try: 
        if v > k:
            aspect_ratios_1.append(np.float(v) / k)
        else:
            aspect_ratios_1.append(np.float(k) / v)
    except ZeroDivisionError:
        aspect_ratios_1.append(0)
    

print len(aspect_ratios_1)


# In[ ]:


major_ellipse_axis1 = [] # length of the major axis of the ellipse 
minor_ellipse_axis1 = []

for i, c in enumerate(contours_1):
    if np.size(c) >5:
        cnt = c[0]
        (x,y),(MA,ma),angle = cv2.fitEllipse(cnt)
        #print x,y, MA, ma, angle
        major_ellipse_axis1.append(ma)
        minor_ellipse_axis1.append(MA)
        
    else:
        major_ellipse_axis1.append(0)
        minor_ellipse_axis1.append(0)


# ## QC for list length against each other and empty channel area

# In[ ]:


print "45) contour Area list:", len(areas_1)
print "46) full Area list:", len(areas1)

print "47) Perimeter list:", len(perimeters_1)
print "48) Circularity list:", len(circularity_1)
print "49) Hull list:", len(hull_1)
print "50) Convexity defects list:", len(convexity_defects_1)
print "51) Hull area list:", len(hull_areas_1)
print "52) Minimum enclosing circle radii list:", len(minc_radii_1)
print "53) Minimum enclosing rotated rectangle list:", len(boxes_1)

print "\n54) Cy list:", len(cy_1)
print "55) Cx list:", len(cx_1)
print "56) Moments list:", len(moments_1)
print "57) Hu Moments list:", len(hu_moments_1)
print "58) Hu-I1:", len(I1_1),',', "59) Hu-I2:", len(I2_1),',', "60) Hu-I3:", len(I3_1),',', "61) Hu-I4:", len(I4_1),',', "62) Hu-I5:", len(I5_1),',',"63) Hu-I6:", len(I6_1),',', "64) Hu-I7:", len(I7_1)
print "65) Log Hu Moments list:", len(log_hu1)
print "66) LogHu-I1:", len(log_i11),',', "67) LogHu-I2:", len(log_i21),',', "68) LogHu-I3:", len(log_i31),',', "69) LogHu-I4:", len(log_i41),',', "70) LogHu-I5:", len(log_i51),',',"71) LogHu-I6:", len(log_i61),',', "72) LogHu-I7:", len(log_i71)

print "\n73) AR list:", len(aspect_ratios_1)
print "74) Convexity / Solidity list:", len(convexity_1)
print "75) Major Ellipse list:", len(major_ellipse_axis1)
print "76) Minor Ellipse list:", len(minor_ellipse_axis1)

if len(areas_1) == len(areas1) == len(perimeters_1) == len(circularity_1) == len(hull_1) == len(convexity_defects_1) == len(hull_areas_1) == len(minc_radii_1) == len(boxes_1) == len(moments_1) == len(cy_1) == len(cx_1)== len(hu_moments_1) == len(I1_1) ==len(I2_1) ==len(I3_1) ==len(I4_1) ==len(I5_1)==len(I6_1)==len(I7_1)==len(log_hu1) ==len(log_i11)==len(log_i21) ==len(log_i31) ==len(log_i41) ==len(log_i51) ==len(log_i61) ==len(log_i71) == len(aspect_ratios_1) == len(convexity_1) ==len(major_ellipse_axis1) == len(minor_ellipse_axis1) == len(areas_0):
    print "\nCheck complete - All lists are of same length"
else: 
    print "\nError: all list lengths are not equal"


# ##### Ch2 (Sox17)

# 1) De-noise

# In[ ]:


#matplotlib.rcParams.update({'figure.max_open_warning': 0})

blurred_2 = []
image_copies_2 = []
blurred_masks_2 = []
blurred_masks_2noisy = []
for i, v in enumerate(files_list[2][0]):
    # readthe original image
    image_2 = np.array(Image.open(files_list[2][0][i]))
    #image_colour0 = np.array(Image.open(files_list_full[0][0][i]))
    
    #Make a copy bc all the transformations change the image, and we don't want to change the original
    image_copy = image_2.copy()
    image_copy2 = image_2.copy()
    image_copies_2.append(image_copy)
    
    #Gaussian Blur 
    blur = cv2.GaussianBlur(image_copy,(15,15),0) 
    blur2 = cv2.medianBlur(image_copy2, 5)
    blurred_2.append(blur2)
    
    blurred_masked_img = cv2.bitwise_and(blur,blur, mask = erosions[i])
    blurred_masks_2.append(blurred_masked_img)
    
    blurred_mask_noisy = cv2.bitwise_and(image_copy,image_copy, mask = threshold_0[i])
    blurred_masks_2noisy.append(blurred_mask_noisy)

    #fig, ax = plt.subplots(1, 3, figsize = (10,3))
    #ax[0].imshow(image_2)
    #ax[1].imshow(blur) 
    #ax[2].imshow(blur2)
    #ax[3].imshow(blurred_masked_img)


# In[ ]:


for i, v in enumerate(blurred_masks_2):

    ## finding the range of the grayscale values so we can do the thresholding properly 
    # https://mmeysenburg.github.io/image-processing/05-creating-histograms/ 
    histogram = cv2.calcHist([v], [0], None, [256], [0, 256])
    plt.plot(histogram)

    plt.title("Grayscale Histogram")
    plt.xlabel("Grayscale Value")
    plt.ylabel("Pixels")
plt.show()


# 2) Thresholding - make sure you fix the thresholds!! 

# In[ ]:


matplotlib.rcParams.update({'figure.max_open_warning': 0})

im_out_2 = []
contours_2 = []
masks_2 = []
th2_in_use = []

for i, v in enumerate(blurred_masks_2):
    ret, th2 = cv2.threshold(v, 150, 255, cv2.THRESH_BINARY)
    mask_inv = cv2.bitwise_not(v)
    th2_in_use.append(th2)
    
    #further removes noise 
    kernel = np.ones((5,5),np.uint8)
    opening = cv2.morphologyEx(th2,cv2.MORPH_OPEN,kernel, iterations = 2)

    im_floodfill_2 = th2.copy()
    
    # Mask used for flood filling; Notice the size needs to be 2 pixels larger than the image.
    h, w, = th2.shape[:2]
    mask = np.zeros((h+2, w+2), np.uint8)

    # Floodfill from point (0, 0)
    cv2.floodFill(im_floodfill_2, mask, (0,0), 255)

    # Invert floodfilled image
    im_floodfill_inv_2 = cv2.bitwise_not(im_floodfill_2)

    # Combine the two images to get the foreground.
    im_out_2s = th2 | im_floodfill_inv_2

    im_out_2.append(im_out_2s)
    
    # try thresholding using the median 
    median_blur = cv2.medianBlur(image_copies_2[i], 9)
    denoised = v - median_blur 
    ret, th = cv2.threshold(median_blur, 245, 254, cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    
    image = th + th2
    
    #create mask 
    masked_img = cv2.bitwise_and(v,v, mask = im_out_1s)
    masks_2.append(masked_img)
    
    #stat = ImageStat.Stat(v, mask=masked_img)
    #print stat.mean, stat.var, stat.median
    
    #plot it 
    fig, ax = plt.subplots(1, 3, figsize = (10,3))
    ax[0].imshow(th2) # probably this?? 
    ax[1].imshow(im_out_2s)
    ax[2].imshow(v)
    #ax[3].imshow(th)
    #ax[3].imshow(median_blur)
#ax[5].imshow(image)


# 3) Feature Extraction

# #### Total area of gastruloid covered by Ch2 fluorescence

# In[ ]:


areas2 = []
for i, v in enumerate(th2_in_use):
    #find areas
    filledareas = cv2.countNonZero(v)
    areas2.append(filledareas)

print areas2



# 4) Measuring fluorescence intensity of Ch2 

# ### QC

# In[ ]:


print "77) Area list:", len(areas2)

if len(areas2) == len(areas_0):
    print "Check complete: all ok"
else: 
    print "List length not correct - check "


# ##### Ch3 (T/Bra)

# 1) De-noise

# In[ ]:


blurred_3 = []
image_copies_3 = []
blurred_masks_3 = []

for i, v in enumerate(files_list[3][0]):
    # readthe original image
    image_3 = np.array(Image.open(files_list[3][0][i]))
    #image_colour0 = np.array(Image.open(files_list_full[0][0][i]))
    
    #Make a copy bc all the transformations change the image, and we don't want to change the original
    image_copy = image_3.copy()
    image_copy2 = image_3.copy()
    image_copies_3.append(image_copy)
    
    #Gaussian Blur 
    blur = cv2.GaussianBlur(image_copy,(15,15),0)
    blur2 = cv2.medianBlur(image_copy2, 5)
    blurred_3.append(blur)
    
    blurred_masked_img = cv2.bitwise_and(blur,blur, mask = erosions[i])
    blurred_masks_3.append(blurred_masked_img)

    #fig, ax = plt.subplots(1, 3, figsize = (10,3))
    #ax[0].imshow(image_3)
    #ax[1].imshow(blur) 
    #ax[2].imshow(blur2)
    #ax[2].imshow(blurred_masked_img)


# In[ ]:


for i, v in enumerate(blurred_masks_3):

    ## finding the range of the grayscale values so we can do the thresholding properly 
    # https://mmeysenburg.github.io/image-processing/05-creating-histograms/ 
    histogram = cv2.calcHist([v], [0], None, [256], [0, 256])
    plt.plot(histogram)

    plt.title("Grayscale Histogram")
    plt.xlabel("Grayscale Value")
    plt.ylabel("Pixels")
plt.show()


# 2) Thresholding - Make sure you fix the thresholds and have a reason for all of them 

# In[ ]:


im_out_3 = []
contours_3 = []
masks_3 = []
th_in_use_3 = []


# In[ ]:


matplotlib.rcParams.update({'figure.max_open_warning': 0})


for i, v in enumerate(blurred_masks_3):
    ret, thch3 = cv2.threshold(v, 150, 255, cv2.THRESH_BINARY)
    mask_inv = cv2.bitwise_not(v)
    th_in_use_3.append(thch3)
    
    #further removes noise - use this if needed 
    #kernel = np.ones((5,5),np.uint8)
    #opening = cv2.morphologyEx(thch3,cv2.MORPH_OPEN,kernel, iterations = 2)

    im_floodfill_3 = thch3.copy()
    
    # Mask used for flood filling; Notice the size needs to be 2 pixels larger than the image.
    h, w, = thch3.shape[:2]
    mask = np.zeros((h+2, w+2), np.uint8)

    # Floodfill from point (0, 0)
    cv2.floodFill(im_floodfill_3, mask, (0,0), 255)

    # Invert floodfilled image
    im_floodfill_inv_3 = cv2.bitwise_not(im_floodfill_3)

    # Combine the two images to get the foreground.
    im_out_3s = thch3 | im_floodfill_inv_3

    im_out_3.append(im_out_3s)
    
    (_, contours, _) = cv2.findContours(thch3, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    contours = sorted(contours, key = cv2.contourArea, reverse = True)[:1]
     
    cv2.drawContours(v, contours, -1, (255, 0, 0), 10)     # draw contours over original image
    if not contours:
        contours = np.array([0])
    contours_3.append(contours)
        
    #create mask 
    masked_img = cv2.bitwise_and(v,v, mask = im_out_1s)
    masks_3.append(masked_img)
    
    #stat = ImageStat.Stat(v, mask=masked_img)
    #print stat.mean, stat.var, stat.median
    
    #plot it 
    fig, ax = plt.subplots(1, 3, figsize = (10,3))
    ax[0].imshow(thch3)
    ax[1].imshow(im_out_3s) # this is the best one here 
    ax[2].imshow(v)


# 3) Feature Extraction

# #### Areas - Proportion of gld covered by Ch3 fluo

# In[ ]:


areas3 = []
for i, v in enumerate(th_in_use_3):
    #find areas
    filledareas = cv2.countNonZero(v)
    areas3.append(filledareas)

areas_3= []
for i, c in enumerate(contours_3):
    if np.size(c)>1:
        cnt = c[0]
        area = cv2.contourArea(cnt)
        areas_3.append(area)
    else:
        areas_3.append(0)


# In[ ]:


matplotlib.rcParams.update({'figure.max_open_warning': 0})

perimeters_3 = []
circularity_3 = []
minc_centres_3 = []
minc_radii_3 = []
min_ellipse_3 = []
boxes_3 = []


# In[ ]:


for (i,c) in enumerate(contours_3): 
    if np.size(c)>1:
        cnt = c[0]

        # perimeters
        perimeters = cv2.arcLength(cnt,True)
        perimeters_3.append(perimeters)

        #circularity
        circularity = ((4*np.pi*areas_3[i])/perimeters**2)
        circularity_3.append(circularity)

        #minimum enclosing circle, radius and centre of the circle are used as parameters 
        (x,y),radius = cv2.minEnclosingCircle(cnt)
        centre = (int(x),int(y))
        radius = int(radius)
        #cv2.circle(image_copies_0[i],centre,radius,(0,255,0),2)
        minc_centres_3.append(centre)
        minc_radii_3.append(radius)

        # fitting minimum rotated rectangle
        rect = cv2.minAreaRect(cnt)
        box = cv2.boxPoints(rect)
        box = np.int0(box)
        boxes_3.append(box)
    
    else: 
        perimeters_3.append(0)
        circularity_3.append(0)
        minc_centres_3.append(0)
        minc_radii_3.append(0)
        boxes_3.append(0)


# #### Calculating moments

# In[ ]:


moments_3 = []
hu_moments_3 = []

cx_3 = [] # the x co-ordinate of the centroid
cy_3 = [] # the y co-ordinate of the cetnroid
for i, v in enumerate(contours_3):
    cnt = v[0]
    M = cv2.moments(cnt)
    
    #Calculate the Hu moments 
    hu = cv2.HuMoments(M).flatten() 
    hu_moments_3.append(hu)
    
    try:
        #Centroid
        cx = int(M['m10']/M['m00'])
        cy = int(M['m01']/M['m00'])

        moments_3.append(M)   
        cx_3.append(cx)
        cy_3.append(cy)
        
    except ZeroDivisionError: 
        cy_3.append(0)
        cx_3.append(0)
        moments_3.append(0)


# #### Separate the 7 Hu moments 

# In[ ]:


I1_3 = []; I2_3 = []; I3_3 = []; I4_3 = []; I5_3 = []; I6_3 = []; I7_3 = []

for i in range(0,len(hu_moments_3)):
    moment = hu_moments_3[i]
    for j in range(0,2): 
        i1 = moment[0]
        i2 = moment[1]
        i3 = moment[2]
        i4 = moment[3]
        i5 = moment[4]
        i6 = moment[5]
        i7 = moment[6]
        #print moment
    
    I1_3.append(i1); I2_3.append(i2); I3_3.append(i3); I4_3.append(i4); I5_3.append(i5); I6_3.append(i6); I7_3.append(i7)


# #### Log-transform the Hu moments - this might not always be necessary but in some cases will show the features better

# In[ ]:


# first the array with all 7 together 
log_hu3 = []
for i, v in enumerate(hu_moments_3):
    log_hu_i = -np.sign(v) * np.log10(np.abs(v))
    log_hu3.append(log_hu_i)


# In[ ]:


# Also for each separate moment 
log_i13 = []; log_i23=[]; log_i33=[]; log_i43=[]; log_i53=[]; log_i63=[]; log_i73=[]
I_3=[I1_3, I2_3, I3_3, I4_3, I5_3, I6_3, I7_3]

for j in range(0,7):
    for i,v in enumerate(I_3[j]):
        log_hu_i = -np.sign(v) * np.log10(np.abs(v))
        if j == 0:
            log_i13.append(log_hu_i)
        if j == 1:
            log_i23.append(log_hu_i)
        if j == 2:
            log_i33.append(log_hu_i)
        if j == 3: 
            log_i43.append(log_hu_i)
        if j == 4:
            log_i53.append(log_hu_i)
        if j == 5:
            log_i63.append(log_hu_i)
        elif j ==6 :
            log_i73.append(log_hu_i)

log_ich3 = [log_i13, log_i23, log_i33, log_i43, log_i53, log_i63, log_i73]


# #### ARs; Post-processing - hulls & hull area

# In[ ]:


hull_areas_3 = []
hull_3 = []
convexity_defects_3 = []
aspect_ratios_3 = []
width3 = []
height3 = []

for (i,c) in enumerate(contours_3): 
    if np.size(c)>1: 
        cnt = c[0]

        #find the hull for for the main contour 
        hull = cv2.convexHull(cnt,returnPoints = False) # this is used for the convexity defects, but is not kept
        hull2 = cv2.convexHull(cnt,returnPoints = True)
        hull_3.append(hull2)

        #Calculating the hull area 
        hullarea = cv2.contourArea(hull2)
        hull_areas_3.append(hullarea)   

        #convexity defects 
        defects = cv2.convexityDefects(cnt,hull)
        convexity_defects_3.append(defects)

        #AR
        x,y,w,h = cv2.boundingRect(cnt)
        width3.append(w)
        height3.append(h)
        
    else: 
        hull_areas_3.append(0)
        aspect_ratios_3.append(0)
        hull_3.append(0)
        convexity_defects_3.append(0)
        
for i, v in enumerate(width3):
    k = height3[i]
    #print k
    if v > k:
        aspect_ratios_3.append(np.float(v) / k)
    else:
        aspect_ratios_3.append(np.float(k) / v)


# #### Convexity / solidity - change the area used for the other ones too

# In[ ]:


convexity_3 = []
for i, v in enumerate(areas_3):
    try: 
        convexity = v / hull_areas_3[i]
        #print convexity
        convexity_3.append(convexity)
    except ZeroDivisionError:
        convexity_3.append(0)


# #### Major & minor ellipse axes - check if this is correct or if they are the other way round

# In[ ]:


major_ellipse_axis3 = [] # length of the major axis of the ellipse 
minor_ellipse_axis3 = []

for i, c in enumerate(contours_3):
    if np.size(c)>1:
        cnt = c[0]
        (x,y),(MA,ma),angle = cv2.fitEllipse(cnt)
        #print x,y, MA, ma, angle
        major_ellipse_axis3.append(ma)
        minor_ellipse_axis3.append(MA)
    else: 
        major_ellipse_axis3.append(0)
        minor_ellipse_axis3.append(0)


# In[ ]:


print "78) contour Area list:", len(areas_3)
print "79) full Area list:", len(areas3)
print "80) Perimeter list:", len(perimeters_3)
print "81) Circularity list:", len(circularity_3)
print "82) Hull list:", len(hull_3)
print "83) Convexity defects list:", len(convexity_defects_3)
print "84) Hull area list:", len(hull_areas_3)
print "85) Minimum enclosing circle radii list:", len(minc_radii_3)
print "86) Minimum enclosing rotated rectangle list:", len(boxes_3)

print "\n87) Cy list:", len(cy_3)
print "88) Cx list:", len(cx_3)
print "89) Moments list:", len(moments_3)
print "90) Hu Moments list:", len(hu_moments_3)
print "91) Hu-I1:", len(I1_3),',', "92) Hu-I2:", len(I2_3),',', "93) Hu-I3:", len(I3_3),',', "94) Hu-I4:", len(I4_3),',', "95) Hu-I5:", len(I5_3),',',"96) Hu-I6:", len(I6_3),',', "97) Hu-I7:", len(I7_3)
print "98) Log Hu Moments list:", len(log_hu3)
print "99) LogHu-I1:", len(log_i13),',', "100) LogHu-I2:", len(log_i23),',', "101) LogHu-I3:", len(log_i33),',', "102) LogHu-I4:", len(log_i43),',', "103) LogHu-I5:", len(log_i53),',',"104) LogHu-I6:", len(log_i63),',', "105) LogHu-I7:", len(log_i73)

print "\n106) AR list:", len(aspect_ratios_3)
print "107) Convexity / Solidity list:", len(convexity_3)
print "108) Major Ellipse list:", len(major_ellipse_axis3)
print "109) Minor Ellipse list:", len(minor_ellipse_axis3)

if len(areas_3) == len(perimeters_3) == len(circularity_3) == len(hull_3) == len(convexity_defects_3) ==len(hull_areas_3) == len(minc_radii_3) == len(boxes_3) == len(moments_3) == len(cy_3) == len(cx_3) ==len(hu_moments_3) == len(I1_3) ==len(I2_3) ==len(I3_3) ==len(I4_3) ==len(I5_3)==len(I6_3)==len(I7_3)==len(log_hu3) ==len(log_i13)==len(log_i23) ==len(log_i33) ==len(log_i43) ==len(log_i53) ==len(log_i63) ==len(log_i73) == len(aspect_ratios_3) == len(convexity_3) ==len(major_ellipse_axis3) == len(minor_ellipse_axis3) == len(areas_0):
    print "\nCheck complete - All lists are of same length"
else: 
    print "\nError: all list lengths are not equal"


# ## Measuring Fluorescence

# #### Performing Illumination and Background Correction

# In[ ]:


# retrieve all the bg images 

bg_1 = []
bg_2 = []
bg_3 = []
for i, v in enumerate(bg_list[1][0]):
    bg1 = np.array(Image.open(bg_list[1][0][i]))
    bg_1.append(bg1)
for i, v in enumerate(bg_list[2][0]):
    bg2 = np.array(Image.open(bg_list[2][0][i]))
    bg_2.append(bg2)

for i, v in enumerate(bg_list[3][0]):    
    bg3 = np.array(Image.open(bg_list[3][0][i]))
    bg_3.append(bg3)  


# In[ ]:


## Correcting for uneven illumination - channel 1

get_ipython().magic(u'matplotlib inline')
fgs_roi1 = []
fgs_subtract1 = []
fgs_norm1 = []

for i, v in enumerate(files_list_full[1][0]):
    image = np.array(Image.open(files_list_full[1][0][i]))
    bg1 = bg_1[0]
    
    img = image.copy()  
    img = cv2.GaussianBlur(img, (5,5),0)
    img = Image.fromarray(img)
    
    #extract pixel values for the 16-bit image 
    WIDTH, HEIGHT = img.size
    pix = img.load()
    data = np.asarray(img.getdata()) # getdata returns the contents of the image as a sequence object containing #pixel values 
    data = data.reshape((HEIGHT,WIDTH)) #reshape the data to the actual dimensions
    
    # extract pixel values for the full image bg 
    bgimg = bg1.copy()
    bgimg = cv2.GaussianBlur(bgimg, (5,5), 0)
    bgimg = Image.fromarray(bgimg)
    WIDTHb, HEIGHTb = bgimg.size
    pixb = bgimg.load()

    datab = np.asarray(bgimg.getdata()) 
    datab = datab.reshape((HEIGHTb,WIDTHb)) 
    
    # extract an ROI to calculate and subtract the background 
    roi = np.array(img.copy())
    roi = roi[0:64, 0:1024]
    roi = Image.fromarray(roi)

    WIDTHr, HEIGHTr = roi.size
    pixr = roi.load()
    datar = np.asarray(roi.getdata())
    datar = datar.reshape((HEIGHTr,WIDTHr))
    
    
    
     
    #returns the column-wise average of the array elements along the axis for each of the types 
    reduced_data = data.mean(axis=0) 
    
    reduced_datar = datar.mean(axis=0) 
    
    reduced_datab = datab.mean(axis=0) 

    fg = reduced_data - reduced_datar
    fgs_roi1.append(fg)
    
    fg6 = reduced_data - reduced_datab
    fgs_subtract1.append(fg6) 
    
    bg = np.matrix(bgimg)
    fgg = np.matrix(img)

    meanbg = np.mean(bgimg)
    meanfg = np.mean(img)
    nbg = bg * meanfg
    nfg = fgg * meanbg 
    result = nfg - nbg
    result = result - np.min(result)
    result = result / np.max(result)
    result2 = result.copy()
    #result2 = cv2.GaussianBlur(result2, (15,15), 0)
    
    result2 = Image.fromarray(result2)
    WIDTHre, HEIGHTre = result2.size
    pixre = result2.load()
    datare = np.asarray(result2.getdata()) 
    datare = datare.reshape((HEIGHTre,WIDTHre))
    reduced_datare = datare.mean(axis=0) 
    fgs_norm1.append(reduced_datare)
    
    
    fig,ax = plt.subplots(1, 2, figsize = (13,3), sharex=False, sharey = False)

    ax[0].plot(fg6)
    ax[0].set_xlabel("Distance / Pixels")
    ax[0].set_ylabel("Grayscale Value")
    ax[0].set_title("Original Image")
    ax[1].plot(fg)
    ax[1].set_xlabel("Distance / Pixels")
    ax[1].set_ylabel("Grayscale Value")
    ax[1].set_title("Normalised from ROI")
    #ax[2].plot(fg6)
    #ax[2].set_xlabel("Distance / Pixels")
    #ax[2].set_ylabel("Grayscale Value")
    #ax[2].set_title("Normalised from Empty Image")
    #ax[3].plot(reduced_datare)
    #ax[3].set_xlabel("Distance / Pixels")
    #ax[3].set_ylabel("Grayscale Value")
    #ax[3].set_title("Normalised from Empty Image")
    
    #plt.tight_layout()
   # plt.show()
    
    #save_locn = "/Users/alexandrabaranowski/Desktop/Images for report /"
    #fig.savefig(save_locn+'Ch1 Background subtraction ROI vs Empty Image.png')


# #### Check to see which method to use; if fine, the preferable method is roi_subtract from the empty image

# In[ ]:


### Checking whether the two methods are similar enough to be able to use the empty channel 
#matplotlib.style.use('ggplot')

# Pearson correlation for the two methods
Pearson = stats.pearsonr(fgs_subtract1[0], fgs_norm1[0])
print 'Correlation Coefficient:', Pearson[0]
print 'p-value:', Pearson[1]


plt.scatter(fgs_subtract1[0], fgs_roi1[0], color='c')
plt.xlabel('De-noised Fluorescence from Empty Image / Pixel Intensity')
plt.ylabel('De-noised & Normalised Fluorescence from \nEmpty Image / Pixel Intensity');
plt.title('Correlation of Background Subtraction Methods')
plt.tight_layout()


#Finding correlation
dfcorr = pd.DataFrame({'Empty': fgs_subtract1[0]})
dfcorr['Normalised'] = dfcorr['Empty'] + fgs_roi1[0] # positively correlated with 'a'

dfcorr.corr()


plt.show()
save_locn = "/Users/alexandrabaranowski/Desktop/Images for report /"
#fig.savefig(save_locn+'2016-03-15 Correlation of Background subtraction Methods.png')


# In[ ]:


x = fgs_subtract1[0]
y = fgs_norm1[0]

def r2(x1, y1):
    return np.float(stats.pearsonr(x, y)[0] ** 2)
sns.jointplot(x, y, kind="reg", stat_func = r2)


# In[ ]:


## Correcting for uneven illumination - channel 2

get_ipython().magic(u'matplotlib inline')
fgs_roi2 = []
fgs_subtract2 = []
fgs_norm2 = []

for i, v in enumerate(files_list_full[2][0]):
    image = np.array(Image.open(files_list_full[2][0][i]))
    bg2 = bg_2[0]
    img = image.copy()  
    img = cv2.GaussianBlur(img, (15,15),0)
    img = Image.fromarray(img)
    
    #extract pixel values for the 16-bit image 
    WIDTH, HEIGHT = img.size
    pix = img.load()
    data = np.asarray(img.getdata()) # getdata returns the contents of the image as a sequence object containing #pixel values 
    data = data.reshape((HEIGHT,WIDTH)) #reshape the data to the actual dimensions
    
    # extract pixel values for the full image bg 
    bgimg = bg2.copy()
    bgimg = cv2.GaussianBlur(bgimg, (15,15), 0)
    bgimg = Image.fromarray(bgimg)
    WIDTHb, HEIGHTb = bgimg.size
    pixb = bgimg.load()

    datab = np.asarray(bgimg.getdata()) 
    datab = datab.reshape((HEIGHTb,WIDTHb)) 
    
    # extract an ROI to calculate and subtract the background 
    roi = np.array(img.copy())
    roi = roi[0:64, 0:1024]
    roi = Image.fromarray(roi)

    WIDTHr, HEIGHTr = roi.size
    pixr = roi.load()
    datar = np.asarray(roi.getdata())
    datar = datar.reshape((HEIGHTr,WIDTHr))
    
    
    
     
    #returns the column-wise average of the array elements along the axis for each of the types 
    reduced_data = data.mean(axis=0) 
    
    reduced_datar = datar.mean(axis=0) 
    
    reduced_datab = datab.mean(axis=0) 

    fg = reduced_data - reduced_datar
    fgs_roi2.append(fg)
    
    fg6 = reduced_data - reduced_datab
    fgs_subtract2.append(fg6) 
    
    bg = np.matrix(bgimg)
    fgg = np.matrix(img)

    meanbg = np.mean(bgimg)
    meanfg = np.mean(img)
    nbg = bg * meanfg
    nfg = fgg * meanbg 
    result = nfg - nbg
    result = result - np.min(result)
    result = result / np.max(result)
    result2 = result.copy()
    #result2 = cv2.GaussianBlur(result2, (15,15), 0)
    
    result2 = Image.fromarray(result2)
    WIDTHre, HEIGHTre = result2.size
    pixre = result2.load()
    datare = np.asarray(result2.getdata()) 
    datare = datare.reshape((HEIGHTre,WIDTHre))
    reduced_datare = datare.mean(axis=0) 
    fgs_norm2.append(reduced_datare)
    
    
    fig,ax = plt.subplots(1, 4, figsize = (13,3), sharex=False, sharey = False)

    ax[0].plot(reduced_data)
    ax[0].set_xlabel("Distance / Pixels")
    ax[0].set_ylabel("Grayscale Value")
    ax[0].set_title("Original Image")
    ax[1].plot(fg)
    ax[1].set_xlabel("Distance / Pixels")
    ax[1].set_ylabel("Grayscale Value")
    ax[1].set_title("Normalised from ROI")
    ax[2].plot(fg6)
    ax[2].set_xlabel("Distance / Pixels")
    ax[2].set_ylabel("Grayscale Value")
    ax[2].set_title("Normalised from Empty Image")
    ax[3].plot(reduced_datare)
    ax[3].set_xlabel("Distance / Pixels")
    ax[3].set_ylabel("Grayscale Value")
    ax[3].set_title("Normalised from Empty Image")
    
    plt.tight_layout()
    plt.show()
    
    #save_locn = "/Users/alexandrabaranowski/Desktop/Images for report /"
    #fig.savefig(save_locn+'Ch1 Background subtraction ROI vs Empty Image.png')


# In[ ]:


## Correcting for uneven illumination - channel 3

get_ipython().magic(u'matplotlib inline')
fgs_roi3 = []
fgs_subtract3 = []
fgs_norm3 = []

for i, v in enumerate(files_list_full[3][0]):
    image = np.array(Image.open(files_list_full[3][0][i]))
    bg3 = bg_3[0]
    img = image.copy()  
    img = cv2.GaussianBlur(img, (15,15),0)
    img = Image.fromarray(img)
    
    #extract pixel values for the 16-bit image 
    WIDTH, HEIGHT = img.size
    pix = img.load()
    data = np.asarray(img.getdata()) # getdata returns the contents of the image as a sequence object containing #pixel values 
    data = data.reshape((HEIGHT,WIDTH)) #reshape the data to the actual dimensions
    
    # extract pixel values for the full image bg 
    bgimg = bg3.copy()
    bgimg = cv2.GaussianBlur(bgimg, (15,15), 0)
    bgimg = Image.fromarray(bgimg)
    WIDTHb, HEIGHTb = bgimg.size
    pixb = bgimg.load()

    datab = np.asarray(bgimg.getdata()) 
    datab = datab.reshape((HEIGHTb,WIDTHb)) 
    
    # extract an ROI to calculate and subtract the background 
    roi = np.array(img.copy())
    roi = roi[0:64, 0:1024]
    roi = Image.fromarray(roi)

    WIDTHr, HEIGHTr = roi.size
    pixr = roi.load()
    datar = np.asarray(roi.getdata())
    datar = datar.reshape((HEIGHTr,WIDTHr))
    
    
    
     
    #returns the column-wise average of the array elements along the axis for each of the types 
    reduced_data = data.mean(axis=0) 
    
    reduced_datar = datar.mean(axis=0) 
    
    reduced_datab = datab.mean(axis=0) 

    fg = reduced_data - reduced_datar
    fgs_roi3.append(fg)
    
    fg6 = reduced_data - reduced_datab
    fgs_subtract3.append(fg6) 
    
    bg = np.matrix(bgimg)
    fgg = np.matrix(img)

    meanbg = np.mean(bgimg)
    meanfg = np.mean(img)
    nbg = bg * meanfg
    nfg = fgg * meanbg 
    result = nfg - nbg
    result = result - np.min(result)
    result = result / np.max(result)
    result2 = result.copy()
    #result2 = cv2.GaussianBlur(result2, (15,15), 0)
    
    result2 = Image.fromarray(result2)
    WIDTHre, HEIGHTre = result2.size
    pixre = result2.load()
    datare = np.asarray(result2.getdata()) 
    datare = datare.reshape((HEIGHTre,WIDTHre))
    reduced_datare = datare.mean(axis=0) 
    fgs_norm3.append(reduced_datare)
    
    
    fig,ax = plt.subplots(1, 2, figsize = (13,3), sharex=False, sharey = False)

    ax[0].plot(fg6)
    ax[0].set_xlabel("Distance / Pixels")
    ax[0].set_ylabel("Grayscale Value")
    ax[0].set_title("Original Image")
    ax[1].plot(fg)
    ax[1].set_xlabel("Distance / Pixels")
    ax[1].set_ylabel("Grayscale Value")
    ax[1].set_title("Normalised from ROI")
   # ax[2].plot(fg6)
   # ax[2].set_xlabel("Distance / Pixels")
   # ax[2].set_ylabel("Grayscale Value")
   # ax[2].set_title("Normalised from Empty Image")
   # ax[3].plot(reduced_datare)
   # ax[3].set_xlabel("Distance / Pixels")
   # ax[3].set_ylabel("Grayscale Value")
   # ax[3].set_title("Normalised from Empty Image")
    
    #plt.tight_layout()
   #plt.show()
    
    #save_locn = "/Users/alexandrabaranowski/Desktop/Images for report /"
    #fig.savefig(save_locn+'Ch1 Background subtraction ROI vs Empty Image.png')


# #### Extracting the fluorescence values and normalising the minimum to zero

# In[ ]:


maxfluo1 = []; maxfluo2 = []; maxfluo3 = [] # this is a list of the max fluorescence value for EACH gastruloid 
meanfluo1 = []; meanfluo2 = []; meanfluo3 = [] # This is a list of the mean fluorescence value for EACH gastruloid
stdevfluo1 = []; stdevfluo2 = []; stdevfluo3 = []
minfluo1 = []; minfluo2 = []; minfluo3 = []

mins1 = []
for i, v in enumerate(fgs_subtract1):
    mins1.append(min(v))

# Make sure that the lowest value for fluorescence is zero; bc the same image is used for background subtraction, it 
# might lead to small deviations / values slightly below zero - this shifts up all values by the minimum 
a = min(mins1)
for i, v in enumerate(fgs_subtract1):
    if a < 0: 
        v = v - a
    else: 
        v = v
    maxfluo1.append(v.max())
    stdevfluo1.append(np.std(v))
    meanfluo1.append(v.mean())
    minfluo1.append(v.min())
    
avg_fluo_1 = np.mean(maxfluo1)
stddev_fluo_1 = np.std(maxfluo1, axis=0)# use the one without the outliers to calculate the avg max fluo

mins2 = []
for i, v in enumerate(fgs_subtract2):
    mins2.append(min(v))

b = min(mins2)
for i, v in enumerate(fgs_subtract2):
    if b < 0: 
        v = v - b
    else: 
        v = v
    maxfluo2.append(v.max())
    stdevfluo2.append(np.std(v))
    meanfluo2.append(v.mean())
    minfluo2.append(v.min())

avg_fluo_2 = np.mean(maxfluo2)
stddev_fluo_2 = np.std(maxfluo2, axis=0)

# Ch3 
mins3 = []
for i, v in enumerate(fgs_subtract3):
    mins3.append(min(v))

c = min(mins3)
for i, v in enumerate(fgs_subtract3):
    if c < 0: 
        v = v - c
    else: 
        v = v
    maxfluo3.append(v.max())
    stdevfluo3.append(np.std(v))
    meanfluo3.append(v.mean())
    minfluo3.append(v.min())

avg_fluo_3 = np.mean(maxfluo3)
stddev_fluo_3 = np.std(maxfluo3, axis=0)

print 'Average Max value across all Ch1 gastruloids is:', avg_fluo_1
print 'Average Max value across all Ch2 gastruloids is:', avg_fluo_2
print 'Average Max value across all Ch3 gastruloids is:', avg_fluo_3
print '\nSt Dev of max fluo across Ch1 gastruloids is', stddev_fluo_1
print 'St Dev of max fluo across Ch2 gastruloids is', stddev_fluo_2
print 'St Dev of max fluo across Ch3 gastruloids is', stddev_fluo_3



# ## QC based on area and circularity of the empty channel 

# #### Field-of-view quality control - removing images with no Gastruloids, debris etc.

# Method: Using area and circularity, but many different methods have been tested below, and can be added or removed as required to get the best option; manual verification against the ch0 contours is required

# In[ ]:


new_area = []

print 'Number of observations before removing outliers is %d' %len(areas_0)
#print areas_0
areas = np.array(areas_0)

mean = np.mean(areas, axis=0)
sd = np.std(areas, axis=0)
bins = 50 
new_area = [x for x in areas_0 if (x > mean - 1.5 * sd)]
new_area = [x for x in new_area if (x < mean + 1.5 * sd)]

#print 'Min area is %d' %min(new_area)
#print 'Max area is %d' %max(new_area)
#print 'Avg area is %d'%(sum(new_area)/len(new_area))

circ = np.array(circularity_0)

mean = np.mean(circ, axis=0)
sd = np.std(circ, axis=0)
bins = 50 
new_circ = [x for x in circularity_0 if (x > mean -1.5 * sd)]
new_circ = [x for x in new_circ if (x < mean + 2 * sd)]


# print new_area
print 'Number of observations within the threshold is %d' %len(new_area)
print 'Number of observations within the threshold is %d' %len(new_circ)

# want to identify which position all the outliers are in so we can remove from all the lists 
# so could do a Boolean of T F whether in new_area or not and print the indices of the ones not? 
outliers_0 = []
outliers_circ_0 = []
print "Outliers are gastruloids at the following indices, with the following areas: "

print "Based on area: "
for (num,item) in enumerate(areas_0):
    if item not in new_area:
        print(num, item)
        outliers_0.append(num)
print "Based on circularity: "
for (num,item) in enumerate(circularity_0):
    if item not in new_circ:
        print(num, item)
        outliers_circ_0.append(num)


# In[ ]:


file_indices = []
for i, v in enumerate(filenames[0]):
    file_indices.append(v)
    
file_indices = np.array(file_indices)
file_indices = np.ndarray.tolist(file_indices)


# In[ ]:


### This is only relevant in time courses 
outliers_24h = []
#print "outliers:" # need to also add the outliers from the previous timepoint, if this is a time-course 
outliers = [10,12,27,55,56,88]
prev_outliers = []
for (num, item) in enumerate(areas_0):
    if num in outliers_24h:
        #print(num, item)
        prev_outliers.append(num)

#### Need change this depending on if the analysis is time-course or not 

# outliers_emptych = np.unique(outliers_0 + outliers_circ_0)
toobig = []
for i, v in enumerate(areas_0):
    max = 1024 * 1024 
    if v > 0.3 * max:
        print i, v 
        toobig.append(i)


outliers_emptych = np.unique(outliers_0 + outliers_circ_0 + toobig)
print "\nCombined outliers are: "
print toobig
print "\nBefore removing area outliers:", len(areas_0)
print "Within based on area:", len(new_area)


# In[ ]:


fig, ax = plt.subplots(2, 2, sharex=False, figsize = (7, 7))

ax[0,0].hist(areas_0, bins=50, density=None, 
         weights=None, cumulative=False, 
         histtype='bar', align='left', orientation='vertical', 
         rwidth=None, log=False, color='black', label= '# Observations', stacked=False, 
         normed=None, data=None)
ax[0,1].hist(new_area, bins=50, density=None, 
             weights=None, cumulative=False, 
             histtype='bar', align='left', orientation='vertical', 
             rwidth=None, log=False, color='black', label= '# Observations', stacked=False, 
             normed=None, data=None)

# add a 'best fit' line

ax[0,0].set_xlabel('Area')
ax[0,1].set_xlabel('Area')
ax[0,0].set_ylabel('# Gastruloids')
ax[0,0].set_title('Area Distribution of Gastruloids \n with outliers')
ax[0,1].set_title('Area Distribution of Gastruloids \n no outliers')
ax[1,0].set_title('Boxplot Area Distribution \n with outliers')
ax[1,1].set_title('Boxplot Area Distribution \n no outliers')
ax[1,0].set_ylabel('Area / # Pixels')
ax[1,0].set_xlabel('Empty Channel')
ax[1,1].set_xlabel('Empty Channel')

new_area = np.array(new_area)
clean_mean = np.mean(new_area, axis=0)
ax[1,0].boxplot(areas_0, manage_xticks=True, autorange=False, zorder=None)
ax[1,1].boxplot(new_area, manage_xticks=True, autorange=False, zorder=None)
plt.tight_layout()

data_locn = '/Users/alexandrabaranowski/Desktop/Images for report /'
#plt.savefig(data_locn+'1 BraGFP 2016-03-16 72h Outlier removal1.png')


# Detecting saturation artifacts - method: % of saturated pixels

# In[ ]:


numbersaturated = []
percentsaturated = []

for k, v in enumerate(files_list[3][0]):  
    # readthe original image
    img = np.array(Image.open(files_list[3][0][k]))  

    saturated = 255
    white = np.sum(img == saturated)
    numbersaturated.append(white)
    width, height = img.shape
    all = width*height
    percent = 100.0*white/all
    percentsaturated.append(percent)

    #print "Total pixels: %d" % all, k
   # print "Saturated pixels: %d (%5.2f%%)" % (white,percent)

for i, v in enumerate(percentsaturated): 
    mean = np.mean(percentsaturated,axis = 0)
    sd = np.std(percentsaturated, axis =0)
    #print sd
    
    if v > mean + 2 * sd or v > 0.2: 
        print i, v


# In[ ]:


###### 2. using the total intensity to find outliers #########
sumofpixels1 = []
channels = [0,1,2,3]
newintensity = []
#for j in range(0, len(channels)):
for i, v in enumerate(blurred_masks_1):  
    # readthe original image
    #img = np.array(Image.open(files_list[3][0][i])) 
    img = v
    x,y,z,a = cv2.sumElems(img)
    sumofpixels1.append(x)


mean = np.mean(sumofpixels1, axis = 0)
sd = np.std(sumofpixels1, axis = 0)

newintensity = [x for x in sumofpixels1 if (x > mean - 1.5 * sd)]
newintensity = [x for x in newintensity if (x < mean + 1.5 * sd)]
print newintensity

outliers_intensity1 = []
print "Based on Sum of Pixels: "
for (num,item) in enumerate(sumofpixels1):
    if item not in newintensity:
        print(num, item)
        outliers_intensity1.append(num)
    
        
        #print sumofpixels


# In[ ]:


########### Using max and mean fluorescence to detect outliers ###########
maxfluo3 
meanfluo3

### Max Fluorescence
mean = np.mean(maxfluo1, axis = 0)
sd = np.std(maxfluo1, axis = 0)

newmaxfluo = [x for x in maxfluo1 if (x > mean - 1.5 * sd)]
newmaxfluo = [x for x in newmaxfluo if (x < mean + 1.5 * sd)]
#print newmaxfluo

outliers_maxfluo1 = []
print "Based on Max Fluorescence: "
for (num,item) in enumerate(maxfluo1):
   if item not in newmaxfluo:
       print(num, item)
       outliers_maxfluo1.append(num)
print num, item

### Mean Fluorescence
mean2 = np.mean(meanfluo1, axis = 0)
sd2 = np.std(meanfluo1, axis = 0)

newmeanfluo = [x for x in meanfluo1 if (x > mean2 - 1.5 * sd2)]
newmeanfluo = [x for x in newmeanfluo if (x < mean2 + 1.5 * sd2)]
#print newmeanfluo

outliers_meanfluo1 = []
print "Based on Mean Fluorescence: "
for (num,item) in enumerate(meanfluo1):
   if item not in newmeanfluo:
       print(num, item)
       outliers_meanfluo1.append(num)
print num, item

### Standard Deviation of Fluorescence
mean3 = np.mean(stdevfluo3, axis = 0)
sd3 = np.std(stdevfluo3, axis = 0)

newsdfluo = [x for x in stdevfluo3 if (x > mean3 - 2 * sd3)]
newsdfluo = [x for x in newsdfluo if (x < mean3 + 2 * sd3)]
#print newmeanfluo

outliers_sdfluo1 = []
print "Based on St Dev of Fluorescence: "
for (num,item) in enumerate(stdevfluo1):
   if item not in newsdfluo:
       print(num, item)
       outliers_sdfluo1.append(num)
print num, item


# # Removing outliers 

# In[ ]:


#total_outliers = outliers_0 + outliers_circ_0

##### When it is a time course and >24h AA, use this: 
total_outliers = outliers_0 + outliers_circ_0 

total_outliers = np.unique(total_outliers)
total_outliers = total_outliers[::-1]
print total_outliers

total_contour_areas = [areas_0, areas_1, areas_3]# areas_3]
total_areas = [areas1, areas2, areas3]
areas = [total_contour_areas, total_areas]

#print total_areas

for d in range(2):
    for j in range(3):
        for i in total_outliers: 
            area = areas[d][j]
            del area[i]


# In[ ]:


for j in range(2): 
    parameter = areas[j]
    for k in parameter:
        #print len(k)
        if len(k) != len(areas_0):
            print "error - check"
        else: 
            continue
print "check complete - all lengths correct "
print len(areas_0)


# Need to now remove all those indices from all the other parameters too 

# In[ ]:


perimeters = [perimeters_0, perimeters_1, perimeters_3]
circularity = [circularity_0, circularity_1, circularity_3]
both = [perimeters, circularity]

for k in range(2):
    for j in range(3):
        for i in total_outliers: 
            perimeter = both[k][j]
            del perimeter[i]


# #### Removing from the remaining parameters

# In[ ]:


hulls = [hull_0, hull_1, hull_3]
hull_areas = [hull_areas_0, hull_areas_1, hull_areas_3]
convexity_defects = [convexity_defects_0, convexity_defects_1, convexity_defects_3]
minc_radii = [minc_radii_0, minc_radii_1, minc_radii_3]
boxes = [boxes_0, boxes_1, boxes_3]

allparameters = [hulls, hull_areas, convexity_defects, minc_radii, boxes]
for k in range(5):
    for j in range(3):
        for i in total_outliers: 
            parameter = allparameters[k][j]
            del parameter[i]


# #### These are arrays so can't delete

# In[ ]:


moments = [moments_0, moments_1, moments_3]
hu = [hu_moments_0, hu_moments_1, hu_moments_3]
loghu = [log_hu, log_hu1, log_hu3]


# #### More parameters to QC 

# In[ ]:


#centroidy = [cy_0, cy_1, cy_3]
#centroidx = [cx_0, cx_1, cx_3]
ARs = [aspect_ratios_0, aspect_ratios_1, aspect_ratios_3]
convexities = [convexity_0, convexity_1, convexity_3]


moreparameters = [ ARs, convexities]#, majors, minors]

for k in range(2): # how many in moreparameters
    for j in range(3): # how many in each list
        for i in total_outliers: 
            parameter1 = moreparameters[k][j]
            del parameter1[i]


# In[ ]:


#min_ellipses = [min_ellipse_0]#, min_ellipse_1, min_ellipse_3]
majors = [major_ellipse_axis]# major_ellipse_axis1]#, major_ellipse_axis3]#, major_ellipse_axis3]
minors = [minor_ellipse_axis]# minor_ellipse_axis1]#, minor_ellipse_axis3]

ellipses_ones = [majors, minors]

for k in range(2):
    for j in range(1):
        for i in total_outliers: 
            e = ellipses_ones[k][j]
            del e[i]


# #### Make sure you do not convert nan to num before, as you cannot delete from the array; so first delete and then remove nan

# In[ ]:


I_total = [I, I_1, I_3]
logI_total = [log_i, log_ich1, log_ich3]

individual_hu = [I_total, logI_total]

for first in range(2):
    for second in range(3): 
        for inside in range(7):
            for i in total_outliers:
                p = individual_hu[first][second][inside]
                del p[i]


# In[ ]:


print len(aspect_ratios_1)


# In[ ]:


log_ich3 = np.nan_to_num(log_ich3)
log_ich1 = np.nan_to_num(log_ich1)
log_i13 = log_ich3[0]; log_i23 = log_ich3[1];log_i33 = log_ich3[2]; log_i43 = log_ich3[3]; log_i53 = log_ich3[4]; log_i63 = log_ich3[5]; log_i73 = log_ich3[6]
log_i11 = log_ich1[0]; log_i21 = log_ich1[1];log_i31 = log_ich1[2]; log_i41 = log_ich1[3]; log_i51 = log_ich1[4]; log_i61 = log_ich1[5]; log_i71 = log_ich1[6]


# ##### Haralick features

# In[ ]:


H_features = [H1ang2ndmoment, H2contrast, H3correlation, H4sumsqvar, H5invdiffmoment,H6sumavg,
              H7sumvar, H8sumentropy, H9entropy, H10diffvar, H11diffentropy, H12imcorr1, H13imcorr2]

hs = [H_features]

for first in range(1):
    for inside in range(13):
        for i in total_outliers:
            p = hs[first][inside]
            del p[i]


# ##### do a check for the individual hu as well

# In[ ]:


normalised_maxfluo1 = normalised_maxfluo1.tolist
normalised_meanfluo1 = normalised_meanfluo1.tolist
normalised_stdev1 = normalised_stdev1.tolist


# ##### Fix fluo intensity lists 

# In[ ]:


max_fluorescence = [maxfluo1, maxfluo2, maxfluo3]
mean_fluorescence = [meanfluo1, meanfluo2, meanfluo3]
st_dev_fluorescence = [stdevfluo1, stdevfluo2, stdevfluo3]

#norm_max_fluorescence = [normalised_maxfluo1]#, maxfluo2, maxfluo3]
#norm_mean_fluorescence = [normalised_meanfluo1]#, meanfluo2, meanfluo3]
#norm_st_dev_fluorescence = [normalised_stdev1]#, stdevfluo2, stdevfluo3]

intensities = [max_fluorescence, mean_fluorescence, st_dev_fluorescence]

for k in range(3):
    for j in range(3):
        for i in total_outliers: 
            pixeli = intensities[k][j]
            del pixeli[i]


# In[ ]:


min_fluorescence = [minfluo1, minfluo2, minfluo3]

for j in range(3):
    for i in total_outliers:
        mini = min_fluorescence[j]
        del mini[i]


# In[ ]:


### Channel 1 
max_max1 = np.max(maxfluo1)
normalised_maxfluo1 = (np.array(maxfluo1) / max_max1)
#print normalised_maxfluo1

max_mean1 = np.max(meanfluo1)
normalised_meanfluo1 = (np.array(meanfluo1) / max_mean1)
#print normalised_meanfluo1

max_stdev1 = np.max(stdevfluo1) 
normalised_stdev1 = (np.array(stdevfluo1) / max_stdev1)
#print normalised_stdev1

### Channel 2 
max_max2 = np.max(maxfluo2)
normalised_maxfluo2 = (np.array(maxfluo2) / max_max2)
#print normalised_maxfluo2

max_mean2 = np.max(meanfluo2)
normalised_meanfluo2 = (np.array(meanfluo2) / max_mean2)
#print normalised_meanfluo2

max_stdev2 = np.max(stdevfluo2) 
normalised_stdev2 = (np.array(stdevfluo2) / max_stdev2)
#print normalised_stdev2

### Channel 3 
max_max3 = np.max(maxfluo3)
normalised_maxfluo3 = (np.array(maxfluo3) / max_max3)
#print normalised_maxfluo3

max_mean3 = np.max(meanfluo3)
normalised_meanfluo3 = (np.array(meanfluo3) / max_mean3)
#print normalised_meanfluo3

max_stdev3 = np.max(stdevfluo3) 
normalised_stdev3 = (np.array(stdevfluo3) / max_stdev3)
#print normalised_stdev3


# #### check all lengths are correct - this is for "allparameters"

# In[ ]:


for j in range(5): 
    parameter = allparameters[j]
    for k in parameter:
        #print len(k)
        if len(k) != len(areas_0):
            print "error - check"
        else: 
            continue
print "check complete - all lengths correct "
print len(circularity_1)


# #### Check all lengths are correct - this is for "ellipses - ones"

# In[ ]:


for j in range(2): 
    e = ellipses_ones[j]
    for k in e:
        #print len(k)
        if len(k) != len(areas_0):
            print "error - check"
        else: 
            continue
print "check complete - all lengths correct"


# #### check all lengths are correct - this is for "moreparameters"

# In[ ]:


for j in range(2): 
    parameter = moreparameters[j]
    for k in parameter:
        #print len(k)
        if len(k) != len(areas_0):
            print "error - check"
        else: 
            continue
print "check complete - all lengths correct "


# #### check all lengths are correct - this is for "intensities"

# In[ ]:


for j in range(3): 
    pixeli = intensities[j]
    for k in pixeli:
        #print len(k)
        if len(k) != len(areas_0):
            print "error - check"
        else: 
            continue
print "check complete - all lengths correct "


# # Finding the areas of each channel and proportions of gastruloid covered by each marker

# ### Looking at the distribution (reproducibility) of the expression area & proportion of each marker across all the gastruloids 

# In[ ]:


areas_0 = np.array(areas_0, dtype=np.float)
areas1 = np.array(areas1, dtype=np.float)
areas2 = np.array(areas2, dtype=np.float)
areas3 = np.array(areas3, dtype=np.float)

prop_ch1 = areas1 / areas_0
prop_ch2 = areas2 / areas_0
prop_ch3 = areas3 / areas_0

print len(prop_ch1)


# In[ ]:


## Some preliminary graphs of the overlap 

sns.set(color_codes=True)

fig, axes = plt.subplots(3,3, figsize = (15,10), sharey=False)

# this has drawn a histogram and plotted a kernel density estimate 

prop_ch1 = pd.Series(prop_ch1, name="Distribution of Sox2 area")
prop_ch2 = pd.Series(prop_ch2, name="Distribution of Sox17 area")
prop_ch3 = pd.Series(prop_ch3, name="Distribution of T/Bra area")

sns.distplot(prop_ch1, color = "y", ax=axes[0,0]) 
sns.distplot(prop_ch2, color = "y", ax=axes[0,1])
sns.distplot(prop_ch3, color = "y", ax=axes[0,2])

# this is a histogram with rugplot instead of kde;the rugplot adds a tick at each observation
sns.distplot(prop_ch1, kde=False, rug=True, ax=axes[2,0])
sns.distplot(prop_ch2, kde=False, rug=True, ax=axes[2,1])
sns.distplot(prop_ch3, kde=False, rug=True, ax=axes[2,2])

sns.kdeplot(prop_ch1, color = "y", ax=axes[1,0], shade=True)
sns.kdeplot(prop_ch2, color = "y", ax=axes[1,1], shade=True)
sns.kdeplot(prop_ch3, color = "y", ax=axes[1,2], shade=True)

plt.tight_layout()

locn = data_locn_mask[0]
plt.savefig(locn+"Wnt3a final 20190429.png")


# ### Finding overlap between the channels 

# #### 1. Looking at just the areas that they overlap 

# In[ ]:


### Ch1 and Ch2 
percent12totalo = []
ch121overlap = []
ch122overlap = []
percent12gld = []
union12 = []

for i, v in enumerate(areas_0):
    ch1_area = areas_1[i]
    ch2_area = areas2[i]
    overlap12 =  np.float(cv2.bitwise_and(ch1_area,ch2_area)[0])
    union_12area = ch1_area + ch2_area - overlap12
    union12.append(union_12area)
    
    print '\nUnion of areas', union_12area
    print 'Ch1 area', ch1_area
    print 'Ch2 area', ch2_area
    print 'Overlap', overlap12
    
    # Overlap as percent of total area covered by those two 
    try:
        percent12total = overlap12 / union_12area
        ch121covered = overlap12 / ch1_area 
        ch122covered = overlap12 / ch2_area
        percent12gldd = overlap12 / v
        
        percent12totalo.append(percent12total*100)
        ch121overlap.append(ch121covered*100)
        ch122overlap.append(ch122covered*100)
        percent12gld.append(percent12gldd*100)
    
    except ZeroDivisionError:
        percent12total = 0
        ch121covered = 0 
        ch122covered = 0 
        percent12gldd = 0
        
        percent12totalo.append(0)
        ch121overlap.append(0)
        ch122overlap.append(0)
        percent12gld.append(0)
    
    print '% Union Area 1 and 2 covered', percent12total *100
    
    # How much of Ch1 area it covers 
    print '% Ch1 covered by overlap', ch121covered*100
    
    # How much of Ch2 area it covers 
    print '% Ch2 covered by overlap', ch122covered*100
    # Overlap as percent of gastruloid area
    print '% Gldcovered by overlap', percent12gldd*100
    


# In[ ]:


### Ch2 and Ch3 
percent23totalo = []
ch232overlap = []
ch233overlap = []
percent23gld = []
union23 = []

for i, v in enumerate(areas_0):
    ch2_area = areas2[i]
    ch3_area = areas_3[i]
    overlap23 =  np.float(cv2.bitwise_and(ch2_area,ch3_area)[0])
    union_23area = ch2_area + ch3_area - overlap23
    union23.append(union_23area)
    
    print '\nUnion of areas', union_23area
    print 'Ch2 area', ch2_area
    print 'Ch3 area', ch3_area
    print 'Overlap', overlap23
  
    try:
        percent23total = overlap23 / union_23area
        ch232covered = overlap23 / ch2_area
        ch233covered = overlap23 / ch3_area
        percent23gldd = overlap23 / v
        
        percent23totalo.append(percent23total*100)
        ch232overlap.append(ch232covered*100)
        ch233overlap.append(ch233covered*100)
        percent23gld.append(percent23gldd*100)
    
    except ZeroDivisionError:
        percent23total = 0
        ch232covered = 0 
        ch233covered = 0 
        percent23gldd = 0
        
        percent23totalo.append(0)
        ch232overlap.append(0)
        ch233overlap.append(0)
        percent23gld.append(0)
    
    print '% of Union of Areas covered', percent23total*100
    
    # How much of Ch1 area it covers 
    print '% Ch2 covered by overlap', ch232covered*100
    
    # How much of Ch2 area it covers 
    print '% Ch3 covered by overlap', ch233covered *100
    
    # Overlap as percent of gastruloid area
    print '% Gld covered by overlap', percent23gldd *100


# In[ ]:


### Ch1 and Ch3 
percent13totalo = []
ch131overlap = []
ch133overlap = []
percent13gld = []
union13 = []

for i, v in enumerate(areas_0):
    ch1_area = areas_1[i]
    ch3_area = areas_3[i]
    overlap13 =  np.float(cv2.bitwise_and(ch1_area,ch3_area)[0])
    union_13area = ch1_area + ch3_area - overlap13
    union13.append(union_13area)
    
    print '\nUnion of areas', union_13area
    print 'Ch1 area', ch1_area
    print 'Ch3 area', ch3_area
    print 'Overlap', overlap13

    try:
        percent13total = overlap13/ union_13area
        ch131covered = overlap13 / ch1_area
        ch133covered = overlap13 / ch3_area
        percent13gldd = overlap13 / v
        
        percent13totalo.append(percent13total)
        ch131overlap.append(ch131covered)
        ch133overlap.append(ch133covered)
        percent13gld.append(percent13gldd)
    
    except ZeroDivisionError:
        percent13total = 0
        ch131covered = 0 
        ch133covered = 0 
        percent13gldd = 0
        
        percent13totalo.append(0)
        ch131overlap.append(0)
        ch133overlap.append(0)
        percent13gld.append(0)
    
    print '% of Union of Areas covered', percent13total*100
    
    # How much of Ch1 area it covers 
    print '% Ch1 covered by overlap', ch131covered*101
    
    # How much of Ch2 area it covers 
    print '% Ch3 covered by overlap', ch133covered *100
    
    # Overlap as percent of gastruloid area
    print '% Gld covered by overlap', percent13gldd *100


# In[ ]:


##### plot pixel intensity values against each other 
max1 = []
ones = []
max2 = []
twos = []
for i, v in enumerate(files_list[1][0][0:2]):
    img1 = Image.open(files_list[1][0][0:2][i])
    ones.append(img1)
    img2 = Image.open(files_list[2][0][0:2][i])
    twos.append(img2)
    
    max1.append(np.max(img1))
    max2.append(np.max(img2))

maxmax1 = np.max(max1)
maxmax2 = np.max(max2)

for i, v in enumerate(ones):
    img1 = np.array(v) / maxmax1
    
    img2 = np.array(twos[i]) / maxmax2
    
    
    plt.scatter(img1, img2, alpha = 0.5)
    plt.show


# In[ ]:


for i, v in enumerate(blurred_2[0:1]):
    #img = np.array(Image.open(files_list[1][0][0:2][i]))
    img = np.array(v)
    img1 = img.copy()

    imgg = np.array(blurred_3[i])
    img2 = imgg.copy()
    #img2 = cv2.GaussianBlur(img2, (15,15), 0)
    
    data1 = np.array(img1)
    data2 = np.array(img2)
    
    #colors = ('data1', 'data2')
    #fig, ax = plt.subplot(1,1, figsize = (10,3))
    
    plt.scatter(data1, data2, alpha = 0.1, cmap='viridis')
    plt.xlim((0,255))
    plt.ylim((0,255))


# In[ ]:


x = blurred_2[0]
y = blurred_3[0]

def r2(x1, y1):
    return np.float(stats.pearsonr(x, y)[0] ** 2)
sns.jointplot(x, y, kind="reg")


# ## Making a dataframe to add all the parameters we have, with each parameter as a column and each gld as a row

# In[ ]:


d = {'TreatmentConditions': 'Drug Experiments','Date': '2019-03-04', 'Hrs Post A': '72', 
      'Cell line': 'RUES2', 'Pre treatment': 'Wnt3a', 'Treatment': '0.5 Chi 3 Chi', 'Cell number': '400',
     'Area Ch0': areas_0, 'Total Area Ch1': areas1, 'Contour Area Ch1': areas_1, 
     'Total Area Ch2': areas2, 'Total Area Ch3': areas3, 'Contour Area Ch3': areas_3,
      'Perimeter Ch0': perimeters_0,'Perimeter Ch1': perimeters_1, 'Perimeter Ch3': perimeters_3, 
    'Circularity Ch0': circularity_0, 'Circularity Ch1': circularity_1,'Circularity Ch3': circularity_3,
     'Min rot rect Ch0': boxes_0, 'Min rot rect Ch1': boxes_1, 'Min rot rect Ch3': boxes_3, 
     'Proportion Ch1': prop_ch1, 'Proportion Ch2': prop_ch2, 'Proportion Ch3': prop_ch3, 
    'MECradius Ch0': minc_radii_0, 'MECradius Ch1': minc_radii_1, 'MECradius Ch3': minc_radii_3,
     'Hu I1 Ch0': I1, 'Hu I2 Ch0': I2,'Hu I3 Ch0': I3,'Hu I4 Ch0': I4,'Hu I5 Ch0': I5,'Hu I6 Ch0': I6,'Hu I7 Ch0': I7,
     'Hu I1 Ch1': I1_1, 'Hu I2 Ch1': I2_1,'Hu I3 Ch1': I3_1,'Hu I4 Ch1': I4_1,'Hu I5 Ch1': I5_1,'Hu I6 Ch1': I6_1,'Hu I7 Ch1': I7_1,
     'Hu I1 Ch3': I1_3, 'Hu I2 Ch3': I2_3,'Hu I3 Ch3': I3_3,'Hu I4 Ch3': I4_3,'Hu I5 Ch3': I5_3,'Hu I6 Ch3': I6_3,'Hu I7 Ch3': I7_3,
     'Log Hu I1 Ch0': log_i1, 'Log Hu I2 Ch0': log_i2,'Log Hu I3 Ch0': log_i3,'Log Hu I4 Ch0': log_i4,'Log Hu I5 Ch0': log_i5,'Log Hu I6 Ch0': log_i6,'Log Hu I7 Ch0': log_i7,
     'Log Hu I1 Ch1': log_i11, 'Log Hu I2 Ch1': log_i21,'Log Hu I3 Ch1': log_i31,'Log Hu I4 Ch1': log_i41,'Log Hu I5 Ch1': log_i51,'Log Hu I6 Ch1': log_i61,'Log Hu I7 Ch1': log_i71,
     'Log Hu I1 Ch3': log_i13, 'Log Hu I2 Ch3': log_i23,'Log Hu I3 Ch3': log_i33,'Log Hu I4 Ch1': log_i43,'Log Hu I5 Ch3': log_i53,'Log Hu I6 Ch3': log_i63,'Log Hu I7 Ch3': log_i73,
     'Hull Area Ch0': hull_areas_0, 'Hull Area Ch1': hull_areas_1,'Hull Area Ch3': hull_areas_3,
     'Convexity Ch0': convexity_0, 'Convexity Ch1': convexity_1,'Convexity Ch3': convexity_3,
     'AR Ch0': aspect_ratios_0, 'AR Ch1': aspect_ratios_1,'AR Ch3': aspect_ratios_3,
      'H1Ang2ndmoment': H1ang2ndmoment, 'H2Contrast': H2contrast, 
     'H3 Correlation': H3correlation, 'H4 Sum Sq Var': H4sumsqvar, 'H5 Ind Diff Moment': H5invdiffmoment,
     'H6 Sum Avg': H6sumavg,'H7 Sum Var': H7sumvar, 'H8 Sum Entropy': H8sumentropy, 'H9 Entropy': H9entropy, 
     'H10 Diff Var': H10diffvar, 'H11 Diff Entropy': H11diffentropy, 'H12 Im Corr 1': H12imcorr1,  
     'H13 Im Corr 2': H13imcorr2,'Norm Max Fluorescence Ch1': normalised_maxfluo1, 'Norm Mean Fluorescence Ch1': normalised_meanfluo1,
     'Norm St Dev Fluorescence Ch1': normalised_stdev1, 'Norm Max Fluorescence Ch2': normalised_maxfluo2, 'Norm Mean Fluorescence Ch2': normalised_meanfluo2,
     'Norm St Dev Fluorescence Ch2': normalised_stdev2, 'Norm Max Fluorescence Ch3': normalised_maxfluo3, 'Norm Mean Fluorescence Ch3': normalised_meanfluo3,
     'Norm St Dev Fluorescence Ch3': normalised_stdev3, 'Overlap of 1 3 %1 3area':percent13totalo,'Ch1 3 Overlap %Ch1 Fluo':ch131overlap,
     'Ch1 3 Overlap % Ch3 Fluo':ch133overlap, 'Ch 1 3 Overlap % Gld': percent13gld,'Overlap of 1 2 %1 2area':percent12totalo,'Ch1 2 Overlap %Ch1 Fluo':ch121overlap,
     'Ch1 2 Overlap % Ch2 Fluo':ch122overlap, 'Ch 1 2 Overlap % Gld': percent12gld,'Overlap of 2 3 %2 3area':percent23totalo,'Ch2 3 Overlap %Ch2 Fluo':ch232overlap,
     'Ch2 3 Overlap % Ch3 Fluo':ch233overlap, 'Ch 2 3 Overlap % Gld': percent23gld, 'St Dev Fluorescence Ch1': stdevfluo1, 'Max Fluorescence Ch2': maxfluo2, 'Mean Fluorescence Ch2': meanfluo2,
     'St Dev Fluorescence Ch2': stdevfluo2, 'Max Fluorescence Ch3': maxfluo3, 'Mean Fluorescence Ch3': meanfluo3, 'St Dev Fluorescence Ch3': stdevfluo3, 'Max Fluorescence Ch1': maxfluo1,
     'Mean Fluorescence Ch1': meanfluo1, 'Union Ch1 Ch2': union12, 'Union Ch2 Ch3': union23, 'Union Ch1 Ch3': union13
    }
df = pd.DataFrame(data=d)
df


# Save as CSV

# In[ ]:


path = '/Users/alexandrabaranowski/Desktop/Project Code/20192904/'

df.to_csv(path+'2019-03-04 Wnt3a 72h 2904.csv') ## this saves it to the Project Code Folder - find a way to 
# save it to a more clear folder and for it to automatically take the name of the original folder

################ N.B. Because conditions are added manually to the df the subgroups are not added at this point, 
#indices that were removed need to then be accounted for when the individual subgroups are added to the csv file manually;
# indices can be found in the first cell of the 'Removing Outliers' section ###################

