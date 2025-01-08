# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:46:40 2024

@author: Courtney
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
import cv2
import pandas as pd
from scipy import ndimage
import warnings
from skimage.feature import peak_local_max
from skimage.segmentation import clear_border, watershed
warnings.filterwarnings("ignore")


def main(): 
    global files_ch0, image_sets, ch1_counts, ch0_counts
    print('\n\n\nThis program will Count the number of cells in an image. \n')
    # User Input Questions 
    wdir=str(input('What is the location of your folder with single channel planes: '))
    os.chdir(wdir)
    for root, dirs, files in os.walk("."):
        files=[FILE for FILE in files if 'ch0' in FILE]
    question=str(input('Do you want to see the masks that are being Counted? (Y/N): ')).upper()
    genotype=sorted(str(input('What Genotypes are you Analyzing? (Comma Separated List ", "): ')).split(', '))
    channel=str(input('What channel are you analyzing (pick 2)? (ch00, ch01, ch02, ch03): ')).split(', ')
    
    # Identify AB Stains and Counter for Analysis
    ab_list = [str(input(f'What Antibody/Stain is in {channel[i]}: ')).upper() for i in range(len(channel))]
    counter_list = [str(input(f'What Compartment is {channel[i]} staining? (nuclear = n, cytoplasm = c): ')).upper() for i in range(len(channel))]
   
    count = len(channel) 
    # Initialize channel_counts
    channel_counts = {
        'ch0_counts': {},
        'ch1_counts': {},
    }
    
    if count >= 3:
        channel_counts['ch2_counts'] = {}
        channel_counts['ch12_counts'] = {}
    
    if count == 4:
        channel_counts['ch3_counts'] = {}
        channel_counts['ch13_counts'] = {}
        channel_counts['ch23_counts'] = {}
    
    # Process genotypes
    print('\n\n\n Counting Cells:  ' + str(ab_list) + '....\n\n')      
    for genotype1 in genotype:
        # Initialize lists for counts for each channel
        ch0_count, ch1_count, ch2_count, ch3_count = [], [], [], []
        if len(channel) >= 3:
            ch12_count = []
        if len(channel) == 4:
            ch13_count, ch23_count = [], []
    
        files_genotype = [FILE for FILE in files if genotype1 in FILE]  # Files for the current genotype
        image_sets = unique_images(files_genotype)  # Get unique image sets
    
        for image_set in image_sets:
            # Create the list of files for the current image set and genotype
            files_ch0 = []
            for j in range(len(channel)):
                # Check if the file contains the corresponding channel info
                files_ch0 += [FILE for FILE in files_genotype if str(image_set + '_' + channel[j]) in FILE]
    
            # Process the first channel (mask_ch0)
            mask_ch0 = filter_mask(files_ch0[0])
            ch0_count.append(object_counter(mask_ch0))  # Count for ch0
    
            if count >= 2 and len(files_ch0) > 1:  # Ensure there are enough files for channel 1
                mask_ch1 = (process_channel(mask_ch0, files_ch0[1], counter_list[1])>0)*mask_ch0
                ch1_count.append(object_counter(mask_ch1))  # Count for ch1

            if count >= 3 and len(files_ch0) > 2:  # Ensure there are enough files for channel 2
                mask_ch2 = (process_channel(mask_ch0, files_ch0[2], counter_list[2])>0)*mask_ch0
                ch2_count.append(object_counter(mask_ch2))  # Count for ch2
                ch12_count.append(object_counter((mask_ch1 > 0) & (mask_ch2 > 0)))  # Count for ch1 * ch2
    
            if count == 4 and len(files_ch0) > 3:  # Ensure there are enough files for channel 3
                mask_ch3 = (process_channel(mask_ch0, files_ch0[3], counter_list[3])>0)*mask_ch0 
                ch3_count.append(object_counter(mask_ch3))  # Count for ch3
                ch13_count.append(object_counter((mask_ch1 > 0) & (mask_ch3 > 0)))  # Count for ch1 * ch3
                ch23_count.append(object_counter((mask_ch2 > 0) & (mask_ch3 > 0)))  # Count for ch2 * ch3
    
            # Display images if needed
            if question == 'Y':
                display_every_10(''.join(files_ch0), mask_ch0, ab_list[0])  # For mask_ch0
                if count >= 2:  # Ensure ch1 exists
                    display_every_10(''.join(files_ch0), mask_ch1, ab_list[1])  # For mask_ch1
                if count >= 3:  # Ensure ch2 exists
                    display_every_10(''.join(files_ch0), mask_ch2, ab_list[2])  # For mask_ch2
                if count == 4:  # Ensure ch3 exists
                    display_every_10(''.join(files_ch0), mask_ch3, ab_list[3])  # For mask_ch3

    
        # Save counts into the channel_counts dictionary for the current genotype
        channel_counts['ch0_counts'][genotype1] = ch0_count
        channel_counts['ch1_counts'][genotype1] = ch1_count
        if count >= 3:
            channel_counts['ch2_counts'][genotype1] = ch2_count
            channel_counts['ch12_counts'][genotype1] = ch12_count
        if count == 4:
            channel_counts['ch3_counts'][genotype1] = ch3_count
            channel_counts['ch13_counts'][genotype1] = ch13_count
            channel_counts['ch23_counts'][genotype1] = ch23_count
    
    # Convert counts to DataFrames
    count_dfs = []
    sheet_names = []  # To keep track of sheet names for Excel
    
    # Export individual channel counts
    for j in range(len(channel)):
        ch_key = f'ch{j}_counts'
        if ch_key in channel_counts:  # Check if the key exists
            count_dfs.append(pd.DataFrame(dict([(k, pd.Series(v)) for k, v in channel_counts[ch_key].items()])))
            sheet_names.append(f"{ab_list[j]} Counts")  # Add corresponding antibody name
    
    # Export combined channel counts (ch12, ch13, ch23)
    if count >= 3:
        if 'ch12_counts' in channel_counts:
            count_dfs.append(pd.DataFrame(dict([(k, pd.Series(v)) for k, v in channel_counts['ch12_counts'].items()])))
            sheet_names.append("Ch1_Ch2 Counts")
    if count == 4:
        if 'ch13_counts' in channel_counts:
            count_dfs.append(pd.DataFrame(dict([(k, pd.Series(v)) for k, v in channel_counts['ch13_counts'].items()])))
            sheet_names.append("Ch1_Ch3 Counts")
        if 'ch23_counts' in channel_counts:
            count_dfs.append(pd.DataFrame(dict([(k, pd.Series(v)) for k, v in channel_counts['ch23_counts'].items()])))
            sheet_names.append("Ch2_Ch3 Counts")
    
    # Export counts to Excel
    print("\n\nExporting to Excel....")
    with pd.ExcelWriter('Cell_Counts.xlsx') as writer:
        for df, sheet_name in zip(count_dfs, sheet_names):
            df.to_excel(writer, sheet_name=sheet_name)

""""
Functions to process images imput into program:
    Filter Mask - Make an mask of alveolar surface with no gaps in walls and 
                    alveolar pores filled in
    Remove Small Objects - Remove any specks that could hinder accuracy
    Invert Mask - invert mask from alveolar tissue to airspace
    Remove Border Objects - Remove any objects touching the border
    Full Alveoli - Return only full alveoli that are not touching the border of 
                    the image
    
Input: image or mask
Output: final mask of alveolar surface or alveolar airspaces to analyze

Use: These functions are used to generate images for analysis to measure 
    different lung features
"""
def filter_mask(image):
    # Load image in grayscale
    im = cv2.imread(image, cv2.IMREAD_GRAYSCALE)
    # Normalize image
    im = im.astype('float64')
    im -= im.min()  # Shift to zero
    im /= im.max()  # Scale to [0, 1]
    im *= 255  # Scale back to [0, 255]  
    im = im.astype('uint8')
    # Apply CLAHE for better contrast (Adaptive Histogram Equalization)
    clahe = cv2.createCLAHE(clipLimit=1.0, tileGridSize=(15, 15))
    im = clahe.apply(im)
    # Use Otsu's method for thresholding
    _, mask = cv2.threshold(im, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    # Morphological operations to remove noise and improve segmentation
    mask = ndimage.binary_erosion(mask, structure=np.ones((5, 5)), iterations=1).astype(mask.dtype)  # Erosion with iteration
    mask = ndimage.binary_dilation(mask, structure=np.ones((5, 5)), iterations=1).astype(mask.dtype)  # Dilation to restore structure
    label_im, nb_labels = ndi.label(mask)
    # Distance transform on the skeleton
    distance = ndi.distance_transform_edt(mask)
    # Local maxima for watershed segmentation
    local_maxi = peak_local_max(distance, indices=False, min_distance=60, labels=label_im)
    markers, _ = ndi.label(local_maxi)
    # Watershed segmentation
    mask = watershed(-distance, markers, mask=mask)
    # Clear objects touching the image border
    mask = clear_border(mask)  # This removes objects touching the borders
    # Remove small and large objects for refinement
    mask = remove_small_objects(mask, 10)
    mask = remove_large_objects(mask, 40000)
    return mask

def remove_small_objects(mask, min_size):
    label_im, nb_labels = ndi.label(mask)
    sizes = ndi.sum(mask, label_im, range(nb_labels + 1))
    # Create a boolean mask for small objects
    small_objects = sizes < min_size
    mask[small_objects[label_im]] = 0  # Remove small objects
    return mask > 0  # Return binary mask

def remove_large_objects(mask, max_size):
    label_im, nb_labels = ndi.label(mask)
    sizes = ndi.sum(mask, label_im, range(nb_labels + 1))
    # Create a boolean mask for large objects
    large_objects = sizes > max_size
    mask[large_objects[label_im]] = 0  # Remove large objects
    return mask > 0  # Return binary mask



"""
Functions to count X in an image

Input: Image mask
Output: total pixels, cells
"""

def object_counter(mask):
    label_im, nb_labels = ndimage.label(mask)
    return nb_labels

def process_channel(mask_ch0, file, counter_type):
    return nuclear_mask(mask_ch0, file) if counter_type == 'N' else cytoplasmic_mask(mask_ch0, file)


def nuclear_mask(dapi, image):
    im=filter_mask(image)
    mask=dapi*im
    mask=mask>.05*mask.max()
    return mask

def cytoplasmic_mask(dapi_mask, image):
    mask = filter_mask(image)
    # Get the dynamic circle radius based on average region size from DAPI
    circle_radius = average_radius(dapi_mask)
    # Find labeled objects in the mask
    labeled_mask, num_features = ndi.label(mask, return_num=True)
    # Create an empty image for circles
    circles_image = np.zeros_like(mask)
    # Loop through each labeled region
    for region in range(1, num_features + 1):
        # Get the coordinates of the region
        coords = np.column_stack(np.where(labeled_mask == region))
        # Calculate the centroid (center of mass) of the region
        centroid_y, centroid_x = np.mean(coords, axis=0).astype(int)
        # Draw a circle at the centroid
        cv2.circle(circles_image, (centroid_x, centroid_y), circle_radius, (255), thickness=-1)
    # Element-wise multiplication with DAPI mask
    combined_mask = (dapi_mask * circles_image.astype(np.float32)) > 0  # Convert to binary (True or False)
    return combined_mask.astype(np.uint8)

def average_radius(dapi_mask):
    # Find labeled objects in the DAPI mask
    labeled_mask, num_features = ndi.label(dapi_mask)
    # Get the sizes (areas) of the labeled regions
    region_sizes = [np.sum(labeled_mask == region) for region in range(1, num_features + 1)]
    # Compute average region size
    if len(region_sizes) == 0:
        return 0  # If no regions are found, return 0
    average_area = np.mean(region_sizes)
    # Estimate radius (circle area = pi * r^2, thus r = sqrt(area / pi))
    radius = np.sqrt(average_area / np.pi)
    return int(radius)

"""
Functions to display masks analysis complete

Input: Mask
Output: Mask displayed in Plots Tab
"""
def show_plane(ax, plane, cmap, title=None):
    ax.imshow(plane, cmap=cmap)
    ax.axis("off")
    if title:
        ax.set_title(title)   
def display_every_10(FILE, mask, identifier, colormap='rainbow'):
    # Check if FILE contains specific substrings
    z_values = ['z10', 'z20', 'z30', 'z40', 'z50', 'z60']
    if any(z in FILE for z in z_values):
        # Display the mask using the provided function
        display_masks(mask, identifier, colormap)
def display_masks(nuclei, ch00, scale):
    # Convert the input mask to float32
    nuclei = np.float32(nuclei)
    # Label the connected components in the nuclei mask
    labeled_nuclei, num_labels = ndimage.label(nuclei)
    # Create a color map for the labeled objects
    # Initialize an empty RGB image for visualization
    colored_image = np.zeros((*labeled_nuclei.shape, 3))  # Create an empty color image
    # Assign a random color to each label
    for i in range(1, num_labels + 1):
        mask = (labeled_nuclei == i)
        # Generate a random color (R, G, B) values between 0 and 1
        color = np.random.rand(3)  # Random color
        colored_image[mask] = color  # Assign random color to the labeled region
    
    # Create the plot
    _, (a) = plt.subplots(ncols=1, figsize=(15, 5))
    a.imshow(colored_image)  # Display the colored image
    a.axis("off")
    a.set_title(str(ch00))

"""
Other Functions to Simplify code
"""
def unique_images(files_genotype):
    image_sets=[]
    for item in range(len(files_genotype)):
        image_sets=np.append(image_sets,str(files_genotype[item].split('_ch')[0]))
        image_sets=list(set(image_sets))    
        # Cycle Through Each Image that Fits each Image Set (ch00, ch01, ch02, ch03)
    return image_sets

def define_method(num, dom):
    return num, dom                              

    
# This provided line is a required at the end of a Python file to call main()
if __name__ =='__main__':
    main()