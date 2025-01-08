# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:46:40 2024

@author: Courtney
"""
#Import Packages:
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
import cv2
import pandas as pd

import warnings
warnings.filterwarnings("ignore")


def main():
    global percentoftot, files_ch0, files_ch1, image_sets
    print('\n\n\nThis program will calculate the % of x pixels in y. \n')

    # User Input Questions
    wdir = input('What is the location of your folder with single channel planes: ').strip()
    try:
        os.chdir(wdir)
    except FileNotFoundError:
        print(f"Error: Directory '{wdir}' not found.")
        return

    files = [FILE for FILE in os.listdir() if 'ch0' in FILE]
    question = input('Do you want to see the masks that are being Counted? (Y/N): ').strip().upper()
    genotype = input('What Genotypes are you Analyzing? (Comma Separated List ", "): ').split(', ')
    channel = input('What channel are you analyzing (pick 2)? (ch00, ch01, ch02, ch03): ').split(', ')

    if len(channel) != 2:
        print("Error: You must select exactly two channels.")
        return

    ch0 = input(f'What Antibody/Stain is in {channel[0]}: ').strip().upper()
    ch1 = input(f'What Antibody/Stain is in {channel[1]}: ').strip().upper()

    # Initialize Dictionaries for Excel Export
    percentoftot = {}
    total_pixels_ch0 = {}
    total_pixels_ch1 = {}

    ratio = input(f'What AB will represent TOTAL: {ch0}, {ch1}. ').strip().upper()
    if ratio not in [ch0, ch1]:
        print("Error: Ratio must be one of the selected antibodies.")
        return

    num, dom = define_method(ch0, ch1) if ratio == ch0 else define_method(ch1, ch0)

    # Counting Cells
    print(f'\n\n\n Calculating {num}/{dom}....\n\n')

    for genotype1 in genotype:
        print(genotype1.upper())
        files_genotype = [FILE for FILE in files if genotype1 in FILE]
        image_sets = unique_images(files_genotype)
        percenttot = []

        # Initialize pixel counts for this genotype
        total_pixels_ch0[genotype1] = 0
        total_pixels_ch1[genotype1] = 0

        for image in image_sets:
            files_ch0 = list({FILE for FILE in files_genotype if f"{image}_{channel[0]}" in FILE})
            files_ch1 = list({FILE for FILE in files_genotype if f"{image}_{channel[1]}" in FILE})

            # Ensure files_ch0 and files_ch1 are matched
            if len(files_ch0) != len(files_ch1):
                print(f"Warning: Mismatch in file count for {image}. Skipping...")
                continue

            # Process Images to Binary
            for file_ch0, file_ch1 in zip(files_ch0, files_ch1):
                if ratio == ch0:
                    mask_dom = filter_mask(file_ch0)
                    mask_num = filter_mask(file_ch1) * mask_dom
                else:
                    mask_dom = filter_mask(file_ch1)
                    mask_num = filter_mask(file_ch0) * mask_dom

                # Generate Pixel Count
                try:
                    dom_count = pixel_counter(mask_dom)
                    num_count = pixel_counter(mask_num)

                    if dom_count == 0:  # Avoid division by zero
                        percenttot.append(0)
                    else:
                        percenttot.append(num_count / dom_count)

                    # Add to total pixel counts
                    total_pixels_ch0[genotype1] += pixel_counter(filter_mask(file_ch0))
                    total_pixels_ch1[genotype1] += pixel_counter(filter_mask(file_ch1))

                except Exception as e:
                    print(f"Error processing {file_ch0} and {file_ch1}: {e}")
                    continue

                if question == 'Y':
                    display_every_10(file_ch0, mask_num, mask_dom, num, dom)

            percentoftot[genotype1] = percenttot

    # Convert to DataFrames
    percentoftot_df = pd.DataFrame({k: pd.Series(v) for k, v in percentoftot.items()})
    total_pixels_df = pd.DataFrame({
        "Genotype": list(total_pixels_ch0.keys()),
        "Total Pixels Channel 0": list(total_pixels_ch0.values()),
        "Total Pixels Channel 1": list(total_pixels_ch1.values()),
    })

    # Export to Excel
    print("\n\nExporting to Excel....")
    output_file = 'Percent_of_Total.xlsx'
    with pd.ExcelWriter(output_file) as writer:
        percentoftot_df.to_excel(writer, sheet_name=f"{num}-{dom}")
        total_pixels_df.to_excel(writer, sheet_name="Total Pixels", index=False)
    print(f"Export completed: {output_file}")



        
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
    # Apply Gaussian Blur to reduce noise
    im = cv2.GaussianBlur(im, (5, 5), sigmaX=0)
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
    # Morphological operations to remove noise and refine segmentation
    # Grow: Dilation to merge nearby pixels
    mask = ndi.binary_dilation(mask, structure=np.ones((5, 5)), iterations=1).astype(mask.dtype)
    # Shrink: Erosion to return to original shape
    mask = ndi.binary_erosion(mask, structure=np.ones((5, 5)), iterations=1).astype(mask.dtype)
    
    return mask


"""
Functions to count the total number of pixels in an image

Input: Image mask
Output: total pixels
"""
def pixel_counter(mask):
    return np.sum(mask*1)
    
    
"""
Functions to display masks analysis complete

Input: Mask
Output: Mask displayed in Plots Tab
"""
def remove_small_objects(mask, min_size):
    label_im, nb_labels = ndi.label(mask)
    sizes = ndi.sum(mask, label_im, range(nb_labels + 1))
    # Create a boolean mask for small objects
    small_objects = sizes < min_size
    mask[small_objects[label_im]] = 0  # Remove small objects
    return mask > 0  # Return binary mask

def show_plane(ax, plane, cmap, title=None):
    ax.imshow(plane, cmap=cmap)
    ax.axis("off")
    if title:
        ax.set_title(title)   
def display_masks(nuclei, ch00, scale):
    nuclei=np.float32(nuclei)
    (n_row, n_col) = nuclei.shape
    _, (a) = plt.subplots(ncols=1, figsize=(15, 5))
    show_plane(a, nuclei, cmap=scale, title=str(ch00))      
    nuclei=0 
def display_every_10(FILE, mask_num,mask_dom, num, dom):
    if 'z10' in FILE:
        display_masks(mask_num, str(num), 'gray')
        display_masks(mask_dom, str(dom), 'gray')
        
    if 'z20' in FILE:
        display_masks(mask_num, str(num), 'gray')
        display_masks(mask_dom, str(dom), 'gray')
            
    if 'z30' in FILE:
        display_masks(mask_num, str(num), 'gray')
        display_masks(mask_dom, str(dom), 'gray')
    if 'z40' in FILE:
        display_masks(mask_num, str(num), 'gray')
        display_masks(mask_dom, str(dom), 'gray')
    if 'z50' in FILE:
        display_masks(mask_num, str(num), 'gray')
        display_masks(mask_dom, str(dom), 'gray')
    if 'z60' in FILE:
            display_masks(mask_num, str(num), 'gray')
            display_masks(mask_dom, str(dom), 'gray')

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