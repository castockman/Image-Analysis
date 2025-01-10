# -*- coding: utf-8 -*-
"""
This program analyzes z-stack images to:
1. Calculate the ratio of pixel counts between two selected binary masks (numerator and denominator) for specified image channels.
2. Process z-stack images into binary masks using image filtering techniques to enhance quality and reduce noise.
3. Compute pixel counts for each binary mask on a per-image basis and store results for further analysis.
4. Export results to an Excel file, including channel ratios and detailed pixel counts for each image.

Requirements:
1. Input images should be stored in a single directory, organized with consistent naming conventions (e.g., genotype_ch00, genotype_ch01, etc.).
2. The program must be run from a directory other than the image folder. Set the working directory to the images folder using user input.
3. Ensure the following Python libraries are installed: numpy, scipy, cv2, pandas, and matplotlib.

Workflow:
1. The user specifies the folder location, genotypes, and channels to analyze, along with the antibody information for the selected channels.
2. For each genotype and z-stack:
   - Images are processed into binary masks using Gaussian blurring, adaptive histogram equalization, and Otsu thresholding.
   - Pixel counts for each binary mask are calculated independently.
   - The ratio of pixel counts between the specified channels is computed for each image.
3. Results, including pixel counts and ratios, are summarized and exported to an Excel file.

Image Processing Functions:
1. `filter_mask`: Processes raw image planes into binary masks using thresholding and morphological operations.
2. `pixel_counter`: Calculates the total pixel count for a given binary mask.
3. `display_every_10`: Provides a visual quality check by displaying masks at specified z-stack intervals.

Data Output:
1. Pixel count ratios between the selected channels, calculated on a per-image basis.
2. Total pixel counts for each channel and image, organized in a detailed Excel file.
@author: Courtney
"""

### Updated Code
# Import Packages:
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
import cv2
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

def main():
    global channel_ratios, files_ch0, files_ch1, image_sets
    print('\n\n\nThis program will calculate the ratio of x:y for all pixels in x and y. \n')

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
    channel_ratios = {}
    total_pixels_image = {}

    ratio = input(f'What AB will represent TOTAL: {ch0}, {ch1}: ').strip().upper()
    if ratio not in [ch0, ch1]:
        print("Error: Ratio must be one of the selected antibodies.")
        return

    num, dom = define_method(ch0, ch1) if ratio == ch0 else define_method(ch1, ch0)

    # Counting Cells
    print(f'\n\n\n Calculating {num}:{dom}....\n\n')

    for genotype1 in genotype:
        print(genotype1.upper())
        files_genotype = [FILE for FILE in files if genotype1 in FILE]
        image_sets = unique_images(files_genotype)
        channel_ratio_per_image = []

        for image in image_sets:
            total_pixels_image[image] = []  # Initialize as a list to store counts for each file
            files_ch0 = list(FILE for FILE in files_genotype if f"{image}_{channel[0]}" in FILE)
            files_ch1 = list(FILE for FILE in files_genotype if f"{image}_{channel[1]}" in FILE)

            # Ensure files_ch0 and files_ch1 are matched
            if len(files_ch0) != len(files_ch1):
                print(f"Warning: Mismatch in file count for {image}. Skipping...")
                continue

            # Process Images to Binary
            for file_ch0, file_ch1 in zip(files_ch0, files_ch1):
                mask_ch02 = filter_mask(file_ch0)
                mask_ch01 = filter_mask(file_ch1)

                # Generate Pixel Count
                try:
                    ch02_count = pixel_counter(mask_ch02)
                    ch01_count = pixel_counter(mask_ch01)

                    total_pixels_image[image].append({
                        'File': file_ch0,
                        'Channel 0': ch02_count,
                        'Channel 1': ch01_count
                    })

                    if ch02_count == 0:  # Avoid division by zero
                        channel_ratio_per_image.append(np.nan)
                    else:
                        channel_ratio_per_image.append(ch01_count / ch02_count)

                except Exception as e:
                    print(f"Error processing {file_ch0} and {file_ch1}: {e}")
                    continue

                if question == 'Y':
                    display_every_10(file_ch0, mask_ch01, mask_ch02, num, dom)

            channel_ratios[genotype1] = channel_ratio_per_image

    # Convert to DataFrames
    channel_ratios_df = pd.DataFrame({k: pd.Series(v) for k, v in channel_ratios.items()})
    total_pixels_df = pd.DataFrame(
        [
            {"Image": img, **counts} for img, counts_list in total_pixels_image.items() for counts in counts_list
        ]
    )

    # Export to Excel
    print("\n\nExporting to Excel....")
    output_file = 'Percent_of_Total.xlsx'
    with pd.ExcelWriter(output_file) as writer:
        channel_ratios_df.to_excel(writer, sheet_name=f"{num}-{dom}")
        total_pixels_df.to_excel(writer, sheet_name="Total Pixels", index=False)
    print(f"Export completed: {output_file}")

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

def pixel_counter(mask):
    return np.sum(mask * 1)

def unique_images(files_genotype):
    image_sets = []
    for item in range(len(files_genotype)):
        image_sets = np.append(image_sets, str(files_genotype[item].split('_ch')[0]))
        image_sets = list(set(image_sets))
    return image_sets

def define_method(num, dom):
    return num, dom

def display_every_10(FILE, mask_ch01, mask_ch02, num, dom):
    if 'z10' in FILE:
        display_masks(mask_ch01, str(num), 'gray')
        display_masks(mask_ch02, str(dom), 'gray')
    if 'z20' in FILE:
        display_masks(mask_ch01, str(num), 'gray')
        display_masks(mask_ch02, str(dom), 'gray')
    if 'z30' in FILE:
        display_masks(mask_ch01, str(num), 'gray')
        display_masks(mask_ch02, str(dom), 'gray')
    if 'z40' in FILE:
        display_masks(mask_ch01, str(num), 'gray')
        display_masks(mask_ch02, str(dom), 'gray')
    if 'z50' in FILE:
        display_masks(mask_ch01, str(num), 'gray')
        display_masks(mask_ch02, str(dom), 'gray')
    if 'z60' in FILE:
        display_masks(mask_ch01, str(num), 'gray')
        display_masks(mask_ch02, str(dom), 'gray')

def display_masks(nuclei, ch00, scale):
    nuclei = np.float32(nuclei)
    (n_row, n_col) = nuclei.shape
    _, (a) = plt.subplots(ncols=1, figsize=(15, 5))
    show_plane(a, nuclei, cmap=scale, title=str(ch00))
    nuclei = 0

def show_plane(ax, plane, cmap, title=None):
    ax.imshow(plane, cmap=cmap)
    ax.axis("off")
    if title:
        ax.set_title(title)

if __name__ == '__main__':
    main()
