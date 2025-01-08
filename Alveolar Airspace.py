# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 13:15:31 2022

This is a program to:
    1. Measure Alveolar features
        - Alveolar Area
        - Alveolar Diameter
        - Gas Exchange Surface Area
        - Volume:Area
        - Mean Linear Intercept
        - Chord Length
        - Septal Width
Measurements for each image will be exported to excel, split up by the treatment groups you have fed into the program

Requirements:
    1. Have all images in the same folder exported as single channels with the label ch00, ch01, ch02 etc
    2. Have this program saved in a folder. You will set a working directory in the code so it does not need to be in the same folder.
    3. Images have a resolution of 1024x1024, if not:
        conversion will be off by a factor of 2 for lengths (Alveolar Diameter, Gas Exchange Surface Area, Chord Length, Septal Width) 
        or 4 for area (Alveolar Area). MLI, Volume:Area not affected. 
@author: stock
"""

#pip install numpy, os, matplotlib, scipy, cv2, imageio, skimage, pandas, xtiff

#Import Packages:
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
import cv2
import pandas as pd
from scipy import ndimage
from skimage.measure import regionprops, label, perimeter
from skimage.segmentation import clear_border
from skimage import measure
from PIL import Image, ImageDraw
import warnings
warnings.filterwarnings("ignore")

def main(): 
    global genotype, alveolar_area, alveolus_diameter, gas_exchange_SA, volume_area, mean_linear_intercept_count, chord_length, septal_width
    print('\n\n\nThis program will Measure Alveolar features: \nAlveolar Area, Alveolar Diameter, Gas Exchange Surface Area, \nVolume:Area, Mean Linear Intercept, Chord Length, Septal Thickness \n')

    # User Input Questions 
    wdir=str(input('What is the location of your folder with single channel planes: '))
    os.chdir(wdir)
    for root, dirs, files in os.walk("."):
        files=files
        files=[FILE for FILE in files if 'ch0' in FILE]   
    question=str(input('Do you want to see the masks that are being Counted? (Y/N): ')).upper()
    genotype=str(input('What Genotypes are you Analyzing? (Comma Separated List ", "): ')).split(', ')
    channel=str(input('What channel are you analyzing (ch00, ch01, ch02, ch03): '))
    if question=='Y':
        marker=str(input('What marker did you use to stain the Alveolar Surface: ')).upper()
    magnification =int(input('What magnification are your images (10x, 20x, 25x): ').replace('x',''))
    res=int(input('What is the resolution of your images (1024, 2048): '))
    cf,sq=pixel_convert(magnification, res)
    # Sorting Files from User Input
    for FILE in files:
        files_ch0=[FILE for FILE in files if str(channel) in FILE] # This will only count images with ch00
      
    # Initialize Dictionaries for Excel Export    
    alveolar_area=dict.fromkeys(genotype,[])
    alveolus_diameter=dict.fromkeys(genotype,[])
    gas_exchange_SA=dict.fromkeys(genotype,[])
    volume_area=dict.fromkeys(genotype,[])
    mean_linear_intercept_count=dict.fromkeys(genotype,[])
    chord_length=dict.fromkeys(genotype,[])
    tissue_area=dict.fromkeys(genotype,[])
    septal_width=dict.fromkeys(genotype,[])
    
    # Measuring Alveolar Features
    print('\n\n\nMeasuring Alveoli....\n')        
    for i in range(len(genotype)):
        print(genotype[i].upper())
        genotype1=str(genotype[i])
        for FILE in files_ch0:
            files_ch00=[FILE for FILE in files_ch0 if str(genotype[i]) in FILE] # This will only count images with ch00  
        count=0       
        for FILE in files_ch00:     
            # Make binary mask for analysis
            mask_ch00=filter_mask(FILE,cf)
            mask_alveoli=full_alveoli(mask_ch00,cf)
            label_im, nb_labels = ndimage.label(mask_alveoli)
            if nb_labels<round((20/magnification)*12):
                continue
            if np.sum(mask_ch00*1)<.6*(res^2):
                continue
            # Display Masks:
            if question == 'Y':
                display_every_10(FILE, mask_ch00,mask_alveoli, marker)
                    
            # Measure Lung Features:             
            if count==0:
                gas_exchange_surface_area=measure_gas_exchange_surface(mask_ch00,cf)
                sizes=measure_alveolar_area(mask_alveoli,cf)
                d=measure_alveolar_diameter(mask_alveoli,cf)
                volume_area_ratio=va_ratio(mask_ch00, ~mask_ch00)
                mli=mean_linear_intercept(mask_ch00,sq)
                tissuea=sum(measure_alveolar_area(mask_ch00,cf))
                chordlength=[chordlen(mask_alveoli,cf,sq)]
                septalwidth=septalwid(mask_ch00,cf,sq)
                count+=1
            else:
                gas_exchange_surface_area=np.append(gas_exchange_surface_area,int(measure_gas_exchange_surface(mask_ch00,cf)))
                sizes.extend(measure_alveolar_area(mask_alveoli,cf))
                d.extend(measure_alveolar_diameter(mask_alveoli,cf))
                volume_area_ratio=np.append(volume_area_ratio,va_ratio(mask_ch00, ~mask_ch00))
                mli=np.append(mli, mean_linear_intercept(mask_ch00,sq))
                tissuea=np.append(tissuea, int(sum(measure_alveolar_area(mask_ch00,cf))))
                chordlength.append(chordlen(mask_alveoli,cf,sq))
                septalwidth.extend(septalwid(mask_ch00,cf,sq))
                
        gas_exchange_SA[genotype1]=gas_exchange_surface_area
        volume_area[genotype1]=volume_area_ratio
        alveolar_area[genotype1]=sizes
        alveolus_diameter[genotype1]=d
        mean_linear_intercept_count[genotype1]=mli
        chord_length[genotype1]=chordlength
        tissue_area[genotype1]=tissuea
        septal_width[genotype1]=septalwidth


    #Convert to DataFrame
    alveolar_area=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in alveolar_area.items()]))
    alveolus_diameter=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in alveolus_diameter.items()]))
    gas_exchange_SA=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in gas_exchange_SA.items()]))
    volume_area=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in volume_area.items()]))
    mean_linear_intercept_count=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in mean_linear_intercept_count.items()]))
    chord_length=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in chord_length.items()]))
    tissue_area=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in tissue_area.items()]))
    septal_width=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in septal_width.items()]))
    
    # Export to Excel:  
    print("\n\nExporting to Excel....")
    with pd.ExcelWriter('Lung_Measurements.xlsx') as writer:
        tissue_area.to_excel(writer, sheet_name='Tissue Area')
        mean_linear_intercept_count.to_excel(writer, sheet_name='Mean Linear Intercept')
        chord_length.to_excel(writer, sheet_name='Chord Lengths')
        gas_exchange_SA.to_excel(writer, sheet_name='Gas Exchange Surface Area')
        volume_area.to_excel(writer, sheet_name='Volume-Area Ratio')
        alveolar_area.to_excel(writer, sheet_name='Alveolar Area')
        alveolus_diameter.to_excel(writer, sheet_name='Alveoalr Diameter')
        septal_width.to_excel(writer, sheet_name='Septal Widths')
        
    # Instructions to User:
    print('\n\nThe following were processed: '+str(genotype))
    print('\nEverything has been exported to the working directory: '+str(wdir))
    if question =='Y':
        print('Preparing to Display Pixel Masks....')
        print('   ')
        print('\nThe masks created are displayed in the "PLOTS" tab above')

""""
Functions to process images imput into program:
    Filter Mask - Make an mask of alveolar surface with no gaps in walls and 
                    alveolar pores filled in
    Remove Small Objects - Remove any specks that could hinder accuracy
    Remove Large Objects - Remove any objects that are too large (duct, airway, artery, etc.)
    Invert Mask - invert mask from alveolar tissue to airspace
    Remove Border Objects - Remove any objects touching the border
    Full Alveoli - Return only full alveoli that are not touching the border of 
                    the image
    
Input: image or mask
Output: final mask of alveolar surface or alveolar airspaces to analyze

Use: These functions are used to generate images for analysis to measure 
    different lung features
"""

def filter_mask(image,cf):
    #LOAD IMAGE
    im=cv2.imread(image,cv2.IMREAD_GRAYSCALE) #"TileScan 003_TileScan_001_Merging_z00_ch00.tif"
    im = im.astype('float64')
    im *= (255/im.mean())
    im = ndi.gaussian_filter(im, 1)
    im *= (255/im.mean())
    mask = im > 255*.80 #WILL CAPTURE POSITIVE PIXELS IN THE TOP 50% OF POSITIVE PIXELS, REMOVES BACKGROUND WITHOUT BEING STRINGENT 
    
    #EXTRAPOLATE TO FIND NEARBY POSITIVE PIXELS    
    mask = ndi.binary_opening(mask)
    mask=ndimage.binary_dilation(mask).astype(mask.dtype)
    mask=ndimage.binary_dilation(mask).astype(mask.dtype)
    mask=ndimage.binary_dilation(mask).astype(mask.dtype)



    #INTERPOLATE TO APROXIMATE CORRECT SIZE OF POSITIVE ZONE
    mask = ndi.binary_closing(mask)
    mask=ndimage.binary_erosion(mask).astype(mask.dtype)
    mask=ndimage.binary_erosion(mask).astype(mask.dtype)
    mask=ndimage.binary_erosion(mask).astype(mask.dtype)
    
    #REMOVE OBJECTS
    mask = remove_large_objects(mask,cf)
    mask = remove_small_objects(mask,cf)
    return mask

def remove_small_objects(mask,cf):
    label_im, nb_labels = ndimage.label(mask)
    #Compute size, mean_value etc of each region
    sizes = ndimage.sum(mask, label_im, range(nb_labels + 1))
    #Clean up small connect components:
    mask_size = sizes < 283.8134765625*cf
    remove_pixel = mask_size[label_im]
    remove_pixel.shape  
    label_im[remove_pixel] = 0
    mask = label_im > label_im.min()
    return mask

def remove_large_objects(mask,cf):
    label_im, nb_labels = ndimage.label(mask)
    #Compute size, mean_value etc of each region
    sizes = ndimage.sum(mask, label_im, range(nb_labels + 1))
    #Clean up small connect components:
    mask_size = sizes > 454101.5625*cf
    remove_pixel = mask_size[label_im]
    remove_pixel.shape  
    label_im[remove_pixel] = 0
    mask = label_im > label_im.min() 
    return mask
    
def invert_mask(im):
    return np.where((im==0)|(im==1), im^1, im)

def remove_border_objects(mask):
    mask=label(mask)
    mask=clear_border(mask, buffer_size=3) 
    return mask>0

def full_alveoli(mask,cf):
    mask_alveoli= invert_mask(mask)*1
    mask_alveoli=remove_small_objects(mask_alveoli,cf)
    return remove_border_objects(mask_alveoli)

""""
Functions to measure properties of distal lung structure:
    Alveolar Area
    Alveolar Diameter
    Gas Exchange Surface Area (Perimeter of alveolar airspace)
    Vollume:Area Ratio
    Mean Linear Intercept - Number of times a grid intersects the alveolar tissue
    Chord Length/Septal Width - lenght of lines that intersect the input mask
    
Input: Mask of alveolar tissue, or airspace
Output: Calculation of lung feature

Use: These functions are used to calculate different features of the lung to 
    quantify differences between the lung structure of samples. 
"""

def measure_alveolar_area(mask,cf):
    label_im, nb_labels = ndimage.label(mask)
    sizes = ndimage.sum(mask, label_im, range(nb_labels + 1))/(cf*cf)
    return list(sizes)
    
def measure_alveolar_diameter(mask,cf):
    properties=regionprops(mask.astype(np.uint8))
    d=[1]
    for p in properties:
        min_row, min_col, max_row, max_col = p.bbox
        d.append(max(max_row - min_row, max_col - min_col)/cf)
    d.pop(0)
    return list(d)

def measure_gas_exchange_surface(mask,cf):
    return perimeter(mask)/cf
   
def va_ratio(mask1, mask2):
    return np.sum(mask2*1)/np.sum(mask1*1) # Mask 1 = Area, Mask 2 = Volume

def mean_linear_intercept(im,sq):
    grid=draw_grid(im*1, sq)
    overlap=grid*im
    label_im, nb_labels = ndimage.label(overlap)
    return nb_labels

def chordlen(im,cf,sq):
    alveolus=invert_mask(im)
    #Vertical Chords
    vert=draw_vertical_lines(im*1,sq)
    chordsv=vert*alveolus
    vert=remove_border_objects(chordsv)
    label_im, nb_labels = ndimage.label(vert)
    sizesv = ndimage.sum(vert, label_im, range(nb_labels + 1))/cf
    
    #Horizontal Chords
    horiz=draw_horizontal_lines(im*1,sq)
    chordsh=horiz*alveolus
    horiz=remove_border_objects(chordsh)
    label_im, nb_labels = ndimage.label(horiz)
    sizesh = ndimage.sum(horiz, label_im, range(nb_labels + 1))/cf
    sizesv=list(sizesv)+list(sizesh)
    avg_chord=sum(sizesv)/len(sizesv)
    return str(avg_chord)

def septalwid(im,cf,sq):
    alveolus1=im
    #Vertical Chords
    grid1=draw_grid(im*1,sq)
    septal_cubes=alveolus1*grid1
    labels = measure.label(septal_cubes)
    props = measure.regionprops(labels)
    d=[1]
    for p in props:
        if p.axis_minor_length < (sq):
            d.append(p.axis_minor_length/cf)
    d.pop(0)
    return list(d)

"""
Functions to generate horizontal, vertical lines or a grid over an image 

Input: image of desired grid size
Output: horizontal, vertical or grid lines every 50 pixels of image

Use: These functions are used to calculate the chord length and the number of times 
    a grid intersects a mask of the alveolar tissue
"""

def draw_vertical_lines(img,sq):
    h, w= img.shape
    image = Image.new(mode='L', size=(h, w), color=255)
    # Draw some lines
    draw = ImageDraw.Draw(image)
    y_start = 0
    y_end = image.height
    step_size = int(image.width/sq*2)//1
    for x in range(0, image.width, step_size):
        line = ((x, y_start), (x, y_end))
        draw.line(line)
    return (np.array(image)>1)*1

def draw_horizontal_lines(img,sq):
    h, w= img.shape
    image = Image.new(mode='L', size=(h, w), color=255)
    # Draw some lines
    draw = ImageDraw.Draw(image)
    x_start = 0
    x_end = image.width
    step_size = int(image.width/sq*2)//1
    for y in range(0, image.height, step_size):
        line = ((x_start, y), (x_end, y))
        draw.line(line)
    return (np.array(image)>1)*1

def draw_grid(img,sq):
    h, w= img.shape
    image = Image.new(mode='L', size=(h, w), color=255)
    # Draw some lines
    draw = ImageDraw.Draw(image)
    y_start = 0
    y_end = image.height
    step_size = int(image.width/sq*2)//1
    for x in range(0, image.width, step_size):
        line = ((x, y_start), (x, y_end))
        draw.line(line)
    x_start = 0
    x_end = image.width

    for y in range(0, image.height, step_size):
        line = ((x_start, y), (x_end, y))
        draw.line(line)
    return (np.array(image)>1)*1


"""
Other Functions for Display:
    Functions to display masks analysis completed
    Pixel to um Conversions for all magnifications on Desai Leica Confocal

Input: Mask
Output: Mask displayed in Plots Tab
"""
def pixel_convert(magnification, res):
    factor=res/1024
    if magnification==10:
        if res==1024:
            dim=1162.5
        elif res==2048:
            dim=1158.73
        cf=res/dim
        sq=32/factor
    elif magnification== 20:
        if res==1024:
            dim=581.25
        elif res==2048:
            dim=579.36
        cf=res/dim
        sq=64/factor
    elif magnification== 25: 
        if res==1024:
            dim=465
        cf=res/dim
        sq=80/factor
    elif magnification== 40:
        if res==1024:
            dim=290.62
        cf=res/dim
        sq=128/factor
    elif magnification== 63:
        if res==1024:
            dim=184.52
        cf=res/dim
        sq=201.6/factor
    return int(cf),int(sq)
    
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

def display_every_10(FILE, mask_ch00,mask_alveoli, marker):
    if 'z10' in FILE:
        display_masks(mask_ch00, str(marker), 'gray')
        display_masks(mask_alveoli, 'Airspace', 'gray')
        
    if 'z20' in FILE:
        display_masks(mask_ch00, str(marker), 'gray')
        display_masks(mask_alveoli, 'Airspace', 'gray')
            
    if 'z30' in FILE:
        display_masks(mask_ch00, str(marker), 'gray')
        display_masks(mask_alveoli, 'Airspace', 'gray')
    if 'z40' in FILE:
        display_masks(mask_ch00, str(marker), 'gray')
        display_masks(mask_alveoli, 'Airspace', 'gray')
    if 'z50' in FILE:
        display_masks(mask_ch00, str(marker), 'gray')
        display_masks(mask_alveoli, 'Airspace', 'gray')
    if 'z60' in FILE:
            display_masks(mask_ch00, str(marker), 'gray')
            display_masks(mask_alveoli, 'Airspace', 'gray')
    
# This provided line is a required at the end of a Python file to call main()
if __name__ =='__main__':
    main()