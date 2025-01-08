# -*- coding: utf-8 -*-
"""
This program processes images to:
1. Remove "dusty" pixels (noise) from individual z-stack images.
2. Create a maximum intensity projection for each z-stack across channels, applying transparency based on pixel intensity.
3. Save the processed z-stack projection as an RGBA TIFF file with transparency, named according to the specified image set.

Requirements:
1. Store all images in a single folder, organized by channel and z-position with a consistent naming convention (e.g., ch00, ch01, ch02, etc.).
2. Run this program from a different directory than the images folder. Set the working directory to the images folder during runtime via user input.
3. Ensure you have the necessary libraries installed: numpy, skimage, and cv2.

Workflow:
1. The program prompts the user to enter the folder location, image set name, and channels to process.
2. For each channel:
    - Images are filtered to remove isolated noisy pixels.
    - A maximum intensity projection is created across z-stack images.
    - Transparency is applied based on intensity (0 is fully transparent, 255 is fully opaque).
3. The final projection is saved as an RGB TIFF file in the working directory, using the specified image set name in the filename.
@author: stock
"""
import re
import numpy as np
import os
from skimage import io
import cv2
import warnings
from tifffile import imsave 
from skimage.morphology import remove_small_objects
import matplotlib.pyplot as plt
from skimage.color import gray2rgb

warnings.filterwarnings("ignore")

def main():
    print("\nThis program will clean images, create a z-stack projection, and export the result as an RGB TIFF.\n")
    
    # User input for folder location
    while True:
        wdir = input('Enter the location of your folder with single-channel planes: ')
        if os.path.isdir(wdir):
            os.chdir(wdir)
            break
        else:
            print("Invalid directory. Please try again.")
    
    # User input for image set names
    image_set_names = input("Enter the names of the image sets (comma separated): ").split(', ')
    
    # User input for channels and colors
    channels = input('Enter the channels to process (e.g., ch00, ch01, ch02, ch03): ').split(', ')
    
    # Input for selecting a range of z-planes
    z_start = int(input("Enter the start z-plane (e.g., 80 for z80): "))
    z_end = int(input("Enter the end z-plane (e.g., 90 for z90): "))
    
    # Get colors for each channel
    channel_colors = {}
    for ch in channels:
        color = input(f"Enter the color for channel {ch} (e.g., green, magenta, red, yellow, cyan): ").lower()
        channel_colors[ch] = color
    
    # Get stain types for each channel (membrane or nuclear)
    stain_types = {}
    for ch in channels:
        stain_type = input(f"Enter the stain type for channel {ch} ('membrane' or 'nuclear'): ").lower()
        while stain_type not in ['membrane', 'nuclear']:
            print("Invalid stain type. Please enter 'membrane' or 'nuclear'.")
            stain_type = input(f"Enter the stain type for channel {ch} ('membrane' or 'nuclear'): ").lower()
        stain_types[ch] = stain_type
    
    # Process each image set individually
    for image_set_name in image_set_names:
        print(f"\nProcessing image set {image_set_name}...\n")
        
        # Gather files for each channel and filter by z-range and image set name
        image_files_per_channel = {}
        for ch in channels:
            image_files_per_channel[ch] = sorted([file for file in os.listdir(wdir)
                                                  if ch in file and image_set_name in file and 
                                                  any(f"z{z}" in file for z in range(z_start, z_end + 1))])
        
        # Initialize a dictionary to store the color projections
        color_projections = {}

        # Process each channel and store the resulting color projection
        for ch, image_files in image_files_per_channel.items():
            if not image_files:
                print(f"No images found for channel {ch} in the range z{z_start}-z{z_end} for {image_set_name}. Skipping.")
                continue
            color_projection = create_z_stack_projection(image_files, image_set_name, ch, channel_colors[ch], stain_types[ch])
            
            if color_projection is not None:
                color_projections[ch] = color_projection

        # Check if any projections were created before merging
        if not color_projections:
            print(f"No valid projections were created for image set {image_set_name}. Skipping...")
            continue

        # Initialize a blank RGB image for the final projection
        first_projection_shape = list(color_projections.values())[0].shape
        final_projection = np.zeros((first_projection_shape[0], first_projection_shape[1], 3), dtype=np.uint8)

        # Assign each channel's projection to the final RGB image based on color
        for ch, projection in color_projections.items():
            # Convert grayscale projection to RGB
            projection_rgb = np.stack([projection] * 3, axis=-1)
            if channel_colors[ch] == 'red':
                final_projection[:, :, 0] = np.maximum(final_projection[:, :, 0], projection_rgb[:, :, 0])
            elif channel_colors[ch] == 'green':
                final_projection[:, :, 1] = np.maximum(final_projection[:, :, 1], projection_rgb[:, :, 1])
            elif channel_colors[ch] == 'blue':
                final_projection[:, :, 2] = np.maximum(final_projection[:, :, 2], projection_rgb[:, :, 2])
            elif channel_colors[ch] == 'yellow':
                final_projection[:, :, 0] = np.maximum(final_projection[:, :, 0], projection_rgb[:, :, 0])
                final_projection[:, :, 1] = np.maximum(final_projection[:, :, 1], projection_rgb[:, :, 1])
            elif channel_colors[ch] == 'cyan':
                final_projection[:, :, 1] = np.maximum(final_projection[:, :, 1], projection_rgb[:, :, 1])
                final_projection[:, :, 2] = np.maximum(final_projection[:, :, 2], projection_rgb[:, :, 2])
            elif channel_colors[ch] == 'magenta':
                final_projection[:, :, 0] = np.maximum(final_projection[:, :, 0], projection_rgb[:, :, 0])
                final_projection[:, :, 2] = np.maximum(final_projection[:, :, 2], projection_rgb[:, :, 2])
            elif channel_colors[ch] == 'white':
                final_projection += projection_rgb

        # Normalize the final projection to [0, 255] range
        final_projection = np.clip(final_projection, 0, 255).astype(np.uint8)

        # Display and save the final merged projection
        plt.figure(figsize=(8, 8))
        plt.imshow(final_projection)
        plt.title(f"Final Merged Projection for {image_set_name}")
        plt.axis('off')
        plt.show()
        
        # Save the merged projection as an RGB TIFF file with a unique name for each image set
        output_file = f"{image_set_name}_merged_projection.tiff"
        imsave(output_file, final_projection)
        print(f"Final merged RGB projection saved as {output_file}")






"""
Cleans the image by following these steps:
1. Normalize the image
2. Apply Gaussian blur
3. Apply intensity threshold to create a binary mask
4. Remove small objects from the mask
5. Erode and dilate the mask
6. Mask the original image with the cleaned mask

:param image: The input image to clean
:param filter_size: The size of the structuring element used in morphological operations
:param percentile: The percentile for the intensity thresholding step
:param sigma: The standard deviation for the Gaussian blur
:param min_size: The minimum size of objects to keep in the binary mask
:return: The cleaned image with noise removed
"""

def clean_image(image, stain_type, filter_size=2, sigma=1, min_size=100):
    """
    Cleans the image by normalizing, smoothing, thresholding, and masking.
    """
    # Normalize the original image to the range 0-255
    normalized_image = cv2.normalize(image, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)

    # Apply Gaussian smoothing to the normalized image
    smoothed_image = cv2.GaussianBlur(normalized_image, (filter_size * 2 + 1, filter_size * 2 + 1), sigma)

    # Set percentile based on stain type
    if stain_type == 'membrane':
        percentile_value = 80  # Lower threshold for membrane stains
    elif stain_type == 'nuclear':
        percentile_value = 98  # Higher threshold for nuclear stains
    else:
        raise ValueError("Invalid stain type. Use 'membrane' or 'nuclear'.")

    # Calculate the automatic threshold based on intensity percentile on the smoothed image
    threshold = np.percentile(smoothed_image, percentile_value)

    # Create a binary mask to isolate bright areas above the threshold
    binary_mask = (smoothed_image > threshold).astype(np.uint8)

    # Remove small objects from the binary mask
    binary_mask = remove_small_objects(binary_mask.astype(bool), min_size=min_size).astype(np.uint8)

    # Structuring element for morphological operations
    kernel = np.ones((filter_size, filter_size), np.uint8)

    # Apply morphological opening (erosion followed by dilation)
    binary_mask = cv2.morphologyEx(binary_mask, cv2.MORPH_OPEN, kernel)

    # Apply morphological closing (dilation followed by erosion) to smooth edges
    binary_mask = cv2.morphologyEx(binary_mask, cv2.MORPH_CLOSE, kernel)

    # Use the cleaned mask to mask the smoothed image
    normalized_image = cv2.normalize(image.astype(np.float32), None, 0, 1, cv2.NORM_MINMAX)
    cleaned_image = normalized_image * binary_mask

    return cleaned_image




"""
Create a maximum intensity projection for a given z-stack of images.

:param image_files: List of image filenames for the z-stack
:param output_name: Output file name (image set)
:param channel: Channel being processed (for naming and color)
:param color: Color to apply to the final projection
:param stain_type: Type of stain ('membrane' or 'nuclear') to determine percentile for thresholding
:return: The max intensity projection image
"""
def create_z_stack_projection(image_files, output_name, channel, color, stain_type):
    # Extract z-positions from the filenames
    z_positions = sorted(set([int(re.search(r'z(\d+)', f).group(1))  # Extracts digits after 'z'
                              for f in image_files if 'z' in f]))  # Ensure 'z' is in the filename

    # Initialize an empty list to store images
    z_stack = []

    # Process images for each z-position
    for z in z_positions:
        z_image_files = [f for f in image_files if f"z{z:02d}" in f]
        if len(z_image_files) != 1:
            print(f"Warning: Expected one image for z{z:02d}, found {len(z_image_files)}. Skipping.")
            continue
        
        # Load the image
        img = io.imread(z_image_files[0])
        
        # Clean the image based on stain type
        img_cleaned = clean_image(img, stain_type=stain_type)

        # Append to the z-stack
        z_stack.append(img_cleaned)

    # Merge all z-planes together (maximum intensity projection)
    if len(z_stack) == 0:
        return None
    
    z_stack = np.array(z_stack)
    max_proj = np.max(z_stack, axis=0)  # Max projection across z-planes

    # If the result is in float format, you can clip values to the range [0, 1] for consistency
    if max_proj.max() > 1.0:  # In case it's in float format
        max_proj = np.clip(max_proj, 0, 1)

    # Ensure the final projection is of type uint8 before saving, if needed
    max_proj_export = (max_proj * 255).astype(np.uint8)  # Convert to uint8 for saving as TIFF
    
    # Plot the final maximum intensity projection
    plt.figure(figsize=(8, 8))
    plt.imshow(max_proj_export)
    plt.title(f"Max Intensity Projection for {channel}")
    plt.axis('off')  # Hide axes for better display
    plt.show()

    # Save the cleaned 3D projection as an RGB TIFF file with high-quality export
    output_file = f"{output_name}_3D_projection_{channel}.tiff"
    save_colored_projection(max_proj, output_file, color)
    return max_proj_export



# COLOR_MAP as defined previously, for convenience.
COLOR_MAP = {
    "green": (0, 255, 0),
    "magenta": (255, 0, 255),
    "cyan": (0, 255, 255),
    "blue": (0, 0, 255),
    "red": (255, 0, 0),
    "yellow": (255, 255, 0),
    "white": (255, 255, 255),
    "black": (0, 0, 0)
}


"""
Save the projection with the specified color as an RGB TIFF.

:param projection: Grayscale projection array.
:param output_name: Base name for the output file.
:param color: Color to apply to the projection ('red', 'green', etc.).
"""
def save_colored_projection(projection, output_name, color):
    try:
        # Ensure the output directory exists
        output_dir = os.path.dirname(output_name)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Normalize the projection to the range [0, 255] for proper display
        projection = np.clip(projection * 255, 0, 255).astype(np.uint8)

        # Convert the grayscale projection to RGB
        projection_rgb = gray2rgb(projection)

        # Apply the color mapping
        if color == "red":
            projection_rgb[:, :, 1:] = 0  # Zero out green and blue
        elif color == "green":
            projection_rgb[:, :, [0, 2]] = 0  # Zero out red and blue
        elif color == "blue":
            projection_rgb[:, :, :2] = 0  # Zero out red and green
        elif color == "yellow":
            projection_rgb[:, :, 2] = 0  # Zero out blue
        elif color == "cyan":
            projection_rgb[:, :, 0] = 0  # Zero out red
        elif color == "magenta":
            projection_rgb[:, :, 1] = 0  # Zero out green
        elif color == "white":
            # Leave all channels as they are
            pass
        else:
            raise ValueError(f"Unsupported color: {color}")

        # Ensure output filename has the proper extension
        output_file = f"{output_name}_{color}.tiff"
        if not output_file.endswith('.tiff'):
            output_file += '.tiff'

        # Save the projection as a colored TIFF
        io.imsave(output_file, projection_rgb)
        print(f"Colored projection saved as: {output_file}")

    except Exception as e:
        print(f"Error saving file {output_name}: {e}")


    
if __name__ == "__main__":
    main()
