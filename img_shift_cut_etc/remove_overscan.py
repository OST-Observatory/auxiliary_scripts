#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
    Automatically detect and remove overscan regions from FITS images
"""

############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

#   Path to the images
image_path: str = '?'

#   Output directory
output_directory: str = '?'

#   Overscan detection parameters
#   Edge detection threshold (fraction of data range)
edge_threshold_fraction: float = 0.1

#   Minimum number of consecutive edges to consider valid
min_edge_group_size: int = 1

#   Maximum gap between edges in a group
max_edge_gap: int = 2

#   Manual overscan bounds (set to None for automatic detection)
#   Format: (start, end) for each axis
overscan_row: tuple[int, int] | None = None
overscan_column: tuple[int, int] | None = None

#   Verbose output
verbose: bool = True

############################################################################
####                            Libraries                               ####
############################################################################

import sys
import numpy as np
from pathlib import Path
import ccdproc as ccdp
from astropy.io import fits
from ost_photometry import checks

############################################################################
####                         Overscan Detection                        ####
############################################################################

def detect_overscan_region(image_data, axis=1):
    """
    Detect overscan region along specified axis using edge detection
    
    Parameters:
    -----------
    image_data : numpy.ndarray
        Image data
    axis : int
        Axis along which to detect overscan (0=rows, 1=columns)
    
    Returns:
    --------
    tuple : (start, end) indices of overscan region, or None if not detected
    """
    if axis == 0:
        # Check rows (top/bottom)
        data = image_data
        size = data.shape[0]
    else:
        # Check columns (left/right)
        data = image_data.T
        size = data.shape[0]
    
    # Calculate mean values for each row/column
    means = np.mean(data, axis=1)
    
    # Calculate gradient (difference between adjacent rows/columns)
    gradient = np.diff(means)
    
    # Find significant edges (where gradient exceeds threshold)
    # Use a threshold based on the overall range of the data
    data_range = np.max(means) - np.min(means)
    edge_threshold = data_range * edge_threshold_fraction
    
    # Find edges that exceed threshold
    significant_edges = np.where(np.abs(gradient) > edge_threshold)[0]
    
    if len(significant_edges) == 0:
        return None
    
    # Look for consistent edges across multiple rows/columns
    # Group nearby edges together
    edge_groups = []
    current_group = [significant_edges[0]]
    
    for edge in significant_edges[1:]:
        if edge - current_group[-1] <= max_edge_gap:  # Allow small gaps
            current_group.append(edge)
        else:
            edge_groups.append(current_group)
            current_group = [edge]
    
    if current_group:
        edge_groups.append(current_group)
    
    # Find the most consistent edge (largest group)
    if not edge_groups:
        return None
    
    # Filter groups by minimum size
    valid_groups = [group for group in edge_groups if len(group) >= min_edge_group_size]
    
    if not valid_groups:
        return None
    
    largest_group = max(valid_groups, key=len)
    edge_position = int(np.median(largest_group))
    
    # Determine if this is the start or end of overscan
    # Check the gradient direction and position
    avg_gradient = np.mean(gradient[largest_group])
    
    if axis == 0:  # Rows (top/bottom)
        if avg_gradient < 0:  # Negative gradient means counts decrease
            # This is likely the start of overscan (end of real data)
            if edge_position < size // 2:
                # Edge in first half, overscan is at the beginning
                return (0, edge_position + 1)
            else:
                # Edge in second half, overscan is at the end
                return (edge_position + 1, size)
        else:  # Positive gradient means counts increase
            # This is likely the end of overscan (start of real data)
            if edge_position < size // 2:
                # Edge in first half, overscan is at the end
                return (edge_position + 1, size)
            else:
                # Edge in second half, overscan is at the beginning
                return (0, edge_position + 1)
    else:  # Columns (left/right)
        if avg_gradient < 0:  # Negative gradient means counts decrease
            # This is likely the start of overscan (end of real data)
            if edge_position < size // 2:
                # Edge in first half, overscan is at the end
                return (edge_position + 1, size)
            else:
                # Edge in second half, overscan is at the beginning
                return (0, edge_position + 1)
        else:  # Positive gradient means counts increase
            # This is likely the end of overscan (start of real data)
            if edge_position < size // 2:
                # Edge in first half, overscan is at the beginning
                return (0, edge_position + 1)
            else:
                # Edge in second half, overscan is at the end
                return (edge_position + 1, size)

def detect_overscan_for_image(image_data):
    """
    Detect overscan regions for a single image
    
    Parameters:
    -----------
    image_data : numpy.ndarray
        Image data
    
    Returns:
    --------
    tuple : (row_overscan, col_overscan) or (None, None) if no overscan detected
    """
    # Detect overscan regions
    row_overscan = detect_overscan_region(image_data, axis=0)
    col_overscan = detect_overscan_region(image_data, axis=1)
    
    return row_overscan, col_overscan

def remove_overscan_with_bounds(image_data, row_bounds, col_bounds):
    """
    Remove overscan regions from image data using predefined bounds
    
    Parameters:
    -----------
    image_data : numpy.ndarray
        Image data
    row_bounds : tuple or None
        (start, end) for row overscan removal
    col_bounds : tuple or None
        (start, end) for column overscan removal
    
    Returns:
    --------
    numpy.ndarray : Image data with overscan removed
    """
    original_shape = image_data.shape
    result = image_data.copy()
    
    if verbose:
        print(f"  Original shape: {original_shape}")
        if row_bounds:
            print(f"  Removing row overscan: {row_bounds}")
        if col_bounds:
            print(f"  Removing column overscan: {col_bounds}")
    
    # Remove overscan regions
    if row_bounds:
        start, end = row_bounds
        # Check if it's at the beginning or end
        if start == 0:
            # Remove from beginning
            result = result[end:, :]
        else:
            # Remove from end
            result = result[:start, :]
    
    if col_bounds:
        start, end = col_bounds
        # Check if it's at the beginning or end
        if start == 0:
            # Remove from beginning
            result = result[:, end:]
        else:
            # Remove from end
            result = result[:, :start]
    
    if verbose:
        print(f"  Final shape: {result.shape}")
    
    return result

############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    #   Check input and output directory
    file_path = checks.check_pathlib_path(image_path)
    checks.check_output_directories(output_directory)
    out_path = Path(output_directory)
    
    #   Create output directory if it doesn't exist
    out_path.mkdir(parents=True, exist_ok=True)

    #   Get image file collection
    ifc = ccdp.ImageFileCollection(file_path)
    total_files = len(ifc.files)
    
    #   Determine overscan bounds (manual or automatic)
    consistent_row_bounds = None
    consistent_col_bounds = None
    
    #   Check if manual bounds are provided
    if overscan_row is not None or overscan_column is not None:
        if verbose:
            print("Using manual overscan bounds:")
            if overscan_row:
                print(f"  Row overscan: {overscan_row}")
            if overscan_column:
                print(f"  Column overscan: {overscan_column}")
        
        consistent_row_bounds = overscan_row
        consistent_col_bounds = overscan_column
        
    else:
        #   Automatic detection mode
        if verbose:
            print(f"Analyzing {total_files} images to determine consistent overscan bounds...")
        
        #   First pass: analyze all images to determine consistent overscan bounds
        all_row_bounds = []
        all_col_bounds = []
        
        for img, file_name in ifc.ccds(
            ccd_kwargs={'unit': 'adu'},
            return_fname=True,
            ):
            
            if verbose:
                print(f"Analyzing {file_name}")
            
            #   Detect overscan for this image
            row_bounds, col_bounds = detect_overscan_for_image(img.data)
            
            if row_bounds:
                all_row_bounds.append(row_bounds)
            if col_bounds:
                all_col_bounds.append(col_bounds)
        
        #   Determine consistent bounds (use maximum range)
        if all_row_bounds:
            # Find the maximum range for row overscan
            max_row_start = max(bounds[0] for bounds in all_row_bounds)
            max_row_end = max(bounds[1] for bounds in all_row_bounds)
            consistent_row_bounds = (max_row_start, max_row_end)
            
            if verbose:
                print(f"Consistent row overscan bounds: {consistent_row_bounds}")
        
        if all_col_bounds:
            # Find the maximum range for column overscan
            max_col_start = max(bounds[0] for bounds in all_col_bounds)
            max_col_end = max(bounds[1] for bounds in all_col_bounds)
            consistent_col_bounds = (max_col_start, max_col_end)
            
            if verbose:
                print(f"Consistent column overscan bounds: {consistent_col_bounds}")
    
    if not consistent_row_bounds and not consistent_col_bounds:
        print("No overscan regions detected or specified!")
        sys.exit(1)
    
    #   Second pass: remove overscan using consistent bounds
    if verbose:
        print(f"\nRemoving overscan from all images using consistent bounds...")
    
    i = 0
    for img, file_name in ifc.ccds(
        ccd_kwargs={'unit': 'adu'},
        return_fname=True,
        ):

        if verbose:
            print(f"Processing {file_name} ({i+1}/{total_files})")
        
        #   Remove overscan using consistent bounds
        img_data = img.data
        cleaned_data = remove_overscan_with_bounds(img_data, consistent_row_bounds, consistent_col_bounds)
        
        #   Create new CCD object with cleaned data
        img.data = cleaned_data
        
        #   Update header information
        img.meta['NAXIS1'] = cleaned_data.shape[1]
        img.meta['NAXIS2'] = cleaned_data.shape[0]
        
        #   Determine detection mode
        if overscan_row is not None or overscan_column is not None:
            img.meta['INFO_0'] = 'Overscan region manually specified and removed'
        else:
            img.meta['INFO_0'] = 'Overscan region automatically detected and removed'
            
        img.meta['INFO_1'] = f'Original size: {img_data.shape[1]}x{img_data.shape[0]}'
        img.meta['INFO_2'] = f'Final size: {cleaned_data.shape[1]}x{cleaned_data.shape[0]}'
        if consistent_row_bounds:
            img.meta['INFO_3'] = f'Row overscan bounds: {consistent_row_bounds}'
        if consistent_col_bounds:
            img.meta['INFO_4'] = f'Column overscan bounds: {consistent_col_bounds}'

        #   Save the result
        img.write(out_path / file_name, overwrite=True)

        #   Write status to console
        i += 1
        if not verbose:
            sys.stdout.write(f"\rProcessed {i}/{total_files} images")
            sys.stdout.flush()
    
    if not verbose:
        sys.stdout.write("\n")
    else:
        print(f"Completed! Processed {i} images.")
        print(f"All images now have consistent dimensions: {cleaned_data.shape[1]}x{cleaned_data.shape[0]}") 