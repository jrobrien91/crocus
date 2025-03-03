# boundary_layer.py

import numpy as np
import warnings

def create_cloud_mask(backscatter, threshold=1e-6):
    """
    Create a simple cloud mask based on a backscatter threshold.
    
    Parameters
    ----------
    backscatter : np.ndarray
        2D array of backscatter data with dimensions [time, height].
    threshold : float
        Backscatter threshold above which data is considered cloud-contaminated.
        
    Returns
    -------
    cloud_mask : np.ndarray (bool)
        2D boolean mask where True indicates 'cloudy' data (above threshold).
    """
    cloud_mask = backscatter >= threshold
    return cloud_mask

def estimate_bl_height_gradient(
    time, height, backscatter,
    cloud_mask=None,
    z_min=0.0, 
    z_max=3000.0,
    smoothing_window=5
):
    """
    Estimate the boundary-layer height using a simple gradient method, 
    ignoring cloud-contaminated data.
    
    Parameters
    ----------
    time : np.ndarray
        1D array of time (e.g., in hours, or datetime objects).
    height : np.ndarray
        1D array of height [m].
    backscatter : np.ndarray
        2D array of backscatter [time, height].
    cloud_mask : np.ndarray (bool), optional
        2D boolean mask marking cloudy (True) vs. non-cloud (False) data.
        If None, no masking is applied.
    z_min : float
        Minimum height [m] to consider for boundary-layer detection.
    z_max : float
        Maximum height [m] to consider for boundary-layer detection.
    smoothing_window : int
        Optional vertical smoothing window (in # of height points) 
        for the backscatter profile before computing the gradient.
        
    Returns
    -------
    bl_height : np.ndarray
        1D array of boundary-layer height estimates at each time step. 
        If detection fails, returns np.nan for that time.
    """
    
    # Ensure height is sorted ascending
    if not np.all(np.diff(height) > 0):
        raise ValueError("Height array must be in ascending order.")
    
    # Restrict to valid height indices (z_min to z_max)
    valid_z = (height >= z_min) & (height <= z_max)
    h_sub = height[valid_z]
    
    # Prepare output array
    bl_height = np.full(len(time), np.nan)

    for i, t in enumerate(time):
        prof = backscatter[i, valid_z].copy()
        
        # Mask clouds if provided
        if cloud_mask is not None:
            prof_mask = cloud_mask[i, valid_z]
            # If large fraction is masked, might skip
            if np.sum(~prof_mask) < 5:  # less than 5 valid points
                continue
            # Set cloud points to some low backscatter or NaN
            prof[prof_mask] = np.nan
        
        # Optionally remove NaNs by interpolation or skip them
        # Simple approach: skip if too many NaNs
        if np.isnan(prof).mean() > 0.5:
            # More than half the profile is cloud or invalid
            continue
        
        # Fill small gaps with interpolation
        prof = fill_nans_with_interpolation(prof)
        
        # Optional smoothing to reduce noise
        if smoothing_window > 1:
            prof = smooth_1d(prof, window_size=smoothing_window)
        
        # Compute vertical derivative
        dprof_dz = np.gradient(prof, h_sub)
        
        # The boundary-layer top is often near the largest negative derivative
        # We'll find the index of the minimum gradient
        idx_bl = np.argmin(dprof_dz)
        
        # The actual BL top estimate
        bl_height[i] = h_sub[idx_bl]
    
    return bl_height

def fill_nans_with_interpolation(profile):
    """
    Simple 1D linear interpolation of NaN values.
    """
    profile_filled = profile.copy()
    nans = np.isnan(profile_filled)
    if np.any(nans):
        valid_x = np.where(~nans)[0]
        valid_y = profile_filled[~nans]
        full_x = np.arange(len(profile_filled))
        # If entire profile is not NaN
        if len(valid_x) >= 2:
            profile_filled[nans] = np.interp(full_x[nans], full_x[valid_x], valid_y)
    return profile_filled

def smooth_1d(x, window_size=5):
    """
    Simple moving-average smoothing.
    """
    if window_size < 2:
        return x
    conv_kernel = np.ones(window_size) / window_size
    # 'same' mode keeps array length the same
    x_smooth = np.convolve(x, conv_kernel, mode='same')
    return x_smooth

