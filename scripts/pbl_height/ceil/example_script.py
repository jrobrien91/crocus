# example_script.py

import numpy as np
import matplotlib.pyplot as plt
from boundary_layer import create_cloud_mask, estimate_bl_height_gradient

# ---------------------------
# 1) Generate or load data
# ---------------------------
# In reality, you would load from a NetCDF or similar:
time = np.arange(0, 24, 0.5)  # 48 half-hour steps, for example
height = np.linspace(0, 3000, 100)  # 0 to 3 km in 100 steps
# Synthetic backscatter: random + gradient for demonstration
backscatter = np.exp(-(height[None, :] - 500)**2 / (2 * 200**2))  # mock boundary layer
backscatter += 0.05 * np.random.rand(*backscatter.shape)  # noise

# Introduce artificial 'cloud' at 1e-4 above some threshold
backscatter[10:20, 30:40] += 1e-3  # 'cloud' block in time steps [10..19], heights [30..39]

# ---------------------------
# 2) Create a cloud mask
# ---------------------------
cloud_threshold = 1e-4
cloud_mask = create_cloud_mask(backscatter, threshold=cloud_threshold)

# ---------------------------
# 3) Estimate BL height
# ---------------------------
bl_height_est = estimate_bl_height_gradient(
    time=time,
    height=height,
    backscatter=backscatter,
    cloud_mask=cloud_mask,
    z_min=0.0,
    z_max=3000.0,
    smoothing_window=5
)

# ---------------------------
# 4) Plot results
# ---------------------------
fig, ax = plt.subplots(1, 2, figsize=(10, 4))

# (a) Backscatter quicklook
pcm = ax[0].pcolormesh(time, height, np.log10(backscatter.T), shading='auto')
ax[0].set_xlabel('Time [hr]')
ax[0].set_ylabel('Height [m]')
ax[0].set_title('Log10(Backscatter)')
fig.colorbar(pcm, ax=ax[0], label='log10(Backscatter)')

# (b) BL height time series
ax[1].plot(time, bl_height_est, marker='o', linestyle='-')
ax[1].set_xlabel('Time [hr]')
ax[1].set_ylabel('BL Height [m]')
ax[1].set_title('Estimated Boundary Layer Height')
ax[1].grid(True)

plt.tight_layout()
plt.show()

