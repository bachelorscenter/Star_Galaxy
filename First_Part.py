from astropy.io import fits
import numpy as np
from astropy.wcs import WCS

# Define the file paths for each spectral band
file_irg = "/Users/batuhank1/PycharmProjects/ArchOfML/frame-irg-003918-3-0213.jpg"
file_r = "/Users/batuhank1/PycharmProjects/ArchOfML/frame-r-003918-3-0213.fits"
file_g = "/Users/batuhank1/PycharmProjects/ArchOfML/frame-g-003918-3-0213.fits"
file_i = "/Users/batuhank1/PycharmProjects/ArchOfML/frame-i-003918-3-0213.fits"
file_u = "/Users/batuhank1/PycharmProjects/ArchOfML/frame-u-003918-3-0213.fits"
file_z = "/Users/batuhank1/PycharmProjects/ArchOfML/frame-z-003918-3-0213.fits"

# Load the R band FITS file using astropy.io.fits
hdulist_r = fits.open(file_r)
data_r = hdulist_r[0].data
hdulist_r.close()

# Load the G band FITS file using astropy.io.fits
hdulist_g = fits.open(file_g)
data_g = hdulist_g[0].data
hdulist_g.close()

# Load the I band FITS file using astropy.io.fits
hdulist_i = fits.open(file_i)
data_i = hdulist_i[0].data
hdulist_i.close()

# Load the U band FITS file using astropy.io.fits
hdulist_u = fits.open(file_u)
data_u = hdulist_u[0].data
hdulist_u.close()

# Load the Z band FITS file using astropy.io.fits
hdulist_z = fits.open(file_z)
data_z = hdulist_z[0].data
hdulist_z.close()

# Calculate the offsets using WCS
wcs_r = WCS(hdulist_r[0].header)

# Convert pixel coordinates to sky coordinates for the reference band (R band)
coords_r = wcs_r.all_pix2world(np.column_stack(np.indices(data_r.shape)).reshape(-1, 2), 0)

# Convert pixel coordinates to sky coordinates for the other bands
coords_g = wcs_r.all_pix2world(np.column_stack(np.indices(data_g.shape)).reshape(-1, 2), 0)
coords_i = wcs_r.all_pix2world(np.column_stack(np.indices(data_i.shape)).reshape(-1, 2), 0)
coords_u = wcs_r.all_pix2world(np.column_stack(np.indices(data_u.shape)).reshape(-1, 2), 0)
coords_z = wcs_r.all_pix2world(np.column_stack(np.indices(data_z.shape)).reshape(-1, 2), 0)

# Calculate the offsets
x_offset_g, y_offset_g = np.split(coords_g - coords_r, 2, axis=1)
x_offset_i, y_offset_i = np.split(coords_i - coords_r, 2, axis=1)
x_offset_u, y_offset_u = np.split(coords_u - coords_r, 2, axis=1)
x_offset_z, y_offset_z = np.split(coords_z - coords_r, 2, axis=1)

# Perform spectral band alignment (assuming aligning with the 'r' band)
aligned_data_g = np.roll(data_g, (int(x_offset_g[0][0]), int(y_offset_g[0][0])), axis=(0, 1))
aligned_data_i = np.roll(data_i, (int(x_offset_i[0][0]), int(y_offset_i[0][0])), axis=(0, 1))
aligned_data_u = np.roll(data_u, (int(x_offset_u[0][0]), int(y_offset_u[0][0])), axis=(0, 1))
aligned_data_z = np.roll(data_z, (int(x_offset_z[0][0]), int(y_offset_z[0][0])), axis=(0, 1))


# Combine the spectral bands into a tensor
tensor_spectral = np.stack([data_r, aligned_data_g, aligned_data_i, aligned_data_u, aligned_data_z], axis=-1)

print("Tensor shape:", tensor_spectral.shape)
