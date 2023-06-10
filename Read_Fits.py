from astropy.io import fits
import numpy as np

fits_file = '/Users/batuhank1/PycharmProjects/ArchOfML/frame-g-003918-3-0213.fits'

# Read the .fits file
hdulist = fits.open(fits_file)
data = hdulist[0].data

# Convert the data to a structured numpy array
structured_array = np.array(data.tolist(), dtype=data.dtype)

# Get the header for the CSV file
header = ','.join(structured_array.dtype.names) if structured_array.dtype.names is not None else ''

# Save the structured array to a CSV file
csv_file = '/Users/batuhank1/PycharmProjects/ArchOfML/g_fits.csv'
np.savetxt(csv_file, structured_array, delimiter=',', header=header, comments='')


from astropy.io import fits

fits_file = '/Users/batuhank1/PycharmProjects/ArchOfML/frame-g-003918-3-0213.fits'  # Replace with the actual path to your FITS file

# Open the FITS file
hdulist = fits.open(fits_file)

# Print the information about the HDU list
hdulist.info()

# Access individual HDUs and their data
for i, hdu in enumerate(hdulist):
    print(f"HDU {i+1}")
    print(hdu.header)  # Print the header information for the HDU
    print(hdu.data)    # Print the data for the HDU

# Close the FITS file
hdulist.close()
