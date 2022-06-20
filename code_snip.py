def WAVE(header):
    '''
    This function takes in the header from a cube and uses the keywords to
    define the transformation from pixel number to wavelength.
    '''
    CDELT=header['CDELT3'] # pixel size
    CRVAL=header['CRVAL3'] # value of pixel at CRPIX
    CRPIX=header['CRPIX3'] # pixel number of CRVAL
    npix=header['NAXIS3'] # number of pixels
    x=CRVAL+CDELT*(np.arange(npix)+1-CRPIX) # creates wavelength array
    return x
