""" Observation class for KCWI"""

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord


class Observation:
    """KCWI observation class

    This class is used for holding all the definitions for a given observation
    and will be used for associations with calibration files and for processing
    steps."""

    date_obs = None  # Date object specifying observation time
    object = None  # Object name
    targname = None  # Target name
    obs_type = None  # Type of observation: obj, cal, test
    img_type = None  # Type of image: obj, sky, tflat, dflat, cflat,
    #  arcflat, cbars, arcbars, test
    image_coords = None  # Coordinate object specifying image center
    targ_coords = None  # Coordinate object specifying target coords
    observer = None
    telescope = None
    instrument = None
    datepclr = None  # Date of pre-clear
    daterend = None  # Date of readout end
    el = None  # Telescope elevation in degrees
    az = None  # Telescope azimuth in degrees
    parang = None  # Parallactic angle in degrees
    rotposn = None  # Rotator position in degrees
    rotrefan = None  # Rotator reference angle in degrees
    rotmode = None  # Rotator mode
    xposure = None  # Shutter open duration (s)
    telapse = None  # Dark current duration (s)
    airmass = None  # Airmass
    obsfname = None  # Original filename for observation
    obsdir = None  # Original directory for observation
    hatpos = None  # hatch position: Open or Closed
    flimagin = None  # dome lamp imaging mode: on or off
    flspectr = None  # dome lamp spectra mode: on or off
    caltype = None  # calibration type: bias, dark, arc, etc.
    frameno = None  # image number
    skyobs = None  # sky observation? 0 - object, 1 - sky
    shuffmod = None  # is this a Nod & Shuffle observation?
    bgratnam = None  # Blue grating id
    bgratnum = None  # Blue graing number
    bgrangle = None  # Blue grating angle (degrees)
    bgrenc = None  # Blue grating rotator encoder steps
    bfiltnam = None  # Blue filter id
    bfiltnum = None  # Blue filter number
    bartang = None  # Blue articulation angle
    bartenc = None  # Blue cmaera articulation encoder steps
    bcwave = None  # Blue central wavelength (Ang)
    bpwave = None  # Blue peak wavelength (Ang)
    bfocpos = None  # Blue focus stage encoder steps
    bfocus = None  # Blue focus (mm)
    bnasnam = None  # Blue mask position name
    bnaspos = None  # Blue mask position
    shufrows = None  # Number of CCD rows shuffled
    calmnam = None  # Cal mirror position name ('Sky', 'Mirror', 'Filter')
    calpnam = None  # Cal position name ('Sky', 'Polar', 'Lens')
    callang = None  # Cal polarizer angle
    ifunum = None  # Slicer number ( 0 -5, - 1 =unknown)
    ifunam = None  # Slicer name ("Small", etc.)
    cwave = None  # central wavelength (Ang)
    gratanom = None  # grating angle anomoly (degrees)
    wave0 = None  # blue end  of wavelength range (Ang)
    wave1 = None  # red end of wavelength range (Ang)
    dwav = None  # average dispersion ( Ang /pix)
    tmpbccd = None  # Blue CCD temperature
    tmpbench = None  # Optical Bench tempaerture
    illum = None  # Illumination source: (sky, FeAr, ThAr, Cont)
    biasrn1 = None  # Amp 1 readnoise in e- from bias
    biasrn2 = None  # Amp 2 readnoise in e- from bias
    biasrn3 = None  # Amp 3 readnoise in e- from bias
    biasrn4 = None  # Amp 4 readnoise in e- form bias
    gain1 = None  # Amp 1 gain in e-/DN
    gain2 = None  # Amp 2 gain in e-/DN
    gain3 = None  # Amp 3 gain in e-/DN
    gain4 = None  # Amp 4 gain in e-/DN
    gainmul = None  # Gain multiplier (1, 2, 5, 10)
    ccdmode = None  # CCD mode (0-slow, 1-fast)
    ampmode = None  # Amp mode name (ALL, TBO, TUP, C, D, E, F)
    nvidinp = None  # Number of amps
    nampsxy = None  # Number of amplifiers in x and y
    xbinsize = None  # Binning in X
    ybinsize = None  # Binning in Y
    xsize = None  # Size of raw image in X
    ysize = None  # Size of raw image in Y
    exptime = None  # Exposure time (s)
    nsskyr0 = None  # Sky region row 0 (bottom, pix)
    nsskyr1 = None  # Sky region row 1 (top, pix)
    nsobjr0 = None  # Object region row 0 (bottom, pix)
    nsobjr1 = None  # Object region row 1 (top, pix)

    def __init__(self, hdr):
        """Initialize observation object using fits header"""

        if 'OBJECT' in hdr:
            self.object = hdr['OBJECT']
        if 'TARGNAME' in hdr:
            self.targname = hdr['TARGNAME']
        if 'INSTRUME' in hdr:
            self.instrument = hdr['INSTRUME']
        else:
            self.instrument = 'KCWI'
        if 'OBSERVER' in hdr:
            self.observer = hdr['OBSERVER']
        if 'TELESCOP' in hdr:
            self.telescope = hdr['TELESCOPE']
        else:
            self.telescope = 'Keck2'
        if 'DATEPCLR' in hdr:
            self.date_obs = Time(hdr['DATEPCLR'], format='isot',
                                 scale='utc')
        if 'TARGRA' in hdr and 'TARGDEC' in hdr:
            self.targ_coords = SkyCoord(hdr['TARGRA'], hdr['TARGDEC'],
                                        unit=(u.hourangle, u.deg))
        if 'RA' in hdr and 'DEC' in hdr:
            self.image_coords = SkyCoord(hdr['RA'], hdr['DEC'],
                                         unit=(u.hourangle, u.deg))
        if 'HATPOS' in hdr:
            self.hatpos = hdr['HATPOS']
        if 'CALMNAM' in hdr:
            self.calmnam = hdr['CALMNAM'].rstrip()
        if 'Mirror' in self.calmnam:
            self.illum = 'Cal'
            if 'LMP0STAT' in hdr and 'LMP0SHST' in hdr:
                if hdr['LMP0STAT'] and hdr['LMP0SHST']:
                    self.illum = hdr['LMP0NAM'].rstrip()
            if 'LMP1STAT' in hdr and 'LMP1SHST' in hdr:
                if hdr['LMP1STAT'] and hdr['LMP1SHST']:
                    self.illum = hdr['LMP1NAM'].rstrip()
            if 'LMP3STAT' in hdr:
                if hdr['LMP3STAT']:
                    self.illum = hdr['LMP3NAM'].rstrip()
        else:
            if 'Open' in self.hatpos:
                self.illum = 'Sky'
            else:
                self.illum = 'Test'

