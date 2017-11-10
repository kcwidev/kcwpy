""" Observation class for KCWI"""

import glob
import sys
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.io.fits as pf


class Observation:
    """KCWI observation class

    This class is used for holding all the definitions for a given observation
    and will be used for associations with calibration files and for processing
    steps."""

    telescope = 'Keck2'

    # Basic observation parameters
    object = None  # Object name
    targname = None  # Target name
    instrument = None
    observer = None
    airmass = None  # Airmass
    imgnum = None  # image number
    ofname = None  # Original filename for observation
    outdir = None  # Original directory for observation
    date_obs = None  # Date object specifying observation time

    # WMKO program/state info
    progname = None  # Program name
    statename = None  # State name
    stateid = None  # State ID

    # Pointing
    image_coords = None  # Coordinate object specifying image center
    targ_coords = None  # Coordinate object specifying target coords
    parang = None  # Parallactic angle in degrees
    rotposn = None  # Rotator position in degrees
    rotrefan = None  # Rotator reference angle in degrees
    rotmode = None  # Rotator mode

    # Detector parameters
    xposure = None  # Shutter open duration (s)
    telapse = None  # Dark current duration (s)
    numopen = None  # Number of shutter opens in exposure

    shufrows = None  # Number of CCD rows shuffled
    ccdmode = None  # CCD mode (0-slow, 1-fast)
    ampmode = None  # Amp mode name (ALL, TBO, TUP, C, D, E, F)
    nvidinp = None  # Number of amps
    nampsxy = None  # Number of amplifiers in x and y
    gainmul = None  # Gain multiplier (1, 2, 5, 10)
    ccdgain = None  # Average CCD gain in e-/DN
    gain1 = None  # Amp 1 gain in e-/DN
    gain2 = None  # Amp 2 gain in e-/DN
    gain3 = None  # Amp 3 gain in e-/DN
    gain4 = None  # Amp 4 gain in e-/DN

    # Instrument configuration
    dome = None  # dome lamp on: True or False
    hatch = None  # hatch position: Open or Closed
    caltype = None  # calibration type: bias, dark, arc, etc.
    calmnam = None  # Cal mirror position name ('Sky', 'Mirror', 'Filter')
    calpnam = None  # Cal position name ('Sky', 'Polar', 'Lens')
    polang = None  # Cal polarizer angle
    slicer = None  # Slicer name ("Small", "Medium", "Large")
    ifunum = None  # Slicer number ( 0 - 5, -1 = unknown)
    filter = None  # Blue filter id
    bfiltnum = None  # Blue filter number
    grating = None  # Blue grating id
    bgratnum = None  # Blue grating number
    grangle = None  # Blue grating angle (degrees)
    artang = None  # Blue articulation angle
    cwave = None  # Blue central wavelength (Ang)
    pwave = None  # Blue peak wavelength (Ang)
    focus = None  # Blue focus (mm)

    # Temperatures
    bccdtemp = None  # Blue CCD temperature
    benchtemp = None  # Optical Bench temperature

    # Derived parameters
    obstype = None  # Type of observation: obj, cal, test
    imgtype = None  # Type of image: obj, sky, tflat, dflat, cflat,
    #  arcflat, cbars, arcbars, test
    skyobs = None  # sky observation? True or False
    shuffmod = None  # is this a Nod & Shuffle observation? True or False
    nasmask = None  # Nod-and-Shuffle mask in? True or False
    nsskyr0 = None  # Sky region row 0 (bottom, pix)
    nsskyr1 = None  # Sky region row 1 (top, pix)
    nsobjr0 = None  # Object region row 0 (bottom, pix)
    nsobjr1 = None  # Object region row 1 (top, pix)
    xbinsize = None  # Binning in X
    ybinsize = None  # Binning in Y
    xsize = None  # Size of raw image in X
    ysize = None  # Size of raw image in Y
    illum = None  # Illumination source: (sky, FeAr, ThAr, Cont)
    biasrn1 = None  # Amp 1 readnoise in e- from bias
    biasrn2 = None  # Amp 2 readnoise in e- from bias
    biasrn3 = None  # Amp 3 readnoise in e- from bias
    biasrn4 = None  # Amp 4 readnoise in e- form bias
    exptime = None  # Exposure time (s)

    def match_bias(self, obs=None):
        """Does the input observation match at CCD bias level?

        Checks the relevant CCD parameters to see if the input observation
        matches such that they could use the same master bias frame.

        Args:
            obs (Observation): observation instance to match against

        Returns:
            bool: True if observation matches otherwise False

        """

        retval = True
        if obs is not None:
            # check relevant CCD parameters
            if self.ampmode not in obs.ampmode:
                retval = False
            if self.ccdmode != obs.ccdmode:
                retval = False
            if self.xbinsize != obs.xbinsize or self.ybinsize != obs.ybinsize:
                retval = False
            if self.gainmul != obs.gainmul:
                retval = False
        # empty observation
        else:
            retval = False

        return retval

    def match_geom(self, obs=None, use_id=True, delang=0.1, delwave=0.1):
        """Does the input observation match at Geometry Configuration level?

        Checks the relevant configuration parameters to see if the input
        observation matches such that they could use the same geometry
        calibration.

        Args:
            obs (Observation): observation instance to match against

        Returns:
            bool: True if observation matches otherwise False

        """

        retval = True
        if obs is not None:
            # just match on state Id
            if use_id:
                if self.stateid not in obs.stateid:
                    retval = False
            # use details of config to match
            else:
                if self.xbinsize != obs.xbinsize or self.ybinsize != obs.ybinsize:
                    retval = False
                if self.grating not in obs.grating:
                    retval = False
                if self.filter not in obs.filter:
                    retval = False
                if self.slicer not in obs.slicer:
                    retval = False
                if abs(self.grangle - obs.grangle) > delang:
                    retval = False
                if abs(self.cwave - obs.cwave) > delwave:
                    retval = False
        # Empty observation
        else:
            retval = False

        return retval

    def __init__(self, hdr):
        """Initialize observation object using fits header"""

        if header_integrity(hdr):

            # Basic observation parameters
            self.object = hdr['OBJECT']
            self.targname = hdr['TARGNAME']
            self.instrument = hdr['INSTRUME']
            self.observer = hdr['OBSERVER']
            self.airmass = hdr['AIRMASS']
            self.imgnum = hdr['FRAMENO']
            self.ofname = hdr['OFNAME']
            self.outdir = hdr['OUTDIR']
            self.date_obs = Time(hdr['DATEPCLR'], format='isot',
                                 scale='utc')
            # WMKO program info
            self.progname = hdr['PROGNAME']
            self.statename = hdr['STATENAM']
            self.stateid = hdr['STATEID']
            # Pointing
            self.image_coords = SkyCoord(hdr['RA'], hdr['DEC'],
                                         unit=(u.hourangle, u.deg))
            self.targ_coords = SkyCoord(hdr['TARGRA'], hdr['TARGDEC'],
                                        unit=(u.hourangle, u.deg))
            self.parang = hdr['PARANG']
            self.rotposn = hdr['ROTPOSN']
            self.rotmode = hdr['ROTMODE']
            # Detector configuration
            self.xposure = hdr['XPOSURE']
            self.telapse = hdr['TELAPSE']
            self.exptime = self.xposure
            self.numopen = hdr['NUMOPEN']
            self.ccdmode = hdr['CCDMODE']
            self.ampmode = hdr['AMPMODE']
            self.nvidinp = hdr['NVIDINP']
            self.nampsxy = hdr['NAMPSXY']
            self.gainmul = hdr['GAINMUL']
            self.ccdgain = hdr['CCDGAIN']
            self.gain1 = hdr['GAIN1']
            if self.nvidinp > 1:
                self.gain2 = hdr['GAIN2']
                if self.nvidinp > 2:
                    self.gain3 = hdr['GAIN3']
                    if self.nvidinp > 3:
                        self.gain4 = hdr['GAIN4']
            self.xsize = hdr['NAXIS1']
            self.ysize = hdr['NAXIS2']
            self.xbinsize = int(hdr['CCDSUM'].split()[0])
            self.ybinsize = int(hdr['CCDSUM'].split()[1])
            self.shufrows = hdr['SHUFROWS']
            if hdr['NSHFUP'] > 0 or hdr['NSHFDN'] > 0:
                self.shuffmod = True
                self.nsskyr0 = 1
                self.nsskyr1 = self.shufrows
                self.nsobjr0 = self.nsskyr1 + 1
                self.nsobjr1 = self.nsobjr0 + self.shufrows - 1
            else:
                self.shuffmod = False
            # Instrument configuration
            self.dome = ('on' in hdr['FLIMAGIN'] or 'on' in hdr['FLSPECTR'])
            self.hatch = hdr['HATPOS']
            self.caltype = hdr['CALTYPE']
            self.calmnam = hdr['CALMNAM']
            self.calpnam = hdr['CALPNAM']
            self.polang = hdr['CALLANG']
            self.slicer = hdr['IFUNAM']
            self.ifunum = hdr['IFUNUM']
            self.filter = hdr['BFILTNAM']
            self.bfiltnum = hdr['BFILTNUM']
            self.grating = hdr['BGRATNAM']
            self.gratnum = hdr['BGRATNUM']
            self.grangle = hdr['BGRANGLE']
            self.artang = hdr['BARTANG']
            self.cwave = hdr['BCWAVE']
            self.pwave = hdr['BPWAVE']
            self.focus = hdr['BFOCMM']
            if 'Mask' in hdr['BNASNAM']:
                self.nasmask = True
            else:
                self.nasmask = False
            # Temperatures
            self.bccdtemp = hdr['BCCDTMP']
            self.benchtemp = hdr['TMPA7']
            # Determine illumination: default is 'Test'
            self.illum = 'Test'
            if 'Mirror' in self.calmnam:
                self.illum = 'IntCal'
                if hdr['LMP0STAT'] and hdr['LMP0SHST']:
                    self.illum = hdr['LMP0NAM']
                if hdr['LMP1STAT'] and hdr['LMP1SHST']:
                    self.illum = hdr['LMP1NAM']
                if hdr['LMP3STAT']:
                    self.illum = hdr['LMP3NAM']
            # External Cal mirror position
            else:
                if 'Open' in self.hatch:
                    self.illum = 'Sky'
                    if self.dome:
                        self.illum = 'Dome'
            # Determine image/observation type
            self.skyobs = False
            self.imgtype = self.caltype
            self.obstype = 'test'
            if 'object' in self.caltype:
                if 'Open' in self.hatch:
                    if self.dome:
                        self.imgtype = 'dflat'
                        self.obstype = 'cal'
                    else:
                        self.imgtype = 'object'
                        self.obstype = 'obj'
                else:
                    self.imgtype = 'test'
                    self.obstype = 'test'
            elif 'arc' in self.caltype:
                self.imgtype = 'arc'
                self.obstype = 'cal'
            elif 'cbars' in self.caltype or 'cflat' in self.caltype:
                self.obstype = 'cal'
            elif 'dark' in self.caltype or 'Dark' in hdr['BNASNAM']:
                self.illum = 'None'
                if self.xposure == 0. and self.telapse > 0.1:
                    self.imgtype = 'dark'
                    self.obstype = 'zero'
                    self.exptime = self.telapse
                elif self.xposure == 0. and self.telapse <= 0.1:
                    self.imgtype = 'bias'
                    self.obstype = 'zero'
            # Check for biases (shutter did not open)
            if self.xposure == 0.:
                self.illum = 'None'
                if self.telapse > 0.1:
                    self.imgtype = 'dark'
                    self.obstype = 'zero'
                    self.exptime = self.telapse
                else:
                    self.imgtype = 'bias'
                    self.obstype = 'zero'
            # Check for direct mode
            if 'None' in self.grating and self.artang < 5. and self.exptime > 0:
                self.obstype = 'dir-' + self.obstype
                if 'test' in self.imgtype:
                    self.imgtype = 'object'
                    self.obstype = 'dir-cal'
        else:
            print("KCWPY Error - bad header")


def header_integrity(hdr):
    """ Test the integrity of the input KCWI image header"""

    # Keywords that are essential
    keys_err = ['OBJECT', 'TARGNAME', 'OBSERVER', 'AIRMASS', 'FRAMENO',
                'OFNAME', 'OUTDIR', 'DATEPCLR', 'PROGNAME', 'STATENAM',
                'STATEID', 'RA', 'DEC', 'TARGRA', 'TARGDEC', 'PARANG',
                'ROTPOSN', 'ROTMODE', 'XPOSURE', 'TELAPSE', 'NUMOPEN',
                'NSHFUP', 'NSHFDN', 'SHUFROWS',
                'CCDMODE', 'AMPMODE', 'NVIDINP', 'NAMPSXY', 'CCDSUM',
                'GAINMUL', 'CCDGAIN', 'GAIN1', 'FLIMAGIN', 'FLSPECTR',
                'HATPOS', 'CALTYPE', 'CALMNAM', 'CALPNAM', 'CALLANG',
                'IFUNAM', 'IFUNUM', 'BFILTNAM', 'BFILTNUM', 'BGRATNAM',
                'BGRATNUM', 'BGRANGLE', 'BGRENC', 'BARTANG', 'BARTENC',
                'BCWAVE', 'BPWAVE', 'BFOCMM', 'BFOCENC', 'BNASNAM',
                'LMP0STAT', 'LMP0SHST', 'LMP1STAT', 'LMP1SHST', 'LMP3STAT']

    # Check keywords
    for k in keys_err:
        if k not in hdr:
            print("KCWPY Error - Missing FITS keyword %s" % k)
            return False

    return True


def read_headers(flist=None, verbose=False):
    """ Read in fits files and return array of Observation structures """

    # check for input list
    if flist is None:
        flist = glob.glob("./kb*.fits")
    # loop over files
    oblist = []
    for f in flist:
        if verbose:
            print(f)
        ff = pf.open(f)
        ob = Observation(ff[0].header)
        oblist.append(ob)

    return oblist


def what(oblist):
    """ Print out a summary of each observation in oblist. """

    # header
    print("# R   = CCD Readout Speed : 0 - slow, 1 - fast")
    print("# G   = Gain multiplier   : 10, 5, 2, 1")
    print("# SSM = Sky, Shuffle, Mask: 0 - no, 1 - yes")
    print("#  #/   N    Imno Bin AMP R  G SSM IFU GRAT FILT     Cwave JDobs"
          "         Expt ObType  Type    Illum     Imno   "
          "RA          Dec            PA     Air   Object")
    # how many?
    nobs = len(oblist)
    # loop over oblist
    for n, ob in enumerate(oblist):
        print("%4d/%4d %7d %1d %1d %3s %1d %2d %1d%1d%1d %3s %-4s %-5s %8.1f"
              "%12.3f%7.1f %-7s %-7s %-6s %7d%13.8f%13.8f %7.2f %7.3f %s" %
              (n, nobs, ob.imgnum, ob.xbinsize, ob.ybinsize, ob.ampmode,
               ob.ccdmode, ob.gainmul, (1 if ob.skyobs else 0),
               (1 if ob.shuffmod else 0), (1 if ob.nasmask else 0),
               ob.slicer[:3], ob.grating[:4], ob.filter, ob.cwave,
               ob.date_obs.jd, ob.exptime,
               ob.obstype, ob.imgtype, ob.illum[:6], ob.imgnum,
               ob.image_coords.ra.degree, ob.image_coords.dec.degree,
               ob.rotposn, ob.airmass,
               (ob.targname if 'object' in ob.imgtype else '-')))


if __name__ == '__main__':

    files = sys.argv[2:]
    obslist = read_headers(files)
    what(obslist)
