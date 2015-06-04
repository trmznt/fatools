#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Single-function module.  

   Function and module names are the same.  See function docstring for 
   description.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: wavelen2rgb.py,v 1.2 2003/12/20 21:26:35 jlin Exp $
#
# Modification History:
# - 4 Dec 2003:  Original by Johnny Lin, Computation Institute,
#   University of Chicago.  Passed passably reasonable tests.
#
# Notes:
# - Written for Python 2.2.
# - Module docstrings can be tested using the doctest module.  To
#   test, execute "python wavelen2rgb.py".
# - No non-"built-in" packages and modules required.
#
# Copyright (c) 2003 by Johnny Lin.  For licensing, distribution 
# conditions, contact information, and additional documentation see
# the URL http://www.johnny-lin.com/pylib.html.
#=======================================================================




#-------------------- Overall Function Declaration ---------------------

def wavelen2rgb(Wavelength, MaxIntensity=100):
    """Calculate RGB values given the wavelength of visible light.

    Arguments:
    * Wavelength:  Wavelength in nm.  Scalar floating.
    * MaxIntensity:  The RGB value for maximum intensity.  Scalar 
      integer.

    Returns:
    * 3-element list of RGB values for the input wavelength.  The
      values are scaled from 0 to MaxIntensity, where 0 is the
      lowest intensity and MaxIntensity is the highest.  Integer
      list.

    Visible light is in the range of 380-780 nm.  Outside of this
    range the returned RGB triple is [0,0,0].

    Based on code by Earl F. Glynn II at:
       http://www.efg2.com/Lab/ScienceAndEngineering/Spectra.htm
    See also:
       http://www.physics.sfasu.edu/astro/color/spectra.html
    whose code is what Glynn's code is based on.

    Example:
    >>> from wavelen2rgb import wavelen2rgb
    >>> waves = [300.0, 400.0, 600.0]
    >>> rgb = [wavelen2rgb(waves[i], MaxIntensity=255) for i in range(3)]
    >>> print rgb
    [[0, 0, 0], [131, 0, 181], [255, 190, 0]]
    """




#--------------------------- Helper Function ---------------------------

    def Adjust_and_Scale(Color, Factor, Highest=100):
        """Gamma adjustment.

        Arguments:
        * Color:  Value of R, G, or B, on a scale from 0 to 1, inclusive,
          with 0 being lowest intensity and 1 being highest.  Floating
          point value.
        * Factor:  Factor obtained to have intensity fall off at limits 
          of human vision.  Floating point value.
        * Highest:  Maximum intensity of output, scaled value.  The 
          lowest intensity is 0.  Scalar integer.

        Returns an adjusted and scaled value of R, G, or B, on a scale 
        from 0 to Highest, inclusive, as an integer, with 0 as the lowest 
        and Highest as highest intensity.

        Since this is a helper function I keep its existence hidden.
        See http://www.efg2.com/Lab/ScienceAndEngineering/Spectra.htm and
        http://www.physics.sfasu.edu/astro/color/spectra.html for details.
        """
        Gamma = 0.80

        if Color == 0.0:
            result = 0
        else:
            result = int( round(pow(Color * Factor, Gamma) * round(Highest)) )
            if result < 0:        result = 0
            if result > Highest:  result = Highest

        return result




#------------------------ Overall Function Code ------------------------

    #- Set RGB values (normalized in range 0 to 1) depending on
    #  heuristic wavelength intervals:

    if (Wavelength >= 380.0) and (Wavelength < 440.0):
        Red   = -(Wavelength - 440.) / (440. - 380.)
        Green = 0.0
        Blue  = 1.0

    elif (Wavelength >= 440.0) and (Wavelength < 490.0):
        Red   = 0.0
        Green = (Wavelength - 440.) / (490. - 440.)
        Blue  = 1.0

    elif (Wavelength >= 490.0) and (Wavelength < 510.0):
        Red   = 0.0
        Green = 1.0
        Blue  = -(Wavelength - 510.) / (510. - 490.)

    elif (Wavelength >= 510.0) and (Wavelength < 580.0):
        Red   = (Wavelength - 510.) / (580. - 510.)
        Green = 1.0
        Blue  = 0.0

    elif (Wavelength >= 580.0) and (Wavelength < 645.0):
        Red   = 1.0
        Green = -(Wavelength - 645.) / (645. - 580.)
        Blue  = 0.0

    elif (Wavelength >= 645.0) and (Wavelength <= 780.0):
        Red   = 1.0
        Green = 0.0
        Blue  = 0.0

    else:
        Red   = 0.0
        Green = 0.0
        Blue  = 0.0


    #- Let the intensity fall off near the vision limits:

    if (Wavelength >= 380.0) and (Wavelength < 420.0):
        Factor = 0.3 + 0.7*(Wavelength - 380.) / (420. - 380.)
    elif (Wavelength >= 420.0) and (Wavelength < 701.0):
        Factor = 1.0
    elif (Wavelength >= 701.0) and (Wavelength <= 780.0):
        Factor = 0.3 + 0.7*(780. - Wavelength) / (780. - 700.)
    else:
        Factor = 0.0


    #- Adjust and scale RGB values to 0 to MaxIntensity integer range:

    R = Adjust_and_Scale(Red,   Factor, MaxIntensity)
    G = Adjust_and_Scale(Green, Factor, MaxIntensity)
    B = Adjust_and_Scale(Blue,  Factor, MaxIntensity)


    #- Return 3-element list value:

    return [R, G, B]




#-------------------------- Main:  Test Module -------------------------

#- Define additional examples for doctest to use:

__test__ = { 'Additional Example 1':
    """
    >>> from wavelen2rgb import wavelen2rgb
    >>> waves = [450.0, 550.0, 800.0]
    >>> rgb = [wavelen2rgb(waves[i], MaxIntensity=255) for i in range(3)]
    >>> print rgb
    [[0, 70, 255], [163, 255, 0], [0, 0, 0]]
    >>> rgb = [wavelen2rgb(waves[i]) for i in range(3)]
    >>> print rgb
    [[0, 28, 100], [64, 100, 0], [0, 0, 0]]
    """ }


#- Execute doctest if module is run from command line:

if __name__ == "__main__":
    """Test the module.

    Tests the examples in all the module documentation strings, plus 
    __test__.

    Note:  In case this is distributed as a package, we add the parent 
    directory to to sys.path for this module testing case to enable 
    doctest to work in more circumstances.
    """
    import doctest, sys, os
    sys.path.append(os.pardir)
    doctest.testmod(sys.modules[__name__])




# ===== end file =====
