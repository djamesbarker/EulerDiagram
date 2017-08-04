# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 14:30:25 2017

@author: leblanckh
"""

import pandas as pd
import numpy as np
from scipy import stats
import math
import os
import sys
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
plt.style.use('ggplot')

#==============================================================================
# make a whole library of small modular functions
# that let you answer given R1, R2, and D what is A
# given R1, R2, and D, what is Beta1 and Beta2 ?
# given R1, R2, and D what is Aisoc1 and Aisoc2 ?
#==============================================================================
RadSmall = 5
RadLg = 10
DistC2C = 4

def Distance_Circle_Center_to_Chord (RadiusBig,RadiusSmall,CCDist):
    """
    given the radii of the two intersecting circles and the distance between their center points, 
    calculate the distance between the center of the large circle and the chord
    
    inputs: RadiusBig, RadiusSmall, CCDist (D)
    return: CenterChordDistance
    """
    
    CenterChordDist = (RadiusBig**2 - RadiusSmall**2 - CCDist**2) / (-2*CCDist)
    return CenterChordDist
    
def Half_Chord_Length_Cusp_to_Cusp (RadiusBig,RadiusSmall,CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points, 
    calculate half the cord length of the Chord line
    """
    Xsolved = Distance_Circle_Center_to_Chord(RadiusBig, RadiusSmall, CCDist)    
    HalfChordLength = math.sqrt( RadiusSmall**2 - Xsolved**2)
    return HalfChordLength
    
def Area_of_Small_Triangle (RadiusBig,RadiusSmall,CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points,
    calculate the area of the triangle from the center point of the small circle to the chord
    """
    Xsolved = Distance_Circle_Center_to_Chord(RadiusBig, RadiusSmall, CCDist)
    Ysolved = Half_Chord_Length_Cusp_to_Cusp(RadiusBig, RadiusSmall, CCDist)
    ASmTri = Xsolved * Ysolved
    return ASmTri
    
def Area_of_Large_Triangle (RadiusBig, RadiusSmall, CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points,
    calculate the area of the triangle from the center point of the large circle to the chord
    """
    Xsolved = Distance_Circle_Center_to_Chord(RadiusBig, RadiusSmall, CCDist)
    Ysolved = Half_Chord_Length_Cusp_to_Cusp(RadiusBig, RadiusSmall, CCDist)
    ALgTri = Ysolved*(CCDist-Xsolved)
    return ALgTri
    
def Angle_of_Small_Triangle(RadiusBig, RadiusSmall, CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points,
    calculate the angle of the triangle from the center point of the small circle to the chord
    """
    Xsolved = Distance_Circle_Center_to_Chord(RadiusBig, RadiusSmall, CCDist)
    Beta1 = 2*(math.acos(Xsolved/RadiusSmall))
    return Beta1
    
def Angle_of_Large_Triangle(RadiusBig, RadiusSmall, CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points,
    calculate the angle of the triangle from the center point of the large circle to the chord
    """
    Xsolved = Distance_Circle_Center_to_Chord(RadiusBig, RadiusSmall, CCDist)
    Beta2 = 2*(math.acos((CCDist-Xsolved)/RadiusBig))
    return Beta2
    
def Area_of_Small_Sector (RadiusBig, RadiusSmall, CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points,
    calculate the area of the circle section from the center point of the small circle to the chord
    """
    AngleSm = Angle_of_Small_Triangle(RadiusBig, RadiusSmall, CCDist)
    AreaSmSect = (AngleSm/(2*math.pi))*math.pi*RadiusSmall**2
    return AreaSmSect
    
def Area_of_Large_Sector(RadiusBig, RadiusSmall, CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points,
    calculate the area of the circle section from the center point of the large circle to the chord
    """
    AngleLg = Angle_of_Large_Triangle(RadiusBig, RadiusSmall, CCDist)
    AreaLgSect = (AngleLg/(2*math.pi))*math.pi*RadiusBig**2
    return AreaLgSect
    
def Area_of_Small_Lens(RadiusBig, RadiusSmall, CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points,
    calculate the area of the small lens
    """
    AreaSmSect = Area_of_Small_Sector(RadiusBig, RadiusSmall, CCDist)
    AreaSmTri = Area_of_Small_Triangle(RadiusBig, RadiusSmall, CCDist)
    AreaSmLens = AreaSmSect-AreaSmTri
    return AreaSmLens
    
def Area_of_Large_Lens(RadiusBig, RadiusSmall, CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points,
    calculate the area of the small lens
    """
    AreaLgSect = Area_of_Large_Sector(RadiusBig, RadiusSmall, CCDist)
    AreaLgTri = Area_of_Large_Triangle(RadiusBig, RadiusSmall, CCDist)
    AreaLgLens = AreaLgSect - AreaLgTri
    return AreaLgLens
    
def Area_of_Overlap(RadiusBig, RadiusSmall, CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points,
    calculate the area of overlap of the two circles
    """
    AreaSmLens = Area_of_Small_Lens(RadiusBig, RadiusSmall, CCDist)
    AreaLgLens = Area_of_Large_Lens(RadiusBig, RadiusSmall, CCDist)
    AreaOverlap = AreaSmLens + AreaLgLens
    return AreaOverlap

if DistC2C>= RadSmall+RadLg:
    print ("There is no overlap between these circles.")

elif DistC2C <= RadLg - RadSmall:
    print ("The smaller circle is completely encircled by the larger circle.")

else:
    Overlap = Area_of_Overlap(RadLg,RadSmall,DistC2C)

print(Overlap)