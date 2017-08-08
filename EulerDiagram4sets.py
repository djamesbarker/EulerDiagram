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
RadLg = 100
DistC2C = 101
AreaOverlap = 28.94795433

def Distance_Circle_Center_to_Chord (RadiusBig,RadiusSmall,CCDist):
    """
    given the radii of the two intersecting circles and the distance between their center points, 
    calculate the distance between the center of the large circle and the chord
    
    inputs: RadiusBig, RadiusSmall, CCDist (D)
    return: CenterChordDistance
    """
    
    CenterChordDist = ( CCDist**2 - RadiusSmall**2 + RadiusBig**2) / (2*CCDist)
    return CenterChordDist
    
def Half_Chord_Length_Cusp_to_Cusp (RadiusBig,RadiusSmall,CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points, 
    calculate half the cord length of the Chord line
    """
    Xsolved = Distance_Circle_Center_to_Chord(RadiusBig, RadiusSmall, CCDist)    
    HalfChordLength = math.sqrt( RadiusBig**2 - Xsolved**2)
    return HalfChordLength
    
def Area_of_Small_Triangle (RadiusBig,RadiusSmall,CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points,
    calculate the area of the triangle from the center point of the small circle to the chord
    """
    Xsolved = Distance_Circle_Center_to_Chord(RadiusBig, RadiusSmall, CCDist)
    Ysolved = Half_Chord_Length_Cusp_to_Cusp(RadiusBig, RadiusSmall, CCDist)
    ASmTri = Ysolved*(CCDist-Xsolved)
    return ASmTri
    
def Area_of_Large_Triangle (RadiusBig, RadiusSmall, CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points,
    calculate the area of the triangle from the center point of the large circle to the chord
    """
    Xsolved = Distance_Circle_Center_to_Chord(RadiusBig, RadiusSmall, CCDist)
    Ysolved = Half_Chord_Length_Cusp_to_Cusp(RadiusBig, RadiusSmall, CCDist)
    ALgTri = Xsolved * Ysolved
    return ALgTri
    
def Angle_of_Small_Triangle(RadiusBig, RadiusSmall, CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points,
    calculate the angle of the triangle from the center point of the small circle to the chord
    """
    Xsolved = Distance_Circle_Center_to_Chord(RadiusBig, RadiusSmall, CCDist)
    X2 = CCDist - Xsolved
    ratio = X2/RadiusSmall
    Beta1 = 2*(math.acos(ratio))
    return Beta1
    
def Angle_of_Large_Triangle(RadiusBig, RadiusSmall, CCDist):
    """
    given the radii of the two intersecting circles, the distance between their center points,
    calculate the angle of the triangle from the center point of the large circle to the chord
    """
    Xsolved = Distance_Circle_Center_to_Chord(RadiusBig, RadiusSmall, CCDist)
    Beta2 = 2*(math.acos(Xsolved/RadiusBig))
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

def Calculate_Distance_for_Given_Overlap(RadiusBig,RadiusSmall,AOverlap):
    """
    guesses a default center to center distance and calculates the area of overlap given that distance
    then compares the calculated area to the desired area of overlap, and adjusts the distance accordingly
    """
    Tolerance = 0.0001
    Dmin = RadiusBig - RadiusSmall
    Dmax = RadiusBig + RadiusSmall
    if AOverlap <= 0:
        print ("There is no overlap between these circles.")
        Dguess = -1
    elif AOverlap >= math.pi*RadiusSmall**2:
        print("The smaller circle is completely encircled by the larger circle.")
        Dguess = 0
    else:
        Dguess = (Dmax+Dmin)/2
        print ('Dguess',Dguess)
        AOCalc = Area_of_Overlap(RadiusBig,RadiusSmall,Dguess)
        print('AoCalc',AOCalc)
        errorVal = AOverlap - AOCalc
        print('errorVal',errorVal)
        loopCounter = 0
        
        while (abs(errorVal)) > Tolerance:
            loopCounter +=1
            #adjust distance
            if errorVal > 0:
                Dmax = Dguess
                Dguess = (Dguess+Dmin)/2
                print('new Dmax = ',Dmax, 'new Dguess =',Dguess)
                AOCalc = Area_of_Overlap(RadiusBig,RadiusSmall,Dguess)
                print('AoCalc',AOCalc)
                errorVal = AOverlap - AOCalc
                print('errorVal',errorVal)
            elif errorVal < 0:
                Dmin = Dguess
                Dguess = (Dmax+Dguess)/2
                print('new Dmin = ',Dmin, 'new Dguess =',Dguess)
                AOCalc = Area_of_Overlap(RadiusBig,RadiusSmall,Dguess)
                print('AoCalc',AOCalc)
                errorVal = AOverlap - AOCalc
                print('errorVal',errorVal)
            else:
                return Dguess
    return Dguess

if DistC2C>= RadSmall+RadLg:
    print ("There is no overlap between these circles.")

elif DistC2C <= RadLg - RadSmall:
    print ("The smaller circle is completely encircled by the larger circle. Overlap = ", math.pi*RadSmall**2)
else:
    Overlap = Area_of_Overlap(RadLg,RadSmall,DistC2C)
    CalcDistC2C = Calculate_Distance_for_Given_Overlap(RadLg,RadSmall,AreaOverlap)

print(Overlap, CalcDistC2C)