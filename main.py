from jplephem.spk import SPK
from izzoLambertSolver import izzo2015
from plottingV2 import Plotting
from fourthRungeKutta import computeRk4

import numpy as np
import pandas as pd

import math

# Constants/Variable Declarations
sunGravParam = 1.327e20

counter = 0

kernel = SPK.open('de440.bsp')
# 1977 8 19 -- voyager 2 launch date

# For loops used to check possible year and month arrangements
for yearIterator in range(2130, 2180, 1):
    for monthIterator in range(1, 13, 2):

        # Keeps count of fly-by's that were done.
        counter = 0
        rawTime = pd.Timestamp(year=yearIterator,  month=monthIterator, day=25,
                      hour=0, second=0, tz='US/Eastern')

        julianTime = rawTime.to_julian_date()

        ignoreThis, velocityEarth = kernel[0, 3].compute_and_differentiate(julianTime)

        # Position of Earth w.r.t. solar system center
        earthPositionKM = kernel[0, 3].compute(julianTime)
        earthPositionM = earthPositionKM*1000

        # Position of Jupiter w.r.t. solar system center
        jupiterPositionKM = kernel[0, 5].compute(julianTime+365)
        jupiterPositionM = jupiterPositionKM * 1000

        # Position of Saturn w.r.t. solar system center
        saturnPositionKM = kernel[0, 6].compute(julianTime)
        saturnPositionM = saturnPositionKM * 1000

        # Time I initially choose
        transferTime_e2j_days = 365
        transferTime_e2j = 86400*transferTime_e2j_days

        #
        v_initial_e2j, v_final_e2j = izzo2015(sunGravParam, earthPositionM, jupiterPositionM, transferTime_e2j)

        velocityEarthAdjusted = (velocityEarth/86400.0)*1000.0


        updatedJulianTime = julianTime + transferTime_e2j_days

        # Used to get the velocity of Jupiter at the time of arrival by the craft;
        # used to ultimately calculate the v_infinity values
        ignoreThis, velocityJupiter1 = kernel[0, 5].compute_and_differentiate(updatedJulianTime)
        velocityJupiter = velocityJupiter1/86400*1000

        v_infinity_jupiter_arrive = v_final_e2j - velocityJupiter

        jupiterPosition1 = (kernel[0, 5].compute(updatedJulianTime)) * 1000

        # +Loop used to iteratively check for what transfer time it would take
        # to get from saturn to jupiter that also has the same v_infinity value
        # that the craft had when reaching jupiter from earth. A v_inf value that is
        # equal in magnitude at arrival and departure of saturn is the definition of a fly-by.
        for i in range(50, 800, 1):
            saturnPosition1 = (kernel[0, 6].compute(updatedJulianTime+i)) * 1000
            v_initial_j2s, v_final_j2s = izzo2015(sunGravParam, jupiterPosition1, saturnPosition1, i*86400)

            ignoreThis, velocityJupiter3 = kernel[0, 5].compute_and_differentiate(updatedJulianTime)
            velocityJupiter2 = velocityJupiter3 / 86400 * 1000

            v_infinity_jupiter_depart = v_initial_j2s - velocityJupiter2
            #print("vinf_a:", v_infinity_jupiter_depart)
            #print("vinf_d:", v_infinity_jupiter_arrive)

            norm_infinity_arrive = np.linalg.norm(v_infinity_jupiter_arrive)
            norm_infinity_depart = np.linalg.norm(v_infinity_jupiter_depart)

            #print("vizzy:", v_initial_j2s)

            #print("n_a:", norm_infinity_arrive)
            #print("n_d:", norm_infinity_depart)



            if abs(norm_infinity_arrive - norm_infinity_depart) < 50:
                #print("DONE")
                #print("vinf_a:", v_infinity_jupiter_depart)
                #print("vinf_d:", v_infinity_jupiter_arrive)
                #print(i)
                counter += 1
                break

        updatedJulianTime1 = updatedJulianTime + 424

        ignoreThis, velocitySaturn3 = kernel[0, 6].compute_and_differentiate(updatedJulianTime1)
        velocitySaturn2 = velocitySaturn3 / 86400 * 1000

        saturnPosition1 = (kernel[0, 6].compute(updatedJulianTime1)) * 1000
        v_infinity_saturn_arrive = v_final_j2s - velocitySaturn2

        # Saving the transfer time for later use
        transferTime_j2s_days = i

        for i in range(50, 1200, 1):
            uranusPosition1 = (kernel[0, 7].compute(updatedJulianTime1+i)) * 1000
            v_initial_s2u, v_final_s2u = izzo2015(sunGravParam, saturnPosition1, uranusPosition1, i*86400)



            v_infinity_saturn_depart = v_initial_s2u - velocitySaturn2
            #print("vinf_a:", v_infinity_jupiter_depart)
            #print("vinf_d:", v_infinity_jupiter_arrive)

            norm_infinity_arrive = np.linalg.norm(v_infinity_saturn_arrive)
            norm_infinity_depart = np.linalg.norm(v_infinity_saturn_depart)



            #print("n_a:", norm_infinity_arrive)
            #print("n_d:", norm_infinity_depart)

            #print("-----")

            if abs(norm_infinity_arrive - norm_infinity_depart) < 50:

                #print("DONE")
                #print("vizzyerino:", v_initial_s2u)
                #print("vinf_a:", v_infinity_saturn_depart)
                #print("vinf_d:", v_infinity_saturn_arrive)
                #print('____')
                #print(v_final_j2s)
                #print(v_initial_s2u)
                #print("_____")
                #print(i)
                counter += 1
                break


        transferTime_s2u_days = i

        updatedJulianTime2 = updatedJulianTime1 + 424+923

        ignoreThis, velocityUranus3 = kernel[0, 7].compute_and_differentiate(updatedJulianTime2)
        velocityUranus2 = velocityUranus3 / 86400 * 1000

        uranusPosition1 = (kernel[0, 7].compute(updatedJulianTime2)) * 1000
        v_infinity_uranus_arrive = v_final_s2u - velocityUranus2


        for i in range(50, 1200, 1):
            neptunePosition1 = (kernel[0, 8].compute(updatedJulianTime2+i)) * 1000
            v_initial_u2n, v_final_u2n = izzo2015(sunGravParam, uranusPosition1, neptunePosition1, i*86400)



            v_infinity_uranus_depart = v_initial_u2n - velocityUranus2
            #print("vinf_a:", v_infinity_jupiter_depart)
            #print("vinf_d:", v_infinity_jupiter_arrive)

            norm_infinity_arrive = np.linalg.norm(v_infinity_uranus_arrive)
            norm_infinity_depart = np.linalg.norm(v_infinity_uranus_depart)



            #print("n_a:", norm_infinity_arrive)
            #print("n_d:", norm_infinity_depart)

            #print("-----")

            if abs(norm_infinity_arrive - norm_infinity_depart) < 50:

                #print("DONE")
                #print("vizzyerino:", v_initial_u2n)
                #print("vinf_a:", v_infinity_saturn_depart)
                #print("vinf_d:", v_infinity_saturn_arrive)
                #print('____')
                #print(v_final_j2s)
                #print(v_initial_u2n)
                #print("_____")
                #print(i)
                counter += 1
                break

        transferTime_u2n_days = i

        timestep = 86400

        # Counter equaling 3 means the craft reached Saturn,Uranus, and Neptune
        if counter == 3:
            print("WINNER!")
            print(yearIterator)

            # Function which handles the plotting of the trajectories
            Plotting(earthPositionM[0], earthPositionM[1], earthPositionM[2],
                 v_initial_e2j[0], v_initial_e2j[1], v_initial_e2j[2], transferTime_e2j_days,
                 v_initial_j2s[0], v_initial_j2s[1], v_initial_j2s[2], transferTime_j2s_days,
                 v_initial_s2u[0], v_initial_s2u[1], v_initial_s2u[2], transferTime_s2u_days,
                 v_initial_u2n[0], v_initial_u2n[1], v_initial_u2n[2], transferTime_u2n_days,
                 timestep, julianTime)
        #else:
        #    print("NO BUENO. GOODBYE!")