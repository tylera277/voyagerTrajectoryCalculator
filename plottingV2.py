import matplotlib.pyplot as plt

from mpl_toolkits import mplot3d
from jplephem.spk import SPK
from fourthRungeKutta import computeRk4

kernel = SPK.open('de440.bsp')

def Plotting(xCraftStart, yCraftStart, zCraftStart,
                xVel_leg1, yVel_leg1, zVel_leg1, tof_leg1,
                 xVel_leg2, yVel_leg2, zVel_leg2, tof_leg2,
                 xVel_leg3, yVel_leg3, zVel_leg3, tof_leg3,
                 xVel_leg4, yVel_leg4, zVel_leg4, tof_leg4,
                 timestep, julianTime):

    mass_sun = 1.989e30

    kernel = SPK.open('de440.bsp')

    xTrajectory, yTrajectory , zTrajectory = [], [], []
    xEarthTrajectory, yEarthTrajectory , zEarthTrajectory = [], [], []
    xJupiterTrajectory, yJupiterTrajectory, zJupiterTrajectory = [], [], []
    xSaturnTrajectory, ySaturnTrajectory, zSaturnTrajectory  = [], [], []
    xUranusTrajectory, yUranusTrajectory, zUranusTrajectory = [], [], []
    xNeptuneTrajectory, yNeptuneTrajectory, zNeptuneTrajectory = [], [], []

    xCraft = xCraftStart
    yCraft = yCraftStart
    zCraft = zCraftStart

    xSun = 0
    ySun = 0
    zSun = 0

    # ###### Leg 1 ########
    tFirstLeg = 0
    tEndFirstLeg = tof_leg1 * 86400

    vxCraft = xVel_leg1
    vyCraft = yVel_leg1
    vzCraft = zVel_leg1

    while tFirstLeg < tEndFirstLeg:
        k = computeRk4(xCraft, yCraft, zCraft, xSun, ySun, zSun, vxCraft, vyCraft, vzCraft, mass_sun, timestep)

        xCraft += (timestep / 6) * (k[0] + 2 * k[6] + 2 * k[12] + k[18])
        yCraft += (timestep / 6) * (k[1] + 2 * k[7] + 2 * k[13] + k[19])
        zCraft += (timestep / 6) * (k[2] + 2 * k[8] + 2 * k[14] + k[20])

        vxCraft += (timestep / 6) * (k[3] + 2 * k[9] + 2 * k[15] + k[21])
        vyCraft += (timestep / 6) * (k[4] + 2 * k[10] + 2 * k[16] + k[22])
        vzCraft += (timestep / 6) * (k[5] + 2 * k[11] + 2 * k[17] + k[23])

        updatedJulianTime = julianTime + (tFirstLeg / 86400.0)
        earthPositionKM = kernel[0, 3].compute(updatedJulianTime)
        earthPositionM = earthPositionKM * 1000

        jupiterPositionKM = kernel[0, 5].compute(updatedJulianTime)
        jupiterPositionM = jupiterPositionKM * 1000

        saturnPositionKM = kernel[0, 6].compute(updatedJulianTime)
        saturnPositionM = saturnPositionKM * 1000

        uranusPositionKM = kernel[0, 7].compute(updatedJulianTime)
        uranusPositionM = uranusPositionKM * 1000

        neptunePositionKM = kernel[0, 8].compute(updatedJulianTime)
        neptunePositionM = neptunePositionKM * 1000

        xEarthTrajectory.append(earthPositionM[0])
        yEarthTrajectory.append(earthPositionM[1])
        zEarthTrajectory.append(earthPositionM[2])

        xJupiterTrajectory.append(jupiterPositionM[0])
        yJupiterTrajectory.append(jupiterPositionM[1])
        zJupiterTrajectory.append(jupiterPositionM[2])

        xSaturnTrajectory.append(saturnPositionM[0])
        ySaturnTrajectory.append(saturnPositionM[1])
        zSaturnTrajectory.append(saturnPositionM[2])

        xUranusTrajectory.append(uranusPositionM[0])
        yUranusTrajectory.append(uranusPositionM[1])
        zUranusTrajectory.append(uranusPositionM[2])

        xNeptuneTrajectory.append(neptunePositionM[0])
        yNeptuneTrajectory.append(neptunePositionM[1])
        zNeptuneTrajectory.append(neptunePositionM[2])

        xTrajectory.append(xCraft[0])
        yTrajectory.append(yCraft[0])
        zTrajectory.append(zCraft[0])

        tFirstLeg += timestep

    # ####### Leg 2 #########
    tSecondLeg = (tof_leg1) * 86400
    tEndSecondLeg = (tof_leg1+tof_leg2)*86400

    vxCraft = xVel_leg2
    vyCraft = yVel_leg2
    vzCraft = zVel_leg2

    while tSecondLeg < tEndSecondLeg:
        k = computeRk4(xCraft, yCraft, zCraft, xSun, ySun, zSun, vxCraft, vyCraft, vzCraft, mass_sun, timestep)

        xCraft += (timestep / 6) * (k[0] + 2 * k[6] + 2 * k[12] + k[18])
        yCraft += (timestep / 6) * (k[1] + 2 * k[7] + 2 * k[13] + k[19])
        zCraft += (timestep / 6) * (k[2] + 2 * k[8] + 2 * k[14] + k[20])

        vxCraft += (timestep / 6) * (k[3] + 2 * k[9] + 2 * k[15] + k[21])
        vyCraft += (timestep / 6) * (k[4] + 2 * k[10] + 2 * k[16] + k[22])
        vzCraft += (timestep / 6) * (k[5] + 2 * k[11] + 2 * k[17] + k[23])

        updatedJulianTime = julianTime + (tSecondLeg / 86400.0)
        earthPositionKM = kernel[0, 3].compute(updatedJulianTime)
        earthPositionM = earthPositionKM * 1000

        jupiterPositionKM = kernel[0, 5].compute(updatedJulianTime)
        jupiterPositionM = jupiterPositionKM * 1000

        saturnPositionKM = kernel[0, 6].compute(updatedJulianTime)
        saturnPositionM = saturnPositionKM * 1000

        uranusPositionKM = kernel[0, 7].compute(updatedJulianTime)
        uranusPositionM = uranusPositionKM * 1000

        neptunePositionKM = kernel[0, 8].compute(updatedJulianTime)
        neptunePositionM = neptunePositionKM * 1000

        xEarthTrajectory.append(earthPositionM[0])
        yEarthTrajectory.append(earthPositionM[1])
        zEarthTrajectory.append(earthPositionM[2])

        xJupiterTrajectory.append(jupiterPositionM[0])
        yJupiterTrajectory.append(jupiterPositionM[1])
        zJupiterTrajectory.append(jupiterPositionM[2])

        xSaturnTrajectory.append(saturnPositionM[0])
        ySaturnTrajectory.append(saturnPositionM[1])
        zSaturnTrajectory.append(saturnPositionM[2])

        xUranusTrajectory.append(uranusPositionM[0])
        yUranusTrajectory.append(uranusPositionM[1])
        zUranusTrajectory.append(uranusPositionM[2])

        xNeptuneTrajectory.append(neptunePositionM[0])
        yNeptuneTrajectory.append(neptunePositionM[1])
        zNeptuneTrajectory.append(neptunePositionM[2])

        xTrajectory.append(xCraft[0])
        yTrajectory.append(yCraft[0])
        zTrajectory.append(zCraft[0])

        tSecondLeg += timestep

    # ###### Leg 3 ########
    tThirdLeg = (tof_leg1+tof_leg2)*86400
    tEndThirdLeg = (tof_leg1+tof_leg2+tof_leg3)*86400

    vxCraft = xVel_leg3
    vyCraft = yVel_leg3
    vzCraft = zVel_leg3

    while tThirdLeg < tEndThirdLeg:
        k = computeRk4(xCraft, yCraft, zCraft, xSun, ySun, zSun, vxCraft, vyCraft, vzCraft, mass_sun, timestep)

        xCraft += (timestep / 6) * (k[0] + 2 * k[6] + 2 * k[12] + k[18])
        yCraft += (timestep / 6) * (k[1] + 2 * k[7] + 2 * k[13] + k[19])
        zCraft += (timestep / 6) * (k[2] + 2 * k[8] + 2 * k[14] + k[20])

        vxCraft += (timestep / 6) * (k[3] + 2 * k[9] + 2 * k[15] + k[21])
        vyCraft += (timestep / 6) * (k[4] + 2 * k[10] + 2 * k[16] + k[22])
        vzCraft += (timestep / 6) * (k[5] + 2 * k[11] + 2 * k[17] + k[23])

        updatedJulianTime = julianTime + (tThirdLeg / 86400.0)
        earthPositionKM = kernel[0, 3].compute(updatedJulianTime)
        earthPositionM = earthPositionKM * 1000

        jupiterPositionKM = kernel[0, 5].compute(updatedJulianTime)
        jupiterPositionM = jupiterPositionKM * 1000

        saturnPositionKM = kernel[0, 6].compute(updatedJulianTime)
        saturnPositionM = saturnPositionKM * 1000

        uranusPositionKM = kernel[0, 7].compute(updatedJulianTime)
        uranusPositionM = uranusPositionKM * 1000

        neptunePositionKM = kernel[0, 8].compute(updatedJulianTime)
        neptunePositionM = neptunePositionKM * 1000

        xEarthTrajectory.append(earthPositionM[0])
        yEarthTrajectory.append(earthPositionM[1])
        zEarthTrajectory.append(earthPositionM[2])

        xJupiterTrajectory.append(jupiterPositionM[0])
        yJupiterTrajectory.append(jupiterPositionM[1])
        zJupiterTrajectory.append(jupiterPositionM[2])

        xSaturnTrajectory.append(saturnPositionM[0])
        ySaturnTrajectory.append(saturnPositionM[1])
        zSaturnTrajectory.append(saturnPositionM[2])

        xUranusTrajectory.append(uranusPositionM[0])
        yUranusTrajectory.append(uranusPositionM[1])
        zUranusTrajectory.append(uranusPositionM[2])

        xNeptuneTrajectory.append(neptunePositionM[0])
        yNeptuneTrajectory.append(neptunePositionM[1])
        zNeptuneTrajectory.append(neptunePositionM[2])

        xTrajectory.append(xCraft[0])
        yTrajectory.append(yCraft[0])
        zTrajectory.append(zCraft[0])

        tThirdLeg += timestep

    # ##### Leg 4 ########

    tFourthLeg = (tof_leg1 + tof_leg2+tof_leg3) * 86400
    tEndFourthLeg = (tof_leg1 + tof_leg2 + tof_leg3+tof_leg4) * 86400

    vxCraft = xVel_leg4
    vyCraft = yVel_leg4
    vzCraft = zVel_leg4

    while tFourthLeg < tEndFourthLeg:
        k = computeRk4(xCraft, yCraft, zCraft, xSun, ySun, zSun, vxCraft, vyCraft, vzCraft, mass_sun, timestep)

        xCraft += (timestep / 6) * (k[0] + 2 * k[6] + 2 * k[12] + k[18])
        yCraft += (timestep / 6) * (k[1] + 2 * k[7] + 2 * k[13] + k[19])
        zCraft += (timestep / 6) * (k[2] + 2 * k[8] + 2 * k[14] + k[20])

        vxCraft += (timestep / 6) * (k[3] + 2 * k[9] + 2 * k[15] + k[21])
        vyCraft += (timestep / 6) * (k[4] + 2 * k[10] + 2 * k[16] + k[22])
        vzCraft += (timestep / 6) * (k[5] + 2 * k[11] + 2 * k[17] + k[23])

        updatedJulianTime = julianTime + (tFourthLeg / 86400.0)
        earthPositionKM = kernel[0, 3].compute(updatedJulianTime)
        earthPositionM = earthPositionKM * 1000

        jupiterPositionKM = kernel[0, 5].compute(updatedJulianTime)
        jupiterPositionM = jupiterPositionKM * 1000

        saturnPositionKM = kernel[0, 6].compute(updatedJulianTime)
        saturnPositionM = saturnPositionKM * 1000

        uranusPositionKM = kernel[0, 7].compute(updatedJulianTime)
        uranusPositionM = uranusPositionKM * 1000

        neptunePositionKM = kernel[0, 8].compute(updatedJulianTime)
        neptunePositionM = neptunePositionKM * 1000

        xEarthTrajectory.append(earthPositionM[0])
        yEarthTrajectory.append(earthPositionM[1])
        zEarthTrajectory.append(earthPositionM[2])

        xJupiterTrajectory.append(jupiterPositionM[0])
        yJupiterTrajectory.append(jupiterPositionM[1])
        zJupiterTrajectory.append(jupiterPositionM[2])

        xSaturnTrajectory.append(saturnPositionM[0])
        ySaturnTrajectory.append(saturnPositionM[1])
        zSaturnTrajectory.append(saturnPositionM[2])

        xUranusTrajectory.append(uranusPositionM[0])
        yUranusTrajectory.append(uranusPositionM[1])
        zUranusTrajectory.append(uranusPositionM[2])

        xNeptuneTrajectory.append(neptunePositionM[0])
        yNeptuneTrajectory.append(neptunePositionM[1])
        zNeptuneTrajectory.append(neptunePositionM[2])

        xTrajectory.append(xCraft[0])
        yTrajectory.append(yCraft[0])
        zTrajectory.append(zCraft[0])

        tFourthLeg += timestep

    plt.ion()
    fig = plt.figure(figsize=(6, 6))


    plt.xlim([-5e12, 5e12])
    plt.ylim([-5e12, 5e12])


    for i in range(0, len(xTrajectory), 30):

        plt.scatter(xTrajectory[i], yTrajectory[i], s=1, c='black')
        plt.scatter(xEarthTrajectory[i], yEarthTrajectory[i], s=1, c='blue')
        plt.scatter(xJupiterTrajectory[i], yJupiterTrajectory[i], s=1, c='orange')
        plt.scatter(xSaturnTrajectory[i], ySaturnTrajectory[i], s=1, c='red')
        plt.scatter(xUranusTrajectory[i], yUranusTrajectory[i], s=1, c='magenta')
        plt.scatter(xNeptuneTrajectory[i], yNeptuneTrajectory[i], s=1, c='purple')

        plt.draw()
        plt.pause(0.00001)