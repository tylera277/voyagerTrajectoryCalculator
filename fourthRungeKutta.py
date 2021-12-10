import numpy as np

G = 6.67e-11

def compute_x_accel(x1,y1,z1,x2,y2,z2,center_mass):
    a_x = (G*center_mass*(x2-x1))/((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**(3/2)
    return a_x

def compute_y_accel(x1,y1,z1,x2,y2,z2,center_mass):
    a_y = (G*center_mass*(y2-y1))/((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**(3/2)
    return a_y

def compute_z_accel(x1,y1,z1,x2,y2,z2,center_mass):
    a_z = (G*center_mass*(z2-z1))/((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**(3/2)
    return a_z

def computeRk4(x1,y1,z1,x2,y2,z2,vx,vy,vz,center_mass,delT):

    k = np.zeros((24, 1))

    k1_x = vx
    k1_y = vy
    k1_z = vz

    k1_vx = compute_x_accel(x1,y1,z1,x2,y2,z2,center_mass)
    k1_vy = compute_y_accel(x1, y1, z1, x2, y2, z2, center_mass)
    k1_vz = compute_z_accel(x1, y1, z1, x2, y2, z2, center_mass)


    k2_x = vx + (delT / 2) * k1_vx
    k2_y = vy + (delT / 2) * k1_vy
    k2_z = vz + (delT / 2) * k1_vz

    k2_vx = compute_x_accel(x1+(delT/2)*k1_x, y1+(delT/2)*k1_y, z1+(delT/2)*k1_z,x2,y2,z2,center_mass)
    k2_vy = compute_y_accel(x1 + (delT / 2) * k1_x, y1 + (delT / 2) * k1_y, z1 + (delT / 2) * k1_z, x2, y2, z2,
                            center_mass)
    k2_vz = compute_z_accel(x1 + (delT / 2) * k1_x, y1 + (delT / 2) * k1_y, z1 + (delT / 2) * k1_z, x2, y2, z2,
                            center_mass)


    k3_x = vx + (delT / 2) * k2_vx
    k3_y = vy + (delT / 2) * k2_vy
    k3_z = vz + (delT / 2) * k2_vz

    k3_vx = compute_x_accel(x1+(delT/2)*k2_x, y1+(delT/2)*k2_y, z1+(delT/2)*k2_z, x2,y2,z2,center_mass)
    k3_vy = compute_y_accel(x1 + (delT / 2) * k2_x, y1 + (delT / 2) * k2_y, z1 + (delT / 2) * k2_z, x2, y2, z2,
                            center_mass)
    k3_vz = compute_z_accel(x1 + (delT / 2) * k2_x, y1 + (delT / 2) * k2_y, z1 + (delT / 2) * k2_z, x2, y2, z2,
                            center_mass)


    k4_x = vx + delT*k3_vx
    k4_y = vy + delT*k3_vy
    k4_z = vz + delT * k3_vz

    k4_vx = compute_x_accel(x1+delT*k3_x, y1+delT*k3_y, z1+delT*k3_z,x2,y2,z2,center_mass)
    k4_vy = compute_y_accel(x1 + delT * k3_x, y1 + delT * k3_y, z1 + delT * k3_z, x2, y2, z2, center_mass)
    k4_vz = compute_z_accel(x1 + delT * k3_x, y1 + delT * k3_y, z1 + delT * k3_z, x2, y2, z2, center_mass)

    k[0] = k1_x
    k[1] = k1_y
    k[2] = k1_z
    k[3] = k1_vx
    k[4] = k1_vy
    k[5] = k1_vz

    k[6] = k2_x
    k[7] = k2_y
    k[8] = k2_z
    k[9] = k2_vx
    k[10] = k2_vy
    k[11] = k2_vz

    k[12] = k3_x
    k[13] = k3_y
    k[14] = k3_z
    k[15] = k3_vx
    k[16] = k3_vy
    k[17] = k3_vz

    k[18] = k4_x
    k[19] = k4_y
    k[20] = k4_z
    k[21] = k4_vx
    k[22] = k4_vy
    k[23] = k4_vz
    return k
