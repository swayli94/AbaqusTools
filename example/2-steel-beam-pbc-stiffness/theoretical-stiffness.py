'''
Calculate the theoretical stiffness of a beam

'''


if __name__ == '__main__':
    
    E  = 2.1E5   # Young's modulus (N/mm^2)
    PR = 0.3     # Poisson ratio
    
    length_x = 100
    length_y = 20
    length_z = 50

    aera = length_x*length_y
    
    Ix = length_x*length_y**3/12.0
    Iy = length_y*length_x**3/12.0
    
    Gz = E/(2*(1+PR))
    
    a = max(length_x, length_y)*0.5
    b = min(length_x, length_y)*0.5
    Ip = a*b**3*(16/3 -3.36*b/a*(1-b**4/a**4/12))
    
    
    C11 = E*aera
    C22 = Gz*Ip
    C33 = E*Ix
    C44 = E*Iy

    print('C11 = %.3E'%(C11))
    print('C22 = %.3E'%(C22))
    print('C33 = %.3E'%(C33))
    print('C44 = %.3E'%(C44))

