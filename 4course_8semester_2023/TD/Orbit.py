import numpy as np

def orbit(Kep, R):
    # Kep[0] -  hp, высота перицентра          -
    # Kep[1] -  ha, высота апоцентра           -
    # Kep[2] -   i, наклонение                rad
    # Kep[3] -   w, аргумент перигелия        rad
    # Kep[4] -  om, долгота восходящего узла  rad
    # Kep[5] - tet, истинная аномалия         rad
    # R      - средний радиус планеты
    # Dec[]  - x,y,x,Vx,Vy,Vz    
    
    Dec = np.zeros(6)
    
    a = (Kep[0] + Kep[1])/2. + R
    e = (Kep[1] - Kep[0]) / (Kep[0] + Kep[1] + 2.*R)

    r = a * (1. - e * e) / (1. + e * np.cos(Kep[5]))
    v = np.sqrt(2. / r - 1. / a)
    
    gam = np.arctan2(e * np.sin(Kep[5]), (1. + e * np.cos(Kep[5])))
    u = Kep[5] + Kep[3]
    ugam = u - gam
    
    Dec[0] = r * (np.cos(u) * np.cos(Kep[4]) - np.sin(u) * np.cos(Kep[2]) * np.sin(Kep[4]))
    Dec[1] = r * (np.cos(u) * np.sin(Kep[4]) + np.sin(u) * np.cos(Kep[2]) * np.cos(Kep[4]))    
    Dec[2] = r * (np.sin(u) * np.sin(Kep[2]))
    
    Dec[3] = v * (-np.sin(ugam) * np.cos(Kep[4]) - np.cos(ugam) * np.cos(Kep[2]) * np.sin(Kep[4]))
    Dec[4] = v * (-np.sin(ugam) * np.sin(Kep[4]) + np.cos(ugam) * np.cos(Kep[2]) * np.cos(Kep[4]))
    Dec[5] = v * (np.cos(ugam) * np.sin(Kep[2]))
    
    return Dec


def orbit_full(Kep, TET, R):

    n = len(TET)
    Dec = np.zeros((n, 6))

    for i in range(n):
        Kep[5] = TET[i]
        Dec[i,:] = orbit(Kep, R)
    
    return Dec


def Topographic_coordinates(R, om_e, dt, lam0, R0):

    r = np.sqrt(R[0]**2 + R[1]**2 + R[2]**2)    
    dlam = np.arctan2(R[0], R[1]) - np.arctan2(R0[0], R0[1])    
    lam = lam0 + dlam - om_e * dt   

    while lam > np.pi:
        lam -= 2 * np.pi
    while lam < -np.pi:
        lam += 2 * np.pi
    
    fi = np.arcsin(R[2] / r)
    
    return lam, fi


def Topographic_coordinates_3D(R, om_e, dt):
    
    fi = om_e * dt
    
    cos_fi = np.cos(fi)
    sin_fi = np.sin(fi)
    r = np.sqrt(R[0]**2 + R[1]**2 + R[2]**2)    
    
    x = (cos_fi*R[0] - sin_fi*R[1]) / r
    y = (sin_fi*R[0] + cos_fi*R[1]) / r
    z = R[2] / r
    
    return x,y,z


def Keplerian_elements (Dec, t, mu=1, par=0):

    # Dec[] - x,y,x,Vx,Vy,Vz UnitR, UnitV
    # t     - момент времени UnitT
    # Kep[ 0] - a, большая полуось				    UnitR
    # Kep[ 1] - p, фокальный параметр			    UnitR
    # Kep[ 2] - e, эксцентриситет				    -
    # Kep[ 3] - rp, радиус перицентра			    UnitR
    # Kep[ 4] - ra, радиус апоцентра				UnitR
    # Kep[ 5] - Incl, наклонение					rad
    # Kep[ 6] - Node, долгота восходящего узла	    rad
    # Kep[ 7] - Peri, аргумент перигелия			rad
    # Kep[ 8] - Nu, истинная аномалия			    rad		
    # Kep[ 9] - E, эксцентрическая аномалия 	    rad		
    # Kep[10] - M, средняя аномалия 				rad		
    # Kep[11] - tp, время прохождения перицентра	UnitT		

    c = np.zeros(3)
    f = np.zeros(3)

    r = np.sqrt(Dec[0]**2 + Dec[1]**2 + Dec[2]**2)
    V = np.sqrt(Dec[3]**2 + Dec[4]**2 + Dec[5]**2)

    h = V**2 - 2.*mu/r
    a = -mu/h

    c[0] = Dec[1]*Dec[5] - Dec[2]*Dec[4]
    c[1] = Dec[2]*Dec[3] - Dec[0]*Dec[5]
    c[2] = Dec[0]*Dec[4] - Dec[1]*Dec[3]

    Cxy = np.sqrt(c[0]**2 + c[1]**2)
    C = np.sqrt(Cxy**2 + c[2]**2)

    p = C**2 / mu

    f[0] = Dec[4]*c[2] - Dec[5]*c[1] - mu/r*Dec[0]
    f[1] = Dec[5]*c[0] - Dec[3]*c[2] - mu/r*Dec[1]
    f[2] = Dec[3]*c[1] - Dec[4]*c[0] - mu/r*Dec[2]

    F = np.sqrt(f[0]**2 + f[1]**2 + f[2]**2)

    e = F / mu

    rp = p / (1.+e)
    ra = p / (1.-e)

    Incl = np.arccos(c[2] / C)

    if np.abs(c[0]) < 1.e-10 and np.abs(c[1]) < 1.e-10:
        Node = 0
    elif c[0] > 0:
        Node = np.arccos(-c[1]/Cxy)
    elif c[0] < 0:
        Node = 2 * np.pi - np.arccos(-c[1]/Cxy)
    elif c[0] == 0 and c[1] < 0:
        Node = 0
    elif c[0] == 0 and c[1] > 0:
        Node = np.pi

    l = np.array([np.cos(Node), np.sin(Node), 0.])

    cosPeri = (l[0]*f[0] + l[1]*f[1] + l[2]*f[2]) / F
    if abs(cosPeri) > 1:
        cosPeri = cosPeri / abs(cosPeri)
    
    Peri = 0.
    if f[2]/F > 1.e-8:
        Peri = np.arccos(cosPeri)
    elif f[2]/F < 1.e-8:
        Peri = 2 * np.pi - np.arccos(cosPeri)
    else:
        if l[0]*f[0] > 0:
            Peri = 0.
        elif l[0]*f[0] < 0:
            Peri = np.pi

    cosNu = (Dec[0]*f[0] + Dec[1]*f[1] + Dec[2]*f[2]) / (r*F)
    if abs(cosNu) > 1:
        cosNu = cosNu / abs(cosNu)
        
    frc = (f[1]*Dec[0] - f[2]*Dec[1]) * c[0] \
        + (f[2]*Dec[0] - f[0]*Dec[2]) * c[1] \
        + (f[0]*Dec[1] - f[1]*Dec[0]) * c[2]

    if frc > 0:
        Nu = np.arccos(cosNu)
    else:
        Nu = 2 * np.pi - np.arccos(cosNu)

    tanNu2 = np.tan(Nu/2.)

    if h < 0:
        E = 2. * np.arctan2(np.sqrt(1.-e) * tanNu2, np.sqrt(1.+e))
        n = np.sqrt(mu / (a**3))
        M = E - e*np.sin(E)
    elif h == 0:
        E = tanNu2
        n = 2. * np.sqrt(mu / (p**3))   # опечатка в оригинале (3.13)
        M = E + 1./3. * E**3
    elif h > 0:                         # Здесь возможны ошибки в переменных E и n
        E = 2. * np.arctanh(complex(np.sqrt((e-1.) / (e+1.)) * tanNu2, 0.)).real
        n = np.sqrt(mu / (-a**3))
        M = e * np.sinh(E) - E

    tp = t - M/n

    if bool(par):
        Kep = np.array([a, p, e, rp, ra, Incl, Node, Peri, Nu, E, M, tp])
    else:
        Kep = np.array([a, e, Incl, Peri, Node, M])
		
    return Kep
