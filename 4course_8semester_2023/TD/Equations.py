import numpy as np
from numpy import sqrt
from math import acos, isnan, sin, cos

from scipy.integrate import ode, solve_ivp
from scipy.optimize import root, minimize
from scipy.interpolate import interp1d    #   https://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html
#import cma

from Orbit import orbit, orbit_full, Topographic_coordinates, Topographic_coordinates_3D, Keplerian_elements

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



ts1 = []
ys1 = []

def fout1(t, y):
    ts1.append(t)
    ys1.append(list(y.copy()))

def ODE(t0, tk, Y0, fun, par):
    ts1.clear()
    ys1.clear()

    r = ode(fun)

    r.set_integrator('dop853', rtol=0, atol=1.e-9, nsteps=1000)

    r.set_solout(fout1)
    r.set_initial_value(Y0, t0).set_f_params(par)
    y1 = r.integrate(tk)

    Y1 = np.array(ys1)

    return ts1, Y1


# модель движения
def fgrav(x,y,par):
        
    f  = np.zeros(6)
    
    x2=y[0]*y[0]
    y2=y[1]*y[1]
    z2=y[2]*y[2]
    r2=x2+y2+z2;
    r=sqrt(r2);
    r1=1/r;
    r3=r1/r2;
    
    f[0]=y[3];
    f[1]=y[4];
    f[2]=y[5];
    
    f[3]=-y[0]*r3;
    f[4]=-y[1]*r3;
    f[5]=-y[2]*r3; 
    
    return f

# считывание таблицы из файла, интерполяция, возврат 2 значений
def Atmosphere(h):
              
    ATM = np.genfromtxt('EARTH.ATM',delimiter = None,skip_header = 1)

    if h < ATM[0,0]:
        return ATM[0,2], ATM[0,3]
    elif h > ATM[-1,0]:
        return ATM[-1,2], ATM[-1,3]
    else:
        fRO = interp1d(ATM[:,0], ATM[:,2], kind='cubic')
        fA = interp1d(ATM[:,0], ATM[:,3], kind='cubic')
    
        ro = float(fRO(h))
        a  = float(fA(h))
    
        return ro, a

# -//- для Cx
def Cx(M):
    
    CXbaz = np.genfromtxt('cx_M.txt',delimiter = None,skip_header = 0)
    fCx = interp1d(CXbaz[:,0], CXbaz[:,1], kind='cubic')      
    
    if M < CXbaz[0,0]:
        return CXbaz[0,1]
    elif M > CXbaz[-1,0]:
        return CXbaz[-1,1]
    else:
        return float(fCx(M))

# функция для мешка
def ALFA(M, M1, M2, alfaM):
    
    if M < M1 or M2 < M:
        ALFA = 0
    else:
        ALFA = -alfaM/2.0*(np.cos(2.0*np.pi*(M-M1)/(M2-M1))-1.0)
    
    return ALFA

# управление, связанное с обкаткой по тангажу
def TET(tet0, dtet, t):
    
    TET = tet0 + dtet*t
    
    return TET

# уравнения, включающие тягу
# в ПАР все необходимые параметры для управления (ступени, ракета)
def Equations_Thirst(x, y, par):
    # par[] = S, P,c, stage, fi, A,  M1, M2, alfaM, tet0, dtet
    #         0  1 2  3      4   5   6   7   8      6     7
    raddeg = 180/np.pi
    om_Earth = 4.17807462*10**(-3)/raddeg        # deg/sec /raddeg = rad/sec
    R_Earth = 6371		     # km
    g0 = 9.80665/1000        # kм/с2

    R = np.zeros(3)
    R[0] = y[0]
    R[1] = y[1]+R_Earth
    R[2] = y[2]
    r = np.sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2])
    h = r - R_Earth
    G = g0*(R_Earth/(R_Earth+h))**2
    
    Vatm = np.zeros(3)
    Votn = np.zeros(3)
    Votn[0] =  om_Earth*np.cos(par[4])*np.cos(par[5])
    Votn[1] =  om_Earth*np.sin(par[4])
    Votn[2] = -om_Earth*np.cos(par[4])*np.sin(par[5])

    Vatm = np.cross(Votn,R)
    Votn[:] = y[3:6] - Vatm[:]

    v = np.sqrt(Votn[0]*Votn[0]+Votn[1]*Votn[1]+Votn[2]*Votn[2])
    ro, a = Atmosphere(h*1000)
    Max = v*1000/a
    X = 0.5*par[0]*Cx(Max)*ro*v*v*1000000
    
    stage = int(par[3])
    if stage == 1:
        tet = np.pi/2 + ALFA(Max, par[6], par[7], par[8])
    else:
        tet = TET(par[6], par[7], x)
    
    Ox = -G*y[0]/r
    Oy = -G*(y[1]+R_Earth)/r
    Oz = -G*y[2]/r 

    f = np.zeros(7)

    PX = (par[1] - X)/y[6]/1000 # km/sec2
    
#    print("time: %8.4f sec, mah: %8.4f, n: %8.4f" %(x, Max, PX/abs(Oy)))

# модель для вычисления дифференциальных уравнений
    f[0] = y[3]                              #   x
    f[1] = y[4]                              #   y
    f[2] = y[5]                              #   z
    
    f[3] = Ox + PX*np.cos(tet)               #   vx
    f[4] = Oy + PX*np.sin(tet)               #   vy
    f[5] = Oz                                #   vz

    f[6] = -par[1]/par[2]                    #   m
    
    return f

# "как это всё полетит с таким управлением"
# получает описание ракеты, описание управления
def Trajectory (u, Nstage, P, c, M, S, lam, fi, A):
# u = [[tact1, M1, M2, alfaM], [tact2, tet0, dtet ], ...]  len(u) = Nstage
# T = [ , , , ]          len(T) = Nstage
# I = [ , , , ]          len(I) = Nstage
# M = [ начальная, сбрасываяемая масса после работы 1 ступени, ..., сбрасываяемая масса после работы N ступени ]
# S = [ , , , ]          len(S) = Nstage
# если tet0 не оптимизируется 

    # константы

    raddeg = 180/np.pi
    
    om_Earth = 4.17807462*10**(-3)/raddeg        # deg/sec /raddeg = rad/sec
    R_Earth = 6371;								 # km 
    tt=[]
    TR=[[],[],[],[],[],[],[]]
    step = [0]
    # Начальные условия
    V0 = np.zeros(3)
    V0[0] = om_Earth*R_Earth*np.cos(fi)*np.sin(A)  # km/sec
    V0[1] = 0.0                                    # km/sec
    V0[2] = om_Earth*R_Earth*np.cos(fi)*np.cos(A)  # km/sec
    t0 = 0
    rvm = [0,0,0,V0[0],V0[1],V0[2],M[0]]
    for i in range(Nstage):
        par = [S[i], P[i], c[i], i+1, fi, A]
        par.extend(u[i][1:])

    # интегрирование
        t, Tr = ODE(t0,t0+u[i][0], rvm, Equations_Thirst, par)

        t0 = t[-1]
        rvm = Tr[-1]
        rvm[6] -= M[i+1]
        step.append(step[-1]+len(t))
        tt.extend(t)
        for j in range(7):
            TR[j].extend(Tr[:,j])
            
    return rvm,tt,np.array(TR),step
    # рвм - 
    # остальное - траектория (тт - время, скорость, масса)
    # степ - 

# rvm - результат интегрирования
# H - высота желаемая
def boundary_condition(rvm, H):
    
    R_Earth = 6371		     # km 
    fM_Earth = 398600.4415   # km^3/sec^2

    r = np.sqrt(rvm[0]*rvm[0] + (rvm[1]+R_Earth)*(rvm[1]+R_Earth) + rvm[2]*rvm[2])
    h = r - R_Earth    
    v  = np.sqrt(rvm[3]*rvm[3]+rvm[4]*rvm[4]+rvm[5]*rvm[5]) 
    rv = np.sqrt(rvm[3]*rvm[0]+rvm[4]*(rvm[1]+R_Earth)+rvm[5]*rvm[2])
    
    ancl = acos(rv/(r*v))
    
#    print('ancl = ', )
    
    V_circ = np.sqrt(fM_Earth/r)
    
    error = np.zeros(4)
    
#    print('h = ', h)
    
    error[0] =  h - H
    error[1] =  ancl - np.pi/2 # много слов...
    error[2] =  v - V_circ    # контроль скорости на круговой обрите
    error[3] = -rvm[6]    # масса
    
    if isnan(error[1]):
        error[1] = 100
    
    print ("%8.4e km %8.4e deg %8.4e km/s Mk = %8.2f kg" %(error[0], ancl*180/np.pi-90, error[2], -error[3]))  
    
    return error

# минимизатор, распределяет управление, собирает начальные данные из исходных
# возвращают подходящую функцию минимизации
def funMin1(x, Nstage, P, c, M, S, lam, fi, A, H, u):

    UX = []
    for i in range(Nstage):
        UX.append(u[i])
            
    k = 0  
    for i in range(Nstage):
        for j in range(len(u[i])):
            UX[i][j] = x[k]*UX[i][j]
            k += 1

    # запрос траектории с управлением
    rvm,tt,TR,step = Trajectory (UX, Nstage, P, c, M, S, lam, fi, A)
    
    error = boundary_condition(rvm, H)
    
    err = np.dot(error[0:3], error[0:3])
    
    if min(TR[6,:]) < 0:
        err = err + 1000
    
    print ("%10.5e %10.4f" %(err, -error[3]))
    
    return error[3] + err

def funMin2(x, Nstage, P, c, M, S, lam, fi, A, H):

    u = []
    for i in range(Nstage):
        if i == 0:
            u.append(x[0:4])
        else:
            u.append(x[4+3*(i-1):4+3*i])       

    rvm,tt,TR,step = Trajectory (u, Nstage, P, c, M, S, lam, fi, A)
    
    error = boundary_condition(rvm, H)
    
    err = np.dot(error[0:3], error[0:3])*100
    
    if min(TR[6,:]) < 0:
        err = err + 1000
    
#    print ("%10.5e %10.4f" %(err, -error[3]))
    
    return error[3] + err

def funMin3(x, Nstage, P, c, M, S, lam, fi, A, H, u):
#   t0, alfa, ti, dfi

    UX = []
    for i in range(Nstage):
        UX.append([0]*len(u[i]))

    for i in range(Nstage):
        if i == 0:
            UX[i][0] = x[0]*u[i][0]
            UX[i][1] = u[i][1]
            UX[i][2] = u[i][2]
            UX[i][3] = x[1]*u[i][3]
        elif i == 1:
            UX[i][0] = x[2]*u[i][0]
            UX[i][1] = np.pi/2
            UX[i][2] = x[3]*u[i][2]
        else:
            UX[i][0] = x[4+2*(i-2)]*u[i][0]
            UX[i][1] = UX[i-1][1] + UX[i-1][2]*UX[i-1][0]
            UX[i][2] = x[5+2*(i-2)]*u[i][2]           
            
#    print(UX)
    rvm,tt,TR,step = Trajectory (UX, Nstage, P, c, M, S, lam, fi, A)
    
    error = boundary_condition(rvm, H)
    
    err = error[0]*error[0] + error[1]*error[1] + error[2]*error[2]*100
                
    return err*100 + error[3]


# решатель
def Solution (u, Nstage, P, c, M, S, lam, fi, A, H):
    
    Nu = 0
    for i in range(Nstage):
        Nu += len(u[i])

    x = [1]*Nu
           
#    res = cma.fmin(funMin1, x, 0.01, options={'maxfevals': 1e2, 'popsize':15}, args=(Nstage, P, c, M, S, lam, fi, A, H, u,));
#    x = res[0]
    
    k = 0  
    for i in range(Nstage):
        for j in range(len(u[i])):
            x[k] = x[k]*u[i][j]
            k += 1
        
    x = minimize(funMin2, x, args=(Nstage, P, c, M, S, lam, fi, A, H), method='SLSQP', options={'maxiter':30})
    print(x)
    
    for i in range(Nstage):
        if i == 0:
            u[0:4] = x.x[0:4]
        else:
            u[4+3*(i-1):4+3*i] = x.x[4+3*(i-1):4+3*i]       

    rvm,tt,TR,step = Trajectory (u, Nstage, P, c, M, S, lam, fi, A)
    error = boundary_condition(rvm, H)
    
    return x, u, rvm,tt,TR,step, error


# 
def Solution2 (u, Nstage, P, c, M, S, lam, fi, A, H):

    rvm,tt,TR,step = Trajectory (u, Nstage, P, c, M, S, lam, fi, A)
    error = boundary_condition(rvm, H)
    
    x = [1]*Nstage*2
     
#    print(x, end = ' ')
                 
#    res = cma.fmin(funMin3, x, 0.001, options={'maxfevals': 1e3, 'popsize':10}, args=(Nstage, P, c, M, S, lam, fi, A, H, u,));
#    x = res[0]    
    
    x = minimize(funMin3, x, args=(Nstage, P, c, M, S, lam, fi, A, H, u), method='Powell', options={'maxiter':5})
    print(x)
 
    for i in range(Nstage):
        if i == 0:
            u[i][0] = x.x[0]*u[i][0]
            u[i][3] = x.x[1]*u[i][3]
        elif i == 1:
            u[i][0] = x.x[2]*u[i][0]
            u[i][1] = np.pi/2
            u[i][2] = x.x[3]*u[i][2]
        else:
            u[i][0] = x.x[4+2*(i-2)]*u[i][0]
            u[i][1] = u[i-1][1] + u[i-1][2]*u[i-1][0]
            u[i][2] = x.x[5+2*(i-2)]*u[i][2]           
     
    rvm,tt,TR,step = Trajectory (u, Nstage, P, c, M, S, lam, fi, A)
    error = boundary_condition(rvm, H)
    
    return x, u, rvm,tt,TR,step, error


def Trajectory_orbit(u, Nstage, P, c, M, S, lam, fi, A, H, inc, W, OM):
 
    # константы
    
    raddeg = 180/np.pi
    
    fM_Earth = 398600.4415 # km^3/sec^2
    fM_Moon = 4902.801 # km^3/sec^2
    fM_Sun = 1.32712440018e11 # km^3/sec^2
    
    R_Earth = 6371;										# km 
    J2_Earth = 1.0826157e-3;
    c_Earth = J2_Earth*fM_Earth*R_Earth*R_Earth/2.0;	# km^5/sec^2
    om_Earth = 4.17807462*10**(-3) *86400/raddeg # deg/sec *86400/raddeg = rad/day
    
    UnitfM = fM_Earth
    UnitR = R_Earth                # km
    UnitV = np.sqrt(UnitfM/UnitR)  # km/sec
    UnitT = (UnitR/UnitV)/86400    # day   
 
    # параметры орбиты

    hp = H
    ha = H
    i = inc/raddeg
    w = W/raddeg
    om = OM/raddeg
    tet = 0/raddeg
    
    rvm,tt,TR,step = Trajectory (u, Nstage, P, c, M, S, lam, fi, A)
    
    error = boundary_condition(rvm, H)
  
    n = len(TR[1,:])
    
    R0 = np.zeros(3)
    R0[0] = cos(fi)*cos(lam)*R_Earth
    R0[1] = cos(fi)*sin(lam)*R_Earth
    R0[2] = sin(fi)*R_Earth
    
    MA = [[ cos(-A), 0 , sin(-A)],
          [ 0      , 1 , 0      ],
          [-sin(-A), 0 , cos(-A)]]
    
    M = [
         [-sin(fi)*cos(lam),-sin(fi)*sin(lam),  cos(fi) ],
         [ cos(fi)*cos(lam), cos(fi)*sin(lam),  sin(fi) ],
         [        -sin(lam),         cos(lam),  0       ]
         ]
    
    MA = np.array(MA)
    M = np.array(M)
    M = np.linalg.inv(M)
    TRG = np.zeros((6,n))
    for i in range(n):
        TRG[0:3,i] = R0[0:3]+np.dot(M,np.dot(MA,TR[0:3,i]))
        TRG[3:6,i] = np.dot(M,np.dot(MA,TR[3:6,i]))
     
    TR[1,:] = TR[1,:] + R_Earth   
    
    # определение высоты
    h = np.zeros(n)
    for i in range(n):
        h[i] = np.sqrt(TR[0,i]*TR[0,i] + TR[1,i]*TR[1,i] + TR[2,i]*TR[2,i]) - R_Earth
    
    # определение дальности
    R0 = [0, R_Earth, 0]
    L = np.zeros(n)
    for i in range(n):
        Rk = TR[0:3,i]
        danc = acos( (R0[0]*Rk[0]+R0[1]*Rk[1]+R0[2]*Rk[2]) / (sqrt(R0[0]*R0[0]+R0[1]*R0[1]+R0[2]*R0[2])*sqrt(Rk[0]*Rk[0]+Rk[1]*Rk[1]+Rk[2]*Rk[2])) )
        L[i] = danc*R_Earth    
    
    
    print('Параметры траектории')
    for i in range(Nstage):
        print('В момент отделения ', str(i+1), 'ступени:')
        print('Положение:', TRG[0:3, step[i+1]-1], 'км')
        print('Скорость:',  TRG[3:6, step[i+1]-1], 'км/с')
        print('Масса:',     TR[6, step[i+1]-1],   'кг')
    
    Kep = Keplerian_elements(TRG[0:6,-1], 0, fM_Earth, 1)

    print('Параметры достигнутой орбиты:')
    print('Высоты перицентра и апоцентра:', [Kep[3]-R_Earth, Kep[4]-R_Earth], 'км')
    print('Наклонение, долгота восх. узла, аргумент перицентра', [Kep[5]*raddeg, Kep[6]*raddeg, Kep[6]*raddeg], 'град.')
 
#    Kep = [hp/UnitR, ha/UnitR, i, w, om, tet]
    Kep = [(Kep[3]-R_Earth)/UnitR, (Kep[4]-R_Earth)/UnitR, Kep[5], Kep[7], Kep[6], tet]

    # по истинной анамалии
    TET = np.linspace(0,2*np.pi,60)
    orb = orbit_full(Kep, TET, 1) 
  
    # графики 
    R = 12000
    
    plt.ion()

    fig = plt.figure(figsize=(8,6))  
    # ax = plt.gca(projection="3d")
    ax = fig.add_subplot(projection='3d')
    
    u, v = np.mgrid[0:2*np.pi:36j, 0:np.pi:18j]
    x = UnitR*np.cos(u)*np.sin(v)
    y = UnitR*np.sin(u)*np.sin(v)
    z = UnitR*np.cos(v)
    ax.plot_wireframe(x, y, z, color="gray")
    
    ax.plot(TRG[0,:],TRG[1,:],TRG[2,:], color='r', linewidth = 4)  
    ax.plot(orb[:,0]*UnitR,orb[:,1]*UnitR,orb[:,2]*UnitR, '--', color='k')  
    
    ax.legend()
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.set_zlim(-R, R)
    ax.set_xlabel('X, km')
    ax.set_ylabel('Y, km')
    ax.set_zlabel('Z, km')#, rotation=90)
    ax.azim = 225 
    plt.draw()
    
    plt.figure(figsize=(8,6))  
    plt.plot(L[:],h[:], color='r', linewidth = 4)
    plt.plot(L[:],[H]*len(TR[0,:]), '--', color='k')  

    plt.ylim(0, H+50)
    plt.xlabel('Дальность полёта (км)')
    plt.ylabel('Высота полёта (км)')    
    plt.draw()
    
    plt.figure(figsize=(8,6))  
    plt.plot(tt[:],TR[6,:], color='r', linewidth = 4)
    plt.xlabel('Время полёта (сек)')
    plt.ylabel('Масса (кг)')    
    plt.draw()    
    
    plt.ioff()
    plt.show()


def boundary_condition_ell(rvm, hp, ha, fi, lam, A):
    
    R_Earth = 6371		     # km 
    fM_Earth = 398600.4415   # km^3/sec^2
    
    R0 = np.zeros(3)
    R0[0] = cos(fi)*cos(lam)*R_Earth
    R0[1] = cos(fi)*sin(lam)*R_Earth
    R0[2] = sin(fi)*R_Earth
    
    MA = [[ cos(-A), 0 , sin(-A)],
          [ 0      , 1 , 0      ],
          [-sin(-A), 0 , cos(-A)]]
    
    M = [
         [-sin(fi)*cos(lam),-sin(fi)*sin(lam),  cos(fi) ],
         [ cos(fi)*cos(lam), cos(fi)*sin(lam),  sin(fi) ],
         [        -sin(lam),         cos(lam),  0       ]
         ]
    
    MA = np.array(MA)
    M = np.array(M)
    M = np.linalg.inv(M)
    RVG = np.zeros(6)
    RVG[0:3] = R0[0:3]+np.dot(M,rvm[0:3])
    RVG[3:6] = np.dot(M,rvm[3:6])   
    
    Kep = Keplerian_elements(RVG, 0, fM_Earth, 1)

    error = np.zeros(3)
    
    error[0] =  Kep[3] - hp - R_Earth
    error[1] =  Kep[4] - ha - R_Earth
    error[2] = -rvm[6]    
       
    print ("%8.4e km %8.4e km Mk = %8.2f kg" %(error[0], error[1], -error[2]))  
    
    return error

def funMin_ell(x, Nstage, P, c, M, S, lam, fi, A, hp, ha, u):
#   t0, alfa, ti, dfi

    UX = []
    for i in range(Nstage):
        UX.append([0]*len(u[i]))

    for i in range(Nstage):
        if i == 0:
            UX[i][0] = x[0]*u[i][0]
            UX[i][1] = u[i][1]
            UX[i][2] = u[i][2]
            UX[i][3] = x[1]*u[i][3]
        elif i == 1:
            UX[i][0] = x[2]*u[i][0]
            UX[i][1] = np.pi/2
            UX[i][2] = x[3]*u[i][2]
        else:
            UX[i][0] = x[4+2*(i-2)]*u[i][0]
            UX[i][1] = UX[i-1][1] + UX[i-1][2]*UX[i-1][0]
            UX[i][2] = x[5+2*(i-2)]*u[i][2]           
            
#    print(UX)
    rvm,tt,TR,step = Trajectory (UX, Nstage, P, c, M, S, lam, fi, A)
    
    error = boundary_condition_ell(rvm, hp, ha, fi, lam, A)    

    err = error[0]*error[0] + error[1]*error[1]
                
    return err*1000 + error[2]


def Solution_ELL (u, Nstage, P, c, M, S, lam, fi, A, hp, ha):

    rvm,tt,TR,step = Trajectory (u, Nstage, P, c, M, S, lam, fi, A)
    error = boundary_condition_ell(rvm, hp, ha, fi, lam, A)
    
    x = [1]*Nstage*2
                    
#    res = cma.fmin(funMin_ell, x, 0.01, options={'maxfevals': 1e3, 'popsize':10}, args=(Nstage, P, c, M, S, lam, fi, A, hp, ha, u,));
#    x = res[0]    
    
    x = minimize(funMin_ell, x, args=(Nstage, P, c, M, S, lam, fi, A, hp, ha, u), method='Powell', options={'maxiter':5})
    print(x)
 
    for i in range(Nstage):
        if i == 0:
            u[i][0] = x.x[0]*u[i][0]
            u[i][3] = x.x[1]*u[i][3]
        elif i == 1:
            u[i][0] = x.x[2]*u[i][0]
            u[i][1] = np.pi/2
            u[i][2] = x.x[3]*u[i][2]
        else:
            u[i][0] = x.x[4+2*(i-2)]*u[i][0]
            u[i][1] = u[i-1][1] + u[i-1][2]*u[i-1][0]
            u[i][2] = x.x[5+2*(i-2)]*u[i][2]           
     
    rvm,tt,TR,step = Trajectory (u, Nstage, P, c, M, S, lam, fi, A)
    error = boundary_condition_ell(rvm, hp, ha, fi, lam, A)
    
    return x, u, rvm,tt,TR,step, error

    
    
    
    
    
    