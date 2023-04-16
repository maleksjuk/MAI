import numpy as np
from math import isnan

from scipy.integrate import ode
from scipy.optimize import minimize
from scipy.interpolate import interp1d
#   https://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html
# import cma

from Orbit import orbit_full, Keplerian_elements

from mpl_toolkits.mplot3d import Axes3D, axes3d
import matplotlib.pyplot as plt

from Base import *


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


def fgrav(x, y, par):
    f = np.zeros(6)

    x2 = y[0] * y[0]
    y2 = y[1] * y[1]
    z2 = y[2] * y[2]
    r2 = x2 + y2 + z2
    r = np.sqrt(r2)
    r1 = 1 / r
    r3 = r1 / r2

    f[0] = y[3]
    f[1] = y[4]
    f[2] = y[5]

    f[3] = -y[0] * r3
    f[4] = -y[1] * r3
    f[5] = -y[2] * r3

    return f


def Atmosphere(h):

    ATM = np.genfromtxt('EARTH.ATM', delimiter=None, skip_header=1)

    if h < ATM[0, 0]:
        return ATM[0, 2], ATM[0, 3]
    elif h > ATM[-1, 0]:
        return ATM[-1, 2], ATM[-1, 3]
    else:
        fRO = interp1d(ATM[:, 0], ATM[:, 2], kind='cubic')
        fA = interp1d(ATM[:, 0], ATM[:, 3], kind='cubic')
        ro = float(fRO(h))
        a = float(fA(h))
        return ro, a


def Cx(M):

    CXbaz = np.genfromtxt('cx_M.txt', delimiter=None, skip_header=0)
    fCx = interp1d(CXbaz[:, 0], CXbaz[:, 1], kind='cubic')

    if M < CXbaz[0, 0]:
        return CXbaz[0, 1]
    elif M > CXbaz[-1, 0]:
        return CXbaz[-1, 1]
    else:
        return float(fCx(M))


def ALFA(M, M1, M2, alfaM):
    """
    Функция для мешка
    """
    if M < M1 or M2 < M:
        ALFA = 0
    else:
        ALFA = -alfaM/2.0 * (np.cos(2.0 * np.pi * (M-M1) / (M2-M1)) - 1.0)
    return ALFA


def TET(tet0, dtet, t):
    """
    Управление, связанное с обкаткой по тангажу
    """
    TET = tet0 + dtet*t
    return TET


def Equations_Thirst(x, y, par):
    """
    Уравнения, включающие тягу. В par все необходимые параметры для управления (ступени, ракета)
    
    par[] = S, P,c, stage, fi, A,  M1, M2, alfaM, tet0, dtet
            0  1 2  3      4   5   6   7   8      6     7
    """
    g0 = 9.80665/1000        # kм/с2

    R = np.zeros(3)
    R[0] = y[0]
    R[1] = y[1] + Constants.R_Earth
    R[2] = y[2]
    r = np.sqrt(np.sum(R ** 2))
    h = r - Constants.R_Earth
    G = g0 * (Constants.R_Earth / (Constants.R_Earth + h))**2

    Vatm = np.zeros(3)
    Votn = np.zeros(3)
    Votn[0] = Constants.om_Earth * np.cos(par[4]) * np.cos(par[5])
    Votn[1] = Constants.om_Earth * np.sin(par[4])
    Votn[2] = -Constants.om_Earth * np.cos(par[4]) * np.sin(par[5])

    Vatm = np.cross(Votn, R)
    Votn[:] = y[3:6] - Vatm[:]

    v = np.sqrt(np.sum(Votn ** 2))
    ro, a = Atmosphere(h * 1000)
    Max = v * 1000 / a
    X = 0.5 * par[0] * Cx(Max) * ro * v * v * 10**6

    stage = int(par[3])
    if stage == 1:
        tet = np.pi/2 + ALFA(Max, par[6], par[7], par[8])
    else:
        tet = TET(par[6], par[7], x)

    Ox = -G * y[0] / r
    Oy = -G * (y[1] + Constants.R_Earth) / r
    Oz = -G * y[2] / r

    f = np.zeros(7)

    PX = (par[1] - X) / y[6] / 1000  # km/sec2

#    print("time: %8.4f sec, mah: %8.4f, n: %8.4f" %(x, Max, PX/abs(Oy)))

    f[0] = y[3]                              # x
    f[1] = y[4]                              # y
    f[2] = y[5]                              # z

    f[3] = Ox + PX*np.cos(tet)               # vx
    f[4] = Oy + PX*np.sin(tet)               # vy
    f[5] = Oz                                # vz

    f[6] = -par[1]/par[2]                    # m

    return f


def Trajectory(u, rocket: Rocket):
    """
    Получает описание ракеты и управления

    # u = [[tact1, M1, M2, alfaM], [tact2, tet0, dtet ], ...]  len(u) = Nstage
    # T = [ , , , ]          len(T) = Nstage
    # I = [ , , , ]          len(I) = Nstage
    # M = [ начальная, сбрасываяемая масса после работы 1 ступени, ...,
    #   сбрасываяемая масса после работы N ступени ]
    # S = [ , , , ]          len(S) = Nstage
    # если tet0 не оптимизируется
    """

    tt = []
    TR = [[], [], [], [], [], [], []]
    step = [0]
    # Начальные условия
    V0 = np.zeros(3)
    V0[0] = Constants.om_Earth * Constants.R_Earth * np.cos(rocket.fi) * np.sin(rocket.A)  # km/sec
    V0[1] = 0.0                                          # km/sec
    V0[2] = Constants.om_Earth * Constants.R_Earth * np.cos(rocket.fi) * np.cos(rocket.A)  # km/sec
    t0 = 0
    rvm = [0, 0, 0, V0[0], V0[1], V0[2], rocket.M[0]]
    for i in range(rocket.Nstage):
        par = [rocket.S[i], rocket.P[i], rocket.c[i], i+1, rocket.fi, rocket.A]
        par.extend(u[i][1:])

        # интегрирование
        t, Tr = ODE(t0, t0 + u[i][0], rvm, Equations_Thirst, par)

        t0 = t[-1]
        rvm = Tr[-1]
        rvm[6] -= rocket.M[i+1]
        step.append(step[-1] + len(t))
        tt.extend(t)
        for j in range(7):
            TR[j].extend(Tr[:, j])

    return rvm, tt, np.array(TR), step


def boundary_condition(rvm, H):
    """
    rvm - результат интегрирования
    H - желаемая высота
    """

    r = np.sqrt(rvm[0]**2 + (rvm[1] + Constants.R_Earth)**2 + rvm[2]**2)
    h = r - Constants.R_Earth
    v = np.sqrt(np.sum(rvm[3:6] ** 2 ))
    rv = np.sqrt(rvm[3] * rvm[0] + rvm[4] * (rvm[1] + Constants.R_Earth) + rvm[5]*rvm[2])

    ancl = np.arccos(rv/(r*v))

    V_circ = np.sqrt(Constants.fM_Earth/r)

    error = np.zeros(4)
    error[0] = h - H
    error[1] = ancl - np.pi/2   # 
    error[2] = v - V_circ       # контроль скорости на круговой орбите
    error[3] = -rvm[6]          # масса

    if isnan(error[1]):
        error[1] = 100

    print("%8.4e km %8.4e deg %8.4e km/s Mk = %8.2f kg" % (error[0], ancl*180/np.pi-90, error[2], -error[3]))

    return error


def funMin1(x, rocket: Rocket, H, u):

    UX = []
    for i in range(rocket.Nstage):
        UX.append(u[i])

    k = 0
    for i in range(rocket.Nstage):
        for j in range(len(u[i])):
            UX[i][j] = x[k] * UX[i][j]
            k += 1

    rvm, tt, TR, step = Trajectory(UX, rocket)

    error = boundary_condition(rvm, H)

    err = np.dot(error[0:3], error[0:3])

    if min(TR[6, :]) < 0:
        err = err + 1000

    print("%10.5e %10.4f" % (err, -error[3]))
    return error[3] + err


def funMin2(x, rocket: Rocket, H):

    u = []
    for i in range(rocket.Nstage):
        if i == 0:
            u.append(x[0:4])
        else:
            u.append(x[4+3*(i-1):4+3*i])

    rvm, tt, TR, step = Trajectory(u, rocket)

    error = boundary_condition(rvm, H)

    err = np.dot(error[0:3], error[0:3])*100

    if min(TR[6, :]) < 0:
        err = err + 1000

    # print ("%10.5e %10.4f" %(err, -error[3]))

    return error[3] + err


def funMin3(x, rocket: Rocket, H, u):
    #   t0, alfa, ti, dfi

    UX = []
    for i in range(rocket.Nstage):
        UX.append([0]*len(u[i]))

    for i in range(rocket.Nstage):
        if i == 0:
            UX[i][0] = x[0] * u[i][0]
            UX[i][1] = u[i][1]
            UX[i][2] = u[i][2]
            UX[i][3] = x[1] * u[i][3]
        elif i == 1:
            UX[i][0] = x[2] * u[i][0]
            UX[i][1] = np.pi / 2
            UX[i][2] = x[3] * u[i][2]
        else:
            UX[i][0] = x[4 + 2 * (i-2)] * u[i][0]
            UX[i][1] = UX[i-1][1] + UX[i-1][2] * UX[i-1][0]
            UX[i][2] = x[5 + 2 * (i-2)] * u[i][2]

    rvm, tt, TR, step = Trajectory(UX, rocket)
    error = boundary_condition(rvm, H)
    err = error[0] * error[0] + error[1] * error[1] + error[2] * error[2] * 100
    return err * 100 + error[3]


def Solution(u, rocket: Rocket, H):
    Nu = 0
    for i in range(rocket.Nstage):
        Nu += len(u[i])
    x = [1] * Nu
    # res = cma.fmin(funMin1, x, 0.01, options={'maxfevals': 1e2, 'popsize':15},
    #                args=(rocket.Nstage, rocket.P, rocket.c, rocket.M, rocket.S, rocket.lam, rocket.fi, rocket.A, H, u,));
    # x = res[0]

    k = 0
    for i in range(rocket.Nstage):
        for j in range(len(u[i])):
            x[k] = x[k] * u[i][j]
            k += 1

    x = minimize(funMin2, x, args=(rocket.Nstage, rocket.P, rocket.c, rocket.M, rocket.S, rocket.lam, rocket.fi, rocket.A, H),
                 method='SLSQP', options={'maxiter': 30})
    print("Solution:", x)

    for i in range(rocket.Nstage):
        if i == 0:
            u[0:4] = x.x[0:4]
        else:
            u[4+3*(i-1):4+3*i] = x.x[4+3*(i-1):4+3*i]

    rvm, tt, TR, step = Trajectory(u, rocket.Nstage, rocket.P, rocket.c, rocket.M, rocket.S, rocket.lam, rocket.fi, rocket.A)
    error = boundary_condition(rvm, H)
    return x, u, rvm, tt, TR, step, error


def Solution2(u, rocket: Rocket, H):

    rvm, tt, TR, step = Trajectory(u, rocket)
    error = boundary_condition(rvm, H)

    x = np.eye(1, rocket.Nstage * 2) # [1] * rocket.Nstage * 2

    x = minimize(funMin3, x, args=(rocket, H, u),
                 method='Powell', options={'maxiter': 5})
    print("Solution2:", x)

    for i in range(rocket.Nstage):
        if i == 0:
            u[i][0] = x.x[0] * u[i][0]
            u[i][3] = x.x[1] * u[i][3]
        elif i == 1:
            u[i][0] = x.x[2] * u[i][0]
            u[i][1] = np.pi / 2
            u[i][2] = x.x[3] * u[i][2]
        else:
            u[i][0] = x.x[4 + 2 * (i-2)] * u[i][0]
            u[i][1] = u[i-1][1] + u[i-1][2] * u[i-1][0]
            u[i][2] = x.x[5+2*(i-2)] * u[i][2]

    rvm, tt, TR, step = Trajectory(u, rocket)
    error = boundary_condition(rvm, H)

    return x, u, rvm, tt, TR, step, error


def Trajectory_orbit(u, rocket: Rocket, orbit: Orbit):

    # параметры орбиты
    tet = 0 / Constants.raddeg

    rvm, tt, TR, step = Trajectory(u, rocket)
    # error = boundary_condition(rvm, orbit.H)

    n = len(TR[1, :])

    R0 = np.zeros(3)
    R0[0] = np.cos(rocket.fi) * np.cos(rocket.lam) * Constants.R_Earth
    R0[1] = np.cos(rocket.fi) * np.sin(rocket.lam) * Constants.R_Earth
    R0[2] = np.sin(rocket.fi) * Constants.R_Earth

    MA = [[ np.cos(-rocket.A), 0 , np.sin(-rocket.A)],
          [ 0      , 1 , 0      ],
          [-np.sin(-rocket.A), 0 , np.cos(-rocket.A)]]

    M = [
         [-np.sin(rocket.fi)*np.cos(rocket.lam), -np.sin(rocket.fi)*np.sin(rocket.lam),  np.cos(rocket.fi) ],
         [ np.cos(rocket.fi)*np.cos(rocket.lam),  np.cos(rocket.fi)*np.sin(rocket.lam),  np.sin(rocket.fi) ],
         [        -np.sin(rocket.lam),          np.cos(rocket.lam),  0       ]
         ]

    MA = np.array(MA)
    M = np.array(M)
    M = np.linalg.inv(M)
    TRG = np.zeros((6, n))
    for i in range(n):
        TRG[0:3, i] = R0[0:3] + np.dot(M,np.dot(MA, TR[0:3, i]))
        TRG[3:6, i] = np.dot(M, np.dot(MA, TR[3:6, i]))

    TR[1, :] = TR[1, :] + Constants.R_Earth

    # определение высоты
    h = np.zeros(n)
    for i in range(n):
        h[i] = np.sqrt(TR[0, i] * TR[0, i] + TR[1, i] * TR[1, i]
                       + TR[2, i] * TR[2, i]) - Constants.R_Earth

    # определение дальности
    R0 = [0, Constants.R_Earth, 0]
    L = np.zeros(n)
    for i in range(n):
        Rk = TR[0:3, i]
        danc = np.arccos((R0[0] * Rk[0] + R0[1] * Rk[1] + R0[2] * Rk[2])
                    / (np.sqrt(R0[0] * R0[0] + R0[1] * R0[1] + R0[2] * R0[2])
                       * np.sqrt(Rk[0] * Rk[0] + Rk[1] * Rk[1] + Rk[2] * Rk[2])))
        L[i] = danc * Constants.R_Earth

    print('Параметры траектории')
    for i in range(rocket.Nstage):
        print('В момент отделения ', str(i+1), 'ступени:')
        print('Положение:', TRG[0:3, step[i+1]-1], 'км')
        print('Скорость:',  TRG[3:6, step[i+1]-1], 'км/с')
        print('Масса:',     TR[6, step[i+1]-1],   'кг')

    Kep = Keplerian_elements(TRG[0:6, -1], 0, Constants.fM_Earth, 1)

    print('Параметры достигнутой орбиты:')
    print('Высоты перицентра и апоцентра:', [Kep[3]-Constants.R_Earth, Kep[4]-Constants.R_Earth], 'км')
    print('Наклонение, долгота восх. узла, аргумент перицентра',
          [Kep[5]*Constants.raddeg, Kep[6]*Constants.raddeg, Kep[6]*Constants.raddeg], 'град.')

    Kep = [(Kep[3]-Constants.R_Earth)/Constants.UnitR, (Kep[4]-Constants.R_Earth)/Constants.UnitR, Kep[5], Kep[7], Kep[6], tet]

    # по истинной анамалии
    TET = np.linspace(0, 2*np.pi, 60)
    orb = orbit_full(Kep, TET, 1)

    # графики
    R = 12000

    plt.ion()

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(projection='3d')  # plt.gca()

    u, v = np.mgrid[0:2*np.pi:36j, 0:np.pi:18j]
    x = Constants.UnitR * np.cos(u) * np.sin(v)
    y = Constants.UnitR * np.sin(u) * np.sin(v)
    z = Constants.UnitR * np.cos(v)
    ax.plot_wireframe(x, y, z, color="gray")

    ax.plot(TRG[0, :], TRG[1, :], TRG[2, :], color='r', linewidth=4)
    ax.plot(orb[:, 0] * Constants.UnitR, orb[:, 1] * Constants.UnitR, orb[:, 2] * Constants.UnitR, '--', color='k')

    ax.legend()
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.set_zlim(-R, R)
    ax.set_xlabel('X, km')
    ax.set_ylabel('Y, km')
    ax.set_zlabel('Z, km')
    ax.azim = 225
    plt.draw()

    plt.figure(figsize=(8, 6))
    plt.plot(L[:], h[:], color='r', linewidth=4)
    plt.plot(L[:], [orbit.H] * len(TR[0, :]), '--', color='k')

    plt.ylim(0, orbit.H+50)
    plt.xlabel('Дальность полёта (км)')
    plt.ylabel('Высота полёта (км)')
    plt.draw()

    plt.figure(figsize=(8, 6))
    plt.plot(tt[:], TR[6, :], color='r', linewidth=4)
    plt.xlabel('Время полёта (сек)')
    plt.ylabel('Масса (кг)')
    plt.draw()

    plt.ioff()

    plt.show()


def boundary_condition_ell(rvm, hp, ha, fi, lam, A):

    R0 = np.zeros(3)
    R0[0] = np.cos(fi) * np.cos(lam) * Constants.R_Earth
    R0[1] = np.cos(fi) * np.sin(lam) * Constants.R_Earth
    R0[2] = np.sin(fi) * Constants.R_Earth

    MA = [[ np.cos(-A), 0 , np.sin(-A)],
          [ 0      , 1 , 0      ],
          [-np.sin(-A), 0 , np.cos(-A)]]

    M = [
         [-np.sin(fi)*np.cos(lam),-np.sin(fi)*np.sin(lam),  np.cos(fi) ],
         [ np.cos(fi)*np.cos(lam), np.cos(fi)*np.sin(lam),  np.sin(fi) ],
         [        -np.sin(lam),         np.cos(lam),  0       ]
         ]

    MA = np.array(MA)
    M = np.array(M)
    M = np.linalg.inv(M)
    RVG = np.zeros(6)
    RVG[0:3] = R0[0:3] + np.dot(M, rvm[0:3])
    RVG[3:6] = np.dot(M, rvm[3:6])

    Kep = Keplerian_elements(RVG, 0, Constants.fM_Earth, 1)

    error = np.zeros(3)

    error[0] = Kep[3] - hp - Constants.R_Earth
    error[1] = Kep[4] - ha - Constants.R_Earth
    error[2] = -rvm[6]

    print("%8.4e km %8.4e km Mk = %8.2f kg" % (error[0], error[1], -error[2]))

    return error


def funMin_ell(x, rocket: Rocket, hp, ha, u):
    #   t0, alfa, ti, dfi

    UX = []
    for i in range(rocket.Nstage):
        UX.append([0] * len(u[i]))

    for i in range(rocket.Nstage):
        if i == 0:
            UX[i][0] = x[0] * u[i][0]
            UX[i][1] = u[i][1]
            UX[i][2] = u[i][2]
            UX[i][3] = x[1] * u[i][3]
        elif i == 1:
            UX[i][0] = x[2] * u[i][0]
            UX[i][1] = np.pi / 2
            UX[i][2] = x[3] * u[i][2]
        else:
            UX[i][0] = x[4 + 2 * (i-2)] * u[i][0]
            UX[i][1] = UX[i-1][1] + UX[i-1][2] * UX[i-1][0]
            UX[i][2] = x[5+2*(i-2)] * u[i][2]
    # print(UX)
    rvm, tt, TR, step = Trajectory(UX, rocket)
    error = boundary_condition_ell(rvm, hp, ha, rocket.fi, rocket.lam, rocket.A)
    err = error[0] * error[0] + error[1] * error[1]
    return err * 1000 + error[2]


def Solution_ELL(u, rocket: Rocket, hp, ha):

    rvm, tt, TR, step = Trajectory(u, rocket)
    error = boundary_condition_ell(rvm, hp, ha, rocket.fi, rocket.lam, rocket.A)
    x = [1] * rocket.Nstage * 2
#    res = cma.fmin(funMin_ell, x, 0.01, options={'maxfevals': 1e3, 'popsize':10}, args=(rocket, hp, ha, u,));
#    x = res[0]
    x = minimize(funMin_ell, x, args=(rocket, hp, ha, u), method='Powell', options={'maxiter':5})
    print(x)

    for i in range(rocket.Nstage):
        if i == 0:
            u[i][0] = x.x[0] * u[i][0]
            u[i][3] = x.x[1] * u[i][3]
        elif i == 1:
            u[i][0] = x.x[2] * u[i][0]
            u[i][1] = np.pi / 2
            u[i][2] = x.x[3] * u[i][2]
        else:
            u[i][0] = x.x[4 + 2 * (i-2)] * u[i][0]
            u[i][1] = u[i-1][1] + u[i-1][2] * u[i-1][0]
            u[i][2] = x.x[5 + 2 * (i-2)] * u[i][2]

    rvm, tt, TR, step = Trajectory(u, rocket)
    error = boundary_condition_ell(rvm, hp, ha, rocket.fi, rocket.lam, rocket.A)

    return x, u, rvm, tt, TR, step, error
