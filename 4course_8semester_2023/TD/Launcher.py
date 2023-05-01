import numpy as np
from Equations import Solution2, Trajectory_orbit, Solution_ELL
from Base import *

if __name__ == "__main__":

    rocket = Rocket()
    control = Control()
    orb = Orbit()

    indata = 0 
    while indata == 1:
        print('Управление:', control.u)
        Trajectory_orbit(control.u, rocket.Nstage, rocket.P, rocket.c, rocket.M,
                         rocket.S, rocket.lam, rocket.fi, rocket.A, 
                         orb.H, orb.inc, orb.W, orb.OM) 
        indata = input("Хотите изменить начальное приближение? (1 - да, 0 - нет) ")
        if indata == 1:
            for i in range(rocket.Nstage):
                if i == 0:
                    control.u[i][0] = input("Время работы %i ступени (сек) " %(i+1))
                    control.u[i][1] = input("Число Маха для начала мешка (-) ")
                    control.u[i][2] = input("Число Маха для конца мешка (-) ")
                    control.u[i][3] = input("Глубина мешка (град) ")
                else:
                    control.u[i][0] = input("Время работы %i ступени (сек) " %(i+1))
                    control.u[i][1] = input("Начальное значение угла тангажа (град) ")
                    control.u[i][2] = input("Угловая скорость по тангажу (град/сек) ")
                    
    #для круговых орбит
    Sol, u, rvm, tt, TR, step, error = Solution2(control.u, rocket.Nstage, rocket.P,
                                                 rocket.c, rocket.M, rocket.S, 
                                                 rocket.lam, rocket.fi, rocket.A, 200)
    Trajectory_orbit(u, rocket.Nstage, rocket.P, rocket.c, rocket.M, rocket.S,
                     rocket.lam, rocket.fi, rocket.A, orb.H, orb.inc, orb.W, orb.OM) 

    #для эллиптических орбит
    hp = 200    # перицентр
    ha = 500    # апоцентр
    control2 = Control()
    control2.u = [[155.3, 0.6, 0.85, -2.75/CONST.raddeg], [427.6, 90/CONST.raddeg, -0.1759/CONST.raddeg ]]
    Trajectory_orbit(control2.u, rocket.Nstage, rocket.P, rocket.c, rocket.M, 
                     rocket.S, rocket.lam, rocket.fi, rocket.A, hp, orb.inc, orb.W, orb.OM) 
    Sol2, u2, rvm2, tt2, TR2, step2, error2 = Solution_ELL(control2.u, rocket.Nstage,
                                                           rocket.P, rocket.c, rocket.M, 
                                                           rocket.S, rocket.lam, 
                                                           rocket.fi, rocket.A, hp, ha)
