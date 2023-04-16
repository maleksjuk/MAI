import numpy as np
from Equations import Solution2, Trajectory_orbit, Solution_ELL
from Base import *


if __name__ == "__main__":

    rocket = Rocket()
    orbit = Orbit()

    ctrl1 = Control()

    indata = 0
    while indata == 1:
        print(ctrl1.u)
        Trajectory_orbit(ctrl1.u, rocket, orbit)
        indata = input("Хотите изменить начальное приближение? (1 - да, 0 - нет) ")
        if indata == 1:
            for i in range(rocket.Nstage):
                if i == 0:
                    ctrl1.u[i][0] = input("Время работы %i ступени (сек) " % (i+1))
                    ctrl1.u[i][1] = input("Число Маха для начала мешка (-) ")
                    ctrl1.u[i][2] = input("Число Маха для конца мешка (-) ")
                    ctrl1.u[i][3] = input("Глубина мешка (град) ")
                else:
                    ctrl1.u[i][0] = input("Время работы %i ступени (сек) " % (i+1))
                    ctrl1.u[i][1] = input("Начальное значение угла тангажа (град) ")
                    ctrl1.u[i][2] = input("Угловая скорость по тангажу (град/сек) ")

    # для круговых орбит
    Sol, ctrl1.u, rvm, tt, TR, step, error = Solution2(ctrl1.u, rocket, 200)

    Trajectory_orbit(Control.u, rocket, orbit)

    # для эллиптических орбит

    hp = 200
    ha = 500
    orbit.H = hp

    ctrl2 = Control()
    ctrl2.u = [[155.3, 0.6, 0.85, -2.75/Constants.raddeg],
               [427.6, 90/Constants.raddeg, -0.1759/Constants.raddeg]]

    Trajectory_orbit(ctrl2.u, rocket, orbit)

    Sol2, u2, rvm2, tt2, TR2, step2, error2 = Solution_ELL(ctrl2.u, rocket, hp, ha)
