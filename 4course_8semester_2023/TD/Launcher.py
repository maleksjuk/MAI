import numpy as np
from Equations import Solution2, Trajectory_orbit, Solution_ELL
from Base import *
from datetime import datetime


def write_input_data(rocket: Rocket, control: Control, od: OrbitData, time_str: str) -> None:
    
    rocket_str = 'Данные ракеты\n'
    rocket_str += f'Количество ступеней: {rocket.Nstage}\n'
    rocket_str += f'Тяга: {rocket.P}\n'
    rocket_str += f'Удельный импульс: {rocket.I}\n'
    rocket_str += f'Скорость истечения: {rocket.c}\n'
    rocket_str += f'Масса: {rocket.M}\n'
    rocket_str += f'Площадь миделя: {rocket.S}\n'
    rocket_str += f'Долгота: {rocket.lam}\n'
    rocket_str += f'Широта: {rocket.fi}\n'
    rocket_str += f'Азимут: {rocket.A}\n'

    control_str = 'Управление\n'
    for i in range(len(control.u)):
        control_str += f'{i + 1} ступень: {control.u[i]}\n'

    od_str = 'Параметры орбиты\n'
    od_str += f'Высота (круговая): {od.H}\n'
    od_str += f'Высота (апоцентр): {od.ha}\n'
    od_str += f'Высота (перицентр): {od.hp}\n'
    od_str += f'Наклонение: {od.inc}\n'
    od_str += f'W: {od.W}\n'
    od_str += f'OM: {od.OM}\n'
    
    file_input_path = './results/' + time_str + '_input.txt'
    file = open(file_input_path, 'w')
    file.write(rocket_str + '\n' + control_str + '\n' + od_str)
    file.close()


def write_solution_data(Sol, u, rvm, tt, TR, step, error):
    
    text = f'Sol: {Sol}\n\n'
    text += f'u: {u}\n\n'
    text += f'rvm: {rvm}\n\n'
    text += f'tt: {tt}\n\n'
    text += f'TR: {TR}\n\n'
    text += f'step: {step}\n\n'
    text += f'error: {error}\n'

    file_solution_path = './results/' + time_str + '_solution.txt'
    file = open(file_solution_path, 'w')
    file.write(text)
    file.close()


if __name__ == "__main__":

    time_str = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')

    rocket = Rocket()
    control = Control()
    od = OrbitData()

    indata = 0 
    while indata == 1:
        print('Управление:', control.u)
        Trajectory_orbit(control.u, rocket.Nstage, rocket.P, rocket.c, rocket.M,
                         rocket.S, rocket.lam, rocket.fi, rocket.A, 
                         od.H, od.inc, od.W, od.OM) 
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

    write_input_data(rocket, control, od, time_str)

    #для круговых орбит
    Sol, u, rvm, tt, TR, step, error = Solution2(control.u, rocket.Nstage, rocket.P,
                                                 rocket.c, rocket.M, rocket.S, 
                                                 rocket.lam, rocket.fi, rocket.A, 200)

    write_solution_data(Sol, u, rvm, tt, TR, step, error)

    Trajectory_orbit(u, rocket.Nstage, rocket.P, rocket.c, rocket.M, rocket.S,
                     rocket.lam, rocket.fi, rocket.A, od.H, od.inc, od.W, od.OM) 
    
    exit()

    #для эллиптических орбит
    hp = 200    # перицентр
    ha = 500    # апоцентр
    control2 = Control()
    control2.u = [[155.3, 0.6, 0.85, -2.75/CONST.raddeg], [427.6, 90/CONST.raddeg, -0.1759/CONST.raddeg ]]
    Trajectory_orbit(control2.u, rocket.Nstage, rocket.P, rocket.c, rocket.M, 
                     rocket.S, rocket.lam, rocket.fi, rocket.A, hp, od.inc, od.W, od.OM) 
    Sol2, u2, rvm2, tt2, TR2, step2, error2 = Solution_ELL(control2.u, rocket.Nstage,
                                                           rocket.P, rocket.c, rocket.M, 
                                                           rocket.S, rocket.lam, 
                                                           rocket.fi, rocket.A, hp, ha)
