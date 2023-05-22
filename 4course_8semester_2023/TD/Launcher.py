import numpy as np
from Equations import Solution2, Trajectory_orbit, Solution_ELL
from Base import *
from datetime import datetime
from os import path, mkdir


def write_input_data(rocket: Rocket, control: Control, od: OrbitData, path_prefix: str) -> None:
    
    rocket_str = 'Данные ракеты\n'
    rocket_str += f'Количество ступеней: {rocket.Nstage}\n'
    rocket_str += f'Тяга: {rocket.P}\n'
    rocket_str += f'Удельный импульс: {rocket.I}\n'
    rocket_str += f'Скорость истечения: {rocket.c}\n'
    rocket_str += f'Масса: {rocket.M}\n'
    rocket_str += f'Площадь миделя: {rocket.S}\n'
    rocket_str += f'Долгота: {rocket.lam * CONST.raddeg: .4f}\n'
    rocket_str += f'Широта: {rocket.fi * CONST.raddeg: .4f}\n'
    rocket_str += f'Азимут: {rocket.A * CONST.raddeg: .4f}\n'

    control_str = '[Управление]\n\n'
    control_str += 'Круговая орбита\n'
    for i in range(len(control.u)):
        control_str += f'{i + 1} ступень:\n'
        control_str += f'  время: {control.u[i][0]}\n'
        if len(control.u[i]) == 4:
            control_str += f'  число Маха для начала мешка: {control.u[i][1]:.4f}\n'
            control_str += f'  число Маха для конца мешка: {control.u[i][2]:.4f}\n'
            control_str += f'  глубина мешка: {control.u[i][3] * CONST.raddeg:.4f}\n'
        elif len(control.u[i]) == 3:
            control_str += f'  начальный угол тангажа: {control.u[i][1] * CONST.raddeg:.4f}\n'
            control_str += f'  угловая скорость по тангажу: {control.u[i][2] * CONST.raddeg:.4f}\n'
    
    control_str += '\nЭллиптическая орбита\n'
    for i in range(len(control.u)):
        control_str += f'{i + 1} ступень:\n'
        control_str += f'  время: {control.u_ell[i][0]}\n'
        if len(control.u_ell[i]) == 4:
            control_str += f'  число Маха для начала мешка: {control.u_ell[i][1]:.4f}\n'
            control_str += f'  число Маха для конца мешка: {control.u_ell[i][2]:.4f}\n'
            control_str += f'  глубина мешка: {control.u_ell[i][3] * CONST.raddeg:.4f}\n'
        elif len(control.u_ell[i]) == 3:
            control_str += f'  начальный угол тангажа: {control.u_ell[i][1] * CONST.raddeg:.4f}\n'
            control_str += f'  угловая скорость по тангажу: {control.u_ell[i][2] * CONST.raddeg:.4f}\n'
    control_str += '---------------\n'

    od_str = 'Параметры орбиты\n'
    od_str += f'Высота (круговая): {od.H}\n'
    od_str += f'Высота (апоцентр): {od.ha}\n'
    od_str += f'Высота (перицентр): {od.hp}\n'
    od_str += f'Наклонение: {od.inc}\n'
    od_str += f'W: {od.W}\n'
    od_str += f'OM: {od.OM}\n'
    

    file_input_path = path_prefix + '_input.txt'
    file = open(file_input_path, 'w')
    file.write(rocket_str + '\n' + control_str + '\n' + od_str)
    file.close()


def write_solution_data(Sol, u, rvm, tt, TR, step, error, path_prefix):
    
    text = f'Sol: {Sol}\n\n'

    control_str = 'Управление\n'
    for i in range(len(u)):
        control_str += f'{i + 1} ступень:\n'
        control_str += f'  время: {u[i][0]}\n'
        if len(u[i]) == 4:
            control_str += f'  число Маха для начала мешка: {u[i][1]:.4f}\n'
            control_str += f'  число Маха для конца мешка: {u[i][2]:.4f}\n'
            control_str += f'  глубина мешка: {u[i][3] * CONST.raddeg:.4f}\n'
        elif len(u[i]) == 3:
            control_str += f'  начальный угол тангажа: {u[i][1] * CONST.raddeg:.4f}\n'
            control_str += f'  угловая скорость по тангажу: {u[i][2] * CONST.raddeg:.4f}\n'
    text += f'u: {control_str}\n\n'

    text += f'rvm: {rvm}\n\n'
    text += f'tt: {tt}\n\n'
    text += f'TR: {TR}\n\n'
    text += f'step: {step}\n\n'
    text += f'error: {error}\n'

    file_solution_path = path_prefix + '_solution.txt'
    file = open(file_solution_path, 'w')
    file.write(text)
    file.close()


def prepare_paths(H, optimize=False) -> tuple[str, str]:
    today_str = datetime.now().strftime('%Y-%m-%d')
    time_str = datetime.now().strftime('%H:%M:%S')

    dir_today = f'{today_str}'
    dir_today_txt = './results/' + dir_today
    dir_today_plt = './plots/' + dir_today
    
    if not path.exists(dir_today_txt):
        mkdir(dir_today_txt)
    if not path.exists(dir_today_plt):
        mkdir(dir_today_plt)

    dir_current_txt = dir_today_txt + f'/{time_str}_{H}km' + ('_opt' if optimize else '')
    dir_current_plt = dir_today_plt + f'/{time_str}_{H}km' + ('_opt' if optimize else '')
    
    if not path.exists(dir_current_txt):
        mkdir(dir_current_txt)
    if not path.exists(dir_current_plt):
        mkdir(dir_current_plt)

    path_prefix_txt = dir_current_txt + f'/{today_str}_{time_str}_{H}km'
    path_prefix_plt = dir_current_plt + f'/{today_str}_{time_str}_{H}km'
    return path_prefix_txt, path_prefix_plt


if __name__ == "__main__":
    
    TYPE = 3 # 0, 1 or 2, 3
    rocket = Rocket()
    control = Control()
    od = OrbitData()

    path_prefix_txt, path_prefix_plt = prepare_paths(od.H if TYPE in [0, 1] else od.hp,
                                                     True if TYPE in [1, 3] else False)
    
    write_input_data(rocket, control, od, path_prefix_txt)

    #для круговых орбит
    if TYPE == 0:
        Trajectory_orbit(control.u, rocket.Nstage, rocket.P, rocket.c, rocket.M, rocket.S,
                         rocket.lam, rocket.fi, rocket.A, od.H, od.inc, od.W, od.OM, path_prefix_txt, path_prefix_plt) 
    elif TYPE == 1:
        Sol, u, rvm, tt, TR, step, error = Solution2(control.u, rocket.Nstage, rocket.P,
                                                    rocket.c, rocket.M, rocket.S, 
                                                    rocket.lam, rocket.fi, rocket.A, od.H)
        write_solution_data(Sol, u, rvm, tt, TR, step, error, path_prefix_txt)
        Trajectory_orbit(u, rocket.Nstage, rocket.P, rocket.c, rocket.M, rocket.S,
                        rocket.lam, rocket.fi, rocket.A, od.H, od.inc, od.W, od.OM, path_prefix_txt, path_prefix_plt) 
    #для эллиптических орбит
    elif TYPE == 2:
        Trajectory_orbit(control.u_ell, rocket.Nstage, rocket.P, rocket.c, rocket.M, 
                        rocket.S, rocket.lam, rocket.fi, rocket.A, od.ha, od.inc, od.W, od.OM, path_prefix_txt, path_prefix_plt)
    if TYPE == 3:
        Sol2, u2, rvm2, tt2, TR2, step2, error2 = Solution_ELL(control.u_ell, rocket.Nstage,
                                                            rocket.P, rocket.c, rocket.M, 
                                                            rocket.S, rocket.lam, 
                                                            rocket.fi, rocket.A, od.hp, od.ha)
        write_solution_data(Sol2, u2, rvm2, tt2, TR2, step2, error2, path_prefix_txt)
        Trajectory_orbit(u2, rocket.Nstage, rocket.P, rocket.c, rocket.M, 
                        rocket.S, rocket.lam, rocket.fi, rocket.A, od.ha, od.inc, od.W, od.OM, path_prefix_txt, path_prefix_plt) 



    # indata = 0
    # while indata == 1:
    #     print('Управление:', control.u)
    #     Trajectory_orbit(control.u, rocket.Nstage, rocket.P, rocket.c, rocket.M,
    #                      rocket.S, rocket.lam, rocket.fi, rocket.A, 
    #                      od.H, od.inc, od.W, od.OM) 
    #     indata = input("Хотите изменить начальное приближение? (1 - да, 0 - нет) ")
    #     if int(indata) == 1:
    #         for i in range(rocket.Nstage):
    #             if i == 0:
    #                 control.u[i][0] = input("Время работы %i ступени (сек) " %(i+1))
    #                 control.u[i][1] = input("Число Маха для начала мешка (-) ")
    #                 control.u[i][2] = input("Число Маха для конца мешка (-) ")
    #                 control.u[i][3] = input("Глубина мешка (град) ")
    #             else:
    #                 control.u[i][0] = input("Время работы %i ступени (сек) " %(i+1))
    #                 control.u[i][1] = input("Начальное значение угла тангажа (град) ")
    #                 control.u[i][2] = input("Угловая скорость по тангажу (град/сек) ")
