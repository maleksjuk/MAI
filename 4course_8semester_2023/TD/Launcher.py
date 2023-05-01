import numpy as np
from Equations import Solution2, Trajectory_orbit, Solution_ELL

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

# параметры ракеты
# брать из википедии
Nstage = 2 
g0=9.80665 # !! константа
P = [4541338,981*1000]              # N     тяга по ступеням
c = [282*g0, 348*g0]                # m/sec скорость истечения по ступеням
M = [415*1000, 25*1000, 3*1000]     # kg 
S = [7, 3.5**2 * np.pi]               # m2    площадь миделева сечения
lam = 80/raddeg                     # rad   долгота
fi  = 51.6/raddeg                   # rad   широта
A   = 80/raddeg                     # rad   азимут старта (вычислить в зависимости от наклонения)

# параметры орбиты
H = 200     # высота
inc = 51.6  # наклонение
W = 0 
OM = 90

# управление
u = [[143.3, 0.2, 0.65, -10.75/raddeg], [430, 90/raddeg, -0.119/raddeg ]]
# 1 1 для времени начальное приближение (честно, не выбираемый параметр. т.к. есть макс заправка): масса топлива умножить на массовый расход
# 1 2.. остальное: по наитию (подобрать не совсем ужасное) - и будет всё хорошо
# 2 1
# 2 2
# 2 3 этот коэф - положить вектор скорости в горизонатль, можно посчитать: за 430 с развернуться до 0: дели, получаем угловую скорость

indata = 0 
while indata == 1:
    print('Управление:', u)
    Trajectory_orbit(u, Nstage, P, c, M, S, lam, fi, A, H, inc, W, OM) 
    indata = input("Хотите изменить начальное приближение? (1 - да, 0 - нет) ")
    if indata == 1:
        for i in range(Nstage):
            if i == 0:
                u[i][0] = input("Время работы %i ступени (сек) " %(i+1))
                u[i][1] = input("Число Маха для начала мешка (-) ")
                u[i][2] = input("Число Маха для конца мешка (-) ")
                u[i][3] = input("Глубина мешка (град) ")
            else:
                u[i][0] = input("Время работы %i ступени (сек) " %(i+1))
                u[i][1] = input("Начальное значение угла тангажа (град) ")
                u[i][2] = input("Угловая скорость по тангажу (град/сек) ")
                
#для круговых орбит
Sol, u, rvm,tt,TR,step, error = Solution2(u, Nstage, P, c, M, S, lam, fi, A, 200)
Trajectory_orbit(u, Nstage, P, c, M, S, lam, fi, A, H, inc, W, OM) 

#для эллиптических орбит
hp = 200    # перицентр
ha = 500    # апоцентр
u2 = [[155.3, 0.6, 0.85, -2.75/raddeg], [427.6, 90/raddeg, -0.1759/raddeg ]]
Trajectory_orbit(u2, Nstage, P, c, M, S, lam, fi, A, hp, inc, W, OM) 
Sol2, u2, rvm2, tt2, TR2, step2, error2 = Solution_ELL(u2, Nstage, P, c, M, S, lam, fi, A, hp, ha)
