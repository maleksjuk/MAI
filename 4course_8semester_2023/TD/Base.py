import numpy as np

class CONST:
    raddeg = 180 / np.pi

    fM_Earth = 398600.4415      # km^3/sec^2
    fM_Moon = 4902.801          # km^3/sec^2
    fM_Sun = 1.32712440018e11   # km^3/sec^2

    R_Earth = 6371				# km
    J2_Earth = 1.0826157e-3
    c_Earth = J2_Earth * fM_Earth * R_Earth * R_Earth / 2.0	 # km^5/sec^2
    # om_Earth = 4.17807462 * 10**(-3) * 86400 / raddeg     # deg/sec*86400/raddeg=rad/day
    om_Earth = 4.17807462 * 10**(-3) / raddeg    # deg/sec /raddeg = rad/sec

    UnitfM = fM_Earth
    UnitR = R_Earth                # km
    UnitV = np.sqrt(UnitfM / UnitR)  # km/sec
    UnitT = (UnitR / UnitV) / 86400    # day

    g0 = 9.80665


class Control:
    # VARIANT 1

    time1_boost = 110
    time1_full = 770
    time2_full = 350
    times = [time1_boost, time1_full - time1_boost, time2_full]
    theta = [90, 90, 5, 0]

    u = [[times[0], #130
          0.4,
          0.65,
          -5. / CONST.raddeg
          ],
         [times[1], # 540 - 130,
          theta[1] / CONST.raddeg,
          (theta[2] - theta[1]) / times[1] / CONST.raddeg
          ],
         [times[2], # 945,
          theta[2] / CONST.raddeg,
          (theta[3] - theta[2]) / times[2] / CONST.raddeg
          ]]
    

    time1_ell_boost = 117
    time1_ell_full = 650
    time2_ell_full = 500
    times_ell = [time1_ell_boost, time1_ell_full - time1_ell_boost, time2_ell_full]
    theta_ell = [90, 90, 20, 0]

    u_ell = [[times_ell[0], #130
              0.4,
              0.65,
              -0. / CONST.raddeg
              ],
             [times_ell[1], # 540 - 130,
              theta_ell[1] / CONST.raddeg,
              (theta_ell[2] - theta_ell[1]) / times_ell[1] / CONST.raddeg
              ],
             [times_ell[2], # 945,
              theta_ell[2] / CONST.raddeg,
              (theta_ell[3] - theta_ell[2]) / times_ell[2] / CONST.raddeg
              ]]
    
    # ORIGIN
    # u = [[143.3, 0.2, 0.65, -10.75/CONST.raddeg], [430, 90/CONST.raddeg, -0.119/CONST.raddeg ]]
    # u = [[155.3, 0.6, 0.85, -2.75/CONST.raddeg], [427.6, 90/CONST.raddeg, -0.1759/CONST.raddeg ]]


class Stage:
    def __init__(self) -> None:
        pass

    def set_control(self, control: Control) -> None:
        pass


def impulse_special(P: list, dm_dt: list) -> float:
    imp = np.sum(P) / np.sum(dm_dt)
    return imp


class Rocket:
    # VARIANT 1
    Nstage = 3

    thrust = {'booster': 5100 * 10**3,
              1: 1400 * 10**3,
              2: 30 * 10**3}
    dm_dt = {'booster': 1551,
             1: 280,
             2: 9}

    P = [(thrust[1] + thrust['booster'] * 2),    # (1371 + 7000 * 2) * 10**3,
         thrust[1],    # 1371 * 10**3,
         thrust[2]    # 27.4 * 10**3
        ]    # N +++      тяга по ступеням
    
    # массовые расходы (кг/с)
    # 1 stage: 
    #   EAP (boosters) * 2 = 3102
    #   EPC = 280
    # 2 stage: EPC = 280
    # 3 stage: ESC-A = 9

    

    I = [284,   # EPC + 2 EAP
         432,                # EPC
         324                 # ESC-A
         ]
    c = [i * CONST.g0 for i in I]            # m/sec +++  скорость истечения по ступеням
    # c = [1 * CONST.g0, 431 * CONST.g0, 324 * CONST.g0]            # m/sec ---  скорость истечения по ступеням

    M = [780 * 10**3,               # start mass
         (33 * 2 + 2.675) * 10**3,  # dry EAP (boosters) + fairing
         14.7 * 10**3,              # dry EPC
         4.540 * 10**3              # dry ESC-A
        ] # kg ++++
    
    S = [((5.4/2)**2 + 2 * (3.05/2)**2 ) * np.pi,
         (5.4/2)**2 * np.pi,
         (5.4/2)**2 * np.pi] # m2 +++     площадь миделева сечения
    
    lam = -52.8 / CONST.raddeg      # rad +     долгота
    fi = 5.2 / CONST.raddeg         # rad +     широта
    A = 90 / CONST.raddeg           # rad +?     азимут старта (вычислить в зависимости от наклонения)

    # ORIGIN
    # Nstage = 2 
    # P = [4541338,981*1000]              # N     тяга по ступеням
    # I = []
    # c = [282 * CONST.g0, 348 * CONST.g0]                # m/sec скорость истечения по ступеням
    # M = [415*1000, 25*1000, 3*1000]     # kg 
    # S = [7, 3.5**2 * np.pi]               # m2    площадь миделева сечения
    # lam = 80/CONST.raddeg                     # rad   долгота
    # fi  = 51.6/CONST.raddeg                   # rad   широта
    # A   = 80/CONST.raddeg                     # rad   азимут старта (вычислить в зависимости от наклонения)



class OrbitData:
    # VARIANT 1
    H = 220     # -
    inc = 5.2   # + or 6
    W = 0       # -
    OM = 90     # -    
    
    # ORIGIN
    # H = 200     # высота
    # inc = 51.6  # наклонение
    # W = 0 
    # OM = 90

    ha = 2429
    hp = 829

    # Орбита
    # Большая полуось [км]
    # Эксцентриситет
    # Долгота восходящего узла [градусы]
    # Наклонение орбиты [градусы]
    # Аргумент перигея [градусы]
    # Период [час]
