import numpy as np

class Constants:
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
    u = [[143.3, 0.2, 0.65, -10.75 / Constants.raddeg],
         [430, 90 / Constants.raddeg, -0.119 / Constants.raddeg]]


class Stage:
    def __init__(self) -> None:
        pass

    def set_control(self, control: Control) -> None:
        pass


class Rocket:
    Nstage = 2
    g0 = Constants.g0
    P = [1115 * 10**3, 27.4 * 10**3]    # N ++      тяга по ступеням
    c = [282*g0, 348*g0]                # m/sec --  скорость истечения по ступеням
    M = [780*1000, 25*1000, 18*1000]    # kg +-+
    S = [7, 3.5 * np.pi**2]             # m2 --     площадь миделева сечения
    lam = -52.8 / Constants.raddeg      # rad +     долгота
    fi = 5.2 / Constants.raddeg         # rad +     широта
    A = 80 / Constants.raddeg           # rad -     азимут старта (вычислить в зависимости от наклонения)


class Orbit:
    H = 200     # -
    inc = 5.2   # + or 6
    W = 0       # -
    OM = 90     # -

    # Орбита
    # Большая полуось [км]
    # Эксцентриситет
    # Долгота восходящего узла [градусы]
    # Наклонение орбиты [градусы]
    # Аргумент перигея [градусы]
    # Период [час]
