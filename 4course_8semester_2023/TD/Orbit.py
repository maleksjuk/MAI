
import math
import numpy

from math import asin, atan2
from numpy import cos, sin, tan, sqrt, abs, pi

def orbit (Kep, R):
    
# Kep[0] -  hp, высота перицентра          -
# Kep[1] -  ha, высота апоцентра           -
# Kep[2] -   i, наклонение                rad
# Kep[3] -   w, аргумент перигелия        rad
# Kep[4] -  om, долгота восходящего узла  rad
# Kep[5] - tet, истинная аномалия         rad
# R      - средний радиус планеты
# Dec[]  - x,y,x,Vx,Vy,Vz    
    
    Dec = numpy.zeros(6)
    
    a = (Kep[0] + Kep[1])/2. + R
    e = (Kep[1] - Kep[0]) / (Kep[0] + Kep[1] + 2.*R)

    r = a * (1. - e * e) / (1. + e * cos(Kep[5]))
    v = sqrt(2. / r - 1. / a)
    
    gam = math.atan2(e * sin(Kep[5]), (1. + e * cos(Kep[5])))
    u = Kep[5] + Kep[3]
    ugam = u - gam
    
    Dec[0] = r * ( cos(u) * cos(Kep[4]) - sin(u) * cos(Kep[2]) * sin(Kep[4]))
    Dec[1] = r * ( cos(u) * sin(Kep[4]) + sin(u) * cos(Kep[2]) * cos(Kep[4]))    
    Dec[2] = r * (                                      sin(u) * sin(Kep[2]))
    
    Dec[3] = v * (-sin(ugam) * cos(Kep[4]) - cos(ugam) * cos(Kep[2]) * sin(Kep[4]))
    Dec[4] = v * (-sin(ugam) * sin(Kep[4]) + cos(ugam) * cos(Kep[2]) * cos(Kep[4]))
    Dec[5] = v * ( cos(ugam) * sin(Kep[2]))
    
    return Dec

def orbit_full(Kep, TET, R):
    
    n = len(TET)
    Dec = numpy.zeros((n,6))

    for i in range(n):
        Kep[5] = TET[i]
        Dec[i,:] = orbit(Kep,R)
    
    return Dec

def Topographic_coordinates(R, om_e, dt, lam0, R0):
        
    r = sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2])    
    dlam = atan2(R[0],R[1]) - atan2(R0[0],R0[1])    
    lam = lam0 + dlam - om_e * dt   

    while lam >  pi: lam -= 2*pi
    while lam < -pi: lam += 2*pi
    
    fi = asin(R[2]/r)
    
    return lam, fi

def Topographic_coordinates_3D(R, om_e, dt):
    
    fi = om_e*dt
    
    cos_fi = cos(fi)
    sin_fi = sin(fi)
    r = sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2])    
    
    x = (cos_fi*R[0] - sin_fi*R[1])/r
    y = (sin_fi*R[0] + cos_fi*R[1])/r
    z = R[2]/r
    
    return x,y,z

def Keplerian_elements (Dec, t, mu=1, par=0):

#		Dec[] - x,y,x,Vx,Vy,Vz UnitR, UnitV
#		t     - момент времени UnitT
#	
#		Kep[ 0] -  a, большая полуось				UnitR
#		Kep[ 1] -  p, фокальный параметр			UnitR
#		Kep[ 2] -  e, эксцентриситет				-
#		Kep[ 3] - rp, радиус перицентра			    UnitR
#		Kep[ 4] - ra, радиус апоцентра				UnitR
#
#		Kep[ 5] - Incl, наклонение					rad
#		Kep[ 6] - Node, долгота восходящего узла	rad
#		Kep[ 7] - Peri, аргумент перигелия			rad
#
#		Kep[ 8] - Nu, истинная аномалия			    rad		
#		Kep[ 9] -  E, эксцентрическая аномалия	    rad		
#		Kep[10] -  M, средняя аномалия				rad		
#
#		Kep[11] - tp, время прохождения перицентра	UnitT		

    c = numpy.zeros(3);
    f = numpy.zeros(3);
    l = numpy.zeros(3);

    r = numpy.sqrt(Dec[0]*Dec[0]+Dec[1]*Dec[1]+Dec[2]*Dec[2]);
    V = numpy.sqrt(Dec[3]*Dec[3]+Dec[4]*Dec[4]+Dec[5]*Dec[5]);

    h = V*V - 2.*mu/r;
    a = -mu/h;

    c[0] = Dec[1]*Dec[5] - Dec[2]*Dec[4];
    c[1] = Dec[2]*Dec[3] - Dec[0]*Dec[5];
    c[2] = Dec[0]*Dec[4] - Dec[1]*Dec[3];

    Cxy = numpy.sqrt(c[0]*c[0]+c[1]*c[1]);
    C = numpy.sqrt(Cxy*Cxy+c[2]*c[2]);

    p = C*C/mu;

    f[0] = Dec[4]*c[2] - Dec[5]*c[1] - mu/r*Dec[0];
    f[1] = Dec[5]*c[0] - Dec[3]*c[2] - mu/r*Dec[1];
    f[2] = Dec[3]*c[1] - Dec[4]*c[0] - mu/r*Dec[2];

    F = numpy.sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);

    e = F/mu;

    rp = p/(1.+e);
    ra = p/(1.-e);

    Incl = math.acos(c[2]/C);

    if	(numpy.abs(c[0])<1.e-10 and numpy.abs(c[1])<1.e-10):
        Node = 0;
    elif (c[0]>0):
        Node = math.acos(- c[1]/Cxy);
    elif (c[0]<0):
        Node = 2*math.pi - math.acos(- c[1]/Cxy);
    elif (c[0]==0 and c[1]<0):
        Node = 0;
    elif (c[0]==0 and c[1]>0):
        Node = math.pi;

    l[0] = numpy.cos(Node);
    l[1] = numpy.sin(Node);
    l[2] = 0.;

    cosPeri = (l[0]*f[0] + l[1]*f[1] + l[2]*f[2])/F;
    if abs(cosPeri)>1:
        cosPeri = cosPeri/abs(cosPeri);
        
    if (f[2]/F>1.e-8):
        Peri = math.acos(cosPeri);
    elif (f[2]/F<1.e-8):
        Peri = 2*math.pi-math.acos(cosPeri);
    elif (-1.e-8<f[2]/F<1.e-8 and l[0]*f[0]>0):
        Peri = 0.;
    elif (-1.e-8<f[2]/F<1.e-8 and l[0]*f[0]<0):
        Peri = math.pi;

    cosNu = (Dec[0]*f[0] + Dec[1]*f[1] + Dec[2]*f[2])/(r*F);
    if abs(cosNu)>1:
        cosNu = cosNu/abs(cosNu);
        
    frc = (f[1]*Dec[0]-f[2]*Dec[1])*c[0] + (f[2]*Dec[0]-f[0]*Dec[2])*c[1] + (f[0]*Dec[1]-f[1]*Dec[0])*c[2];

    if (frc>0):
        Nu = math.acos(cosNu);
    else:
        Nu = 2*math.pi - math.acos(cosNu);

    tanNu2 = numpy.tan(Nu/2.);

    if (h<0):
        E = 2.*math.atan2(numpy.sqrt(1.-e)*tanNu2,numpy.sqrt(1.+e));
        n = numpy.sqrt(mu/(a*a*a));
        M = E - e*numpy.sin(E);

    elif (h==0):

        E = tanNu2;
        n = 2.*numpy.sqrt(mu/p*p*p);
        M = E + 1./3.*E*E*E;

    elif (h>0):

        E = 2.*numpy.arctanh(complex (numpy.sqrt((e-1.)/(e+1.))*tanNu2, 0.)).real;
        n = numpy.sqrt(mu/(-a*a*a));
        M = e*numpy.sinh(E) - E;

    tp = t - M/n;

    if bool(par) is True:
       
        Kep = numpy.zeros(12);

        Kep[0]= a; Kep[1]= p; Kep[2]=e; 
        Kep[3]=rp; Kep[4]=ra;

        Kep[5]=Incl; Kep[6]=Node; Kep[7]=Peri;

        Kep[8]=Nu;   Kep[9]=E;    Kep[10]=M;

        Kep[11]=tp; 

    else:
        Kep = numpy.zeros(6);

        Kep[0] = a;
        Kep[1] = e;
        Kep[2] = Incl;
        Kep[3] = Peri;
        Kep[4] = Node;
        Kep[5] = M;
		
    return Kep