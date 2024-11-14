#Imports
import math 
import numpy as np
from scipy.integrate import quad

#Konstanten + Astronomische Werte(Buch Sonnensystem Quelle [43])
G=6.67430e-11 #Gravitationskonstante G
c=299792458 #Lichtgeschwindigkeit c
r_Son=695500 #Sonnenradius
M_Son= 1.989*10**30 #Sonnenmasse
r_E= 6.378*10**6 #Erdradius
M_E=5.976 *10**24 #Erdmasse

#=======================================================Theoretische Grundlagen=========================================================#
#Berechnung von Schwarzschildradien für verschiedene Objekte. Abschnitt 2.4.3 
print("Theoretische Grundlagen Berechnungen Tabelle 2.1")
#Erde
E_s=(2*G*M_E)/c**2
print("Schwarzschildradius der Erde=",E_s,"m","=",E_s*1000,"mm")
#Sonne
Son_s=(2*G*M_Son)/c**2
print("Schwarzschildradius der Sonne=",Son_s,"m","=",Son_s/1000,"km")
print()
#=====================================================Resultate=======================================================================#
#Gegebene Werte (siehe Arbeit) 
M= 1.989e+38 #Masse M von Gargantua in kg
chi=0.99999999999999 #Dimensionsloser Spin-Parameter chi
a=(chi*G*M)/c**2 #Dimensionierter Spin-Parameter a

#===========Abschnitt 4.1.2 Ereignishorizont===============#
print("Resultate Abschnitt 4.1.2 Ereignishorizont:")
#Berechnung Ereignishorizont mit Schwarzschild-Metrik
R_s= 2*G*M/c**2
print("Ereignishorizont Schwarzschild-Metrik =",R_s,"m","=",R_s/1000,"km")
#Berechnung Ereignishorizont mit Kerr-Metrik
R_k= (G*M/c**2)*(1+ math.sqrt(1-chi**2))
print("Äusserer Ereignishorizont Kerr-Metrik =",R_k,"m","=", R_k/1000, "km")

#Vergleich von R_s und R_k
Vergleich= R_k/R_s
print("Der Kerr-Ereignishorizont ist",Vergleich, "mal grösser als der Schwarzschild-Ereignishorizont")

#Vergleich der Radien mit Sonnenradius
S_Vergleich_s= (R_s/1000)/r_Son 
print("R_s ist", S_Vergleich_s,"mal grösser als der Sonnenradius")

S_Vergleich_k=(R_k/1000)/r_Son 
print("R_k ist", S_Vergleich_k,"mal grösser als der Sonnenradius")

print()

#============Abschnitt 4.1.3 ISCO-Orbit=======================#
print("Resultate Abschnitt 4.1.3 ISCO-Orbit:")
#ISCO-Orbit mit Schwarzschild-Metrik
R_s_ISCO= (6*G*M)/c**2
print("Schwarzschild: Radiale Koordinate des ISCO=", R_s_ISCO,"m","=" ,R_s_ISCO/1000,"km") 

#ISCO-Orbit mit Kerr Metrik
Z_1=1 + (1 - chi**2)**(1/3)*((1+chi)**(1/3)+(1-chi)**(1/3))
Z_2=math.sqrt(3*chi**2+Z_1**2)
R_ISCO_retrograd=(G*M/c**2)*(3+Z_2+ math.sqrt((3-Z_1)*(3+Z_1+(2*Z_2))))
print("Kerr: Retrograde radiale Koordinate ISCO=", R_ISCO_retrograd,"m","=", R_ISCO_retrograd/1000,"km")
R_ISCO_prograd=(G*M/c**2)*(3+Z_2 - math.sqrt((3-Z_1)*(3+Z_1+(2*Z_2))))
print("Kerr: Prograde radiale Koordinate ISCO=", R_ISCO_prograd,"m","=",R_ISCO_prograd/1000,"km")

print()

#=================Abschnitt 4.2.1 Eigenschaften=================#
print("Resultate Abschnitt 4.2.1 Eigenschaften")
r_M=1.3*r_E
print("Radius Miller-Planet=", r_M)
M_M= (1.3**3)*M_E
print("Masse Miller-Planet=",M_M,"kg")

F_GM=(G*(M_M)/(r_M**2))
print("Bindende Kraft=",F_GM,"m' m/s^2")
print()

#==================Abschnitt 4.2.2 Zeitdilatation===================#
print("Resultate Abschnitt 4.2.2 Zeitdilatation")
#Berechnung Gamma
Delta_t_M= 3600 #Vergangene Zeit auf Miller Planet in s
Delta_t_E= 7*365*24*60*60 #Vergangene Zeit auf Erde in s
Gamma= Delta_t_E/Delta_t_M
print("Relativistischer Faktor=", Gamma)

#Schwarzschild-Metrik radiale Koordinate
r_s= (2*G*M)/(c**2*(1-((Delta_t_M/Delta_t_E)**2)))
print("Schwarzschild: Radiale Koordinate des Planeten=",r_s,"m","=",r_s/1000,"km")

#Kerr-Metrik radiale Koordinate
theta=(3*math.pi/8)
A=((c**2)*((1-(Delta_t_M/Delta_t_E)**2)))/(2*G*M)
b=-1
B=(a**2)*(math.cos(theta)**2)
z=(b**2-4*(A**2)*B) #Diskrimante


r_prograd=(-b-math.sqrt((z)))/(2*A)
print("Prograde radiale Koordinate für Planet=", r_prograd,"m","=",r_prograd/1000,"km")
r_retrograd=(-b+math.sqrt((z)))/(2*A)
print("Retrograder radiale Koordinate für Planet=", r_retrograd,"m","=", r_retrograd/1000,"km")

print()
#==================Abschnitt 4.2.3 Gezeitenkräfte===================#
#Schwarzschild-Metrik Beschleunigung
print("Resultate Abschnit 4.2.3 Gezeitenkräfte")
a_s=((2*G*M)/r_s**3)*2*r_M
print("Schwarzschild Gezeitenbeschleunigung =",a_s, "m/s^2")

if a_s>F_GM:
    print("Schwarzschild: Planet überlebt nicht")
if a_s<F_GM:
    print("Schwarzschild: Planet überlebt")

#Kerr-Metrik Beschleunigung
a_k=(2*G*M*r_retrograd*(2*r_M)*(r_retrograd**2-3*a**2))/((r_retrograd**2 + a**2)**3)
print("Kerr Gezeitenbeschleunigung=", a_k, "m/s^2")

if a_k>F_GM:
    print("Kerr: Planet überlebt nicht")
if a_k<F_GM:
    print("Kerr: Planet überlebt")