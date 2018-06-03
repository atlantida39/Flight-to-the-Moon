import copy
import math
M = 4670
F = 15600
u = 3050
rm = 1738000 #м,Радиус Луны
gm = 1.62
Rm = 385000000 #м,Радиус орбиты луны
w0 = 0.00031 #рад/с
Mm = 7.35 * (10**22) #кг,масса луны
x = int(input())#это будут координаты лунного модуля в момент посадки в СО, связанной с луной
y = int(input())
X = int(input())#это будут координаты лунного корабля на орбите в момент посадки ЛМ в СО, связанной с луной
Y = int(input())
V = int(input())
t0 = int(input())
H = math.hypot(X,math.hypot(Y,Z))
mu = F/u
phi0 = math.acos(x/rm)
Phi0 = math.acos(x/H)
phi = copy.copy(phi0)
Phi = copy.copy(Phi0)
W = V/H
w = copy.copy(w0)
v = 0
t = 0
h = rm
g = copy.copy(gm)
dt = 0.01
while h < H:
    h += v * dt + F/M * (dt**2)/2 * math.sin(math.radians(45) - g * (dt**2)/2)
    phi += w * dt + F/M * (dt**2) * math.cos(math.radians(45)) / h
    Phi += W * dt
    M -= mu * dt
    g /= (rm/h)**2
    v += F/M * dt * math.sin(math.radians(45))
    w += F/M * (dt**2) * math.cos(math.radians(45)) / h
    phi = phi % (2*math.pi)
    Phi = Phi % (2 * math.pi)
    if M <= 2335:
        mu = 0
        F = 0
    #t += dt
    #print(h*math.cos(phi),h*math.cos(phi))
if (phi - Phi) > 0:
    Waitingtime = (phi - Phi)/(W - w0)
else :
    Waitingtime = (2*math.pi - (phi - Phi)) / (W - w0)
phi = copy.copy(phi0)
Phi = copy.copy(Phi0)
W = V/H
w = copy.copy(w0)
v = 0
t = Waitingtime
h = rm
dt = 0.01
while h < H:
    h += v * dt + F / M * (dt ** 2) / 2 * math.sin(math.radians(45) - g * (dt ** 2) / 2)
    phi += w * dt + F / M * (dt ** 2) * math.cos(math.radians(45)) / h
    Phi += W * dt
    M -= mu * dt
    g /= (rm / h) ** 2
    v += F / M * dt * math.sin(math.radians(45))
    w += F / M * (dt ** 2) * math.cos(math.radians(45)) / h
    phi = phi % (2 * math.pi)
    Phi = Phi % (2 * math.pi)
    if M <= 2335:
        mu = 0
        F = 0
    t += dt
    print(h*math.cos(phi),h*math.cos(phi))
print (t,h*math.cos(phi),h*math.cos(phi),math.hypot(v,w * h)*math.cos(phi),math.hypot(v,w * h)*math.sin(phi),Phi-phi)

