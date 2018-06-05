import copy
import math
M0 = 4670
G = 6.67 * 10**11
F = 15600
u = 3050
rm = 1738000 #м,Радиус Луны
gm = 1.62
Rm = 385000000 #м,Радиус орбиты луны
w0 = 0.00031 #рад/с
Mm = 7.35 * (10**22) #кг,масса луны
mu = F/u
dt = 0.01
def Waitingtime(x,y,X,Y,Z,V,alpha):
    global M0,mu,dt
    H = math.hypot(X, math.hypot(Y, Z))
    phi = math.acos(x / rm)
    Phi = math.acos(x / H)
    W = V / H
    v = 0
    M = M0
    t = 0
    w = w0
    h = rm
    g = gm

    while h < H:
        h += v * dt + F/M * (dt**2)/2 * math.sin(math.radians(alpha) - g * (dt**2)/2)
        phi += w * dt + F/M * (dt**2) * math.cos(math.radians(alpha)) / h
        Phi += W * dt
        M -= mu * dt
        alpha -= 1 * dt
        if alpha <= 0:
            alpha = 0
        g = G*Mm/(h)**2
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
        return (phi - Phi)/(W - w0)
    else :
        return (2*math.pi - (phi - Phi)) / (W - w0)
def output(x,y,X,Y,Z,V,t0,alpha):
    global M0, mu,dt
    H = math.hypot(X, math.hypot(Y, Z))
    phi = math.acos(x / rm)
    Phi = math.acos(x / H)
    W = V / H
    v = 0

    M = M0
    t = Waitingtime(x,y,X,Y,Z,V,alpha)
    w = w0
    h = rm
    g = gm

    while h < H:
        h += v * dt + F / M * (dt ** 2) / 2 * math.sin(math.radians(alpha) - g * (dt ** 2) / 2)
        phi += w * dt + F / M * (dt ** 2) * math.cos(math.radians(alpha)) / h
        Phi += W * dt
        M -= mu * dt
        alpha -= 1*dt
        if alpha <= 0:
            alpha = 0
        g = G * Mm / (h) ** 2
        v += F / M * dt * math.sin(math.radians(45))
        w += F / M * (dt ** 2) * math.cos(math.radians(45)) / h
        phi = phi % (2 * math.pi)
        Phi = Phi % (2 * math.pi)
        if M <= 2335:
            mu = 0
            F = 0
        t += dt
        print(h*math.cos(phi),h*math.sin(phi),alpha)
    if alpha == 0:
        return [t,h*math.cos(phi),h*math.cos(phi),math.hypot(v,w * h)*math.cos(phi),math.hypot(v,w * h)*math.sin(phi),Phi-phi]
    else: return output(x,y,X,Y,Z,V,t0,alpha-1)

print (output(int(input()),int(input()),int(input()),int(input()),int(input()),int(input()),int(input()),45))