Vk = 15.873016;
Fk = 3000;
Vk_1 = sqrt((2*ds/m)*Fk + (1 - 2*ds*CdA/m)*(Vk*Vk) - 2*ds*g*(sin(0)+crr*cos(0)))

dt = 2 * ds / (21.541950 + 21.541950);
Tk = 25.102041;
Qk = -1020.408163;
Tk_1 = Tk + (dt/Cth)*(Qk + Qsun + Qpas + (Tamb - Tk)/Rth)