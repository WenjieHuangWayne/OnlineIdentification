%Normal critic value 

vs = 50;
vl = 500;
vt = 5;

Ms = ((2*1.645)^2)*vs/ds;
Ml = ((2*1.645)^2)*vl/dl;
Mt = ((2*1.645)^2)*vt/dt;

Mtrue = ceil(max([Ms,Ml,Mt]));