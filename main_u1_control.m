% This code computes the optimal dose in estrogen monotherapy that inhibits ovulation.

tic

[J, vmin]      = obj_fcn();

optimal cost  = sprintf('%0.7f', J)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N               = 140;
u1(1:N)         = vmin;
u1(N+1)         = u1(1);	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t0              = 0;
tf              = 28;
step            = tf/N;
time            = t0:step:tf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RPLH0           = 167.57;                   
LH0             = 11.81;                     
RPFSH0          = 14.48;                  
FSH0            = 11.41;                    
RcF0            = 2.10;                     
GrF0            = 4.12;                     
DomF0           = 0.46;                    
Sc10            = 1.06;                      
Sc20            = 1.67;                    
Lut10           = 4.16;                    
Lut20           = 13.03;                   
Lut30           = 16.48;                   
Lut40           = 10.29;   
Hist            = [RPLH0 LH0 RPFSH0 FSH0 RcF0 GrF0 DomF0...
                    Sc10 Sc20 Lut10 Lut20 Lut30 Lut40 0];                  % history for the state and control
options         = ddeset('RelTol', 1e-7, 'AbsTol', 1e-7);

solfiner        = dde23(@delayDE, [1.5], Hist, [t0 tf], options, vmin); 
timespline      = solfiner.x;                             
RPLHspline      = solfiner.y(1,:)';
LHspline        = solfiner.y(2,:)';
RPFSHspline     = solfiner.y(3,:)';
FSHspline       = solfiner.y(4,:)';
RcFspline       = solfiner.y(5,:)';
GrFspline       = solfiner.y(6,:)';
DomFspline      = solfiner.y(7,:)';
Sc1spline       = solfiner.y(8,:)';
Sc2spline       = solfiner.y(9,:)';
Lut1spline      = solfiner.y(10,:)';
Lut2spline      = solfiner.y(11,:)';
Lut3spline      = solfiner.y(12,:)';
Lut4spline      = solfiner.y(13,:)';
e0              = 57.60;              
e1              = 0.0269;             
e2              = 0.4196;             
e3              = 0.4923;             
p1              = 0.0032;             
p2              = 0.1188;             
h0              = 0.6606;             
h1              = 0.0193;             
h2              = 0.0159;             
h3              = 0.0119; 

u1spline        = interpn(time, u1, timespline, 'makima')';     % u1spline is the estrogen optimal control

P4spline        =  p1.*Lut3spline + p2.*Lut4spline;  
E2spline        =  u1spline + e0 + e1.*GrFspline + e2.*DomFspline + e3.*Lut4spline;
Inhspline       =  h0 + h1.*DomFspline + h2.*Lut2spline + h3.*Lut3spline;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc



function [J, vmin]         = obf_fcn()
N               = 140;                                                     % number of mesh intervals
IG              = zeros(N,1);                                              % initial guess for the components of the control vector u 
lowerb          = zeros(N,1);                                              % lower bound for the control

options         = optimoptions('fmincon', 'Maxiter', 9000, 'MaxFunEvals', 50000, 'TolCon', 1e-7, 'TolFun', 1e-7);

[vmin, J, exitflag]        = fmincon(@objective, IG, [], [], [], [], lowerb, [], [], options);

end


function [J]    = objective(v)
t0              = 0;
tf              = 28;
RPLH0           = 167.57;                   
LH0             = 11.81;                     
RPFSH0          = 14.48;                  
FSH0            = 11.41;                    
RcF0            = 2.10;                     
GrF0            = 4.12;                     
DomF0           = 0.46;                    
Sc10            = 1.06;                      
Sc20            = 1.67;                    
Lut10           = 4.16;                    
Lut20           = 13.03;                   
Lut30           = 16.48;                   
Lut40           = 10.29;   
Hist            = [RPLH0 LH0 RPFSH0 FSH0 RcF0 GrF0 DomF0...
                    Sc10 Sc20 Lut10 Lut20 Lut30 Lut40 0];                  % history for the state and control
options         = ddeset('RelTol', 1e-7, 'AbsTol', 1e-7);

solution        = dde23(@delayDE, [1.5], Hist, [t0 tf], options, v); 
time            = solution.x;
value           = solution.y;
J               = value(14, end);
end


function dstate = delayDE(t, state, delay, v)
N               = 140;
t0              = 0;
tf              = 28;
step            = tf/N;
time            = t0:step:tf;


kLH             =  0.9661;           
V0LH            =  550.03;          
V1LH            = 3329.19;          
KmLH            =  136.05;                                 
KiLHP           =  6.78;           
cLHE            =  0.0060;          
cLHP            =  1.98;           
VFSH            =  294.90;          
kFSH            = 14.59;            
cFSHE           = 0.0151;          
KiFSHInh        = 16.83;               
cFSHP           = 52.31;           
b               = 0.0453;             
c1              = 0.1036;             
c2              = 0.0577;             
c3              = 0.0170;             
c4              = 1.14;               
d1              = 0.7537;             
d2              = 0.6866;             
k1              = 0.6699;             
k2              = 0.6388;             
k3              = 0.9191;            
k4              = 1.88;               
alpha           = 0.9505;          
gamma           = 0.1615;          
e0              = 57.60;              
e1              = 0.0269;             
e2              = 0.4196;             
e3              = 0.4923;             
p1              = 0.0032;             
p2              = 0.1188;             
h0              = 0.6606;             
h1              = 0.0193;             
h2              = 0.0159;             
h3              = 0.0119;             
w               = 9.21;


RPLH            =  state(1); 
LH              =  state(2); 
RPFSH           =  state(3); 
FSH             =  state(4); 
RcF             =  state(5);   
GrF             =  state(6);
DomF            =  state(7);
Sc1             =  state(8);
Sc2             =  state(9);
Lut1            =  state(10);
Lut2            =  state(11);
Lut3            =  state(12);
Lut4            =  state(13);
OF              =  state(14);


u(1:N)          =   v;
u(N+1)          =   u(1);
E2              =   interpn(time, u, t, 'makima') + e0 + e1.*GrF + e2.*DomF + e3.*Lut4;
P4              =   p1.*Lut3 + p2.*Lut4;

statelag1       = delay(:,1);

DomFdelayInh    = statelag1(7);
Lut2delayInh    = statelag1(11);
Lut3delayInh    = statelag1(12);

InhdelayInh     = h0 + h1.*DomFdelayInh + h2.*Lut2delayInh + h3.*Lut3delayInh;



dRPLHdt         =   ((V0LH + (V1LH.*E2.^8./(KmLH.^8 + E2.^8)))./(1 + (P4./KiLHP))) - ((kLH.*(1 + cLHP.*P4).*RPLH)./(1 + cLHE.*E2));
dLHdt           =   (1./2.5).*((kLH.*(1 + cLHP.*P4).*RPLH)./(1 + cLHE.*E2)) - 14.*LH;
dRPFSHdt        =   (VFSH./(1 + (InhdelayInh./KiFSHInh) + (P4./w))) - ((kFSH.*(1 + cFSHP.*P4).*RPFSH)/(1 + cFSHE.*E2.^2));
dFSHdt          =   (1./2.5).*((kFSH.*(1 + cFSHP.*P4).*RPFSH)./(1 + cFSHE.*E2.^2)) -  8.21.*FSH;  
dRcFdt          =   (b + c1.*RcF).*(FSH./(1 + (P4./5.11))) - c2.*(LH.^alpha).*RcF;
dGrFdt          =    c2.*(LH.^alpha).*RcF - c3.*LH.*GrF; 
dDomFdt         =    c3.*LH.*GrF - c4.*(LH.^gamma).*DomF;
dSc1dt          =    c4.*(LH.^gamma).*DomF - d1.*Sc1;
dSc2dt          =    d1.*Sc1 - d2.*Sc2;
dLut1dt         =    d2.*Sc2 - k1.*Lut1;
dLut2dt         =    k1.*Lut1 - k2.*Lut2;
dLut3dt         =    k2.*Lut2 - k3.*Lut3;
dLut4dt         =    k3.*Lut3 - k4.*Lut4;
dOFdt           =    LH + (5e-06).*(interpn(time, u, t, 'makima'));




dstate          = zeros(14,1);

dstate(1)       = dRPLHdt;
dstate(2)       = dLHdt;
dstate(3)       = dRPFSHdt;
dstate(4)       = dFSHdt;
dstate(5)       = dRcFdt;
dstate(6)       = dGrFdt;
dstate(7)       = dDomFdt;
dstate(8)       = dSc1dt;
dstate(9)       = dSc2dt;
dstate(10)      = dLut1dt;
dstate(11)      = dLut2dt;
dstate(12)      = dLut3dt;
dstate(13)      = dLut4dt;
dstate(14)      = dOFdt;

end