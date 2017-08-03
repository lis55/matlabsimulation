d=3;                  %%%% Dimensiones espaciales
n=3;                  %%%% Número de partículas

a=100*10^(-9);          %%%% Radio de las partículas [m]
masa=8.29*10^(-12);     %%%% Masa de las partícula [kg]
%masa=8.29*10^(-18);
I=(2/5)*masa*a^2;       %%%% Momento de inercia de las partículas


E=10^9;                 %%%% Módulo de Young
nu=0.5;                 %%%% Coeficiente de Poisson
mu=0.5;                 %%%% Coeficiente de roce
AA=9*10^-10;
AA=8*10^-7;
kn= 2/3*sqrt(a/2)*(E/(1-nu^2));     %%%%Constante elástica
bn= AA*2/3*sqrt(a/2)*(E/(1-nu^2));                      %%%%Constante viscosa 


u0=4*pi*10^(-7);      %%%% Permeabilidad en el vacío
ur=12.3;              %%%% Permeabilidad de las partículas
c2=4*pi*(ur-1)*(2*a)^3/(u0*(ur+2)*8);  %%%% Constantes c1 y c2
c1=u0/(4*pi);

bo=[0 1 0];           %%%% Campo magnético [T]	
dxbo=[0 0 0];         %%%% Derivadas espaciales de bo [T/m]
dybo=[0 0 0];
dzbo=[0 0 0];

ham=32;               %%%% Constante de Hamaker

vis=10^-3;            %%%% Viscosidad de la sangre
dv=3*10^-5;           %%%% Diámetro del canal
vf=10^-2;             %%%% Velocidad del fluido

%%%% tiempo
steps=1000;
t=1*10^(-4);
dt=t/steps;
time(1)=0;

%%%% Posiciones y velocidades iniciales 
%x=0.3*[1 3 2.5 3 3 8 7 4 1.6 6]*10^-6;
%y=0.3*[2 5 9 7 10 9 3 7 2 4]*10^-6;
%z=0.3*[0 0 2 4 1 2 8 2 4 7]*10^-6;
%x=0.3*[1 1 2 2 3 4 5 6 6 6 1 1 2 2 3 4 5 6 6 6]*10^-6;
%y=0.3*[0 3 1 6 0 3 1 1 2 1 5 2 6 2 7 0 4 1 3 8]*10^-6;
%z=0.3*[0 1 3 4 1 2 5 6 8 1 2 4 2 1 0 0 1 2 0 2]*10^-6;
%x=3*rand(1,n)'*10^-6;  
%y=3*rand(1,n)'*10^-6;
%z=3*rand(1,n)'*10^-6;
%x=[0 1 1.9]'*10^(-7);  
%y=[0 1.5 0]'*10^(-7);
%vx=[1 -1]*10^(-2);
%vy=[0 0];
x=[0 2 4]'*10^(-7);
y=[0 2*sqrt(3) 0]'*10^(-7);
z=[0 0 0];
%x=[4 3*cos(pi/4)]*10^-7;
%y=[0 3*sin(pi/4)]*10^-7;
%vx=rand(1,n)*10^-7;
%vy=rand(1,n)*10^-7;
%vz=rand(1,n)*10^-7;
vx=zeros(1,n);
vy=zeros(1,n);
vz=zeros(1,n);
%x=[0 0]*10^-7;
%y=[0 3]*10^-7;
%z=[0 0]*10^-7;
%x=[0 1 2]'*10^(-7);
%y=[0 3.45*sin(pi/6) 0]'*10^(-7);
%vy=[sin(pi/8) -2 sin(pi/8)];
%vx=[cos(pi/8) 0 -cos(pi/8)];
omega=zeros(1,n);
wy=zeros(1,n);
wx=zeros(1,n);
wz=zeros(1,n);


