clc
clear all
close all
format long

parametros
condiciones_iniciales

%%%% tiempo
steps=1000;
t=1*10^(-4);
dt=t/steps;
time(1)=0;

tt=0:pi/30:2*pi;
circle=[sin(tt) ; cos(tt)]*100*10^(-9);

for i=1:n

    hold on
    plot(circle(1,:)+x(i),circle(2,:)+y(i))
   

end
hold off

for step=1:steps


    Lx=0;
    Ly=0;
    Lz=0;
    Fm=0;
    Ftotal=0;
    DXB=0;
    DYB=0;
    DZB=0;
    

    fuerzas_mgneticas
    fuerzas_contacto
    Fm=F_dip;
    Fm=Fp;


%Ftotal=Fm+F_contact+Fh;
Ftotal=Fm+F_contact+Fh;
Mtotal=0*M_contact;
%Ftotal=F_contact;
%Ftotal=Fm;
step
for q=1:n %Integración por el método de euler
       
       ax(q)=Ftotal(3*q-2)/masa;
       ay(q)=Ftotal(3*q-1)/masa;
       az(q)=Ftotal(3*q)/masa;
       vx(q)= Vx(q,step)+dt*ax(q);  
       vy(q)= Vy(q,step)+dt*ay(q);
       vz(q)= Vz(q,step)+dt*az(q);
       x(q)= X(q,step)+dt*Vx(q,step);
       y(q)= Y(q,step)+dt*Vy(q,step);
       z(q)= Z(q,step)+dt*Vz(q,step);
       alfax(q)=Mtotal(3*q-2)/I;
       alfay(q)=Mtotal(3*q-1)/I;
       alfaz(q)=Mtotal(3*q)/I;
       wx(q)=WX(q,step)+dt*alfax(q);
       wy(q)=WY(q,step)+dt*alfay(q);
       wz(q)=WZ(q,step)+dt*alfaz(q);
       
       X(q,step+1)=x(q);
       Y(q,step+1)=y(q);
       Z(q,step+1)=z(q);
       Vx(q,step+1)=vx(q);
       Vy(q,step+1)=vy(q);
       Vz(q,step+1)=vz(q);
       Ax(q,step+1)=ax(q);
       Ay(q,step+1)=ay(q);
       Az(q,step+1)=az(q);
       WX(q,step+1)=wx(q);
       WY(q,step+1)=wy(q);
       WZ(q,step+1)=wz(q);
       MX(q,step+1)=Mtotal(3*q-2);
       MY(q,step+1)=Mtotal(3*q-1);
       MZ(q,step+1)=Mtotal(3*q);
       FTx(q,step+1)=Ftotal(3*q-2);
       FTy(q,step+1)=Ftotal(3*q-1);
       FTz(q,step+1)=Ftotal(3*q);
       FMX(q,step+1)=Fm(3*q-2);
       FMY(q,step+1)=Fm(3*q-1);
       FMZ(q,step+1)=Fm(3*q);
       FCX(q,step+1)=F_contact(3*q-2);
       FCY(q,step+1)=F_contact(3*q-1);
       FCZ(q,step+1)=F_contact(3*q);
       FCTX(q,step+1)=F_T(3*q-2);
       FCTY(q,step+1)=F_T(3*q-1);
       FCTZ(q,step+1)=F_T(3*q);
       FCNX(q,step+1)=F_N(3*q-2);
       FCNY(q,step+1)=F_N(3*q-1);
       FCNZ(q,step+1)=F_N(3*q);
       FHX(q,step+1)=Fh(3*q-2);
       FHY(q,step+1)=Fh(3*q-1);
       FHZ(q,step+1)=Fh(3*q);
       
       sumatoriafx=sum(FTx(:,step+1));
       sumatoriafy=sum(FTy(:,step+1));
       sumatoriafz=sum(FTz(:,step+1));
       sumatoriafxm=sum(FMX(:,step+1));
       sumatoriafym=sum(FMY(:,step+1));
       sumatoriafzm=sum(FMZ(:,step+1));
       sumatoriaMx=sum(MX(:,step+1));
       sumatoriaMy=sum(MY(:,step+1));
       sumatoriaMz=sum(MZ(:,step+1));

       plot(circle(1,:)+x(q),circle(2,:)+y(q))
       hold on
end
hold off  
F(step) = getframe;
time(step+1) = time(step) + dt;
step

end

graficas
