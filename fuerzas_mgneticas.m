S=0;
R_ij=0;
m=0;
F_dip=0;
for i=1:n
r_ij=0;

    for j=1:n
    
        if (i~=j)
            dx=x(i)-x(j);  
            dy=y(i)-y(j);
            dz=z(i)-z(j);
            r_ij(j)=sqrt(dx^2+dy^2+dz^2);

        else    
       r_ij(j) = 0;
        end
        
        R_ij(i,j)=r_ij(j);
        
       
    end  
    
end

for i=1:n
  
    for k=1:d
           fil=d*(i-1)+k;
           
           B(fil)=bo(k);
           dxB(fil)=dxbo(k);
           dyB(fil)=dybo(k);
           dzB(fil)=dzbo(k);
    end

    for j=1:n      
        
        A=0; 
        
       if i==j
           
           A(1,1)=1/c2;
           A(2,2)=1/c2;
           A(3,3)=1/c2; 
           
       else
           if(R_ij(i,j)<5*a)
           dx=x(i)-x(j);
           dy=y(i)-y(j);
           dz=z(i)-z(j);
           dr=R_ij(i,j);
           
           A(1,1)=-c1*(3*dx^2/dr^5-1/dr^3);
           A(1,2)=-c1*(3*dx*dy/dr^5);
           A(1,3)=-c1*(3*dx*dz/dr^5);
           A(2,1)=-c1*(3*dx*dy/dr^5);
           A(2,2)=-c1*(3*dy^2/dr^5-1/dr^3);
           A(2,3)=-c1*(3*dy*dz/dr^5);
           A(3,1)=-c1*(3*dx*dz/dr^5);
           A(3,2)=-c1*(3*dy*dz/dr^5);
           A(3,3)=-c1*(3*dz^2/dr^5-1/dr^3);
       else
           A(1,1)=0;
           A(2,2)=0;
           A(3,3)=0;
       end
       end
       
   for k=1:d

       for u=1:d
           col=d*(j-1)+u;
           fil=d*(i-1)+k;
           
           S(col,fil)=A(k,u);
       end
       
   end
          
    end
    
end

%%%% Calculo del vector m
m=(S\(B'));

for i=1:n
    
    fdip(1)=0;
    fdip(2)=0;
    fdip(3)=0;
    F_dmtrix(1,:)=0;
    F_dmtrix(2,:)=0;
    F_dmtrix(3,:)=0;
    for j=1:n
        if i~=j
            dx=x(i)-x(j);
            dy=y(i)-y(j);
            dz=z(i)-z(j);
            rij=sqrt(dx^2+dy^2+dz^2);
            f1=[dx;dy;dz]*(m(3*i-2)*m(3*j-2)+m(3*i-1)*m(3*j-1)+m(3*i)*m(3*j))/rij;
            f2=[m(3*i-2);m(3*i-1);m(3*i)]*(dx*m(3*j-2)+dy*m(3*j-1)+dz*m(3*j))/rij;
            f3=[m(3*j-2);m(3*j-1);m(3*j)]*(dx*m(3*i-2)+dy*m(3*i-1)+dz*m(3*i))/rij;
            f4=-5/rij^3*[dx;dy;dz]*(dx*m(3*i-2)+dy*m(3*i-1)+dz*m(3*i))*(dx*m(3*j-2)+dy*m(3*j-1)+dz*m(3*j));
            fdip(:)=3*u0/(4*pi*rij^4)*(f1+f2+f3+f4);
            F_dmtrix(1,j)=fdip(1);
            F_dmtrix(2,j)=fdip(2);
            F_dmtrix(3,j)=fdip(3);

        end
   
    end

   F_dip(3*i-2)=sum(F_dmtrix(1,:)); 
   F_dip(3*i-1)=sum(F_dmtrix(2,:));
   F_dip(3*i)=sum(F_dmtrix(3,:));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Segundo metodo

for i=1:n


    for j=1:n
lx=0;
ly=0;
lz=0;
        if i==j
           lx(1,1)=0;
           lx(2,2)=0;
           lx(3,3)=0;
           ly(1,1)=0;
           ly(2,2)=0;
           ly(3,3)=0;
           lz(1,1)=0;
           lz(2,2)=0;
           lz(3,3)=0;
           
        else
           if(R_ij(i,j)<5*a) %%ojo
           dx=x(i)-x(j);
           dy=y(i)-y(j);
           dz=z(i)-z(j);
           dr=R_ij(i,j);
          
           lx(1,1)=c1*(9*dx/dr^5-15*dx^3/dr^7);
           lx(2,1)=c1*(3*dy/dr^5-15*dx^2*dy/dr^7);
           lx(3,1)=c1*(3*dz/dr^5-15*dx^2*dz/dr^7);
           lx(1,2)=c1*(3*dy/dr^5-15*dx^2*dy/dr^7);
           lx(2,2)=c1*(-15*dy^2*dx/dr^7+3*dx/dr^5);
           lx(3,2)=c1*(-15*dx*dy*dz/dr^7);
           lx(1,3)=c1*(3*dz/dr^5-15*dx^2*dz/dr^7);
           lx(2,3)=c1*(-15*dx*dy*dz/dr^7);
           lx(3,3)=c1*(3*dx/dr^5-15*dx*dz^2/dr^7);
           
           ly(1,1)=c1*(-15*dx^2*dy/dr^7+3*dy/dr^5);
           ly(2,1)=c1*(-15*dy^2*dx/dr^7+3*dx/dr^5);
           ly(3,1)=c1*(-15*dy*dx*dz/dr^7);
           ly(1,2)=c1*(-15*dy^2*dx/dr^7+3*dx/dr^5);
           ly(2,2)=c1*(9*dy/dr^5-15*dy^3/dr^7);
           ly(3,2)=c1*(3*dz/dr^5-15*dy^2*dz/dr^7);
           ly(1,3)=c1*(-15*dy*dx*dz/dr^7);
           ly(2,3)=c1*(3*dz/dr^5-15*dy^2*dz/dr^7);
           ly(3,3)=c1*(3*dy/dr^5-(15*dy*dz^2/dr^7));
           
           lz(1,1)=c1*(-15*dx^2*dz/dr^7+3*dz/dr^5);
           lz(2,1)=c1*(-15*dy*dz*dx/dr^7);
           lz(3,1)=c1*(3*dx/dr^5-15*dz^2*dx/dr^7);
           lz(1,2)=c1*(-15*dy*dz*dx/dr^7);
           lz(2,2)=c1*(3*dz/dr^5-15*dz*dy^2/dr^7);
           lz(3,2)=c1*(3*dy/dr^5-15*dz^2*dy/dr^7);
           lz(1,3)=c1*(3*dx/dr^5-15*dz^2*dx/dr^7);
           lz(2,3)=c1*(3*dy/dr^5-15*dz^2*dy/dr^7);
           lz(3,3)=c1*(9*dz/dr^5-15*dz^3/dr^7);
           else
           lx(1,1)=0;
           lx(2,2)=0;
           lx(3,3)=0;
           ly(1,1)=0;
           ly(2,2)=0;
           ly(3,3)=0;
           lz(1,1)=0;
           lz(2,2)=0;
           lz(3,3)=0;     
           end 
       end
       
   for k=1:d
       for u=1:d
   
           col=d*(j-1)+u;
           fil=d*(i-1)+k;
           
           Lx(col,fil)=-lx(k,u);
           Ly(col,fil)=-ly(k,u);
           Lz(col,fil)=-lz(k,u);
           
       end
   end
    end
end

DXB=Lx*m+dxB';
DYB=Ly*m+dyB';
DZB=Lz*m+dzB';

%%%% Cálculo de la fuerza mágnetica 
for i=1:n

    fm=0;
    fm(1)=(m(3*i-2)*DXB(3*i-2))+(m(3*i-1)*DYB(3*i-2))+(m(3*i)*DZB(3*i-2));
    fm(2)=(m(3*i-2)*DXB(3*i-1))+(m(3*i-1)*DYB(3*i-1))+(m(3*i)*DZB(3*i-1));
    fm(3)=(m(3*i-2)*DXB(3*i))+(m(3*i-1)*DYB(3*i))+(m(3*i)*DZB(3*i));


     for k=1:d
  
           fil=d*(i-1)+k;
           
           Fp(fil)=fm(k);%ojo
           
       end
end