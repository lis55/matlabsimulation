
    number_of_contacts= 0;
    particles_in_contact=0;
    F_contact(1:3*n)=0;
    M_contact(1:3*n)=0;
    F_T(1:3*n)=0;
    F_N(1:3*n)=0;
    
for i=1:n
    % Solapamiento
    solap = zeros(1,n);

    for j=1:n
    
        if (i~=j)
            dx=x(i)-x(j);  
            dy=y(i)-y(j);
            dz=z(i)-z(j);
            rij=sqrt(dx^2+dy^2+dz^2);
            solap(j)=(2*a)-rij;
                  nij=[dx dy dz]*1/rij; 

        else    
            % Una partícula no puede estar en contacto consigo misma
        solap(j) = -1;
        end
        
        SOLAP(i,j)=solap(j);
        
       
    end

    temp = find(solap>=0);  % store into array temp the identities of the particles 
    % that overlap with particle i, i.e. those whose delta is >0. The
    % identities of each particle are equal to the corresponding indices
    % in the delta array.

    % store the identity number of the particles that are in contact with
    % particle i
    particles_in_contact(i,1:length(temp)) = temp;
    
    % temporary variable
    no_of_cont = length(find( particles_in_contact(i,:)~=0 ));
    
    % store the number of particles that are in contact with the i
    % particles
    number_of_contacts(i) = no_of_cont;
    
    
end

for i=1:n 
    F_nmatrix=zeros(3,n);
    F_tmatrix=zeros(3,n);
    tor=zeros(3,n);

    if(number_of_contacts(i) ~= 0) 
    temp = find(particles_in_contact(i,:)~=0);
    c_particle = particles_in_contact(i,temp);

    
   
             for k=1:length(c_particle)  % (j loop)
                 F_norm=zeros(1,n);
                 
            j = c_particle(k); % partícula que se encuentra en contacto con la partícula i
                        %% ============ Fuerza normal =============
                        
                  %% término elástico
                              
                  dx=x(i)-x(j);  
                  dy=y(i)-y(j);
                  dz=z(i)-z(j);
                  rij=sqrt(dx^2+dy^2+dz^2);%ojo
                  sol=(2*a)-rij;
                  nij=[dx dy dz]*1/rij; 
                  vn=dot([vx(i)-vx(j) vy(i)-vy(j) vz(i)-vz(j)],nij)*nij; 
                  dsol=-dot(vn,nij); 
                  fcn_ela=kn*sol^(3/2);       
            
                  %% término viscoso
                  
                  fcn_vis=(bn*sqrt(sol)*dsol); 
                
               if(fcn_ela+fcn_vis>=0)   
                   
                  F_norm(1)=(fcn_ela+fcn_vis)*nij(1);
                  F_norm(2)=(fcn_ela+fcn_vis)*nij(2);
                  F_norm(3)=(fcn_ela+fcn_vis)*nij(3);
                  %F_norm(:,j)=fcn1_ela*nij+(fcn1_vis*vn*(1/norm(vn)))
                  F_nmatrix(1,j)=F_norm(1); % Almacenaje de Fnij
                  F_nmatrix(2,j)=F_norm(2);
                  F_nmatrix(3,j)=F_norm(3);
                  
               else
                   F_norm(1)=0; 
                   F_norm(2)=0;
                   F_norm(3)=0; 
                   F_nmatrix(1,j)=F_norm(1); % Almacenaje de Fnij
                   F_nmatrix(2,j)=F_norm(2);
                   F_nmatrix(3,j)=F_norm(3);
               end
                  
                  %% ============ Fuerza tangencial =============

wi=0;
wj=0;
wi=[wx(i) wy(i) wz(i)];
wj=[wx(j) wy(j) wz(i)];
vtij=[vx(i)-vx(j) vy(i)-vy(j) vz(i)-vz(j)]-[vn(1) vn(2) vn(3)]-a*(cross(nij,wi)+cross(nij,wj));

if(norm(vtij)~=0)
tij=[vtij(1) vtij(2) vtij(3)]*1/norm(vtij);

            % Fuerza elástica
            zeta=norm(vtij*dt);
            Ft_el=kn*zeta;
            
            %% Fuerza viscosa
            %bt_visc = 2*sqrt(kt_stiff*mm);
            bt=bn;

            ftvis=bt*norm(vtij)-0*Ft_el;
            fr=mu*norm(F_norm);              
           
            if(ftvis<fr)
                  F_tmatrix(1,j)=-ftvis*tij(1);
                  F_tmatrix(2,j)=-ftvis*tij(2);
                  F_tmatrix(3,j)=-ftvis*tij(3);
            else
                F_tmatrix(1,j)=-fr*tij(1);
                F_tmatrix(2,j)=-fr*tij(2);   
                F_tmatrix(3,j)=-fr*tij(3);
            end
else
    F_tmatrix(1,j)=0;
    F_tmatrix(2,j)=0;
    F_tmatrix(3,j)=0;
    
end                      
     tor(:,j)=cross(rij*nij,F_tmatrix(:,j)+F_nmatrix(:,j));   
     
        end %%%% Fin del ciclo de partículas en contacto (ciclo j )
        
        %% ====================== Fuerza de contacto total ======================
        
        F_T(3*i-2)=sum(F_tmatrix(1,:));
        F_T(3*i-1)=sum(F_tmatrix(2,:));
        F_T(3*i)=sum(F_tmatrix(3,:));       
        
        F_N(3*i-2)=sum(F_nmatrix(1,:));
        F_N(3*i-1)=sum(F_nmatrix(2,:));
        F_N(3*i)=sum(F_nmatrix(3,:));
        
       F_contact(3*i-2)= sum(F_nmatrix(1,:))+sum(F_tmatrix(1,:));
       F_contact(3*i-1)=sum(F_nmatrix(2,:))+sum(F_tmatrix(2,:));
       F_contact(3*i)= sum(F_nmatrix(3,:))+sum(F_tmatrix(3,:));
        
        %% ====================== Torque de contacto =====================
        M_contact(3*i-2) = sum(tor(1,:));
        M_contact(3*i-1) = sum(tor(2,:));
        M_contact(3*i) = sum(tor(3,:));
        
    else
        % Las fuerzas y torques son nulos cuando no hay contacto
        F_contact(3*i-2) =0;
        F_contact(3*i-1)=0;
        F_contact(3*i)=0;
        
        F_T(3*i-2)=0;
        F_T(3*i-1)=0;
        F_T(3*i)=0;       
        
        F_N(3*i-2)=0;
        F_N(3*i-1)=0;
        F_N(3*i)=0;
        
        M_contact(3*i-2) = 0;
        M_contact(3*i-1) = 0;
        M_contact(3*i) = 0;
        
    end     
    
end %%%% Fin del ciclo en i