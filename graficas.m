
%%%% Desplazamientos

subplot 211
hold on
plot(time,X(1,:))
plot(time,X(2,:))
plot(time,X(3,:))
title('Desplazamientos X ')

subplot 212
hold on
plot(time,Y(1,:))
plot(time,Y(2,:))
plot(time,Y(3,:))
title('Desplazamientos Y ')


%%%% fuerzas Magneticas
figure
subplot 221
hold on
plot(time,FMX(1,:))
title('Fuerza magnetica X p1')

subplot 222
hold on
 plot(time,FMX(2,:))
title('Fuerza magnetica X p2')

subplot 223
hold on
plot(time,FMY(1,:))
title('Fuerza magnetica Y p1')

subplot 224
hold on
 plot(time,FMY(2,:))
title('Fuerza magnetica Y p2')

figure
subplot 211
hold on
plot(time,FMX(3,:))
title('Fuerza magnetica X p3')

subplot 212
hold on
plot(time,FMY(1,:))
title('Fuerza magnetica Y p3 ')

%%%% fuerzas contacto
figure
subplot 221
hold on
plot(time,FCX(1,:))
title('Fuerza magnetica X p1')

subplot 222
hold on
 plot(time,FCX(2,:))
title('Fuerza magnetica X p2')

subplot 223
hold on
plot(time,FCY(1,:))
title('Fuerza magnetica Y p1')

subplot 224
hold on
 plot(time,FCY(2,:))
title('Fuerza magnetica Y p2')

figure
subplot 211
hold on
plot(time,FCX(3,:))
title('Fuerza magnetica X p3')

subplot 212
hold on
plot(time,FCY(1,:))
title('Fuerza magnetica Y p3 ')