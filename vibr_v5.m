function PT
global a m k0 k eta dt Vt
m=0.1; % weight
k0=1000; % модуль Юнга  (что-то типа упругости)
k=k0/3.; % модуль Юнга  (что-то типа упругости) но тут для другой пружины (333.3333)
V_=1.; % Начальная скорость стенда

nT=2000; % количество отрезков, на которое мы разбиваем время
nOm=100; % - количество рассматриваемых частот

% Для абсолютно упругой системы можно проверить по более простым формулам (крайний случай)
% k_=k*k0/(k+k0); % for the case of eta=0
% om_=sqrt(k_/m);
% nu_=om_/2./pi % ------- resonant friquency in Hz
% T_=1./nu_


om=zeros(nOm,1); % - омега - угловая частота
del=zeros(3, nOm); % - дельта - относительная погрешность (три массива для положения, скорости, силы)
om0=20; % - частота от (2*pi*10)
om1=200; % - частота до (2.*pi*100.)
dom=(om1-om0)/(nOm-1); % - шаг частот
for j=0:nOm-1
  om(j+1)=om0+dom*j; % - массив частот для проверки
end

eta=1000; % вязкость демфера (то насколько гасятся колебания)
T=.5; % сколько всего времени мы смотрим
dt=T/nT; % шаг времени

t=(1:nT)*dt; % массив рассматриваемых точек времени

a=4.*m/dt/dt/k0+(k+k0)/eta/k0*2.*m/dt; % константа часто повторяющаяся в вычислениях

amplituda = zeros(nOm,1);
amplituda_ = zeros(nOm,1);

% цикл по частотам
for jj=1:nOm
  % массивы численных результатов
  z=zeros(nT,1); % положение
  v=zeros(nT,1); % скорость
  F=zeros(nT,1); % сила

  % массивы точных результатов
  zt=zeros(nT,1); % положение
  vt=zeros(nT,1); % скорость
  Ft=zeros(nT,1); % сила

  % точное решение для нулевого времени
  [z_ v_ F_]=sol(om(jj),V_);

  % для точного решения задаем крайние условия
  zt(1)=z_;
  vt(1)=v_;
  Ft(1)=F_;

  % для численного решения задаем крайние условия
  z(1)=z_;
  v(1)=v_;
  F(1)=F_;

  for j=1:nT-1
    % Точное решение для нужного времени получаем по простым формулам
    Ft(j+1)=F_*exp(i*om(jj)*t(j));
    vt(j+1)=v_*exp(i*om(jj)*t(j));
    zt(j+1)=z_*exp(i*om(jj)*t(j));

    Vt=V_*exp(i*om(jj)*t(j)); % текущая скорость стенда

    % численное решение
    [z(j+1) v(j+1) F(j+1)]=step(z(j),v(j),F(j));
  end

  % figure % ------- test
  % [po]=plot(t,real(vt));
  % set(po,'linewidth',2);
  % colormap hsv;
  % grid on;
  % hold on;

  % figure
  % [pf]=plot(t,real(v),'b-');
  % set(pf,'linewidth',2);
  % colormap hsv;
  % grid on;
  % hold on;

  del(1, jj)=norm(z - zt) / norm(zt); % относительная погрешность положения
  del(2, jj)=norm(v - vt) / norm(vt); % относительная погрешность скорости
  del(3, jj)=norm(F - Ft) / norm(Ft); % относительная погрешность силы

  % считаем амплитуду
  % amplituda(jj) = abs(max(zt) - min(zt));
  % amplituda_(jj) = abs(max(z) - min(z));
end

disp('Max delta z:')
disp(max(del(1)))
figure('Name','delta z');
plotDelta=plot(om, del(1, :));
set(plotDelta,'linewidth',2);
colormap hsv;
grid on;
hold on;

disp('Max delta v:')
disp(max(del(2)))
figure('Name','delta v');
[plotDelta]=plot(om, del(2, :));
set(plotDelta,'linewidth',2);
colormap hsv;
grid on;
hold on;

disp('Max delta F:')
disp(max(del(3)))
figure('Name','delta F');
[plotDelta]=plot(om, del(3, :));
set(plotDelta,'linewidth',2);
colormap hsv;
grid on;
hold on;

% figure
% [pom_]=plot(om,amplituda_);
% set(pom_,'linewidth',2);
% colormap hsv;
%                  
% grid on;
% hold on;
% 
% figure
% [pom]=plot(om,amplituda);
% set(pom,'linewidth',2);
% colormap hsv;
%                  
% grid on;
% hold on;


% figure
% [pf]=plot(t,real(z),'b-');
% set(pf,'linewidth',2);
% colormap hsv;
%                    
% grid on;
% hold on;

function [z_ v_ F_]=step(z,v,F);
  global a m k0 k eta dt Vt
  z_=z*((1.+a)/dt-k/2./eta)+a*(v-Vt)-2.*F/dt/k0;
  z_=z_/((1.+a)/dt+k/2./eta);
  v_=2.*(z_-z)/dt-v+2*Vt;
  F_=2.*m*(v-v_)/dt-F;

%-----------------------------
function [u v]=dudt(dt,k2,u0,v0,f)
  v=(1.-k2*dt*dt/4.)*v0-k2*dt*u0+dt*f;
  v=v/(1.+k2*dt*dt/4.);
  u=u0+dt*(v+v0)/2.;
%------------------
function [z_ v_ F_]=sol(om,V_);
  global a m k0 k eta dt Vt
  a11=om*om*m;
  a12=-1;
  b1=V_*i*om*m;
  a22=(k+k0)/k0+i*om*eta/k0;
  a21=-k-i*om*eta;
  b2=0;
  d=a11*a22-a12*a21;
  z_=(b1*a22-b2*a12)/d;
  F_=(-b1*a21+b2*a11)/d;
  v_=-F_/(i*om*m);
% a11*z_+a12*F_-b1
% a21*z_+a22*F_-b2
%------------------

% function d=customNorm(z,v,F);
%   global a m k0 k eta dt Vt
%   d=m*v*v+F*F/k0+k*(z-F/k0)*(z-F/k0);
%   d=sqrt(d/2.);