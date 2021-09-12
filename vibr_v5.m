function PT
global a m k0 k eta dt Vt
m=0.1;
k0=1.e3;
k=k0/3.;
V_=1.; % ------- amplitude

k_=k*k0/(k+k0); % for the case of eta=0
om_=sqrt(k_/m); 
nu_=om_/2./pi % ------- resonant friquency in Hz
T_=1./nu_

n=2;
om=zeros(n,1);
del=zeros(n,1);
om0=2*pi*10;
om1=2.*pi*100.;
dom=(om1-om0)/(n-1);
for j=0:n-1
   om(j+1)=om0+dom*j;
end

eta=1.e3; % ------- !!!!!!!
T=.5;
nT=2000; % steps
dt=T/nT;
t=(1:nT)*dt;
a=4.*m/dt/dt/k0+(k+k0)/eta/k0*2.*m/dt;

amplituda_ = zeros(n,1);
amplituda = zeros(n,1);
    
z=zeros(nT,1);
v=zeros(nT,1);
F=zeros(nT,1);

zt=zeros(nT,1);
vt=zeros(nT,1);
Ft=zeros(nT,1);

for jj=1:n
   [z_ v_ F_]=sol(om(jj),V_);
   z(1)=z_;
   v(1)=v_;
   F(1)=F_;
   
   zt(1)=z_;
   vt(1)=v_;
   real(v_)
   Ft(1)=F_;
   for j=1:nT-1
       Vt=V_*exp(i*om(jj)*(t(j)-dt/2.));
       [z(j+1) v(j+1) F(j+1)]=step(z(j),v(j),F(j));
       % ------- precise solution      
       zt(j+1)=z_*exp(i*om(jj)*t(j)); 
       vt(j+1)=v_*exp(i*om(jj)*t(j));
       Ft(j+1)=F_*exp(i*om(jj)*t(j));
   end
 
   
   
   
figure % ------- test
[po]=plot(t,real(vt));
set(po,'linewidth',2);
colormap hsv;                
grid on;
hold on;

 figure
 [pf]=plot(t,real(v),'b-');
 set(pf,'linewidth',2);
 colormap hsv;
 grid on;
 hold on;

 %  del(jj)=norm(z(j)-zt(j),v(j)-vt(j),F(j)-Ft(j)); % ????????
   d1=norm(zt(j),v(j),F(j));
   del(jj)=del(jj)/d1;
   
   rz = zeros(nT, 1);
   for j=1, nT
        rz(j)=real(z_(j));
   end

   amplituda_(jj) = abs(max(rz) - min(rz)) ;
   amplituda(jj) = abs(max(zt) - min(zt)) ;
   
end

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
function d=norm(z,v,F);
global a m k0 k eta dt Vt
d=m*v*v+F*F/k0+k*(z-F/k0)*(z-F/k0);
d=sqrt(d/2.);