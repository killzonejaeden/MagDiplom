function vibr_v2
   clear all;close all;
   global dt
mk=[1. 2.5 ]; % ------- adjacent masses [kg]
Jk=[1. 1.5 ]*1.e-5; %inertia moments [kg*m^2]
alk=[4. 2.5 ]*1.e6; % springs stiffness [N/m]
alk0=[4. 2.5 ]*1.e6; 
btk0=[4000. 2500. ];
btk=[4000. 2500. ]/3;
etak=[1.e-3, 1.e-3];
gmk=[.4 .25 ]*1.e-3;
nm=numel(mk);

l=1.; % ------- length of a rod [m]
h=.005; % ------- thikness [m]
w=.01; % ------- width [m]
xk=[.25 .75];%points of masses connection
ro=7850.; % ------- density [kg/m^3]
S=w*h; % ------- cross section [m^2]
E=2.1e11; % ------- Young modulus [Pa]
Iz=h*h*h*w/12.; % moment of inertia [m^4]
D=E*Iz/S; % ------- bending stiffness [Pa*m^2]
c=sqrt(E/ro);
%nu=225 res gor1   5 vert
nu=4; %frequency of loading
%nu=400;
om=2.*pi*nu;

N_time=200000; % ------- steps in time variable
n=100; % ------- nodes of a grid
dx=l/n; % ------- step in space variable
pK=1.; %Courant parameter
dt=pK*dx/c; % ------- time step

u=zeros(1,n+1); % ------- displacements
u_=zeros(1,n+1);%for next time step
w=zeros(1,n+1);
w_=zeros(1,n+1);
wwx=zeros(N_time,n+1); % for graph ainmation
uux=zeros(N_time,n+1); % for graph ainmation
fu=zeros(1,n+1); % ------- right hand sides
fw=zeros(1,n+1);
ffi=zeros(1,n+1);

uk=zeros(1,nm); % displacements of masses
ut=zeros(1,nm); % ------- velocities
wk=zeros(1,nm);
wt=zeros(1,nm);
fik=zeros(1,nm); % ------- angles of rotation
fit=zeros(1,nm);

a=zeros(1,n+1);
b=zeros(1,n+1);
c=zeros(1,n+1);
d=zeros(1,n+1);
e=zeros(1,n+1);
f=zeros(1,n+1);

x=(0:n)*dx;

for j=3:n-1 % ------- coefficients of w_xxxx
   a(j)=dt*dt*D/ro/dx^4;
   b(j)=-4.*dt*dt*D/ro/dx^4;
   c(j)=1.+6.*dt*dt*D/ro/dx^4;
   d(j)=b(j);
   e(j)=a(j);
end

for k=1:nm % ------- nearby nodes for masses
   [xm(k), km(k)]=min(abs(x(:)-xk(k)));
end

t=dt;
%u_(1)=sin(om*dt); % ------- given displacement
u(1)=0.; % ------ fixed point

f_t=fopen('t.dat','w');
f_u=fopen('u.dat','w');
f_w=fopen('w.dat','w');
for k=1:nm % open files for each point mass 
  f_uk(k)=fopen(['uk',num2str(k),'.dat'],'w');
  f_wk(k)=fopen(['wk',num2str(k),'.dat'],'w');
end

k_in=10; % ------- node for output into files

wmax=0; % max amplitude for w in graph
umax=0;

for it=1:N_time  % ------- solution of equations
   fprintf(f_t,'%e ',t);
   t=t+dt;

   for k=1:nm
      [uk(k) ut(k) fu(k)] = step(uk(k),ut(k),fu(k), mk(k), alk0(k), alk(k), etak(k),(u(km(k))-u_(km(k)))/dt );
      [wk(k) wt(k) fw(k)] = step(wk(k),wt(k),fw(k), mk(k), btk0(k), btk(k), etak(k),(w(km(k))-w_(km(k)))/dt);
   %   [fik(k) fit(k) ffi(k)] = step(fik(k),fit(k),ffi(k),mk(k), btk0(k), btk(k), etak(k),ffi(k));
    %  [uk(k) ut(k)]=dudt(dt,alk(k)/mk(k), uk(k),ut(k),alk(k)*u_(km(k))/mk(k));
    %  [wk(k) wt(k)]=dudt(dt,btk(k)/mk(k), wk(k),wt(k),btk(k)*w_(km(k))/mk(k));
    [fik(k) fit(k)]=dudt(dt,gmk(k)/Jk(k), fik(k),fit(k),...
     gmk(k)*(w_(km(k)+1)-w_(km(k)-1))/Jk(k)/2./dx);

      fprintf(f_uk(k),'%e ',uk(k)); %  mass 
      fprintf(f_wk(k),'%e ',wk(k)); %  mass 
   end

 %if it==1
 u(1)=sin(om*t)*(10e-6);
% end    
%   u(1)=0.;
   for j=2:n %
      df=0.;
      for k=1:nm
         if j==km(k)
            df=df+alk(k)*(uk(k)-u_(j))/S;
            cf=1.;
         end
      end
      df=E*(u_(j+1)-2.*u_(j)+u_(j-1))/dx/dx+df;
      u(j)=2.*u_(j)-u(j)+dt*dt*df/ro;
   end       
   u(n+1)=0.; %

   fprintf(f_u,'%e ',u_(k_in));
   
   c(1)=1.; % 
   c(2)=1.;
   d(2)=-1.;
%   if it==1
f(1)=sin(om*t)*10e-6;
%   end

   % a(1)=1.; % ------- M=0
   % b(1)=-2.;
   % c(1)=1.;
   % f(1)=1.;

   % a(2)=1.; % ------- Q=0
   % b(2)=-3.;
   % c(2)=3.;
   % e(2)=-1.;
   % f(2)=1.;

   c(n+1)=1.;  %
   c(n)=-1.;   % ------- w_x=0
   b(n)=1.;

   %c(n+1)=.1; % ------- M=0
   %d(n+1)=-2.;
   %e(n+1)=1.;
   %f(n+1)=-1.;

   %b(n)=1.; % ------- Q=0
   %c(n)=-3.;
   %d(n)=3.;
   %e(n)=-1.;
   %f(n)=-10.e5;

   f(3:n-1)=0.;
   for j=3:n-1 %
      df=0.;
      dfp=0.; 
      dfm=0.; 
      for k=1:nm
         if j==km(k)
            df=df+btk(k)*(wk(k)-w_(j))/S; % было alk
            dfik=fik(k)-(w_(j+1)-w_(j-1))/dx/2.;
            dfp=dfp-gmk(k)*dfik/S/2.;
            dfm=dfm+gmk(k)*dfik/S/2.;
         end
      end
      f(j)=f(j)+2.*w_(j)-w(j);
      f(j)=f(j)+dt*dt*df/ro;
      f(j+1)=dt*dt*dfp/ro;
      f(j-1)=f(j-1)+dt*dt*dfm/ro;
   end       

   w=progon5(a,b,c,d,e,f,n+1);
   
   fprintf(f_w,'%e ',w_(k_in)); %);,wk(1)
   
  % fprintf(f_w,'%e ',w_); %);,wk(1)



   for j=1:n+1 %
      uu=u_(j);
      u_(j)=u(j);
      u(j)=uu;
      ww=w_(j);
      
      wwx(it,j) = w_(j); % for graph w
      if abs(w_(j))>wmax      % search max amplitude for graph
        wmax = abs(w_(j));
     end
     uux(it,j) = u_(j)+(j-1)/n;
      if uux(it,j)>umax      % search max amplitude for graph
        umax = uux(it,j);
     end

     w_(j)=w(j);
     w(j)=ww;
  end
end



fclose(f_w);
fclose(f_u);
fclose(f_t);

for k=1:nm % close array of files for mass point
  fclose(f_uk(k));
  fclose(f_wk(k));
end;

tt=load('t.dat');
ywk1 =load('wk1.dat');
ywk2 = load('wk2.dat');
yuk1 =load('uk1.dat');
yuk2 = load('uk2.dat');

%---------------------Graf animations

F=figure('name',' Колебания стержня с двумя массами ');
A=axes(F);

%zxczxczcxzc
filename = 'testAnimated.gif';
%zxczxczxczxc


[P]=plot(A,uux(1,1:n+1),wwx(1,1:n+1) );
hold on;
[P1]=plot(A,[yuk1(1)*50+km(1)/n, yuk1(2)*50+km(2)/n ], [wwx(1,km(1))+ ywk1(1) wwx(1,km(2))+ ywk2(1)],'r*' );
hold on;
%[P2]=plot(A,yuk1(2)+km(2)/n,wwx(1,km(2))+ ywk2(1),'gx' );
% set(P1,'linewidth',4);
hold on;
ylim ([-200*wmax, 200*wmax]); % y ases limits
xlim([0, umax]);

grid on;
xlabel('u(x,t)');ylabel('w(x,t)');

drawnow;

delay = 0.01
for j=2:1000:N_time

  P.YData = wwx(j,1:n+1);
  P.XData = uux(j,1:n+1);

  P1.XData(1) = yuk1(j)*1.e3+km(1)/n;
  P1.XData(2) = yuk2(j)*1.e3+km(2)/n;
  P1.YData(1) = wwx(j,km(1))+(wwx(j,km(1))-ywk1(j))*50;
  P1.YData(2) =  wwx(j,km(2))+(wwx(j,km(2))-ywk2(j))*50;
  title(['t=',num2str(j)]); 
  drawnow;

  %Capture the plot as an image
  frame = getframe(F);
  im = frame2im(frame);
  [imind,cm] = rgb2ind(im,256);
  if j == 2
    imwrite(
      imind,
      cm,
      filename,
      'gif',
      'Loopcount',
      inf
    );
  else
    imwrite(
      imind,
      cm,
      filename,
      'gif',
      'WriteMode',
      'append',
      'DelayTime',
      delay
    );
  end

end



%figure('name',' Резонанс');
%for j=1:length(tt)
%figure('name',['t= ',num2str(t(tt(j)))]);
%plot(x,u(1:Nx,tt(j)));hold on;
%end%j
%grid on;title('\omega=1');
%xlabel('x');ylabel('u(x,t)');


figure('name',' Вертикальные колебания масс относительно стержня');
[pf]=plot(tt,minus(ywk1, wwx(km(1)) ) ,'b');
hold on;
[pf]=plot(tt,minus(ywk2, wwx(km(2)) ) ,'r');
%[pf1]=plot(tt,ywk2-wwx(j,km(2)),'r');
set(pf,'linewidth',1);
colormap hsv;

grid on;
hold on;

figure('name',' Горизонтальные колебания масс относительно стержня');
[pf]=plot(tt, minus(yuk1, uux(km(1)) )  ,'b');
hold on;
[pf]=plot(tt, minus(yuk2, uux(km(2)) ) ,'r');
%[pf1]=plot(tt,ywk2-wwx(j,km(2)),'r');
set(pf,'linewidth',1);
colormap hsv;

grid on;
hold on;

function x=progon5(a,b,c,d,e,f,n);
   p=zeros(1,n);
   q=zeros(1,n);
   r=zeros(1,n);
   p(1)=-b(1)/c(1);
   q(1)=-a(1)/c(1);
   r(1)=f(1)/c(1);
   dl=c(2)+d(2)*p(1);
   p(2)=-(b(2)+d(2)*q(1))/dl;
   q(2)=-a(2)/dl;
   r(2)=f(2)/dl;
   for j=3:n
      pr=d(j)+e(j)*p(j-2);
      dl=c(j)+e(j)*q(j-2)+pr*p(j-1);
      p(j)=-(b(j)+pr*q(j-1))/dl;
      q(j)=-a(j)/dl;
      r(j)=(f(j)-e(j)*r(j-2)-pr*r(j-1))/dl;
   end
   x(n)=r(n);
   x(n-1)=p(n-1)*x(n)+r(n-1);
   for j1=2:n-1
      j=n-j1;
      x(j)=p(j)*x(j+1)+q(j)*x(j+2)+r(j);
   end

%-----------------------------
function [u v]=dudt(dt,k2,u0,v0,f)
   v=(1.-k2*dt*dt/4.)*v0-k2*dt*u0+dt*f;
   v=v/(1.+k2*dt*dt/4.);
   u=u0+dt*(v+v0)/2.;
%------------------
function [z_ v_ F_]=step(z,v,F, m , k0, k ,eta, Vt );
%global a m k0 k eta dt Vt
global dt
a=4.*m/dt/dt/k0+(k+k0)/eta/k0*2.*m/dt;

ztt_=z*((1.+a)/dt-k/2./eta)+a*(v-Vt)-2.*F/dt/k0;
z_=ztt_/((1.+a)/dt+k/2./eta);
v_=2.*(z_-z)/dt-v+2*Vt;
F_=2.*m*(v-v_)/dt-F;