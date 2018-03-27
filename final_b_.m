clear;
clc;
close all;
%Define the parameters given
 
lx=4;
ly=4;
alpha=0.01;
lambda=10;
k=4;
Nx=41;
Ny=51;
q0=1;
dx=lx/(Nx-1);
dy=ly/(Ny-1);
x=0:dx:lx;
y=0:dy:ly;
tol=1E-5;

%Define T matrix T(y,x)
T=zeros(Ny,Nx);
T(1,:)=4*sin(2*pi*k*x(1:end)/lx);
%heat source q(y,t)
q=zeros(Ny,1);
for i=1:Ny;
    q(i)=q0*exp(-lambda*(y(i)-ly/2).^2);
end
q1=-0.4*q;
q2=0*q;
q3=0.4*q;

%Initial Condition f(y,t) 
f=zeros(size(q));            %f(y,0)=0
f(1)=4;              %f(0,0)=4
f(end)=0;            %f(ly,0)=0

%time step
dt=0.01;
tf=100;
t=0:dt:tf;
%Define f for different steady state
T1 = T; T2 = T; T3 = T;
fs1=f;
fs2=f;
fs3=f;
cri1=0;
cri2=0;
cri3=0;

for i=1:length(t)-1;
     df = zeros(size(fs1(:,1)));
    % every column stands for one time step
    for j=2:Ny-1;
        % every row in each column stands for different position y
        df(j)=alpha*(-4*pi.^2*k.^2/(lx.^2)*fs1(j,i)+(fs1(j+1,i)-2*fs1(j,i)+fs1(j-1,i))/(dy.^2))+q1(j); 
    end
    fs1(:,i+1)=dt*df+fs1(:,i);
   T1(:,:,i+1) = zeros(size(T1(:,:,1)));
    
    for j = 1:Ny;
        for m = 1:Nx;
            
            T1(j,m,i+1) = sin(2*pi*k*x(m)/lx)*fs1(j,i+1);
            
        end
    end
    
          
    if abs(max(T1(:,:,i+1)-T1(:,:,i))) < tol;
        break
    end
   
end
cri1=i;
for i=1:length(t)-1;
     df = zeros(size(fs2(:,1)));
    % every column stands for one time step
    for j=2:Ny-1;
        % every row in each column stands for different position
        df(j)=alpha*(-4*pi.^2*k.^2/(lx.^2)*fs2(j,i)+(fs2(j+1,i)-2*fs2(j,i)+fs2(j-1,i))/(dy.^2))+q2(j); 
    end
    fs2(:,i+1)=dt*df+fs2(:,i);
    T2(:,:,i+1) = zeros(size(T2(:,:,1)));
    
    for j = 1:Ny;
        for m = 1:Nx;
            
            T2(j,m,i+1) = sin(2*pi*k*x(m)/lx)*fs2(j,i+1);
            
        end
    end
    
          
    if abs(max(T2(:,:,i+1)-T2(:,:,i))) < tol;
        break
    end
   
end
cri2=i;

for i=1:length(t)-1;
     df = zeros(size(fs3(:,1)));
    % every column stands for one time step
    for j=2:Ny-1;
        % every row in each column stands for different position
        df(j)=alpha*(-4*pi.^2*k.^2/(lx.^2)*fs3(j,i)+(fs3(j+1,i)-2*fs3(j,i)+fs3(j-1,i))/(dy.^2))+q3(j); 
    end
    fs3(:,i+1)=dt*df+fs3(:,i);
    T3(:,:,i+1) = zeros(size(T3(:,:,1)));
    
    for j = 1:Ny;
        for m = 1:Nx;
            
            T3(j,m,i+1) = sin(2*pi*k*x(m)/lx)*fs3(j,i+1);
            
        end
    end
    
          
    if abs(max(T3(:,:,i+1)-T3(:,:,i))) < tol;
        break
    end
   
end
cri3=i;




subplot(3,1,1);
plot(x,(T1(13,:,end)+T1(14,:,end))/2,'o-');
xlabel('X','fontsize',14);
ylabel('T','fontsize',14);
title('Steady Temperature at y=ly/4 by EE Method','fontsize',14,'fontweight','bold');
grid on;
legend('q0=-0.4');

subplot(3,1,2);
plot(x,(T2(13,:,end)+T2(14,:,end))/2,'o-');
xlabel('X','fontsize',14);
ylabel('T','fontsize',14);
grid on;
legend('q0=0');

subplot(3,1,3);
plot(x,(T3(13,:,end)+T3(14,:,end))/2,'o-');
xlabel('X','fontsize',14);
ylabel('T','fontsize',14);
grid on;
legend('q0=0.4');

