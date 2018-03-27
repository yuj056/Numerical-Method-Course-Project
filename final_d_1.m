clear;
clc;
close all;
%Define the parameters given
 
lx=4;
ly=4;
alpha=0.0001;
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
difnum=2;
dt=difnum*dy^2/alpha;


%Define T matrix T(y,x)
T=zeros(Ny,Nx);
T(1,:)=4*sin(2*pi*k*x(1:end)/lx);
%heat source q(y,t)
q=zeros(Ny,1);
for i=1:Ny;
    q(i)=q0*exp(-lambda*(y(i)-ly/2).^2);
end
q1=0.4*q;

%Initial Condition f(y,t) 
f=zeros(size(q));            %f(y,0)=0
f(1)=4;              %f(0,0)=4
f(end)=0;            %f(ly,0)=0

%%time step

tfCN=1280000000;
t=0:dt:tfCN;
%Define f for different steady state
T1 = T; T2 = T; T3 = T;
fscn1=f;

cricn1=0;

beta=1/2*dt*alpha/dy^2;
gamma=1/2*dt*alpha*4*pi^2*k^2/lx^2;

a=ones(Ny,1);
b=ones(Ny,1);
c=ones(Ny,1);

a=-beta*a;
a(1)=0;
a(end)=0;

b=(1+gamma+2*beta)*b;
b(1)=1;
b(end)=1;

c=-beta*c;
c(1)=0;
c(end)=0;

for i=1:length(t)-1;
    fscn1(:,i+1)=thomas1(a,b,c,fscn1(:,i),q1,beta,gamma,dt);
 T1(:,:,i+1) = zeros(size(T1(:,:,1)));
    
    for j = 1:Ny;
        for m = 1:Nx;
            
            T1(j,m,i+1) = sin(2*pi*k*x(m)/lx)*fscn1(j,i+1);
            
        end
    end
    
          
    if abs(max(T1(:,:,i+1)-T1(:,:,i))) < tol;
        break
    end
    
end
cricn1=i;


figure;
plot(x,(T1(13,:,end)+T1(14,:,end))/2,'o-');
xlabel('X','fontsize',14);
ylabel('T','fontsize',14);
title('Steady Temperature at y=ly/4 by CN Method','fontsize',14,'fontweight','bold');
grid on;
legend('alpha=0.0001');

T1_s = 0.5*(T1(13,9,end)+T1(14,9,end));
