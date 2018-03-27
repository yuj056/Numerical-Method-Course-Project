function [B]=thomas1(a,b,c,f,q,beta,gamma,dt)

n=length(b);

% construct B for every single time step
B=zeros(n,1);
B(1)=4;
B(end)=0;
B(2)=beta*f(1)+(1-2*beta-gamma)*f(2)+beta*f(3)+dt*q(2)+beta*4;
for j=3:n-1;
B(j)=beta*f(j-1)+(1-2*beta-gamma)*f(j)+beta*f(j+1)+dt*q(j);
end

for i=2:n-1;                               % FORWARD SWEEP
   a(i+1)=-a(i+1)/b(i);               
   b(i+1)=b(i+1)+a(i+1)*c(i);        
   B(i+1)=B(i+1)+a(i+1)*B(i);      
end 

B(n)=B(n)/b(n);                             % BACK SUBSTITUTION
for i=n-1:-1:2,                             
   B(i)=(B(i)-c(i)* B(i+1))/b(i);
end
end
