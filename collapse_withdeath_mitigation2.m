
clc
clear all
% global a b k  h mu p h0 B m n y0 x1 x2 p0  q0 dxx1 dyy1 c1 dxx2 c2 sigma
% a=-0.3;k=0;h=0.2;sigma=0.1;e=0.6;
mu=0.0001;p=0.5;
q=[];q1=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mutualistic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%strength%%%%%%%%%%%%%%%%%%%%%%
% T0=293*ones(1,41);sigma=6.5;
% s=2*(sigma)^2*ones(1,41);
% g1=1;
% % g1=g1/s;
% T=273:1:313;
% % T=T';
% % for i=1:length(T)
%  gamma =g1*exp((-(T-T0).^(2))./s);
% %  gamma1=[gamma1 gamma];
% end
% plot(T,gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%birth rate%%%%%%%%%%%%%%%%%

T0=293*ones(1,41);sigma=5;
s=2*(sigma)^2*ones(1,41);
g1=0.35;
% g1=g1/s;
T=273:1:313;
% T=T';
% for i=1:length(T)
a =g1*exp((-(T-T0).^(2))./s);
% %  gamma1=[gamma1 gamma];
% end
% plot(T,a,'Linewidth',1.8)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%handling time%%%%%%%%%%%%%%%%%

T0=293*ones(1,41);sigma=15;
s=2*(sigma)^2*ones(1,41);
g1=0.15;
% g1=g1/s;
T=273:1:313;
% T=T';
% for i=1:length(T)
h =g1*exp(((T-T0).^(2))./s);
% %  gamma1=[gamma1 gamma];
% end
% plot(T,h,'Linewidth',1.8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%death rate%%%%%%%%%%%%%%%%%
T0=293*ones(1,41);sigma=5;
s=2*(sigma)^2*ones(1,41);
g1=0.1;
% g1=g1/s;
T=273:1:313;
% T=T';
% for i=1:length(T)
k =g1*exp((10000*(1./T0-1./T)));

% %  gamma1=[gamma1 gamma];
% end
% plot(T,k,'Linewidth',1.8)

c1=[];c2=[];
% load A1.dat
% B=A1;

d=dir('*.csv')
n1=length(d)
data=cell(1,n1);

%  i=116 %%%low nested
% i=97 %%%%%%%nested0.25
i=59 %%%%%%%nested0.84
data{1,i}=csvread(d(i).name);
B=data{1,i};



[n m]=size(B);
for i=1:n
    for j=1:m
if B(i,j)>0
    B(i,j)=1;
else B(i,j)=0;
end
    end
end
% % b2=0.01+0.04*rand(m^2,1);
% % b3=0.01+0.04*rand(n^2,1);
% % 
% for i=1:m
%     for j=1:m
% if i==j
%     b(i,j)=1;
% else  
% %     b(i,j)=b2(j+(i-1)*m);
% b(i,j)=0;
%     end
%     end
% end
% for i=1:n
%     for j=1:n
% if i==j
%     b1(i,j)=1;
% else  
% %     b1(i,j)=b3(j+(i-1)*n);
% b1(i,j)=0;
%     end
%     end
% end
b=eye(m);
b1=eye(n);
k1=sum(B,1);
k2=sum(B,2);
% k=0.1;
% ts=0:0.05:3;
y0=[];
 y0=[4*rand(m,1); 4*rand(n,1)];
 y0=y0';
 y0=reshape(y0,[1,m+n]);
 p0=y0(1:m);
q0=y0(m+1:m+n);
x=p0;
y=q0;
 t=0;t_max =400; dt=0.05;
 m1=t_max/dt;
 
%  k=0; % vary epsilon1/epsilon2 in the range [0.5 2.6] and y variable is v
% k_max=2;
% n_max=m1;
% dk=(k_max-k)/n_max;
% 
% rng('shuffle')
% gamma1=[];
% gamma=[];
% % T0=293;s=2*(11.5)^2;
% g1=0.27;
% g1=g1/s;
% T=273:1:313;
% T=T';
% for i=1:length(T)
%  gamma(i) ={g1*exp(T(i)-T0)^2};
% %  gamma1=[gamma1 gamma];
% end
% T0=293*ones(1,41);s=2*(5)^2;
% g1=3;
% g1=g1/s;
% T=273:1:313;
% gamma=[];
% % T=T';
% for i=1:length(T)
%  gamma =g1.*exp(-(T-T0).^2/s);
%  gamma1=[gamma1 gamma];
  t=0.0;
%   gama1=gamma;
  gama1=0:0.1:3;
  for   jjj=33:41
%      while t <  t+dt 
         for jjjj=1:length(gama1)
             jjjj
%        jjj=41;      
    g=gama1(jjjj); 
    B1=[];
    B2=[];
for ii=1:m
    B1(:,ii)=(B(:,ii)./(k1(ii)^p))*g;
end
for ii=1:n
    B2(ii,:)=(B(ii,:)./(k2(ii)^p))*g;
end
          
   for i=1:m1  
%        m2=4000;
%        for j=1:m
%            
% %          for i=1:n
%     c1(j)=B1(:,j)'*y';
%     c1(j)=c1(j)/(1+h*c1(j));
%  end
c1=B1'*y';
c1=c1./(1+h(jjj)*c1);
% for j=1:n
%     c2(j)=B2(j,:)*x';
%     c2(j)=c2(j)/(1+h*c2(j));
% end
c2=B2*x';
c2=c2./(1+h(jjj)*c2);
         B3=b*(x'.*x);
         B4=b1*(y'.*y);
%           if i<m2
       
         k4=k(jjj)*ones(1,m);
         k4(1,1)=0.01;   %%%mitigation2
         %%%%%%%%%%%%%%%%%%without noise%%%%%%%%%%%%%%%%%%
        x=x+ (a(jjj)*x-k4.*x-diag(B3)'+mu+c1'.*x)*dt; 
        y=y+(a(jjj)*y-diag(B4)'+mu+c2'.*y)*dt;
%           else 
%                [x1 i1]=min(x);
%                 k4=k(jjj)*ones(1,m);
%          k4(1,i1)=0.1;   %%%mitigation2
%          %%%%%%%%%%%%%%%%%%%without noise%%%%%%%%%%%%%%%%%%
%         x=x+ (a(jjj)*x-k4.*x-diag(B3)'+mu+c1'.*x)*dt; 
%         y=y+(a(jjj)*y-diag(B4)'+mu+c2'.*y)*dt;
%           end
% xppaut
% x(:,1)=0.2;  %%%%%mitigation1
        x=abs(x);
        y=abs(y);
%%%%%%%%%%%%%%%%%%%%%%%%Demographic noise%%%%%%%%%%%%%%%%%%%%%%%
%   x=x+ (a*x-k*x-diag(B3)'+mu+c1.*x)*dt+e.*randn(1,m).*sqrt(x)*sqrt(dt);
%   y=y+(a*y-diag(B4)'+mu+c2.*y)*dt+e.*randn(1,n).*sqrt(y)*sqrt(dt);
% %%%%%%%%%%%%%%%%%%%%%%%%Environmental noise%%%%%%%%%%%%%%%%%%%%%%%
%    x=x+ (a*x-k*x-diag(B3)'+mu+c1.*x)*dt+sigma.*randn(1,m)*sqrt(dt);
%   y=y+(a*y-diag(B4)'+mu+c2.*y)*dt+sigma.*randn(1,n)*sqrt(dt);
  
%   
%  for  i=1:m
% if x(:,i)<0.00001
% x(:,i)=0.05;
% end
% end
% for i=1:n
% if y(:,i)<0.00001
% y(:,i)=0.05;
% end
% end
 
   end
     q=[q;g T(jjj) sum(x)/m sum(y)/n];
%      q1=[q1;k y];
         end
  end
%       plot(q(:,1),q(:,2:m),'r') 
 hold on
%     plot(q(:,1),q(:,63:end),'b') 
     set(findobj(gcf,'type','axes'),'FontSize',20,'TickDir','out', 'LineWidth', 1.8);
         