clc
clear all
% global a b k  h mu p h0 B m n y0 x1 x2 p0  q0 dxx1 dyy1 c1 dxx2 c2 sigma
% a=-0.3;k=0;h=0.2;sigma=0.1;e=0.6;
mu=0.0001;p=0.5;
q1=[];q3=[];q=[];

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

% d=dir('*.csv')
% n1=length(d)
% data=cell(1,n1);
% i=100;
% data{1,i}=csvread(d(i).name);
% B=data{1,i};
% [nodf,qb,Nm] = cal_structure(B);








%%%%%%%%%deletion Plant%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % pa1=[0.05:0.05:0.95];
% % 
% % for l=1:length(pa1)
% % 
% %     load A1.dat
% % B=A1;
% % [n m]=size(B);
% % [n2 m2]=size(B);
% % for i=1:n
% %     for j=1:m
% % if B(i,j)>0
% %     B(i,j)=1;
% % else B(i,j)=0;
% % end
% %     end
% % end
% % 
% % b=eye(m);
% % b1=eye(n);
% % k1=sum(B,2);
% % k2=sum(B,1);
% % 
% % if length(k1)>round(pa1(l)*n2)
% % [iii,d]=mink(k1,round(pa1(l)*n2));
% % % d=randperm(length(k1),round(pa1(l)*n2))
% % % fid = fopen(sprintf('frac_temp_1000.dat'),'a' )
% % % [iii,d]=maxk(k1,7);
% % % 
% % % % d=[];
% %     B([d],:)=[];
% % 
% % k1=sum(B,2); %% degree plant
% % 
% % k2=sum(B,1)'; %% degree Animal
% % [ii1 d]=find(k2(:,1)==0);
% % B1=B;
% % % for i=1: length(ii1)
% %     B(:,[ii1])=[];
% % % end
% % k1=sum(B,2);
% % k2=sum(B,1);
% % [n m]=size(B);
% % % % 
% % % [n m]=size(B);
% % b=eye(m);
% % b1=eye(n);
%%%%%%%%%%%deletion Pollinator%%%%%%%%%%%%%%%%%%%%%%%%%%%

pa1=[0.05:0.05:0.95];

for l=1:length(pa1)
    
        load A1.dat
B=A1;
[n m]=size(B);
[n2 m2]=size(B);
for i=1:n
    for j=1:m
if B(i,j)>0
    B(i,j)=1;
else B(i,j)=0;
end
    end
end

b=eye(m);
b1=eye(n);
k1=sum(B,2);
k2=sum(B,1);
    
    
    
    
    
if length(k2)>round(pa1(l)*m2)
[iii,d]=mink(k2,round(pa1(l)*m2)); %%Specialist
% [iii,d]=maxk(k2,round(pa1(l)*m2)); %%Generalist
% d=randperm(length(k2),round(pa1(l)*m2)); %%Random

    B(:,[d])=[];
k1=sum(B,2); %% degree plant

k2=sum(B,1)'; %% degree Animal
[ii1 d]=find(k1(:,1)==0);
B1=B;
% for i=1: length(ii1)
    B([ii1],:)=[];
% end
k1=sum(B,2);
k2=sum(B,1);
[n m]=size(B);
% % 
% [n m]=size(B);
b=eye(m);
b1=eye(n);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k=0.1;
% ts=0:0.05:3;

 t=0;t_max =400; dt=0.05;
 m1=t_max/dt;
 
  t=0.0;
%   gama1=gamma;
  gama1=3:-0.1:0;
%     gama1=3:-0.1:0;
%      while t  t+dt 
   for jjj=[1 6 11 31 36 41]
       %         q=[];
         for jjjj=1:length(gama1)
             jjjj
            
        
    g=gama1(jjjj); 
    B1=[];
    B2=[];
for ii=1:m
    B1(:,ii)=(B(:,ii)./(k2(ii)^p))*g;
end
for ii=1:n
    B2(ii,:)=(B(ii,:)./(k1(ii)^p))*g;
end
   y0=[];
 y0=[4*rand(m,1); 4*rand(n,1)];
 y0=y0';
 y0=reshape(y0,[1,m+n]);
 p0=y0(1:m);
q0=y0(m+1:m+n);
x=p0;
y=q0; 

   for i=1:m1  
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
         
         
         %%%%%%%%%%%%%%%%%%%without noise%%%%%%%%%%%%%%%%%%
        x=x+ (a(jjj)*x-k(jjj)*x-diag(B3)'+mu+c1'.*x)*dt; 
        y=y+(a(jjj)*y-diag(B4)'+mu+c2'.*y)*dt;
%         x(:,1)=0.2;
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
     q=[q; pa1(l) T(jjj) gama1(jjjj) sum(x)/m n m];
%      if max(q)<0.1 
%      q3=[q3;pa1(l),T(jjj),g,m,n,nodf ];
%      fprintf(fid,'%f %f %f %f\n',pa1(l),T(jjj),g);
%      break
%      else 
%   
% continue
%      end
     
     
     
%      q1=[q1;k y];
         end
end       
end

%          figure
%       plot(q(:,1),q(:,2:m),'r') 
%  hold on
%     plot(q(:,1),q(:,m+1:m+n),'b') 
%     x=q3(:,1);
% y=q3(:,2)-273;
% z=q3(:,3);
% figure;
% scatter(x,y,380,z,'filled')
end      