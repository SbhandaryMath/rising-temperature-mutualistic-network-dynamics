
clc
clear all

mu=0.0001;p=0.5;
q=[];q1=[];

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

d=dir('*.csv')
n1=length(d)
data=cell(1,n1);

%  i=116 %%%low nested
% i=97 %%%%%%%nested0.25
% i=59 %%%%%%%nested0.84
% data{1,i}=csvread(d(i).name);
% B=data{1,i};

c=[1 5 11 17 38 52 76 99 118 128 134 140 142 145 147 148 ];
load network_148_data.dat
A=network_148_data;
[B,sortIdx] = sort(A(:,1),'ascend');
d=dir('*.csv')
n1=length(d)
data=cell(1,n1);
for iii=1:16
    ii=sortIdx(c(iii));
data{1,ii}=csvread(d(ii).name);
B=data{1,ii};
[nodf,qb,Nm] = cal_structure(B);


% load A1.dat
% B=A1;
[n m]=size(B);


c1=[];c2=[];
% load A1.dat
% B=A1;
[n m]=size(B);
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
k1=sum(B,1);
k2=sum(B,2);
[k3 ind]=max(k1);
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
 

  t=0.0;
%   gama1=gamma;
  gama1=0:0.1:3;
%   u=[0.1:0.2:1.4]
%   for iii=1:5
 for jjj=33:2:41; 
%      while t <  t+dt 
         for jjjj=1:length(gama1)
           
%        jjj=31;      
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
         x(:,ind)=0.2;
           x=x+ (a(jjj)*x-k(jjj)*x-diag(B3)'+mu+c1'.*x)*dt; 
        y=y+(a(jjj)*y-diag(B4)'+mu+c2'.*y)*dt;
%         x(:,1)=u(iii);

%         x=abs(x);
%         y=abs(y);

 
   end
     q=[q; nodf g T(jjj)  (sum(x))/m (sum(y))/n];
%      q1=[q1;k y];
         end

 hold on

 end
end
% for i=1:5 
%     plot(q(31*(i-1)+1:31*(i),2),q(31*(i-1)+1:31*(i),4))
%     hold on 
% end
save recovery_jul18.dat q -ascii