clc
clear all
% The below code simulates the change in abundance for change in temperature at different mutualistic strengths.
% param T: temperature
% param a: intrinsic growth rate
% param h: handling time
% param k1: degree of pollinator
% param k2: degree of plant
% param p: mutualisitc trade-off
% param mu: immigration term
% param gama1: interaction strength

mu=0.0001;p=0.5;
q=[];q1=[];q2=[];a=[];a1=[];h1=[];h=[];kk=[];k=[];

load A1.dat   % Real interaction matrix of dim 17*61
B=A1;
[n m]=size(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%functional response function for birth rate%%%%%%%%%%%%%%%%%
for i=1:m+n
T0=293*ones(1,41);sigma=5;sigma1=0.1;
s=2*(sigma)^2*ones(1,41);
g1=0.35;
T=273:1:313;
a1 =g1*exp((-(T-T0).^(2))./s);
a=[a a1'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%functional response function for handling time%%%%%%%%%%%%%%%%%

T0=293*ones(1,41);sigma=15;
s=2*(sigma)^2*ones(1,41);
g1=0.15;
T=273:1:313;
h1=g1*exp(((T-T0).^(2))./s);
h=[h h1'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%functional response function for death rate%%%%%%%%%%%%%%%%%
T0=293*ones(1,41);sigma=5;
s=2*(sigma)^2*ones(1,41);
g1=0.1;
T=273:1:313;
kk =g1*exp((10000*(1./T0-1./T)));
k=[k kk'];
end



c1=[];c2=[];
M2=[];A=[];
C=[];D=[];E=[];


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


 t=0;t_max =400; dt=0.01;
 m1=t_max/dt;
 
  t=0.0;
  q=[];
gama1=0:0.05:3;
for jjjj=[11 16 21 26 31 41]
    g=gama1(jjjj);
    q=[];
 y0=[];
 y0=[4*rand(m,1); 4*rand(n,1)];   
 y0=y0';
 y0=reshape(y0,[1,m+n]);
 p0=y0(1:m);       % initial pollinator abundance
q0=y0(m+1:m+n);    % initial plant abundance

x=p0;    
y=q0;      
for jjj=1:1:41
  B1=[];
    B2=[];
for ii=1:m
    B1(:,ii)=(B(:,ii)./(k1(ii)^p))*g;
end
for ii=1:n
    B2(ii,:)=(B(ii,:)./(k2(ii)^p))*g;
end
         
   for i=1:m1  

c1=B1'*y';
c1=c1'./(1*ones(1,m)+h(jjj,1:m).*c1');        % growth due to mutualism for pollinators   
c2=B2*x';
c2=c2'./(1*ones(1,n)+h(jjj,m+1:end).*c2');     % growth due to mutualism for pollinators   

         B3=b*(x'.*x);    % competition faced by pollinators 
         B4=b1*(y'.*y);   % competition faced by plants
      % ================================================================================================= 
                    %                      Model Equation 
      % ================================================================================================= 
        x=x+ (a(jjj,1:m).*x-k(jjj,1:m).*x-diag(B3)'+mu+c1.*x)*dt; 
        y=y+(a(jjj,m+1:end).*y-diag(B4)'+mu+c2.*y)*dt;
   end
     q=[q;T(jjj) g x y];
end
         save(['g_' num2str(jjjj-1) '.mat'],'q')
end
