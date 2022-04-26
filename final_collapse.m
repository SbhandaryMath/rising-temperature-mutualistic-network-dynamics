clc
clear all
% The below code calculates the first and final point collapse for all the real networks  available in the Web of Life: Ecological Networks Database~``http://www.web-of-life.es/".
% In the process the code deploys a function calculate_NODF to calculate nestedness for each interaction matrix as in Lever et al, 2014.
% param T: temperature
% param a: intrinsic growth rate
% param h: handling time
% param k1: degree of pollinator
% param k2: degree of plant
% param p: mutualisitc trade-off
% param mu: immigration term
% param gama1: interaction strength

mu=0.0001;p=0.5;
q=[];q1=[];q2=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%functional response function for birth rate%%%%%%%%%%%%%%%%%

T0=293*ones(1,41);sigma=5;
s=2*(sigma)^2*ones(1,41);
g1=0.35;
T=273:1:313;
a =g1*exp((-(T-T0).^(2))./s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%functional response function for handling time%%%%%%%%%%%%%%%%%

T0=293*ones(1,41);sigma=15;
s=2*(sigma)^2*ones(1,41);
g1=0.15;
T=273:1:313;
h =g1*exp(((T-T0).^(2))./s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%functional response function for death rate%%%%%%%%%%%%%%%%%
T0=293*ones(1,41);sigma=5;
s=2*(sigma)^2*ones(1,41);
g1=0.1;
T=273:1:313;
k =g1*exp((10000*(1./T0-1./T)));


c1=[];c2=[];
   d=dir('*.csv')    % this line loads all the networks present in the folder (in .csv format)
n1=length(d)
data=cell(1,n1);

M2=[];A=[];
C=[];D=[];E=[];
for i=1:n1
   
    fid = fopen(sprintf('name of network.dat'),'a' );
  
data{1,i}=csvread(d(i).name);
B=data{1,i};
[nodf,qb,Nm] = calculate_NODF(B);   



[n m]=size(B);
for i=1:n
    for j=1:m
if B(i,j)>0
    B(i,j)=1;
else B(i,j)=0;
end
    end
end
connectance=sum(sum(B))/(m*n);
b=eye(m);
b1=eye(n);
k1=sum(B,1);
k2=sum(B,2);


 t=0;t_max =400; dt=0.05;
 m1=t_max/dt;
 

  t=0.0;
  q=[];

  gama1=3:-0.1:0;

for jjj=1:41
    y0=[];
 y0=[4*rand(m,1); 4*rand(n,1)];
 y0=y0';
 y0=reshape(y0,[1,m+n]);
 p0=y0(1:m);                        % initial pollinator abundance
q0=y0(m+1:m+n);                     % initial plant abundance
x=p0;
y=q0;
    
         for jjjj=1:length(gama1)
    
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

c1=B1'*y';
c1=c1./(1+h(jjj)*c1);           % growth due to mutualism for pollinators
c2=B2*x';
c2=c2./(1+h(jjj)*c2);           % growth due to mutualism for plants 
         B3=b*(x'.*x);          % competition faced by pollinators
         B4=b1*(y'.*y);         % competition faced by plants
      % ================================================================================================= 
                    %                      Model Equation 
      % ================================================================================================= 
   
        x=x+ (a(jjj)*x-k(jjj)*x-diag(B3)'+mu+c1'.*x)*dt; 
        y=y+(a(jjj)*y-diag(B4)'+mu+c2'.*y)*dt;

     
 
   end
   q1=[x ];
     q=[q;g x y];
     if max(q1)<0.01     % Final pt collapse 
%       if any(q1<0.01 )   % First pt collapse   
     q2=[q2;nodf T(jjj) g ];
     fprintf(fid,'%f %f %f %f\n',connectance,nodf,T(jjj),g);
     break
     else 
  
continue
     end

         end

end

end
fclose(fid);
