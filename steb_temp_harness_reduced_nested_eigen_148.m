clc
clear all
% The below code calculates the stable and unstable steady states for the reduced two dimensional mutualistic network model and the dominant eigenvalue corresponding to the stable steady state
% Please see the supplementary material to find the analytical expression for the steady states used in the code for further calculation
% param T: temperature
% param a: intrinsic growth rate
% param h: handling time
% param k1: degree of plant
% param k2: degree of pollinator
% param p: mutualisitc trade-off
% param mu: immigration term
% param gama_p and gama_A: average mutualistic strength for plant and pollinator calculated using any of the three methods: unweighted, degree weoighted or eigenvectector weightedmethod

%%%%%%%%%%%%%%%%%%%%%%%%%%%functional response function for birth rate%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%functional response function for handling time%%%%%%%%%%%%%%%%%

T0=293*ones(1,41);sigma=15;
s=2*(sigma)^2*ones(1,41);
g1=0.15;
% g1=g1/s;
T=273:1:313;
% T=T';
% for i=1:length(T)
h1 =g1*exp(((T-T0).^(2))./s);
% %  gamma1=[gamma1 gamma];
% end
% plot(T,h,'Linewidth',1.8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%functional response function for death rate%%%%%%%%%%%%%%%%%
T0=293*ones(1,41);sigma=5;
s=2*(sigma)^2*ones(1,41);
g1=0.1;
% g1=g1/s;
T=273:1:313;
% T=T';
% for i=1:length(T)
k3 =g1*exp((10000*(1./T0-1./T)));
% %  gamma1=[gamma1 gamma];
% end
% plot(T,k,'Linewidth',1.8)
d=dir('*.csv')               % this line loads all the networks present in the folder (in .csv format)
n3=length(d)
data=cell(1,n3);

A2=[];
A4=[];
A6=[];


   
 for  iiii=1:n3
     iiii
data{1,iiii}=csvread(d(iiii).name);
B=data{1,iiii};
[n1 , n2]=size(B);  %%%n1: plant, n2: Animal
for i1=1:n1
    for j1=1:n2
if B(i1,j1)>0
    B(i1,j1)=1;
else B(i1,j1)=0;
end
    end
end

[nodf,qb,Nm] = calculate_NODF(B);

k1=sum(B,2); %% degree plant
k2=sum(B,1)'; %% degree Animal

gama1=3:-0.1:0;
t=0.5; %%%trade off
 for j=21
    
for i=1:5:length(T)

A1= [];
P1= [];

P2=[];
jacob2=[];

% for i=1
%     for j=1
% alpha=a(i);
h=h1(i);
alpha=a(i);
beta=1;
k=k3(i);
g=gama1(j);
% %%%%%%%unweighted%%%%%%%%%%%%

% gama_p=sum(g*k1.^(1-t))/n1;
% gama_a=sum(g*k2.^(1-t))/n2;

% %%%%%%%%degree  weighted%%%%%%%%%%%%

gama_p=sum(g*k1.^(1-t).*k1)/sum(k1);
gama_a=sum(g*k2.^(1-t).*k2)/sum(k2);

%  
% 
% % %%%%%%%% eigenvector weighted%%%%%%%%%%%%
% [VP DP]=eigs(B*B');  
% [VA DA]=eigs(B'*B); 
% 
% gama_p=sum(g*k1.^(1-t).*VP(:,1))/sum(VP(:,1));
% gama_a=sum(g*k2.^(1-t).*VA(:,1))/sum(VA(:,1));
 % ================================================================================================= 
                                   % Calculation of stable and unstable steady states 
% ================================================================================================= 
q1=-(((beta^2)*h*gama_p)+(beta*h*gama_a*gama_p)+(beta*h^2*alpha*gama_a*gama_p));
q2=-(beta^2)-(h*alpha*beta*gama_a)+(h*alpha*beta*gama_p)+(gama_a*gama_p)+(2*h*alpha*gama_a*gama_p)+(h^2*alpha^2*gama_a*gama_p)-(k*(h*beta*gama_p)+(h*gama_a*gama_p)+(h^2*alpha*gama_a*gama_p));
q3=(alpha*beta)+(alpha*gama_a)+(h*alpha^2*gama_a)-(k*(beta+(h*alpha*gama_a)));
A3=(-q2+sqrt(q2^2-4*q1*q3))/(2*q1);    % unstable state for pollinator nodes
A4=(-q2-sqrt(q2^2-4*q1*q3))/(2*q1);    % stable state for pollinator nodes

P3=(alpha+(gama_p*A3./(1+h*gama_p*A3)));  % unstable state for plant nodes
P4=(alpha+(gama_p*A4./(1+h*gama_p*A4)));  % stable state for plant nodes


if imag(A3)==0 && A3>0
    A1=A3;
    P1=P3;
 a11=(alpha-2*P1*beta+((gama_p*A1)./(1+h*gama_p*A1)));
 a12=(-h*gama_p.^2*A1.*P1./(1+h*gama_p*A1).^2)+(gama_p*P1./(1+h*gama_p*A1));
 a13=(-h*gama_a.^2*A1.*P1./(1+h*gama_a*P1).^2)+(gama_a*A1./(1+h*gama_a*P1));
 a14=(alpha-k-2*A1*beta+((gama_a*P1)./(1+h*gama_a*P1)));
%  for iii=1:2
 jacob=[a11 a12; a13 a14];
 jacob1=eig(jacob);
%  if real(jacob1)<0
     A2=[A2; nodf T(i) g real(jacob1(1)) real(jacob1(2)) min(real(jacob1))];
     u1(i,j)=real(jacob1(1));
     u2(i,j)=real(jacob1(2));
else
     A2=[A2; nodf T(i) g NaN NaN NaN];
      u1(i,j)=NaN;
     u2(i,j)=NaN;
   v2(i,j)=real(jacob1(2));
end
 

if imag(A4)==0 && A4>0
    A1=A4;
    P1=P4;
 a11=(alpha-2*P1*beta+((gama_p*A1)./(1+h*gama_p*A1)));
 a12=(-h*gama_p.^2*A1.*P1./(1+h*gama_p*A1).^2)+(gama_p*P1./(1+h*gama_p*A1));
 a13=(-h*gama_a.^2*A1.*P1./(1+h*gama_a*P1).^2)+(gama_a*A1./(1+h*gama_a*P1));
 a14=(alpha-k-2*A1*beta+((gama_a*P1)./(1+h*gama_a*P1)));
 jacob=[a11 a12; a13 a14];
 jacob1=eig(jacob);
     A6=[A6; nodf T(i) g real(jacob1(1)) real(jacob1(2)) min(real(jacob1))];
     v1(i,j)=real(jacob1(1));
     v2(i,j)=real(jacob1(2));     %%%%%eigenvalue of the jacobian matrix
else
    
     A6=[A6; nodf T(i) g NaN NaN NaN];
      v1(i,j)=NaN;
     v2(i,j)=NaN;

end








end

    end
end
 
% save g05_steb.dat A6 -ascii
% save g05_unsteb.dat A2 -ascii
