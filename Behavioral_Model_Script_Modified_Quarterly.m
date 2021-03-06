%Behavioral Model
%pstar=0;       %Central Bank's Inflation Target 
%a1=.5;         %Coefficient of expected output in output equation 
%a2=-.2;        %a is the intereste elasticty of output demand
%b1=.5;         %b1 is the coefficient of of expected inflation 
%b2=.05;        %b2 is coefficient of output in inflation equation
%c1=1.5;        %c1 is coefficient of inflation in Taylor equation
%c2=.5;         %c2 is coefficient of output in Taylor equation 
%c3=.5;         %interest smoothing parameter in Taylor equation
%beta=1;        %fixed divergence of beliefs
%delta=2;       %variable component in divergence of beliefs
%gamma=1;       %intensity of choice parameter
%sigma1=.5;     %standard deviation shocks output
%sigma2=.5;     %standard deviation shocks inflation 
%sigma3=.5;     %standard deviation shocks Taylor
%rho=.5;        % rho measures the speed of declining weights in mean square errors (memory parameter)                %errors (memory parameter) 


%Rational Model
%pstar=0;       %Central Bank's Inflation Target 
%a1=.5;         %Coefficient of expected output in output equation 
%a2=-.2;        %a is the intereste elasticty of output demand
%b1=.5;         %b1 is the coefficient of of expected inflation 
%b2=.05;        %b2 is coefficient of output in inflation equation
%c1=1.5;        %c1 is coefficient of inflation in Taylor equation
%c2=.5;         %c2 is coefficient of output in Taylor equation 
%c3=.5;         %interest smoothing parameter in Taylor equation
%sigma1=.5;     %standard deviation shocks output
%sigma2=.5;     %standard deviation shocks inflation 
%sigma3=.5;     %standard deviation shocks Taylor


%% Parameters of the model  
mm = 1;    %switching parameter gamma in Brock Hommes
pstar = 0.02; % the central bank's inflation target   
eprational=0;   % if all agents have rational forecast of inflationthis parameter is 1%
epextrapol=0;   % if all agents use inflation extrapolation this parameter is 1%
a1 = .5;     %coefficient of expected output in output equation
a2 = -0.2;   %a is the interest elasticity of output demand
b1 = .5;     %b1 is coefficient of expected inflation in inflation equation
b2 = 0.05;   %b2 is coefficient of output in inflation equation
c1 = 1.5;   %c1 is coefficient of inflation in Taylor equation
c2 = 0.5;   %c2 is coefficient of output in Taylor equation
c3 = 0.5;   %interest smoothing parameter in Taylor equation
A = [1 -b2;-a2*c1 1-a2*c2]; 
B = [b1 0;-a2 a1]; 
C = [1-b1 0;0 1-a1]; 
T = 2000; 
TI = 250; 
K = 50;        %length of period to compute divergence    
sigma1 = 0.5;  %standard deviation shocks output    
sigma2 = 0.5;  %standard deviation shocks inflation    
sigma3 = 0.5;  %standard deviation shocks Taylor    
rho=0.5;       %rho in mean squares errors    
rhoout=0.0;    %rho in shocks output       
rhoinf=0.0;    %rho in shocks inflation       
rhotayl=0.0;   %rho in shocks Taylor   
rhoBH=0.0; 
epfs=pstar;    %forecast inflation targeters       
p = zeros(T,1); 
y = zeros(T,1); 
plagt = zeros(T,1); 
ylagt = zeros(T,1); 
r = zeros(T,1); 
epf = zeros(T,1); 
epc = zeros(T,1); 
ep = zeros(T,1); 
ey = zeros(T,1); 
CRp = zeros(T,1); 
FRp = zeros(T,1); 
alfapt = zeros(T,1); 
eyfunt = zeros(T,1); 
CRy = zeros(T,1); 
FRy = zeros(T,1); 
alfayt = zeros(T,1); 
anspirits = zeros(T,1); 
epsilont = zeros(T,1); 
etat = zeros(T,1); 
ut = zeros(T,1); 

%%%%%%%%%%%%%%%
%Model
%%%%%%%%%%%%%%%%
    alfap=0.5; 
    alfay=0.5; 
    K1=K+1; 
for t=2:3;
    epsilont(t) = rhoout*epsilont(t-1) + sigma1*randn;  %shocks in output equation (demand shock)  
    etat(t)= rhoinf*etat(t-1) + sigma2*randn;    %shocks in inflation equation (supply shock)    
    ut(t) = rhotayl*ut(t-1) + sigma3*randn;       %shocks in Taylor rule (interest rate shock)  
    epsilon = epsilont(t); 
    eta = etat(t); 
    u = ut(t); 
    shocks = [eta;a2*u+epsilon]; 
    epcs=p(t-1); 
if eprational==1 ;epcs=pstar; 
end
    eps=alfap*epcs+(1-alfap)*epfs; 
if epextrapol==1;eps=p(t-1); 
end
    eychar=y(t-1); 
    eyfun=0+randn/2; 
    eyfunt(t)=eyfun; 
    eys=alfay*eychar+(1-alfay)*eyfun; 
    forecast = [eps;eys]; 
    plag=p(t-1); 
    ylag=y(t-1); 
    rlag=r(t-1); 
    lag = [plag;ylag]; 
    smooth = [0;a2*c3]; 
    D = B*forecast + C*lag + smooth*rlag + shocks; 
    X = A\D;           
    p(t)= X(1,1); 
    y(t)= X(2,1); 
    r(t)= c1*p(t)+c2*y(t)+c3*r(t-1)+u; 
if  (r(t))^2== 1; r(t)= c1*(p(t))^2+c2*y(t)+c3*r(t-1)+u;  %it says squared =1 in book code?
end
    plagt(t)=p(t-1); 
    ylagt(t)=y(t-1);       
    CRp(t) = rho*CRp(t-1) - (1-rho)*(epcs-p(t))^2; 
    FRp(t) = rho*FRp(t-1) - (1-rho)*(epfs-p(t))^2; 
    CRy(t) = rho*CRy(t-1) - (1-rho)*(eychar-y(t))^2;
    FRy(t) = rho*FRy(t-1) - (1-rho)*(eyfun-y(t))^2;
    alfap = rhoBH*alfapt(t-1)+(1-rhoBH)*exp(mm*CRp(t))/(exp(mm * CRp(t)) + exp(mm * FRp(t))); 
    alfay = rhoBH*alfayt(t-1)+(1-rhoBH)*exp(mm*CRy(t))/(exp(mm * CRy(t)) + exp(mm * FRy(t))); 
    alfapt(t) = alfap; 
    alfayt(t) = alfay; 
if eychar>0;anspirits(t)=alfay; 
end
if eychar<0;anspirits(t)=1-alfay; 
end
end %Line 123 is broken

for t=4:T;
    epsilont(t) = 0.6*rhoout*epsilont(t-1) + 0.3*rhoout*epsilont(t-2) + 0.1*rhoout*epsilont(t-3) + sigma1*randn;  %shocks in output equation (demand shock)  
    etat(t)= 0.6*rhoinf*etat(t-1) + 0.3*rhoinf*etat(t-2) + 0.1*rhoinf*etat(t-3) + sigma2*randn;    %shocks in inflation equation (supply shock)    
    ut(t) = 0.6*rhotayl*ut(t-1) + 0.3*rhotayl*ut(t-2) + 0.1*rhotayl*ut(t-3) + sigma3*randn;       %shocks in Taylor rule (interest rate shock)  
    epsilon = epsilont(t); 
    eta = etat(t); 
    u = ut(t); 
    shocks = [eta;a2*u+epsilon]; 
    epcs=0.6*p(t-1) + 0.3*p(t-1) + 0.1*p(t-3); 
if eprational==1 ;epcs=pstar; 
end
    eps=alfap*epcs+(1-alfap)*epfs; 
if epextrapol==1;eps=0.6*p(t-1) + 0.3*p(t-2) + 0.1*p(t-3); 
end
    eychar=0.6*y(t-1) + 0.3*y(t-2) + 0.1*y(t-3); 
    eyfun=0+randn/2; 
    eyfunt(t)=eyfun; 
    eys=alfay*eychar+(1-alfay)*eyfun; 
    forecast = [eps;eys]; 
    plag=0.6*p(t-1) + 0.3*p(t-2) + 0.1*p(t-3); 
    ylag=0.6*y(t-1) + 0.3*y(t-2) + 0.1*y(t-3); 
    rlag=0.6*r(t-1) + 0.3*r(t-2) + 0.1*r(t-3); 
    lag = [plag;ylag]; 
    smooth = [0;a2*c3]; 
    D = B*forecast + C*lag + smooth*rlag + shocks; 
    X = A\D;           
    p(t)= X(1,1); 
    y(t)= X(2,1); 
    r(t)= c1*p(t)+c2*y(t)+0.6*c3*r(t-1)+0.3*c3*r(t-2)+0.1*c3*r(t-3)+u; 
if  (r(t))^2== 1; r(t)= c1*(p(t))^2+c2*y(t)+0.6*c3*r(t-1)+0.3*c3*r(t-2)+0.1*c3*r(t-3)+u;  %it says squared =1 in book code?
end
    plagt(t)=0.6*p(t-1) + 0.3*p(t-2) + 0.1*p(t-3); 
    ylagt(t)=0.6*y(t-1) + 0.3*y(t-2) + 0.1*y(t-3);       
    CRp(t) = 0.6*rho*CRp(t-1) + 0.3*rho*CRp(t-2) + 0.1*rho*CRp(t-3) - (1-rho)*(epcs-p(t))^2; 
    FRp(t) = 0.6*rho*FRp(t-1) + 0.3*rho*FRp(t-2) + 0.1*rho*FRp(t-3) - (1-rho)*(epfs-p(t))^2; 
    CRy(t) = 0.6*rho*CRy(t-1) + 0.3*rho*CRy(t-2) + 0.1*rho*CRy(t-3) - (1-rho)*(eychar-y(t))^2;
    FRy(t) = 0.6*rho*FRy(t-1) + 0.3*rho*FRy(t-2) + 0.1*rho*FRy(t-3) - (1-rho)*(eyfun-y(t))^2;
    alfap = 0.6*rhoBH*alfapt(t-1)+0.3*rhoBH*alfapt(t-2)+0.1*rhoBH*alfapt(t-3)+(1-rhoBH)*exp(mm*CRp(t))/(exp(mm * CRp(t)) + exp(mm * FRp(t))); 
    alfay = 0.6*rhoBH*alfayt(t-1)+0.3*rhoBH*alfayt(t-2)+0.1*rhoBH*alfayt(t-3)+(1-rhoBH)*exp(mm*CRy(t))/(exp(mm * CRy(t)) + exp(mm * FRy(t))); 
    alfapt(t) = alfap; 
    alfayt(t) = alfay; 
if eychar>0;anspirits(t)=alfay; 
end
if eychar<0;anspirits(t)=1-alfay; 
end
end
autocory = corrcoef(y,ylagt); 
autocorp = corrcoef(p,plagt); 
coroutputanimal = corr(y,anspirits); %%  mean, median, max, min, standard deviation, kurtosis 
Kurt    = kurtosis(y); %% jarque-bera test
[jb,pvalue,jbstat] = jbtest(y,0.05); 

