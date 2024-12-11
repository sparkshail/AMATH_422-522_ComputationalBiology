clear all

% model parameters in 1/hours: 
% r is stem division rate, 
% l stands for lambda, TA division rate
% g stands for gamma, FD apoptosis rate 

l=1/30;
g=1/(3.5*24);
r=1/(2.5*24);

% number of stem cells
n0=18;
% average number of cells per crypt
ntot=2392.10;

% calculation of TA differentiation rate d
syms y
d=solve((1+r/(y-l)+r*y/(g*(y-l)))*n0==ntot,y);
d=double(d);

% calculation of auxiliary quantities 
a=(d-l)/g;
B=d*l/g^2;
rr=r*n0/l;
b=l/g;

h = @(x)(B.*(1 - x)).^(a./2).*besseli(-a, 2.*(B.*(1 - x)).^(1./2));
q1=@(x)a+2*(B*(1-x)).^(1/2).*besseli(a+1,2*(B*(1-x)).^(1/2))./besseli(a,2*(B*(1-x)).^(1/2));
q2=@(x)-a+2*(B*(1-x)).^(1/2).*besseli(-a+1,2*(B*(1-x)).^(1/2))./besseli(-a,2*(B*(1-x)).^(1/2));

% probability-generating function for TA cell population
F_TA=@(x)((d-l)./(d-l.*x)).^rr;
% probability-generating function for FD cell population
F_FD = @(x)(gamma(1 - a).*h(x).*(1 - (a + q2(x))./(a + q1(x)))).^rr;
% probability-generating function for total cell population
F_tot = @(x)x.^(n0).*(gamma(1 - a).*h(x).*(1 - (a +2*b*(1-x)+ q2(x))./(a+2*b*(1-x) + q1(x)))).^rr;


% contour on complex plane
C = [1+1i -1+1i -1-1i 1-1i];

% calculation of probability of TA cell population
X1=zeros(1601,1);
P_TA=zeros(1601,1);
parfor k=0:1600
    X1(k+1)=k;
    fun_TA=@(x)F_TA(x)./x.^(k+1);
    % calculation via Cauchy integration
    P_TA(k+1)=real((1/(2*pi*1i))*integral(fun_TA,(1+1i),(1+1i),'Waypoints',C));
end

% calculation of probability of FD and total cell population
X2=zeros(5001,1);
P_FD=zeros(5001,1);
P_tot=zeros(5001,1);
parfor k=0:5000
    X2(k+1)=k;
    fun_FD=@(x)F_FD(x)./x.^(k+1);
    fun_tot=@(x)F_tot(x)./x.^(k+1);
    P_FD(k+1)=real((1/(2*pi*1i))*integral(fun_FD,(1+1i),(1+1i),'Waypoints',C));
    P_tot(k+1)=real((1/(2*pi*1i))*integral(fun_tot,(1+1i),(1+1i),'Waypoints',C));
end

% calculation of explicit formula for TA cell population distribution
x_TA = 0:4000;
p_TA = nbinpdf(x_TA,rr,(d-l)/d);

% calculation of approximation for FD cell population distribution
x_FD = 0:6000;
dd=1;
for n_TA=0:4000
    p(dd,:)=poisspdf(x_FD,d*n_TA/g);
    dd=dd+1;
end 
p_FD=p_TA*p;

% data from Bravo and Axelrod 2013
TA=[880.00
360.00
360.00
760.00
612.00
880.00
936.00
792.00
552.00
400.00
552.00
782.00
690.00
546.00
840.00
520.00
782.00
660.00
432.00
192.00
420.00
400.00
720.00
440.00
912.00
384.00
580.00
456.00
440.00
320.00
414.00
756.00
1364.00
924.00
520.00
1080.00
966.00
600.00
640.00
880.00
600.00
624.00
414.00
728.00
546.00
462.00
320.00
352.00
780.00]-18;

FD=[2552.00
2000.00
1155.00
2400.00
2196.00
1672.00
1872.00
1892.00
1403.00
2075.00
1909.00
1863.00
1495.00
1785.00
2100.00
1716.00
1495.00
2112.00
3051.00
1200.00
1407.00
1320.00
1600.00
1640.00
1680.00
1344.00
2349.00
1691.00
2320.00
2560.00
2139.00
2025.00
1829.00
1617.00
1248.00
1431.00
897.00
1825.00
1280.00
1584.00
1675.00
1920.00
2415.00
2028.00
1449.00
1407.00
1400.00
1320.00
1300.00];


figure
hold on
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 2;
histogram(TA,'Normalization','pdf','LineWidth',1,'FaceColor','#FEC471','EdgeColor','#BB8041','FaceAlpha',1)
plot(X1,P_TA,'color','#a600ff','linewidth',4)
plot(x_TA(1:70:end),p_TA(1:70:end),'ko','MarkerFaceColor','#5ce1e6','MarkerSize',8,'linewidth',1.5)
legend('experimental data','Cauchy integration','explicit formula','fontsize',14)
xlabel('Number of TA cells')
ylabel('Probability mass function')
xlim([0 1600])
ylim([0 2.2*10^(-3)])

figure
hold on
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 2;
histogram(FD,'Normalization','pdf','LineWidth',1,'FaceColor','#E4CCBB','EdgeColor','#9D8678','FaceAlpha',1)
plot(X2,P_FD,'color','#a600ff','linewidth',4)
plot(x_FD(1:180:end),p_FD(1:180:end),'ko','MarkerFaceColor','#ffac50','MarkerSize',8,'linewidth',1.5)
legend('experimental data','Cauchy integration','approximation formula','fontsize',14)
xlabel('Number of FD cells')
ylabel('Probability mass function')
 xlim([0 5000])
 ylim([0 8*10^(-4)])

figure
hold on
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 2;
histogram(TA+FD+18,'Normalization','pdf','LineWidth',1,'FaceColor','#faa0a0','EdgeColor','#c88080','FaceAlpha',1)
plot(X2,P_tot,'color','#a600ff','linewidth',4)
legend('experimental data','Cauchy integration','fontsize',14)
xlabel('Total number of cell per crypt')
ylabel('Probability mass function')
xlim([0 5000])
 ylim([0 8*10^(-4)])
