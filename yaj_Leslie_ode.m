clear; clc;
rr=1;
K=10;
m=0.54;
a=1.25;
b=-1.65;
h=1;
s=0.4;

uo = [0.6307];
vo = [1.089];
col=['b'];
tfinal=[1000];
for j=1:length(uo)
Trophic=@(t,x) [rr*x(1)*(1-x(1)/K)-(m*x(1)^(2)*x(2))/(a*x(1)^2+b*x(1)+1);...
                s*x(2)*(1-x(2)/(h*x(1)))];
uo_tem=[uo(j) vo(j)]; 
tspan=[0 tfinal(j)];
options = odeset('InitialStep',1,'MaxStep',1);
[t, A]=ode45(Trophic,tspan,uo_tem,options);

u=A(:,1);
v=A(:,2);

% subplot(1,3,1)
plot(u,v,'g','LineWidth',1.); hold on;
axis([0 7 0 7]) 
set(gcf,'color','w')
xlabel('x','FontSize',13);
ylabel('y','FontSize',13);
box on;
% plot(1,1,'k.','MarkerSize',15); hold on;
% txt = {'E_{1}'}; text(1.1,1.1,txt);
% plot(4.3,5.,'k.','MarkerSize',5); hold on;
plot(4.1,4.1,'ko','MarkerSize',15); hold on;
% plot(2.2,4.4,'k.','MarkerSize',15); hold on;
% plot(1.5,4.2,'k.','MarkerSize',15); hold on;
% plot(0.5,3.1,'k.','MarkerSize',15); hold on;
% plot(0.4,2.2,'k.','MarkerSize',15); hold on;
% plot(0.6,1.3,'k.','MarkerSize',15); hold on;

% txt = {'E_{3}'}; text(4.,3.5,txt);
% txt = {'E_{2}'}; text(2.2,1.8,txt);
% txt = {'\Gamma_{2}'}; text(0.6,0.9,txt);
% txt = {'\Gamma_{1}'}; text(0.75,1.3,txt);
% txt = {'\Gamma_{3}'}; text(4.1,4.7,txt);
% 

% subplot(1,3,3)
figure(2)
plot(t,u,col(j),'LineWidth',1.); hold on;
% axis([3 4.7 3.3 4.7]) 
% set(gcf,'color','w')
% xlabel('x','FontSize',13);
% ylabel('y','FontSize',13);
% box on;
% plot(4,4,'k.','MarkerSize',15); hold on;
% % txt = {'E_{3}'}; text(4.,3.5,txt);


end
