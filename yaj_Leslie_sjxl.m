clc;clear;
rr=1;
K=10;
m=0.54;
a=1.25;
% a = 1.37;
b=-1.65;
h=1;
s=0.4;
sigma=0.01;

T=1000; X=50*2^10; dt=T/X;
R=4; Dt=R*dt; L=X/R;
r=randn(1,X);
dW=sqrt(dt)*r;
W=cumsum(dW);
Xzero=0.6307;
Yzero=1.089;
Xtemps=Xzero;
Ytemps=Yzero;
Xstroges=zeros(1,L);Xstroges(1)=Xtemps;
Ystroges=zeros(1,L);Ystroges(1)=Ytemps;


for j=1:L
    Winc=sum(dW(R*(j-1)+1:R*j));
    
    Xtemps=Xtemps+(rr*Xtemps*(1-Xtemps/K) - m*Xtemps^(2)*Ytemps/(a*Xtemps^2+b*Xtemps+1))*Dt+sigma*Xtemps*Winc;
    Ytemps=Ytemps+(s*Ytemps*(1-Ytemps/(h*Xtemps)))*Dt+sigma*Ytemps*Winc;
    
    Xstroges(j)=Xtemps;
    Ystroges(j)=Ytemps;
end
%% 

figure(1)
% plot([0:Dt:(T-Dt)],Xstroges,'r-'); hold on
% axis([0,T,0,7]);
% xlabel('t','FontSize',13);
% ylabel('x','FontSize',13);

plot([0:Dt:(T-Dt)],Ystroges,'g-');
axis([0,T,0,6]);
xlabel('t','FontSize',13);
ylabel('y','FontSize',13);

figure(2)
plot(Xstroges,Ystroges,'b'); hold on;
axis([0,8,0,8]);
xlabel('x','FontSize',13);
ylabel('y','FontSize',13);
plot(4,4,'k.','MarkerSize',15); hold on;
txt = {'E_{3}'}; text(4.,3.5,txt);
plot(1,1,'k.','MarkerSize',15); 
txt = {'E_{1}'}; text(1.,0.6,txt);
hold off;

figure(3)
[p1,q1]=ksdensity(Ystroges);
plot(q1,p1,'*b');
xlabel('y','FontSize',13);
ylabel('p(x,y)','FontSize',13);
xlim([0,8]);