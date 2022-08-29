clear; clc;
X  = 8.;          
dx = 0.1;       
nx = 1+X/dx;    
Y  = 8.;         
dy = 0.1;       
ny = 1+Y/dy;    
T  = 2;        
dt = 0.01;      
Nt = T/dt;      
r=1;s=0.4;m=0.54;K=10;a=1.25;b=-1.65;h=1;D1=0.1;D2=0.1;
a1=@(x,y ) r*x*(1-x/K) - (m*x^(2)*y)/(a*x^2+b*x+1);
a2=@(x,y ) s*y*(1-y/(h*x));
b11=@(x,y) D1^2*x^2;
b22=@(x,y) D2^2*y^2;
alpha_1=@(x,y) -(dt/dx)*( r*x*(1-x/K) - (m*x^(2)*y)/(a*x^2+b*x+1) )/4 -...
    (dt/dx^2)*( D1^2*x^2 )/4;
alpha_2=@(x,y) -(dt/dy)*( s*y*(1-y/(h*x)) )/4 - (dt/dy^2)*( D2^2*y^2 )/4;
beta_1 =@(x,y) 1 + (dt/dx^2)*( D1^2*x^2 )/2;
beta_2 =@(x,y) 1 + (dt/dy^2)*( D2^2*y^2 )/2;
eta_1  =@(x,y) 1 - (dt/dx^2)*( D1^2*x^2 )/2;
eta_2  =@(x,y) 1 - (dt/dy^2)*( D2^2*y^2 )/2;
gamma_1=@(x,y) (dt/dx)*( r*x*(1-x/K) - (m*x^(2)*y)/(a*x^2+b*x+1) )/4 - ...
    (dt/dx^2)*( D1^2*x^2 )/4;
gamma_2=@(x,y) (dt/dy)*( s*y*(1-y/(h*x)) )/4 - (dt/dy^2)*( D2^2*y^2 )/4;

%%
for i=1:nx
    x(i)=(i-1)*dx;
end

for j=1:ny
    y(j)=(j-1)*dy;
end

for l=1:Nt 
    tk(l) = (l-1)*dy/2;
end
%%
for k = 1:Nt
    for i = 1:nx
        for j = 1:ny
            U(i,j,k) = uexact(x(i),y(j),3);
        end
    end
end
%%
p = zeros(ny,ny);
for k=1:Nt
    for i = 1:nx
        for j = 1:ny
            if j==1
                U(i,j,k)    = 0;   
            elseif i==nx
                U(i,j,k)    = 0;   
            elseif i==1
                U(i,j,k)    = 0;   
            elseif j==ny
                U(i,j,k)    = 0;   
            end
        end
    end
end
%%
A1=zeros(nx,nx);
A2=zeros(ny,ny);
for k=2:2:Nt-1
    for j=2:ny-1
        for i=2:nx-1
            p(i,j) = -alpha_1(x(i-1),y(j)).*U(i-1,j,k-1) + eta_1(x(i),y(j)).*U(i,j,k-1) - gamma_1(x(i+1),y(j)).*U(i+1,j,k-1);
            A1(i,i)   = beta_1(x(i),y(j)); 
            A1(i,i-1) = alpha_1(x(i-1),y(j)); 
            A1(i,i+1) = gamma_1(x(i+1),y(j)); 
        end  
            A1(1,1:2)=[beta_1(x(1),y(j)), gamma_1(x(1+1),y(j))];
            A1(nx,nx-1:nx)=[alpha_1(x(nx-1),y(j)), beta_1(x(nx),y(j))];
            U(:,j,k) = A1\p(:,j);        
    end
    
    for  i=2:nx-1
        for j=2:ny-1
            p(i,j) = -alpha_2(x(i),y(j-1)).*U(i,j-1,k) + eta_2(x(i),y(j)).*U(i,j,k) - gamma_2(x(i),y(j+1)).*U(i,j+1,k);
            A2(j,j)   = beta_2(x(i),y(j)); 
            A2(j,j-1) = alpha_2(x(i),y(j-1)); 
            A2(j,j+1) = gamma_2(x(i),y(j+1)); 
        end   
            A2(1,1:2)=[beta_2(x(i),y(1)), gamma_2(x(i),y(1+1))];
            A2(nx,ny-1:ny)=[alpha_2(x(i),y(1)), beta_2(x(i),y(1))];
            U(i,:,(k+1)) = A2\(p(i,:)');        
    end
end
%% 
[Y,X] = meshgrid(0:dx:X);
Z = U(:,:,Nt-1);

figure(1)
s = surface(X,Y,Z,'FaceAlpha',0.5);
xlabel('x','FontSize',13);
ylabel('y','FontSize',13);
zlabel('p(x,y)','FontSize',13);
colorbar;
box on;
view(43,40); 

figure(2)
pcolor(X,Y,Z); 
shading interp;  
colorbar; colormap(jet);
xlabel('x');ylabel('y');
% title('(d) T=24');
title('(a) \sigma=0.1','FontSize',13);
