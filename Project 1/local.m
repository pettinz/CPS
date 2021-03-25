%% project modeling and control - localization
clear all
close all
clc

%% map

l_room=10; 
l_p=1;
p=100;

xg = 0:1:10;
yg = 0:1:10;
%xs = linspace(1, 9, 5);
%ys = linspace(1, 9, 5);

%% sensors uniformly distribuited 
n=25;
x=l_room*rand(n, 1);
y=l_room*rand(n, 1);
figure(1)
make_grid(xg,yg,x,y);
hold on

Q=eye(n);
eps=0.1;
r=4;

for i=1:n
    for j=n:-1:i
        if norm([x(i), y(i)]- [x(j), y(j)])<=r && i~=j
            Q(i ,j)=eps;
            Q(i, i)=Q(i, i)-eps;
            Q(j, i)=eps;
            Q(j, j)=Q(j, j)-eps;
        end
    end
end

G=graph(Q);
figure(2)
plot(G)
eigenvalue=abs(eig(Q));

%% sensors grid topology ...

%% build A

Pt=25;
dev_stand=0.5;
var=0.5^2;

for k=1:p
    k_n=k-1; % matlab fa partire i cicli da 1 mannaggia a lui
    x_broad=fix(k_n/10)+l_p/2;
    y_broad=mod(k_n, 10)+l_p/2;
    figure(1)
    pl=plot(x_broad, y_broad,'.g');
    
    for i=1:n
        d=norm([x_broad, y_broad]- [x(j), y(j)]);
        if d<=8
            Rss=Pt-40.2-20*log10(d)+var*randn();
        else
            Rss=Pt-58.5-33*log10(d)+var*randn();
        end
        A(i, k)=Rss;
    end
    
    pause()
    delete(pl);
end
