%% project modeling and control - localization
clear all
close all
clc

%% map
l_room=10; 
l_p=1;
p=100;
n=25; %%number of sensors

xg = 0:1:10;
yg = 0:1:10;

uniform=0; %%1->usiamo i sensori distribuiti in maniera uniforme, 0-> usiamo la grid topology

%% sensors uniformly distribuited 
if uniform==1
x_sens=l_room*rand(n, 1);
y_sens=l_room*rand(n, 1);
figure(1)
make_grid(xg,yg,x_sens,y_sens);
hold on

r=4;
Q=make_Q_rand(n,r,x_sens,y_sens);
end

%% sensors grid topology ...
if uniform==0
tmp=linspace(1, 9, 5)';
x_sens=[tmp(1)*ones(5,1);tmp(2)*ones(5,1);tmp(3)*ones(5,1);tmp(4)*ones(5,1);tmp(5)*ones(5,1);];
y_sens=[tmp; tmp; tmp; tmp; tmp;];
figure(1)
make_grid(xg,yg,x_sens,y_sens);
hold on

Q=make_Q_grid(n);
end

%% check connettivity
G=graph(Q);
figure(2)
plot(G)
eigenvalue=sort(abs(eig(Q)))

%% build A
Pt=25;
dev_stand=0.5;
var=0.5^2;
x_ref=zeros(100,1);
y_ref=zeros(100,1);

for k=1:p
    k_n=k-1; % matlab fa partire i cicli da 1 mannaggia a lui
    x_ref(k)=fix(k_n/10)+l_p/2;
    y_ref(k)=mod(k_n, 10)+l_p/2;
    
    for i=1:n
        d=norm([x_ref(k), y_ref(k)]-[x_sens(i), y_sens(i)]);
        if d<=8
            Rss=Pt-40.2-20*log10(d)+dev_stand*randn();
        else
            Rss=Pt-58.5-33*log10(d)+dev_stand*randn();
        end
        A(i, k)=Rss;
    end
    
end

figure(1)
plot(x_ref, y_ref,'.g');

%% Runtime Phase
figure(3)
make_grid(xg, yg, x_sens, y_sens);
hold on

ni=50; %%number of iterations
lam = 1e-4;
tau = 0.7;
th = 0.5;
max_iter = 1e5;
min_eps = 1e-6;
success = 0;

u_A=mutual_coherence(A); %coherence before orth
Om=(orth(A'))'; %orthogonalization of A
Apseudo=A'*inv(A*A');
u_Om=mutual_coherence(Om); %coherence after orth
cell=randperm(p, ni);

for i=1:ni
    x_measured=x_ref(cell(i));
    y_measured=y_ref(cell(i));
    figure(3), p1=plot(x_measured, y_measured, 'sb', 'MarkerSize', 10);
    pause(1.5)
    y=zeros(n, 1);
    for j=1:n
        d=norm([x_measured, y_measured]-[x_sens(j), y_sens(j)]);
        if d<=8
            y(j)=Pt-40.2-20*log10(d)+dev_stand*randn();
        else
            y(j)=Pt-58.5-33*log10(d)+dev_stand*randn();
        end
    end
    
    if u_Om<u_A
        yp=Om*Apseudo*y;
        Ap=Om;
    else
        yp=y;
        Ap=A;
    end
    
    xt_1 = zeros(p,1);
    
    for j=1:max_iter    
        xt = soft_tresh(xt_1+tau.*(Ap'*(yp-Ap*xt_1)), lam);
        eps = norm(xt-xt_1,2)^2;
        
        if eps <= min_eps && sum((abs(xt)>th))==1
            break
        end
        xt_1 = xt;
    end
  
    [~, p_cell] = max(abs(xt));
    p_cell=p_cell-1;
    x_estimated=fix(p_cell/10)+l_p/2;
    y_estimated=mod(p_cell, 10)+l_p/2;

    if cell(i)==p_cell+1
        fprintf('Success\nnum iter: %d\n',j);
        success=success+1;
    else
        fprintf('Fail\nnum iter: %d, Distance: %f\n', j, norm([x_estimated y_estimated]...
        - [x_measured y_measured]));
    end
    fprintf('Cella: %d, Predetta: %d\n', cell(i), p_cell+1);
    figure(3), p2=plot(x_estimated, y_estimated, '.k', 'MarkerSize', 10);
    pause(0.5)
    delete(p1)
    delete(p2)
end
fprintf('Success rate: %2.0f%%\n',(success/ni*100));
