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

uniform=1; %%1->usiamo i sensori distribuiti in maniera uniforme, 0-> usiamo la grid topology

%% sensors uniformly distribuited 
if uniform==1
n=25; %%number of sensors
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
%%...
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
figure()
make_grid(xg, yg, x_sens, y_sens);
hold on

ni=50; %%number of iterations
lam = 1e-4;
tau = 0.7;
th = 0.5;
max_iter = 1e3;
min_eps = 1e-6;
success = 0;

u_A=mutual_coherence(A); %coherence before orth
Om=(orth(A'))'; %orthogonalization of A
Ap=A'*inv(A*A');
u_Om=mutual_coherence(Om); %coherence after orth

for i=1:ni
    x_measured=x_ref(ceil(p*rand()));
    y_measured=x_ref(ceil(p*rand()));
    p1=plot(x_measured, y_measured, 'sb', 'MarkerSize', 10);
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
        y=Om*Ap*y;
        A=Om;
    end
    
    xt_1 = zeros(p,1);
    
    for j=1:max_iter    
        xt = soft_tresh(xt_1+tau.*(A'*(y-A*xt_1)), lam);
        eps = norm(xt-xt_1,2)^2;
        
        if eps <= min_eps
            break
        end
        xt_1 = xt;
    end
%     
%     if sum((xt>th)) ~= 1
%         continue
%     end
    [~, p_cell] = max(abs(xt));
    p_cell=p_cell-1;
    x_estimated=fix(p_cell/10)+l_p/2;
    y_estimated=mod(p_cell, 10)+l_p/2;

    pause()
    if x_estimated==x_measured && y_estimated==y_measured
        fprintf('Success\nnum iter: %d, eps: %f\n',j,eps);
        success=success+1;
    end
    p2=plot(x_estimated, y_estimated, '.k', 'MarkerSize', 10);
    pause()
    delete(p1)
    delete(p2)
end

fprintf('Success rate: %2.2f\n',(success/ni));