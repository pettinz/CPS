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
A=zeros(n,p);
B=A;
for k=1:p
    k_n=k-1; % matlab fa partire i cicli da 1 mannaggia a lui
    x_ref(k)=fix(k_n/10)+l_p/2;
    y_ref(k)=mod(k_n, 10)+l_p/2;
    
    %Calcolo RSS
    A(:,k) = get_RSS([x_ref(k), y_ref(k)],[x_sens, y_sens],Pt,dev_stand);
end

figure(1)
plot(x_ref, y_ref,'.g');

%% Runtime Phase - IST
figure(3)
make_grid(xg, yg, x_sens, y_sens);
title('IST')
hold on

ni=50; %%number of iterations
lam = 1e-4;
tau = 0.7;
th = 0.5;
max_iter = 1e4;
min_eps = 1e-5;
success = 0;
iter=zeros(ni, 1);
dist=zeros(ni, 1);

u_A=mutual_coherence(A); %coherence before orth
Om=(orth(A'))'; %orthogonalization of A
Apseudo=A'*inv(A*A');
u_Om=mutual_coherence(Om); %coherence after orth
cell=randperm(p, ni);

for i=1:ni
    x_measured=x_ref(cell(i));
    y_measured=y_ref(cell(i));
    figure(3), p1=plot(x_measured, y_measured, 'sb', 'MarkerSize', 10);
    %pause(1.5)
    y=zeros(n, 1);
    
    %Calcolo RSS
    y(:) = get_RSS([x_measured, y_measured],[x_sens, y_sens],Pt,dev_stand);
    
    if u_Om<u_A
        yp=Om*Apseudo*y;
        Ap=Om;
    else
        yp=y;
        Ap=A;
    end
    
    [xt, iter(i)]=ist(max_iter, tau, Ap, yp, min_eps, lam);
    
    [~, p_cell] = max(abs(xt));
    p_cell=p_cell-1;
    x_estimated=fix(p_cell/10)+l_p/2;
    y_estimated=mod(p_cell, 10)+l_p/2;
    dist(i)= norm([x_estimated y_estimated] - [x_measured y_measured]);
 
    if cell(i)==p_cell+1
        fprintf('Success\n');
        success=success+1;
    else
        fprintf('Fail\n');
    end
    fprintf('Cella: %d, Predetta: %d\n', cell(i), p_cell+1);
    figure(3), p2=plot(x_estimated, y_estimated, '.k', 'MarkerSize', 10);
    %pause(0.5)
    pause()
    delete(p1)
    delete(p2)
end
fprintf('\n\nSuccess rate: %2.0f%%\nAverage number of iterations: %d\n',...
    (success/ni*100), round(mean(iter)));
figure()
plot([1:ni], dist, '--*')
xlabel('iteration')
ylabel('distance(m)')
title('IST')
xlim([1 50])
pause()

figure()
plot([1:ni], iter, '--*')
hold on 
plot([1 ni], [mean(iter) mean(iter)], '--r')
xlabel('iteration')
ylabel('number of iterations')
title('IST')
ylim([min(iter) max(iter)])
xlim([1 50])
legend('number of iterations', 'average number of iterations', 'Location', 'southwest')
legend('boxoff')
pause()

%% DIST
figure(6)
make_grid(xg, yg, x_sens, y_sens);
title('DIST')
hold on

lam = 1e-4;
tau = 0.7;
max_iter = 5e4;
min_eps = 1e-6;
success = 0;
iter=zeros(ni,1);

for it=1:ni
    x_measured=x_ref(cell(it));
    y_measured=y_ref(cell(it));
    figure(6), p1=plot(x_measured, y_measured, 'sb', 'MarkerSize', 10);
    %pause(1.5)
    y=zeros(n, 1);
    xt_0 = rand(p,n);
    
    %Calcolo RSS
    y(:) = get_RSS([x_measured, y_measured],[x_sens, y_sens],Pt,dev_stand);
    
    if u_Om<u_A
        yp=Om*Apseudo*y;
        Ap=Om;
    else
        yp=y;
        Ap=A;
    end
    
    [xt, iter(it)]=distt(Ap, yp, xt_0, max_iter, Q, tau, lam, min_eps); 
    %round(xt,2) %% see if consensus
    
    [~, p_cell] = max(abs(xt));
    p_cell=p_cell-1;
    x_estimated=fix(p_cell/10)+l_p/2;
    y_estimated=mod(p_cell, 10)+l_p/2;
    dist(it)= norm([mode(x_estimated) mode(y_estimated)] - [x_measured y_measured]);
 
    if cell(it)==mode(p_cell+1)
        fprintf('Success\nnum iter: %d\n', iter(it));
        success=success+1;
    else
        fprintf('Fail\nnum iter: %d\n', iter(it));
    end
    fprintf('Cella: %d, Predetta: %d\n', cell(it), mode(p_cell+1));
    figure(6), p2=plot(x_estimated, y_estimated, '.k', 'MarkerSize', 10);
    %pause(0.5)
    pause()
    delete(p1)
    delete(p2)
end
fprintf('\n\nSuccess rate: %2.0f%%\nAverage number of iterations: %d\n',...
    (success/ni*100), round(mean(iter)));
figure()
plot([1:ni], dist, '--*')
xlabel('iteration')
ylabel('distance(m)')
title('DIST')
pause()

figure()
plot([1:ni], iter, '--*')
hold on 
plot([1 ni], [mean(iter) mean(iter)], '--r')
xlabel('iteration')
ylabel('number of iterations')
title('DIST')
ylim([min(iter) max(iter)])
xlim([1 50])
legend('number of iterations', 'average number of iterations', 'Location', 'southwest')
legend('boxoff')

%% O-DIST
figure(7)
make_grid(xg, yg, x_sens, y_sens);
title('O-DIST')
hold on

path = [1 12 23 34 44 43 33 34 45 56 67 77 76 66 67 78 89 100];
i=1;
y=zeros(n, 1);
xt_0 = zeros(p,n);
xt_0(path(1),:) = 1;
max_iter = 100;
cumulative_dist = zeros(length(path),1);
coord_measured = zeros(length(path),2);
coord_estimated = zeros(length(path),2);
while(true)
    if i > length(path) %fine del percorso
        break;
    end
    x_measured=x_ref(path(i)); %il target si è mosso
    y_measured=y_ref(path(i));
    coord_measured(i,:) = [x_measured, y_measured];
    
    
    %Calcolo RSS
    y(:) = get_RSS([x_measured, y_measured],[x_sens, y_sens],Pt,dev_stand);
    
    [xt, ~]=distt(Ap, yp, xt_0, max_iter, Q, tau, lam, min_eps);
    [~, p_cell] = max(abs(xt));
    p_cell=p_cell-1;
    x_estimated=fix(p_cell/10)+l_p/2;
    y_estimated=mod(p_cell, 10)+l_p/2;
    coord_estimated(i,:) = [mean(x_estimated), mean(y_estimated)];
    
    cumulative_dist(i) = norm(coord_estimated(i,:) - coord_measured(i,:));
    
    xt_0 = xt;
    i = i+1;
end
figure(7), plot(coord_measured(:,1), coord_measured(:,2), '-o', 'MarkerSize', 10);
figure(7), plot(coord_estimated(:,1), coord_estimated(:,2), '-o', 'MarkerSize', 10, 'Color', 'r');
fprintf('\nEnd of path. Cumulative distance: %2.2f\n',sum(cumulative_dist));