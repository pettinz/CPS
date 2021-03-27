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
        d=norm([x_broad, y_broad]- [x(i), y(i)]);
        if d<=8
            Rss=Pt-40.2-20*log10(d)+dev_stand*randn();
        else
            Rss=Pt-58.5-33*log10(d)+dev_stand*randn();
        end
        A(i, k)=Rss;
    end
    
    pause()
    delete(pl);
end


%% runtime : centralized localization
reduce_coh = true;
if reduce_coh
    coh = mutual_coherence(A) %coherence before orth
    Om = (orth(A'))'; %orthogonalization of A
    Ap = A'*inv(A*A');
    K = Om*Ap;
    coh = mutual_coherence(Om) %coherence after orth
else
    Om = A;
end

lam = 1e-4;
tau = 0.7;
th = 0.5;
max_iter = 1e3;
min_eps = 1e-6;
num_exp = 50;
success = 0;
for k=1:num_exp
    Y = zeros(n,1);
    
    cell = randi([0 p],1,1);
    x_broad=fix(cell/10)+l_p/2;
    y_broad=mod(cell, 10)+l_p/2;
    
    for i=1:n
        d=norm([x_broad, y_broad]- [x(i), y(i)]);
        if d<=8
            Rss=Pt-40.2-20*log10(d)+dev_stand*randn();
        else
            Rss=Pt-58.5-33*log10(d)+dev_stand*randn();
        end
        Y(i)=Rss;
    end
    
    if reduce_coh
        Y = K*Y;
    end
    
    xt_1 = zeros(p,1);
    for i=1:max_iter
        
        xt = soft_tresh(xt_1+tau.*(Om'*(Y-Om*xt_1)), lam);
        eps = norm(xt-xt_1,2)^2;
        
        if eps <= min_eps
            break
        end
        xt_1 = xt;
    end
    
    if sum((xt>th)) ~= 1
        continue
    end
    [~, p_cell] = max((xt>th));
    fprintf('Ex %d:\nActual cell: %d\nPredicted cell: %d\nnum iter: %d, eps: %f\n',k,cell,p_cell-1,i,eps);
    if cell == p_cell-1
        success=success+1;
    end
    %pause()
end
fprintf('Success rate: %2.2f\n',(success/num_exp));