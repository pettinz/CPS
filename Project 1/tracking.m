clear all
close all
clc

%% Build map
l_room = 10; 
l = 1;
p = 100;
n = 25; % number of sensors

% get cell coordinates using a meshgrid
tmpc = .5:10;
[xc,yc] = meshgrid(tmpc,tmpc);

% get sensor coordinates using a meshgrid
tmps = linspace(1,9,5);
[xs,ys] = meshgrid(tmps,tmps);

figure(1)
make_grid_v2(xc,yc,xs,ys,l_room)

Q = make_Q_grid(n);

%% Build A
Pt = 25;
dev_std = 0.5;
var = 0.5^2;

A = zeros(n,p);
for k = 1:p
    [xm,ym] = get_ref(k,l,p);
    d = vecnorm(([xm,ym]-[xs(:),ys(:)])')';
    
    A(:,k) = get_rss_v2(Pt,dev_std,d);
end

%% O-DIST
figure(1), hold on

lam = 1e-4;
tau = 0.7;
max_iter = 5e2;
min_eps = 1e-6;
success = 0;

ni = l_room;
dist=zeros(l_room,1);

[c_is_lower,Om,Apseudo] = reduce_coherence(A);

for it = 1:ni
    c = (it-1)*l_room+it;
    [xm,ym] = get_ref(c,l,p);  % position from measured cell
    
    p1 = plot(xm,ym,'sb','MarkerSize',10);
    
    % RSS computation
    d = vecnorm(([xm,ym]-[xs(:),ys(:)])')';
    y = get_rss_v2(Pt,dev_std,d);
    xt_0 = zeros(p,n);
    
    if c_is_lower
        yp=Om*Apseudo*y;
        Ap=Om;
    else
        yp=y;
        Ap=A;
    end
    
    [xt, iter]=distt(Ap, yp, xt_0, max_iter, Q, tau, lam, min_eps); 
    
    [~, ce] = max(abs(xt)); % estimated cell
    [xe,ye] = get_ref(ce,l,p);  % position from estimated cell
    p2 = scatter(xe, ye,'filled','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5);
    
    dist(it)= norm([mean(xe),mean(ye)] - [xm,ym]);
 
    if sum(ce == c) > n/2
        fprintf('Success\nnum iter: %d\n', iter);
        success = success+1;
    else
        fprintf('Fail\nnum iter: %d\n', iter);
    end
    fprintf('Position: %d, Estimation: %d\n', c, mode(ce));
    
    pause()
    delete(p1), delete(p2)
end

dist = cumsum(dist);
fprintf('\n\nSuccess rate: %2.0f%%\nCumulative sum: %f\n',...
    (success/ni*100), dist(end))