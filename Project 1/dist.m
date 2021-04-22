function [xt] = dist(y, A, Q, lambda, tau, iter, min_eps)
%DIST Summary of this function goes here
%   Detailed explanation goes here
%     arguments
%         y (1,1) double {mustBeInteger}
%         A (:,:) double {mustBeMatrix}
%         Q (:,:) double {mustBeMatrix}
%         lambda (1,1) double
%         tau (1,1) double
%         iter (1,1) double {mustBeInteger}
%         min_eps (1,1) double
%     end
    
    [n, p] = size(A);
    
    xt_0 = zeros(p,n);
    xt = zeros(p,n);
    
    done = zeros(n,1);
    
    for k=1:iter
        for i=1:n
            if done(i)==1
                continue;
            end
            
            xt(:,i) = soft_tresh(xt_0*Q(i,:)'+tau.*(A(i,:)'*(y(i)-A(i,:)*xt_0(:,i))), lambda);
            eps = norm(xt(:,i)-xt_0(:,i),2);
            if eps <= min_eps
                done(i)=1;
            end
        end
        if(min(done)==1)
            break;
        end
        xt_0=xt;
    end

end

