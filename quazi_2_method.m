function [p_1, q_1, tau_1] = quazi_2_method(n,m,u_min,u_max, Tau_max)
    hold on;
    tau = -rand()*Tau_max;
    p = lipatov_1(n,u_min,u_max);
    %p = BRT(n,u_min,u_max);
    
    q_ = rand(1,m)*(min(p(1,n-m+1:n+1))-u_min)+u_min;
    q__ = rand()*(p(1,n+1)-u_min)+u_min;
    q = [q_ q__]; 
    
    while ~check_godograph(p,q,tau)
        tau = -rand()*Tau_max;
        p=lipatov_1(n,u_min,u_max);
        %p = BRT(n,u_min,u_max);
        q_ = rand(1,m)*(min(p(1,n-m+1:n+1))-u_min)+u_min;
        q__ = rand()*(p(1,n+1)-u_min)+u_min;
        q = [q_ q__];
    end
    
    p_1 = p;
    q_1 = q;
    tau_1 = tau;
end
