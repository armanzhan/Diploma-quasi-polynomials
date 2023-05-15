function [a_ret] = lipatov_1(n,u_min,u_max, M)
    m_real = 0;% количество релально построенных полиномов
    r = roots([1 2 1 -1]);
    lambda_max = r(3,1);
    a_0 = rand()*(u_max-u_min) + u_min;
    a_1 = rand()*(u_max-u_min) + u_min;
    a_2 = rand()*(u_max-u_min) + u_min;
    a = [a_0 a_1 a_2 0 0 0];
    
    while m_real < M
        flag = false;
        for i = 3:n
            if ~flag
                if lambda_max > a(1,i-2)/a(1,i-1)/a(1,i)*u_min
                    a_k = a(1,i-2)/a(1,i-1)/a(1,i);
                    lambda = rand()*a_k*(u_max-u_min)+u_min*a_k;
                    while lambda > lambda_max
                        lambda = rand()*a_k*(u_max-u_min)+u_min*a_k;
                    end
                    a(1,i+1) = lambda /a(1,i-2)*a(1,i-1)*a(1,i);
                else
                    flag = true;
                    a_0 = rand()*(u_max-u_min) + u_min;
                    a_1 = rand()*(u_max-u_min) + u_min;
                    a_2 = rand()*(u_max-u_min) + u_min;
                    a = [a_0 a_1 a_2 0 0 0];
                end
            end
        end
        if ~flag
            m_real = m_real+1;
        end
    end
    a_ret = a;
end