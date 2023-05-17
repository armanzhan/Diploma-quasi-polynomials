%u_1  нижняя граница 
%u_2  верхняя граница  
% n - степень полиномов.   
%M - количество полиномов

function [polynomial] = BRT(n, u_1, u_2,M) 
m_real = 0; %   количество реально построенных полиномов
while m_real < M
    flag = false; 
    p = (u_2-u_1)*[rand() rand() rand()] + u_1;
    for i = 3:n
        if ~flag
            a = [0 p];
            p_1 = p;
            p_1(1,2:2:end) = 0;
            b = [p_1 0];
            
            a_1 = u_1;
            b_nz = b(3:2:end);
            a_   = a(3:2:end);
            b_1 = min([u_2 (u_2 - a_).* p(1,1)./b_nz]);

            if (b_1 > a_1)
                h_1 = rand()*max(0, b_1 - a_1) + a_1;
                p = a + b*h_1/p(1,1);
            else 
                flag = true;
            end
        end
    end
    if ~flag
        m_real = m_real+1;
    end
end
polynomial = p;
end
