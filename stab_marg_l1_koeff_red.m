%p,q,tau - задаем квазиполином
%delta - величина такая, что внутри любого промежутка величины дельта будет не больее одного корня, подбирается эвристически
%err - точность выисления корня
%a - левый конец (ноль),
%b - правый конец, после достижения аргумента годографом которого годограф не пересекает оси

function [m] = stab_marg_l1_koeff_red(p,q,tau,delta,err,a,b)
    n = length(p);
    m = length(q);
    marg_var = zeros(1,n+m);
    for i=1:n+m
        marg_var(1,i) = marg_koef_i(p,q,tau,i,delta,err,a,b);
    end
%m = marg_var;
m = min(abs(marg_var),abs([p q]));
end

%эта функция неверна, надо искать минимум с помощью оптимизации
function [mm_] = find_zeros_trig(p,q,tau,i, delta, err,a,b)
    c = a:delta:b;
    var_zeros = 0;
    for k = 1:length(c)-1
        f_left = imag(quazi_val(p,q,tau,1i*c(1,k)  )/(1i*c(1,k)  )^i/exp(tau*1i*c(1,k)));
        f_rigt = imag(quazi_val(p,q,tau,1i*c(1,k+1))/(1i*c(1,k+1))^i/exp(tau*1i*c(1,k+1)));
        if f_left*f_rigt<=0
            %ищем сам корень
            x_curr = c(1,k);
            x_prev = c(1,k+1);
%---------------первая итерация--------------------------------------
            
            f_prev = imag(quazi_val(p,q,tau,1i*x_prev)/(1i*x_prev)^i/exp(tau*1i*x_prev));
            f_curr = imag(quazi_val(p,q,tau,1i*x_curr)/(1i*x_curr)^i/exp(tau*1i*x_curr));
            
            x_next = x_curr - f_curr * (x_prev - x_curr) / (f_prev - f_curr);
%---------------------------------------------------------------------
            while abs(x_curr-x_next)> err
                f_prev = imag(quazi_val(p,q,tau,1i*x_prev)/(1i*x_prev)^i/exp(tau*1i*x_prev));
                f_curr = imag(quazi_val(p,q,tau,1i*x_curr)/(1i*x_curr)^i/exp(tau*1i*x_curr));

                tmp = x_next;
                x_next = x_curr - f_curr * (x_prev - x_curr) / (f_prev - f_curr);
                x_prev = x_curr;
                x_curr = tmp;
            end
            if var_zeros ==0
                var_zeros = x_next;
            else
                var_zeros = [var_zeros x_next];
            end
        end
    end
    
    mm_ = var_zeros;
    
end

function [z] = find_zeros(p,q,tau,i,delta,err,a,b)
    c = a:delta:b;
    var_zeros = 0;
    for k = 2:length(c)-1
        f_left = imag(quazi_val(p,q,tau,1i*c(1,k)  )/(1i*c(1,k))^i);
        f_rigt = imag(quazi_val(p,q,tau,1i*c(1,k+1))/(1i*c(1,k+1))^i);
        if f_left.*f_rigt<=0
            x_curr = c(1,k);
            x_prev = c(1,k+1);
%---------------первая итерация--------------------------------------
            
            f_prev = imag(quazi_val(p,q,tau,1i*x_prev)/(1i*x_prev)^i);
            f_curr = imag(quazi_val(p,q,tau,1i*x_curr)/(1i*x_curr)^i);
            
            x_next = x_curr - f_curr * (x_prev - x_curr) / (f_prev - f_curr);
%---------------------------------------------------------------------
            while abs(x_curr-x_next)> err
                f_prev = imag(quazi_val(p,q,tau,1i*x_prev)/(1i*x_prev)^i);
                f_curr = imag(quazi_val(p,q,tau,1i*x_curr)/(1i*x_curr)^i);

                tmp = x_next;
                x_next = x_curr - f_curr * (x_prev - x_curr) / (f_prev - f_curr);
                x_prev = x_curr;
                x_curr = tmp;
            end
            if var_zeros ==0
                var_zeros = x_next;
            else
                var_zeros = [var_zeros x_next];
            end
        end
    end
    if mod(i,2)
        z = [0 var_zeros];
    else
        z = var_zeros;
    end
end
function [mm] = marg_koef_i(p,q,tau,i,delta,err,a,b)
% находим минимальное значение w при котором достигается граница устойчивости
n = length(p);
m = length(q);
    if i<=n %то есть это коэффициенты p
        k = n - i;
        w_zer = find_zeros(p,q,tau,k,delta,err,a,b);
        w_min = min(abs(real(quazi_val(p,q,tau,1i*w_zer)./(1i*w_zer).^k)));
    else % то есть это коэффициенты q
        k = m + n - i;
        w_zer = find_zeros_trig(p,q,tau,k,delta,err,a,b);
        w_min = min(abs(real(quazi_val(p,q,tau,1i*w_zer)./(1i*w_zer).^k./exp(tau*1i*w_zer))));
    end
    mm = w_min;
end
