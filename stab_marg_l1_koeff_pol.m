%p - задаем полином
%delta - величина такая, что внутри любого промежутка величины дельта будет не больее одного корня, подбирается эвристически
%err - точность выисления корня
%a - левый конец (ноль),
%b - правый конец, после достижения аргумента годографом которого годограф не пересекает оси

function [m] = stab_marg_l1_koeff_pol(p,delta,err,a,b)
    n = length(p);
    marg_var = zeros(1,n);
    for i=1:n
        marg_var(1,i) = marg_koef_i(p,i,delta,err,a,b);
    end
m = min(abs(marg_var), p);
end

function [z] = find_zeros(p, i, delta, err, a, b)
    c = a:delta:b;
    var_zeros = 0;
    for k = 2:length(c)-1
        f_left = imag(polyval(p,1i*c(1,k)  )/(1i*c(1,k))^i);
        f_rigt = imag(polyval(p,1i*c(1,k+1))/(1i*c(1,k+1))^i);
        if f_left*f_rigt<=0
            %ищем сам корень
            x_curr = c(1,k);
            x_prev = c(1,k+1);
%---------------первая итерация--------------------------------------
            
            f_prev = imag(polyval(p, 1i*x_prev)/(1i*x_prev)^i);
            f_curr = imag(polyval(p, 1i*x_curr)/(1i*x_curr)^i);
            
            x_next = x_curr - f_curr * (x_prev - x_curr) / (f_prev - f_curr);
%---------------------------------------------------------------------
            while abs(x_curr-x_next)> err
                f_prev = imag(polyval(p, 1i*x_prev)/(1i*x_prev)^i);
                f_curr = imag(polyval(p, 1i*x_curr)/(1i*x_curr)^i);

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
function [mm] = marg_koef_i(p,i,delta,err,a,b)
% находим минимальное значение w при котором достигается граница устойчивости
n = length(p);
k = n - i;
w_zer = find_zeros(p,k,delta,err,a,b);
w_min = min(abs(real(polyval(p, 1i*w_zer)./(1i*w_zer).^k)));
    
mm = w_min;
end
