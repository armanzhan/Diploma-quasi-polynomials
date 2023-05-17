function [p_1 q_1 tau_1 ] = quazi_1_method(n,m,u_1,u_2, tau_max)

p_q = BRT(n, u_1, u_2);%
%p_q = lipatov_1(n,u_1, u_2);

%строим рандомный q
lambda = rand(1,m+1); 
q = lambda*(u_2-u_1)+u_1;  
p = polysum(p_q,-q)     ; %p будет разностью p_q и q

p_i = p.*1i.^matrix_of_deg(length(p));
p_i_ = p.*(-1i).^matrix_of_deg(length(p));
q_i = q.*1i.^matrix_of_deg(length(q));
q_i_ = q.*(-1i).^matrix_of_deg(length(q));

PP = conv(p_i,p_i_);
QQ = conv(q_i,q_i_);

pp_qq = real(polysum(PP,-QQ));%
PP_QQ = pp_qq(pp_qq~=0);%

%вычисляем корни
roots_ = roots(PP_QQ);
roots_ = roots_(imag(roots_) == 0 & real(roots_) >=0 ); %выбираем неотрицательные корни
roots_ = unique(roots_);
%w_i = unique([sqrt(roots_);-sqrt(roots_)]) %убрал так как все равно корни все симметричным получается, в итоге tau дублируются
w_i = unique(sqrt(roots_));

%подстановка.
phi = 1i.*log(-polyval(p_i,w_i)./polyval(q_i,w_i));

t_max = minimum_tau(phi,w_i); % граничное значение \tau

tau = -t_max;%-rand()*t_max;
%godograph(p,q,tau);

p_1 = p;
q_1 = q;
tau_1 = tau;%max(tau,-tau_max);
end

function [tau_] = minimum_tau(phi,roots_)
l = length(phi);
res = zeros(1,l);
for i = 1:l
    a = phi(i,1)/roots_(i,1);
    b = abs(2*pi/roots_(i,1));
    if a > 0
        while a - b > 0
            a = a - b;
        end
    else
        while a + b < 0
            a = a + b;
        end
    end
    res(1,i) = a;
end
%res
tau_ = real(min(res));
end
