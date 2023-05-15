%код строит рандомизированные квазиполиномы вида p(z)+q(z) * exp(tau * z)
% условие задачи:
N_ = 1;         % количество построенных квазиполиномов
n = 5;          % степень полинома p
m = 3;          % степень полинома q
tau_max = 1;    % ограничение для запаздывания
vec_massive = 0;
u_min = 1;      % ограничения для коэффициентов снизу
u_max = 10;     % ограничения для коэффициентов сверху
%--------------------------------------------------------------------------
table_ = 0;
flag = true;

w=0:0.1:10;
for i = 1:N_
   %[p, q, tau] = quazi_1_method(n,m,u_min,u_max, tau_max);
  %pol = polysum(p,q);

   [p,q,tau] = quazi_2_method(n,m,u_min,u_max, tau_max);
   pol = p;
   %q
   
   vec_ = [p q tau];
   if vec_massive == 0
       p_massive = p;
       q_massive = q;
       vec_massive = vec_;
   else
       vec_massive =  [vec_massive; vec_];
       p_massive = [p_massive; p];
       q_massive = [q_massive; q];
   end
   % построение графиков--------------------------------------
hold on
%    y = check_dinamics(pol,p,q,tau);
%    if y
%        display(mod(y,2));
%        subplot(2,1,2);
%        hold on;
%        plot(1:n+1,vec_(1:n+1),'r-');
%        plot(n+2:n+m+2,vec_(n+2:n+m+2),'b-');
%        plot(n+m+3,vec_(n+m+3),'k.');
%    end
   %display(y);
   
%           fi = imag(quazi_val(p,q,tau,1i*w));
%           fr = real(quazi_val(p,q,tau,1i*w));
%           plot(w,fi);
%           plot(w,fr);
 
% %--------------------------------------------------------------------------
%    if length(vec_) == n+m+3   
        for i = 0:vec_(n+m+3)/5:vec_(n+m+3)
           hold on
           l_gp = l_godograph_p(pol);
           l_g = l_godograph(vec_(1:n+1),vec_(n+2:n+m+2),i,'b');
        end
        l_g = l_godograph(vec_(1:n+1),vec_(n+2:n+m+2),vec_(n+m+3),'b');
       %l_gp = l_godograph_p(pol);
% %        if l_g ~= -1 
%            subplot(2,1,1);
%            hold on;
%            plot(1:n+1,vec_(1:n+1),'r-');
%            plot(n+2:n+m+2,vec_(n+2:n+m+2),'b-');
%            plot(n+m+3,vec_(n+m+3),'k.');
%            if flag
%                 table_= [vec_ l_g];
%                 flag = false;
%            else
%                 table_ = [table_; vec_ l_g ];
%                 delta_l = l_g - l_gp;         
%                 %display([l_g l_gp delta_l]);
%            end
%--------------------------------------------------------------------------
           
%vec_(n+m+3)
%             l_g = l_godograph(vec_(1:n+1),vec_(n+2:n+m+2),vec_(n+m+3),'b');
%            l2_koeff = stab_marg_l1_koeff_red(vec_(1:n+1),vec_(n+2:n+m+2),vec_(n+m+3),0.1,0.01,0,15);
%            display(l2_koeff);
%            vec_m = vec_;
%            for i = 1:1%n+1%+m+2
%                vec_m = vec_;
%                vec_m(1,i) = vec_(1,i) + l2_koeff(1,i);%/10*9;
%                godograph(vec_m(1:n+1),vec_m(n+2:n+m+2),vec_m(n+m+3));
% 
%                vec_m = vec_;
%                vec_m(1,i) = vec_(1,i) - l2_koeff(1,i);%/10*9;
%                godograph(vec_m(1:n+1),vec_m(n+2:n+m+2),vec_m(n+m+3));
%            end
%        end
    end
   % ----------------------------------------------------------------------
   %amplitude_phaze_har(vec_(1:n+1),vec_(n+2:n+m+2),vec_(n+m+3))
   
   % расчет ограничения для запаса устойчивости.  
   %l2 = stab_marg_l2(vec_(1:n+1),vec_(n+2:n+m+2),vec_(n+m+3));
   
% end 
% [k,vol] = convhulln(p_massive);
% vol/(u_max-u_min)^(n-m)/(2*(u_max-u_min))^(m+1)
% [k,vol] = convhulln(q_massive);
% vol/(u_max-u_min)^(m+1)
%[k,vol] = convhulln(vec_massive)
%лучше критерий найквеста использовать, он нагляднее в общем и целом

%теперь надо посмотреть корреляцию между коэффициентами и запасами
% k_ = zeros(1,n+m+3);
% for i = 1:n+m+3
%     A = table_(:,i);
%     B = table_(:,n+m+4);
%     k = corrcoef(A,B);
%     k_(1,i) = k(1,2);
% end
% k_
% %table_