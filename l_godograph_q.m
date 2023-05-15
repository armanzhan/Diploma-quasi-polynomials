function [l_g] = l_godograph_q(p,q,tau)
l_var = 0;
t = 50;
W = 0:0.01:t;
iW = 1i * W;

K = size(iW,2);

Y_im = zeros(1,K);
Y_re = zeros(1,K);

n = length(p)-1;
ang_ = 0;
delta_ang = 0;
for i = 1:K
    iw = iW(1,i);
    
    sigma_ = polyval(p,iw) + polyval(q,iw) * exp(tau*iw);
    jw1_n =(iw + 1)^n;
    gdgrph = sigma_/jw1_n;
    
    delta_ang = delta_ang + my_angle(angle(gdgrph), ang_);
    ang_ = angle(gdgrph);

    Y_re(1,i) = real(gdgrph);
    Y_im(1,i) = imag(gdgrph);

    if i==1
        l_var = norm(gdgrph);
    else
        l_var = min(l_var,norm(gdgrph));
    end
end

if abs(delta_ang) >= pi
    %display('non stable');
    l_g = -1;

%   hold on
%   subplot(2,1,2);
%   hold on;
%   plot([min(Y_re), max(Y_re)],[0,0]);
%   plot([0,0],[min(Y_im), max(Y_im)]);
%   plot(Y_re,Y_im);
   %delay = 0.001;  % seconds
   %comet3(Y_re,Y_im,zeros(1,K));
else
    %display('stable');
    l_g = l_var;

%   hold on
%   subplot(2,1,2);
%   hold on;
%   plot([min(Y_re), max(Y_re)],[0,0]);
%   plot([0,0],[min(Y_im), max(Y_im)]);
%   plot(Y_re,Y_im,color);
   %delay = 0.001;  % seconds
   %comet3(Y_re,Y_im,zeros(1,K));
end


end
function [ang] = my_angle(ang_1, ang_0)
    if ang_1-ang_0>pi
        ang = -2*pi+ang_1-ang_0;
    elseif ang_1-ang_0 < -pi
        ang = 2*pi-ang_0+ang_1;
    else
        ang = ang_1 - ang_0;
    end
end