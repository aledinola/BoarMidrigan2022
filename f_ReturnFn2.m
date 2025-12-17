function F = f_ReturnFn2(h,aprime,cash,tau_s,theta,gamma)

resources = cash - aprime;
% consumption
c = resources/(1+tau_s);

F = -Inf;
if c>0
    if theta==1
        utility_c = log(c);
    else
        utility_c = c^(1-theta)/(1-theta);
    end

    disutility_h = (h^(1+gamma))/(1+gamma);
    F = utility_c - disutility_h;
end

end %end function