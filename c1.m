function y = c1(x,a,b)
    if strcmp(b,'upper')
        y = -log(1 + exp(x)) + a;
    elseif strcmp(b,'lower')
        y = log(1 + exp(x)) + a;
    else
        error('choose either upper or lower for variable b')
    end
end