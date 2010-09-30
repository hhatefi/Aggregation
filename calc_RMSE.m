function [rmse, t, f] = calc_RMSE(t1, f1, t2, f2)

t1_len = length(t1);
t2_len = length(t2);
t_len = t1_len + t2_len;

t = zeros(t_len, 1);
f = zeros(t_len, 1);

i = 1;
j = 1;
k = 1;

while i <= t1_len && j <= t2_len,
    if t1(i) == t2(j),
        t(k) = t1(i);
        f(k) = f1(i) - f2(j);
        i = i + 1;
        j = j + 1;
        k = k + 1;
    elseif t1(i) < t2(j),
        % interpolate f2
        f2_intp = intp([t2(j - 1), t2(j)], [f2(j - 1), f2(j)], t1(i));
        
        t(k) = t1(i);
        f(k) = f1(i) - f2_intp;
        
        i = i + 1;
        k = k + 1;
    else
        % interpolate f1
        f1_intp = intp([t1(i - 1), t1(i)], [f1(i - 1), f1(i)], t2(j));
        
        t(k) = t2(j);
        f(k) = f1_intp - f2(j);
        
        j = j + 1;
        k = k + 1;
    end
end

while i <= t1_len,
    % interpolate f2
    f2_intp = intp([t2(t2_len - 1), t2(t2_len)], [f2(t2_len - 1), f2(t2_len)], t1(i));
    
    t1(k) = t1(i);
    f(k) = f1(i) - f2_intp;

    i = i + 1;
    k = k + 1;
end


while j <= t2_len,
    % interpolate f1
    f1_intp = intp([t1(t1_len - 1), t1(t1_len)], [f1(t1_len - 1), f1(t1_len)], t2(j));
    
    t(k) = t2(j);
    f(k) = f1_intp - f2(j);

    j = j + 1;
    k = k + 1;
end

t_len = k - 1;
t = t(1:t_len);
f = f(1:t_len);
mse = sum(f.*f)/t_len;
rmse = sqrt(mse);

    function y_i = intp(x, y, x_i)
        y_i = (y(2) - y(1)) / (x(2) - x(1)) * (x_i - x(1)) + y(1);
        
    end

end