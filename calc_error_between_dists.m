function [t, mse] = calc_error_between_dists(t_app, app, t_ref, ref)

t_app_len = size(t_app, 1);
t_ref_len = size(t_ref, 1);
t_len = t_app_len + t_ref_len;

nStates = size(app, 2);

t = zeros(t_len, 1);
f = zeros(t_len, nStates);

i = 1;
j = 1;
k = 1;

while i <= t_app_len && j <= t_ref_len,
    if t_app(i) == t_ref(j),
        t(k) = t_app(i);
        f(k,:) = app(i,:) - ref(j,:);
        i = i + 1;
        j = j + 1;
        k = k + 1;
    elseif t_app(i) < t_ref(j),
        % interpolate ref
        ref_intp = intp([t_ref(j - 1); t_ref(j)], [ref(j - 1,:); ref(j,:)], t_app(i));
        
        t(k) = t_app(i);
        f(k,:) = app(i,:) - ref_intp;
        
        i = i + 1;
        k = k + 1;
    else
        % interpolate app
        app_intp = intp([t_app(i - 1); t_app(i)], [app(i - 1, :); app(i,:)], t_ref(j));
        
        t(k) = t_ref(j);
        f(k,:) = app_intp - ref(j,:);
        
        j = j + 1;
        k = k + 1;
    end
end

while i <= t_app_len,
    % interpolate ref
    ref_intp = intp([t_ref(t_ref_len - 1); t_ref(t_ref_len)], [ref(t_ref_len - 1,:); ref(t_ref_len,:)], t_app(i));
    
    t_app(k) = t_app(i);
    f(k,:) = app(i,:) - ref_intp;

    i = i + 1;
    k = k + 1;
end


while j <= t_ref_len,
    % interpolate app
    app_intp = intp([t_app(t_app_len - 1); t_app(t_app_len)], [app(t_app_len - 1,:); app(t_app_len,:)], t_ref(j));
    
    t(k) = t_ref(j);
    f(k,:) = app_intp - ref(j,:);

    j = j + 1;
    k = k + 1;
end

t_len = k - 1;
t = t(1:t_len);
f = f(1:t_len,:);
mse = sum(f.*f, 2)/nStates;

    function y_i = intp(x, y, x_i)
        y_i = (y(2,:) - y(1,:)) / (x(2) - x(1)) * (x_i - x(1)) + y(1,:);
        
    end

end