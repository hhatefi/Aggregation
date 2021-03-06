function [t, mse, mre] = calc_error_between_dists(t_app, app, t_ref, ref, THRESHOLD)

t_app_len = size(t_app, 1);
t_ref_len = size(t_ref, 1);
t_len = t_app_len + t_ref_len;

nStates = size(app, 2);

t = zeros(t_len, 1);
diffrence = zeros(t_len, nStates);
relative_error = zeros(t_len, nStates);

i = 1;
j = 1;
k = 1;

while i <= t_app_len && j <= t_ref_len,
    if t_app(i) == t_ref(j),
        t(k) = t_app(i);
        diffrence(k,:) = app(i,:) - ref(j,:);
        relative_error(k,:) = (ref(j,:) >= THRESHOLD) .* abs(diffrence(k,:)) / ((ref(j,:) >= THRESHOLD) .* ref(j,:) + (ref(j,:) < THRESHOLD));

        i = i + 1;
        j = j + 1;
        k = k + 1;
    elseif t_app(i) < t_ref(j),
        % interpolate ref
        ref_intp = intp([t_ref(j - 1); t_ref(j)], [ref(j - 1,:); ref(j,:)], t_app(i));
        
        t(k) = t_app(i);
        diffrence(k,:) = app(i,:) - ref_intp;
        relative_error(k,:) = (ref_intp >= THRESHOLD) .* abs(diffrence(k,:)) / ((ref_intp >= THRESHOLD) .* ref_intp + (ref_intp < THRESHOLD));
        
        i = i + 1;
        k = k + 1;
    else
        % interpolate app
        app_intp = intp([t_app(i - 1); t_app(i)], [app(i - 1, :); app(i,:)], t_ref(j));
        
        t(k) = t_ref(j);
        diffrence(k,:) = app_intp - ref(j,:);
        relative_error(k,:) = (ref(j,:) >= THRESHOLD) .* abs(diffrence(k,:)) / ((ref(j,:) >= THRESHOLD) .* ref(j,:) + (ref(j,:) < THRESHOLD));
        
        j = j + 1;
        k = k + 1;
    end
end

while i <= t_app_len,
    % interpolate ref
    ref_intp = intp([t_ref(t_ref_len - 1); t_ref(t_ref_len)], [ref(t_ref_len - 1,:); ref(t_ref_len,:)], t_app(i));
    
    t_app(k) = t_app(i);
    diffrence(k,:) = app(i,:) - ref_intp;
    relative_error(k,:) = (ref_intp >= THRESHOLD) .* abs(diffrence(k,:)) / ((ref_intp >= THRESHOLD) .* ref_intp + (ref_intp < THRESHOLD));

    i = i + 1;
    k = k + 1;
end


while j <= t_ref_len,
    % interpolate app
    app_intp = intp([t_app(t_app_len - 1); t_app(t_app_len)], [app(t_app_len - 1,:); app(t_app_len,:)], t_ref(j));
    
    t(k) = t_ref(j);
    diffrence(k,:) = app_intp - ref(j,:);
    relative_error(k,:) = (ref(j,:) >= THRESHOLD) .* abs(diffrence(k,:)) / ((ref(j,:) >= THRESHOLD) .* ref(j,:) + (ref(j,:) < THRESHOLD));

    
    j = j + 1;
    k = k + 1;
end

t_len = k - 1;
t = t(1:t_len);
diffrence = diffrence(1:t_len,:);
relative_error = relative_error(1:t_len,:);
mse = sum(diffrence .* diffrence, 2)/nStates;
mre = max(relative_error,[],2);

    function y_i = intp(x, y, x_i)
        y_i = (y(2,:) - y(1,:)) / (x(2) - x(1)) * (x_i - x(1)) + y(1,:);
    end

end