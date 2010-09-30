% set parameters

% the name of temporary file
filename = 'tmp.csv';

% Number of reactions
nReact = 200;

% Enzyme
minEnz = 10;
maxEnz = 100;

% Substrate
minSub = 10;
maxSub = 100;

% C1
minC1 = 0.001;
maxC1 = 0.1;

% C1 * (1 - maxAbsDisC1toC2) <= C2 <= C1 * (1 + maxAbsDisC1toC2)
maxAbsDisC1toC2 = 0.2;

% Rate factor  = C1 / C3
minRF = 4;
maxRF = 10;

% min time
minTime = 50;
maxTime = 150;

% Threshold for calculation Relative Error
Threshold = 0.05;

clear('nMethods', 'statistics');

for i = 1:nReact,
        
    % open file
    fid = fopen(filename, 'w');
    % Generating constants
    c1 = rand() * (maxC1 - minC1) + minC1;
    c2 = (1 + 2 * rand() * maxAbsDisC1toC2 - maxAbsDisC1toC2) * c1;
    c3 = c1 / (rand() * (maxRF - minRF) + minRF);
    % Generating number of molecules
    x1_0 = floor ( rand() * (maxEnz - minEnz) + minEnz );
    x2_0 = floor ( rand() * (maxSub - minSub) + minSub );
    
    % Generating Interval
%    interval = [0, rand() * (maxTime - minTime) + minTime];
    interval = [0, rand() * maxTime];

    fprintf('# %d\nK1 = %f, K2 = %f, K3 = %f, X1(0) = %d, X2(0) = %d, Interval = [0, %f]\n\n', i, c1, c2, c3, x1_0, x2_0, interval(2));
    
    % Run compare
    compare(c1, c2, c3, x1_0, x2_0, interval, 0, fid, Threshold);
        
    % close file
    fclose(fid);
    
    % analyse the results
    results = csvread('tmp.csv');
    
    if ~exist('nMethods','var'),
        nMethods = size(results, 1);
    end;

    if ~exist('statistics','var'),
        statistics = zeros (nMethods, nReact, 2);
    end;
    
    for j = 1:nMethods,
        statistics(j, i, :) = results(j, [1 3]);
    end

end

for j = 1:nMethods,
    figure;
    plot((1:nReact)', statistics(j,:,1), (1:nReact)', statistics(j,:,2));
    legend('RMSE', 'MRE');
    mean_rmse = sum(statistics(j,:,1))/nReact;
    sd_rmse = sqrt(sum(statistics(j,:,1) .* statistics(j,:,1) ) / nReact - mean_rmse^2);
    fprintf('# %d:\nRMSE: Mean = %f,   SD = %f\n', j, mean_rmse, sd_rmse);
    mean_mre = sum(statistics(j,:,1))/nReact;
    sd_mre = sqrt(sum(statistics(j,:,2) .* statistics(j,:,2) ) / nReact - mean_mre^2);
    fprintf('MRE:  Mean = %f,   SD = %f\n\n\n', mean_mre, sd_mre);
end