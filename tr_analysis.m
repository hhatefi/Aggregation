function [time_samples, transient_prob]= tr_analysis(Q, pi0, interval, solver)
Q_T = Q';
[time_samples, transient_prob]=solver (@deriva, interval, pi0');

    function pdot = deriva(t, p)
        pdot = Q_T * p;
    end

end
