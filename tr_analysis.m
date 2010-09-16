function [transient_prob,istate,msg]= tr_analysis(Q, pi0, interval, solver)
Q_T = Q';
[transient_prob,istate,msg]=solver (@deriva, interval, pi0');

    function pdot = deriva(t, p)
        pdot = Q_T * p;
    end

end
