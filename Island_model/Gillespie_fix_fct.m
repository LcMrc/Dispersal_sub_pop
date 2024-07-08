% This code was written by Loïc Marrec

function [XBlist, nfixlist, tfixlist] = Gillespie_fix_fct(Nit, D, K, fA, gA, fB, gB, mA, mB)

    NA = K*(1-gA/fA);       % Wild-type equilibrium size
    NB = K*(1-gB/fB);       % Mutant equilibrium size

    if fB == fA
   
        phiA = 1/NB;        % Fixation of a wild type in a mutant deme (neutral drift)
        phiB = 1/NA;        % Fixation of a mutant in a wild-type deme (neutral drift)
        
    else
        
        phiA = (1-(fB/fA))/(1-(fB/fA)^NB);      % Fixation of a wild type in a mutant deme 
        phiB = (1-(fA/fB))/(1-(fA/fB)^NA);      % Fixation of a mutant in a wild-type deme 
        
    end
    
    XBlist = NaN(1, Nit);       % List of final number of mutants 
    nfixlist = NaN(1, Nit);     % List of number of local fixations
    tfixlist = NaN(1, Nit);     % List of fixation times
    
    i = 1;
    
    while i <= Nit
    
        XB = 1;     % Initial number of mutant demes
        nfix = 0;   % Initialization of the number of local fixations
        t = 0;      % Initialization of time

        while sum(XB) ~= 0 && sum(XB) ~= D

            Tsup = XB*(D-XB)*NB*phiB*mB;    % Transition rate associated with an increase in the number of mutant demes
            Tinf = (D-XB)*XB*NA*phiA*mA;    % Transition rate associated with a decrease in the number of mutant demes

            r = rand;

            if r*(Tsup+Tinf) < Tsup

                XB = XB+1;      % Increase the number of mutant demes

            else

                XB = XB-1;      % Decrease the number of mutant demes

            end
            
            t = Sampling_time(t, Tsup+Tinf);    % Update time
            
            nfix = nfix+1;        % Update of the number of local fixations

        end
        
        XBlist(i) = XB;
        nfixlist(i) = nfix;
        tfixlist(i) = t;
        
        if XB ~= 0
            
            i = i+1;
            
        end
        
    end
    
end   