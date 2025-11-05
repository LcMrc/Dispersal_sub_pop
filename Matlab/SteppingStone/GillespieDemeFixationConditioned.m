% ========================================================================
% Gillespie Simulation of Deme Fixation Dynamics (Conditioned on Fixation)
% in the Stepping Stone Model
% Author: Loïc Marrec
%
% Description:
%   Runs a Gillespie simulation of mutant demes spreading across a
%   metapopulation, but only retains simulations where the mutant type
%   ultimately becomes fixed in all demes.
%
% Inputs:
%   Nit - Number of successful (mutant-fixation) runs desired
%   D   - Total number of demes
%   K   - Carrying capacity per deme
%   bW, dW - Birth and death rates of wild-type individuals
%   bM, dM - Birth and death rates of mutant individuals
%   mW, mM - Migration rates of wild-type and mutant individuals
%
% Outputs:
%   DMlist   - Final number of mutant demes per successful simulation
%   nfixlist - Number of local fixation events
%   tfixlist - Fixation times
%
% ========================================================================

function [DMlist, nfixlist, tfixlist] = GillespieDemeFixationConditioned(Nit, D, K, bW, dW, bM, dM, mW, mM)

    % Compute equilibrium deme sizes
    NW = K * (1 - dW / bW);    % Wild-type deme equilibrium size
    NM = K * (1 - dM / bM);    % Mutant deme equilibrium size

    % Compute fixation probabilities within a deme
    if bM / dM == bW / dW

        % Neutral drift
        phiW = 1 / NM;  % Fixation of wild type in mutant deme
        phiM = 1 / NW;  % Fixation of mutant in wild-type deme

    else

        % Selection
        phiW = (1 - (bM / bW)) / (1 - (bM / bW)^NM);
        phiM = (1 - (bW / bM)) / (1 - (bW / bM)^NW);

    end

    % Preallocate result arrays
    DMlist   = NaN(1, Nit);   % Final number of mutant demes
    nfixlist = NaN(1, Nit);   % Number of local fixations
    tfixlist = NaN(1, Nit);   % Fixation times

    % Main simulation loop (conditioned on fixation)
    i = 1;
    while i <= Nit
        
        % Initialization
        DM = 1;     % Initial number of mutant demes
        nfix = 0;   % Local fixation counter
        t = 0;      % Simulation time
        
        % Gillespie simulation
        while DM > 0 && DM < D
            
            % Transition rates
            Tsup = NM * phiM * mM;   % Mutant deme expansion
            Tinf = NW * phiW * mW;   % Mutant deme loss
            
            % Event selection
            r = rand;
            if r * (Tsup + Tinf) < Tsup

                DM = DM + 1;  % Increase mutant demes

            else

                DM = DM - 1;  % Decrease mutant demes

            end
            
            % Time update via sampling function
            t = t + log(1 / rand) / (Tsup + Tinf);
            
            % Local fixation counter update
            nfix = nfix + 1;

        end
        
        % Only record if mutant fixed (DM = D)
        if DM == D

            DMlist(i) = DM;         % Final number of mutant demes
            nfixlist(i) = nfix;     % Number of local fixations
            tfixlist(i) = t;        % Fixation time
            i = i + 1;

        end
        
    end

end
