% ========================================================================
% Gillespie Simulation of Deme Fixation Dynamics (Conditioned on Fixation)
% in the Island Model
% Author: Loïc Marrec
%
% Description:
%   Runs a Gillespie simulation of mutant demes spreading across a
%   metapopulation, but only retains simulations where the mutant type
%   ultimately becomes fixed in all demes.
%
% Inputs:
%   Nit - Number of stochastic replicates (with mutant fixation)
%   D   - Total number of demes
%   K   - Carrying capacity of each deme
%   bW, dW - Birth and death rates of wild-type individuals
%   bM, dM - Birth and death rates of mutant individuals
%   mW, mM - Migration rates of wild-type and mutant individuals
%
% Outputs:
%   DMlist   - Final number of mutant demes (should be D in all cases)
%   nfixlist - Number of local fixation events
%   tfixlist - Fixation times for each replicate
% ========================================================================

function [DMlist, nfixlist, tfixlist] = GillespieDemeFixationConditioned(Nit, D, K, bW, dW, bM, dM, mW, mM)

    % Compute equilibrium deme sizes
    NW = K * (1 - dW / bW);   % Wild-type equilibrium size
    NM = K * (1 - dM / bM);   % Mutant equilibrium size

    % Compute fixation probabilities within a deme
    if bM / dM == bW / dW

        % Neutral drift 
        phiW = 1 / NM;  % Fixation of wild-type in mutant deme
        phiM = 1 / NW;  % Fixation of mutant in wild-type deme

    else

        % Evolution
        phiW = (1 - (bM / bW)) / (1 - (bM / bW)^NM);
        phiM = (1 - (bW / bM)) / (1 - (bW / bM)^NW);

    end

    % Preallocate result arrays
    DMlist   = NaN(1, Nit);
    nfixlist = NaN(1, Nit);
    tfixlist = NaN(1, Nit);

    % Main simulation loop (conditioned on fixation)
    i = 1;
    while i <= Nit
        
        % Initialization
        DM   = 1;  % Start with one mutant deme
        nfix = 0;  % Number of local fixations
        t    = 0;  % Time
        
        % Gillespie simulation
        while DM > 0 && DM < D

            % Transition rates
            Tsup = DM * (D - DM) * NM * phiM * mM / (D - 1);  % Mutant spread
            Tinf = DM * (D - DM) * NW * phiW * mW / (D - 1);  % Mutant loss
            
            % Event selection
            r = rand;

            if r * (Tsup + Tinf) < Tsup

                DM = DM + 1;    % Mutant deme gained

            else

                DM = DM - 1;    % Mutant deme lost

            end
            
            % Update time and counters
            t = t + log(1 / rand) / (Tsup + Tinf);
            nfix = nfix + 1;

        end
        
        % Only record if mutant fixed (DM = D)
        if DM == D

            DMlist(i) = DM;
            nfixlist(i) = nfix;
            tfixlist(i) = t;
            i = i + 1;

        end

    end

end
