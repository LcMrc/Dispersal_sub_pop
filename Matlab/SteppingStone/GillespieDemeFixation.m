% ========================================================================
% Gillespie Simulation of Deme Fixation Dynamics in the Stepping Stone Model
% Author: Loïc Marrec
%
% Description:
%   Simulates the stochastic dynamics of mutant demes spreading (or going
%   extinct) across a metapopulation using a Gillespie algorithm.
%
% Inputs:
%   Nit - Number of stochastic replicates
%   D   - Total number of demes
%   K   - Carrying capacity per deme
%   fW, gW - Birth and death rates of wild-type individuals
%   fM, gM - Birth and death rates of mutant individuals
%   mW, mM - Migration rates of wild-type and mutant individuals
%
% Outputs:
%   DMlist   - Final number of mutant demes per simulation
%   nfixlist - Number of local fixation events
%   tfixlist - Fixation times
%
% ========================================================================

function [DMlist, nfixlist, tfixlist] = GillespieDemeFixation(Nit, D, K, bW, dW, bM, dM, mW, mM)

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
        phiW = (1 - (bM / bW)) / (1 - (bM / bW)^NM);    % Fixation of wild type in mutant deme
        phiM = (1 - (bW / bM)) / (1 - (bW / bM)^NW);    % Fixation of mutant in wild-type deme

    end

    % Preallocate result arrays
    DMlist   = NaN(1, Nit);   % Final number of mutant demes
    nfixlist = NaN(1, Nit);   % Number of local fixations
    tfixlist = NaN(1, Nit);   % Fixation times

    % Main simulation loop
    for i = 1 : Nit
        
        % Initialization
        DM = 1;     % Initial number of mutant demes
        nfix = 0;   % Number of local fixations
        t = 0;      % Time
        
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
            
            % Time update 
            t = t + log(1 / rand) / (Tsup + Tinf);
            
            % Local fixation counter update
            nfix = nfix + 1;

        end
        
        % Store results
        DMlist(i)   = DM;       % Number of mutant demes
        nfixlist(i) = nfix;     % Number of local fixations
        tfixlist(i) = t;        % Fixation time
        
    end

end
