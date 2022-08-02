


function result_array = calcium_model_simulation_v1(A, D, P, lambda, Kcdpk)

% simulates flagellar length regulation via calcium activation of CDPK
% which then inactivates kinesin
% assumptions:
% 1.  calcium concentration is a linear function of length (over a
% physiological range)
%  2.  CDPK activity depends on the calcium concentration according to a
%  Michaelis-Menten saturable binding relation with some Km for calcium
%  binding
%  3.  kinesin activity depends on the CDPK activity, again through a
%  saturable, Michaelis-type relation.  This assumption includes the
%  assumption that kinesin phosphoryation and dephosphorylation is taking
%  place at such a rapid timescale compared to length variation that it is essentially at a quasi-steady
%  state
%  4.  the number of IFT particles injected per unit time is proportional
%  to the proportion of non-phosphorylatd kinesins 
%  5.  the precursor loading to each IFT particle is proportional to the
%  total size of the available precursor pool, i.e. binding is far from
%  saturation.   this is the same assumption made in the original balance
%  point formulation in 2005.

%   parameters:
%  lambda - the KM for calcium binding to CDPK divided by the
%       proportionality constant relating calcium concentration to length.  this
%       allows the equation for CDPK activity to be expressed in terms of length
%       based on the pan lab data the phosphorylation of kinesin seems to
%       track the length pretty closely, so it seems like the system tends
%       not to get into a saturated regime, so this number should
%       relatively high, probably close to the steady state length
%   Kcdpk - the KM for kinesin phosphorylation by CDPK so that as CDPK
%       activity increases, the kinesin phosphorylation saturates at 100%
%   P  -  pool size (in  units of microns)
%   D -   disassembly rate of axonemal microtubules (microns per second,
%          default value is 0.01)
%   A  -  proportionality constant that includes the binding affinity of
%          tubulin for IFT particles and also the elongate in microns per added
%          tubulin at the tip, as well as the numbe rof potential tubulin binding
%          sites per IFT particle.

%  for interpreating the plots remember that output values are only stored
%  once every 50 time steps, which corresponds to 1 second.
%   so 1 hour is 3,600 seconds

%  example function call
%  calcium_model_simulation_v1(0.0004, 0.005, 45, 12, 1)


num_iterations = 360000;  % total simulation time is two hours 


time=[];


flagella_length =[];

length_measurements = []; % arrays to store output values for graphing




t = 0; % initial time
time_step = 0.02; % in seconds - previously used 0.02 for Elisa's paper
                  % so for regeneration to take say 120 minutes, that would
                  % require 7200 seconds which is 360,000 time steps




length_fl = 0.1; % initial flagella
kinesin_active = 1;  % initial fraction of dephosphorylated kinesin

length_fl_0 = 0.1;  % this is a floor length to prevent it getting to zero



length_measurements(1) = 0;

current_time_counter = -1;




for j = 1:(num_iterations)
    
    if mod(j,50) == 0  % store data every 50 time steps which is one second
        current_time=j
        current_time_counter = current_time_counter + 1;
        if current_time_counter > 0
            store_flag_length(current_time_counter) = length_fl;
            store_phosphorylated_kinesin(current_time_counter) = (1-kinesin_active);
        end
    end
    
    
    CDPKactive = length_fl/(lambda + length_fl);   % fraction of active kinases
    
    kinesin_active = 1 - (CDPKactive/(CDPKactive + Kcdpk)); % fracton of active kinesin
    
    length_fl = length_fl + (A*kinesin_active*(P - 2*length_fl) - D)*time_step;
    
    if length_fl < length_fl_0
        length_fl =length_fl_0;
    end
    

    
    t = t + time_step;
    
    
    %     microtubule_ultimate(j) = microtubule(end);
    %
    %     microtubule_penultimate(j) = microtubule(end-1);
    
    
end




flagella_length = length_fl




phosphorylated_kinesin_array = store_phosphorylated_kinesin;


flagellar_length_array = store_flag_length;


figure;
plot(phosphorylated_kinesin_array, 'r');  % compare to Figure 5E from Liang et al 2014

figure;
plot(flagellar_length_array, 'b');

end




