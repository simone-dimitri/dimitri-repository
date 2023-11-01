%% MODELING: based on Huang-Kim formalism for the bunching
%% Beam size and energy vary with length. Normalized emittances are constant. Betatron functions are average values along linac sections.
%% 3-D LSC averaged (note: sign of bunching is consistent with R56<0 in a chicane). Staged LSC-microbunching.
%% Z_LSC is averaged over the trasverse distribution, whose average area is constant and = pi*r^2, in free space [M.Venturini+J.Qiang]
%% Effective beam radius is 1.747(sigx+sigy) [M.Venturini]
%% Z_LSC_rr is like Z_LSC but for a round beam in a round vacuum chamber (only important at long wavelengths) [K.Wang,Y.Li]
%% Laser Heater effect is for a beam initial Gaussian energy and transverse distribution -> hypergeometric function [Z.Huang et al.]
%% CSR and transverse Landau damping is for a 4-dipoles chicane [Z.Huang,K.Kim]
%% IBS in all straight sections (constant and varying energy) and chicanes (half+half) [Bane's approx. of B-M formalism]
%% Blog includes a correction factor of 2 w.r.t. Bane, and matches Z. Huang. Raubenheimer's cut with damping time is replaced by time of traveling.
%% Low gain terms are optional.

clc
clear variables
%close all

%% Layout:
% Drift (D0) -> Linac (L1) -> Disp.Section (BC1) -> Linac (L2+L3) -> Disp.Section (BC2) -> Linac (L4)


%% Constant
A_csr = 1.63-1i*0.9454;   % CSR free vacuum impedance
Z0 = 120*pi;              % vacuum impedance n [Ohm]
e = 1.6e-19;              % electron charge in [C]
me = 0.511;               % electron rest mass in [MeV]
re = 2.818e-15;           % electron classical radius in [m]
IA = 17045;               % Alfven current in [A]
c = 2.998e8;              % light speed in vacuum in [m/s]

%%%%%%%%%%%%%%
%% SWITCHES %%
%%%%%%%%%%%%%%

switch_lh = 0;            % 1 means LH ON, 0 means LH OFF
switch_R561_spec = 1;     % 1 means user-specified R561, 0 means R561 of a 4-dipole chicane (theta1)
switch_R562_spec = 1;     % 1 means user-specified R561, 0 means R562 of a 4-dipole chicane (theta2)
switch_csr = 1;           % 1 means CSR ON, 0 means CSR OFF
switch_lowgain = 1;       % 1 means bunching in low gain included, 0 means high gain regime
switch_bane1 = 0;         % 1 means IBS effects ON, 0 means IBS effects OFF
switch_shield = 0;        % 1 means LSC impedance with shielding, 0 means LSC impedance in vacuum

flagLH = {'noLH','LH'};
flagIBS = {'noIBS','IBS'};


%%%%%%%%%%%%%%%%%%%%%%
%% Input Parameters %%
%%%%%%%%%%%%%%%%%%%%%%

%% General
Q = 50e-12;                        % total bunch charge in [C]
I0 = 7;                            % initial peak current in [A]
E0 = 100;                            % initial mean energy in [MeV]
enx = 1e-6;                         % normalized horizontal emittance rms in [m rad]
eny = 1e-6;                         % normalized vertical emittance rms in [m rad]
sigdE = 4e-3;                       % natural initial uncorrelated energy spread rms in [MeV]

%% Laser Heater (Hypergeometric Function)
Blh = 1.2;            % ratio of laser over electron beam radius (>=1) in the LH
sigdE_lh = 6e-3;    % RMS energy spread induced by the LH in [MeV]

%% LH Area
L0 = 5;            % drift path length in [m]
Ei0 = E0;           % initial mean energy in [MeV]
Ef0 = E0;           % final mean energy in [MeV]
betax0 = 8;        % average horizontal betatron function in [m]
betay0 = 8;        % average vertical betatron function in [m]
rw0 = 10e-3;        % average inner radius of round vacuum chamber in [m]

%% LINAC1
L1 = 95;            % Linac1 path length in [m]
Ei1 = Ef0;          % initial mean energy in [MeV]
Ef1 = 580;          % final mean energy in [MeV]
betax1 = 10;        % average horizontal betatron function in [m]
betay1 = 10;        % average vertical betatron function in [m]
rw1 = 10e-3;        % average inner radius of round vacuum chamber in [m]

%% BC1
C1 = 1.25;            % linear compression factor
theta1 = 0.39;     % dipole bending angle in [rad]
Lb1 = 0.67;        % dipole length in [m]
DL1 = 1;         % drift length between outer and inner dipole of chicane, in [m]
alpha1_in = 1;      % horizontal alpha-function at entrance of chicane
alpha1_w = 0;
beta1_in = 5;      % horizontal betatron function at entrance of chicane, in [m]
beta1_w = 5;        % horizontal betatron function at waist in the second half of chicane, in [m]
R561_spec = -0.13;      % R561 user-defined in [m], replaces R561 calculated for a chicane from theta1

%% LINAC2+LINAC3
L2 = 1;            % Linac(2+3) path length in [m]
Ef2 = 580;          % final mean energy in [MeV]
betax2 = 5;        % average horizontal betatron function in [m]
betay2 = 5;        % average vertical betatron function in [m]
rw2 = 10e-3;        % average inner radius of round vacuum chamber in [m]

%% BC2
C2 = 0.8;             % linear compression factor
theta2 = 0.039;         % dipole bending angle in [rad]
Lb2 = 0.67;        % dipole length in [m]
DL2 = 2.65;         % drift length between outer and inner dipole of chicane, in [m]
alpha2_in = 1;      % horizontal alpha-function at entrance of chicane
alpha2_w = 0;
beta2_in = 5;      % horizontal betatron function at entrance of chicane, in [m]
beta2_w = 5;        % horizontal betatron function at waist in the second half of chicane, in [m]
R562_spec = -0.13;      % R562 user-defined in [m], replaces R561 calculated for a chicane from theta2

%% LINAC4
L3 = 0.01;            % Linac4 path length in [m]
Ef3 = 580;         % final mean energy in [MeV]
betax3 = 15;        % average horizontal betatron function in [m]
betay3 = 15;        % average vertical betatron function in [m]
rw3 = 10e-3;        % average inner radius of round vacuum chamber in [m]


%%%%%%%%%%%%%
%% Derived %%
%%%%%%%%%%%%%

sigd0 = sigdE/E0;                                            % initial fractional uncorrelated energy spread rms   
DeltaE_lh = 2*sigdE_lh;                                      % Half-Amplitude of the LH-induced energy modulation, in [MeV]
lb = c*Q/(5.5*I0);                                           % bunch length
halfL = DL1+Lb1;                                             % Half length of the compressor
halfL2 = DL2+Lb2;
k = @(lambda) 2*pi./lambda;                                  % uncompressed wave number in [1/m]
b0 = @(lambda) sqrt(2*e*c./(I0*lambda));                       % initial bunch factor form shot noise                  
eta1_max = theta1*(Lb1+DL1);                                 % maximum dispersion function in the chicane, in [m]
eta2_max = theta2*(Lb2+DL2);                                 % maximum dispersion function in the chicane, in [m]
% R56 of DS1-chicane in [m]
if switch_R561_spec == 0
    R561 = -2*theta1^2*((2/3)*Lb1+DL1);                          % R56 of DS1-chicane in [m]
else
    R561 = R561_spec;
end
if switch_R562_spec == 0
    R562 = -2*theta2^2*((2/3)*Lb2+DL2);                          % R56 of DS1-chicane in [m]
else
    R562 = R562_spec;
end
% R56 of DS2-chicane in [m]
if R561 ~= 0
    h1 = abs((1-1/C1)/R561);                                     % linear energy chirp at DS1
elseif R561 == 0
    h1 = 0;
end
if R562 ~= 0
    h2 = abs((1-1/C2)/R562);                                     % linear energy chirp at DS2
elseif R562 == 0
    h2 = 0;
end
sigd_BC1_cor = abs(h1*Q*c/(I0*sqrt(2*pi)));                    % correlated fractional energy RMS spread at BC1
sigd_BC2_cor = abs(h2*Q*c/(C1*I0*sqrt(2*pi)));                 % correlated fractional energy RMS spread at BC2

H11 = (theta1^2*beta1_in);                                   % horizontal H-function at the exit of first dipole of BC1-DS1, in [m]
H12 = (theta1^2*beta1_w+2*eta1_max^2/beta1_w);               % horizontal H-function at the exit of second+third dipole of BC1-DS1, in [m]
H21 = (theta2^2*beta2_in);                                   % horizontal H-function at the exit of first dipole of BC2-DS2, in [m]
H22 = (theta2^2*beta2_w+2*eta2_max^2/beta2_w);               % horizontal H-function at the exit of second+third dipole of BC2-DS2, in [m]

D = re^2*(Q/e)/(8*lb*enx);                                    % constant factor in Bane's approximation                                 

options = optimset('MaxFunEvals',10000,'TolFun',1e-15,'TolX',1e-15);


%%%%%%%%%%%%%
%% D0: IBS %%
%%%%%%%%%%%%%

% Derived
gamma0 = E0/me;  % Lorentz factor for mean energy along D0

% Q max
qmax_L0 = sqrt(L0*(Q/e)*re^2/(2*gamma0^(3/2)*enx^(3/2)*lb*sqrt(betax0)));

% IBS-induced rms relative energy spread with Bane's approx. of the B-M expression
if switch_bane1 == 1
    sigdL0_ibs = sqrt(2*L0*D*log(qmax_L0*enx/(re*2*sqrt(2)))/(gamma0^2*sqrt(enx*betax0./gamma0)));
elseif switch_bane1 == 0
    sigdL0_ibs = 0;
end

% IBS-induced absolute energy spread in [MeV] cumulated through L0
sigEL0_ibs = sigdL0_ibs*Ef0;


%%%%%%%%%%%%%
%% D0: LSC %%
%%%%%%%%%%%%%

% Average beam radius in [m] and related quantities
rb0 = 0.8735*(sqrt(enx*betax0/gamma0)+sqrt(eny*betay0/gamma0));
a0 =  @(lambda) k(lambda)*rb0/gamma0;  
I01 = @(lambda) besseli(1,a0(lambda));
K00 = @(lambda) besselk(0,a0(lambda));
K01 = @(lambda) besselk(1,a0(lambda));
aw0 =  @(lambda) k(lambda)*rw0/gamma0;
I00_w = @(lambda) besseli(0,aw0(lambda));
K00_w = @(lambda) besselk(0,aw0(lambda));

% 3-D LSC impedance averaged over transverse dimension (_rr stays for
% chamber shielding)
if switch_shield == 0
Z0_LSC = @(lambda) 1i*(Z0./(pi*k(lambda)*rb0^2)).*(1-2*I01(lambda).*K01(lambda));
elseif switch_shield ~= 0
Z0_LSC_rr = @(lambda) 1i*(Z0./(pi*k(lambda)*rb0^2))*(1-2*(I01(lambda)/I00_w(lambda))*(I00_w(lambda)*K01(lambda)+I01(lambda)*K00_w(lambda)));
Z0_LSC = @(lambda) Z0_LSC_rr(lambda);
end

% LSC-induced energy modulation amplitude in me unit
Dgamma0 = @(lambda) -(4*pi/Z0)*(I0/IA)*b0(lambda).*Z0_LSC(lambda)*L0;
Dg0_rel = @(lambda) Dgamma0(lambda)/gamma0;


%%%%%%%%%%%%%%%%%
%% LINAC1: IBS %%
%%%%%%%%%%%%%%%%%

% Derived
G1 = (Ef1-Ei1)/L1;                       % L1 average accelerating gradient in [MeV/m]
gamma1 = @(s1) (Ei1+G1*s1)/me;           % Lorentz factor for mean energy along L1
gammaf1 = Ef1/me;                        % Lorentz factor for mean energy at the END of L1
gammam1 = (gammaf1 - Ei1/me)/2;          % Mean Lorentz factor of L1

sigd_L1 = sqrt(sigd0^2 + sigdL0_ibs^2);  % Uncorrelated energy spread at the entrance of L1

% Max angle
qmax_L1 = sqrt(L1*(Q/e)*re^2/(2*enx^(3/2)*lb*sqrt(betax1)));

% Coulomb Logarithm
ln1 = log(qmax_L1*enx/(re*2*sqrt(2)));

% IBS-induced rms absolute energy spread in [MeV] with Bane's approximation
% of B-M expression, cumulated through L1
if switch_bane1 == 1
        % Di Mitri-Perosa closed form for energy-dependent Clog 
        sigdL1_DMP1 = sqrt((2*me*D/(3*gammaf1^2*G1*sqrt(enx*betax1)))*(gammaf1^(3/2)-gamma0^(3/2)+2*ln1*(gammaf1^(3/2)-gamma0^(3/2))+gamma0^(3/2)*log(gamma0^(3/4))-gammaf1^(3/2)*log(gammaf1^(3/4))));
        % Di Mitri-Perosa closed form for energy-INdependent Clog 
        sigdL1_DMP2 = sqrt((4/3)*(re^2*(Q/e)*ln1*me/(8*enx^(3/2)*sqrt(betax1)*lb*G1))*((gammaf1^(3/2)-gamma0^(3/2))/gammaf1^2));

     sigdL1_ibs = sigdL1_DMP1;
     sigEL1_ibs = sigdL1_ibs*Ef1;

elseif switch_bane1 == 0
    sigdL1_ibs = 0;
    sigEL1_ibs = 0;
end


%%%%%%%%%%%%%%%%%
%% LINAC1: LSC %%
%%%%%%%%%%%%%%%%%

% Average beam radius in [m] and related quantities
rb1 = @(s1) 0.8735*(sqrt(enx*betax1./gamma1(s1))+sqrt(eny*betay1./gamma1(s1)));
a1 =  @(s1,lambda) k(lambda).*rb1(s1)./gamma1(s1);                                                 
I11 = @(s1,lambda) besseli(1,a1(s1,lambda));
K10 = @(s1,lambda) besselk(0,a1(s1,lambda));
K11 = @(s1,lambda) besselk(1,a1(s1,lambda));
aw1 =  @(s1,lambda) k(lambda)*rw1./gamma1(s1);
I10_w = @(s1,lambda) besseli(0,aw1(s1,lambda));
K10_w = @(s1,lambda) besselk(0,aw1(s1,lambda));

% 3-D LSC impedance averaged over transverse dimension (_rr stays for
% chamber shielding)
if switch_shield == 0
Z1_LSC = @(s1,lambda) 1i*(Z0./(pi*k(lambda).*rb1(s1).^2)).*(1-2*I11(s1,lambda).*K11(s1,lambda));
elseif switch_shield ~= 0
Z1_LSC_rr = @(s1,lambda) 1i*(Z0./(pi*k(lambda).*rb1(s1).^2)).*(1-2*(I11(s1,lambda)./I10_w(s1,lambda)).*(I10_w(s1,lambda).*K11(s1,lambda)+I11(s1,lambda).*K10_w(s1,lambda)));
Z1_LSC = @(s1,lambda) Z1_LSC_rr(s1,lambda);
end
% Approximated LSC impedance
%Z1_LSC_1 = @(s1,lambda) 1i*(Z0./(pi*k(lambda)*rb1(s1).^2))*(1-a1(s1,lambda)*K11(s1,lambda));
%Z1_LSC_2 = @(s1,lambda) 1i*(Z0./(2*pi*gamma1(s1).*rb1(s1)))*sqrt(1+a1(s1,lambda).^2*(K10(s1,lambda).^2-K11(s1,lambda).^2));

% LSC-induced energy modulation amplitude in me unit
Z1_LSC_real = @(s1,lambda) real(Z1_LSC(s1,lambda));
Z1_LSC_imag = @(s1,lambda) imag(Z1_LSC(s1,lambda));

Dgamma1_real = @(lambda) -(4*pi/Z0)*(I0/IA)*b0(lambda).*integral(@(s1)Z1_LSC_real(s1,lambda),0,L1,'ArrayValued',true);
Dgamma1_imag = @(lambda) -(4*pi/Z0)*(I0/IA)*b0(lambda).*integral(@(s1)Z1_LSC_imag(s1,lambda),0,L1,'ArrayValued',true);
Dgamma1 = @(lambda) Dgamma1_real(lambda)+1i*Dgamma1_imag(lambda);
Dg1_rel = @(lambda) -(4*pi/Z0)*(I0/IA)*b0(lambda).*integral(@(s1)Z1_LSC(s1,lambda)./gamma1(s1),0,L1,'ArrayValued',true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BC1: Compression factor %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigd01 = sqrt((sigd0*Ef0)^2+sigEL0_ibs^2+sigEL1_ibs^2)/Ef1;                                         % fractional energy spread rms, normalized at the DS1 energy
I1 = I0*C1;                                                                                         % peak current compressed by C1 in DS1 in [A]
ex1 = enx/gammaf1;                                                                                  % geometric horizontal emittance at chicane in [m rad]
lambda_co_long_1 = 2*pi*abs(R561)*C1.*sqrt(sigdE^2+switch_lh*sigdE_lh^2+sigEL0_ibs^2+sigEL1_ibs^2)/Ef1;   % Energy Landau damping cutoff wavelength in [m]
lambda_co_tran_1 = 2*pi*sqrt(ex1*H11);                                                              % Energy Landau damping cutoff wavelength in [m]


%%%%%%%%%%%%%%
%% BC1: IBS %%
%%%%%%%%%%%%%%

% Derived
sigd_BC1_i = sqrt(sigd_L1^2*Ef0^2 + sigEL1_ibs^2)/Ef1;      % Uncorrelated energy spread at the entrance of BC1
qmax_BC1_i = sqrt(halfL*(Q/e)*re^2/(2*gammaf1^(3/2)*enx^(3/2)*lb*sqrt(beta1_in)));
qmax_BC1_w = sqrt(C1*halfL*(Q/e)*re^2/(2*gammaf1^(3/2)*enx^(3/2)*lb*sqrt(beta1_w)));

% optics functions in the first half (_i) and second half (_w) of BC1
etap1 = theta1;
eta1  = theta1*(DL1+Lb1);
gamma1_in = (1+alpha1_in^2)/beta1_in;
gamma1_w  = (1+alpha1_w^2)/beta1_w;

if switch_bane1 == 1
    if theta1 ~= 0
    % IBS-induced fractional rms energy spread with Bane's formula for the first half of BC1
        dispH1_i = beta1_in*etap1^2+2*alpha1_in*eta1*etap1+gamma1_in*eta1^2;
        
        aBC11 = vpa(2*D/(gammaf1^(3/2)*sqrt(enx*beta1_in)));
        bBC11 = vpa(qmax_BC1_i*enx/(2*sqrt(2)*re));
        hBC11 = vpa(dispH1_i*gammaf1/enx);
        solBC11 = @(y) vpa(abs( -ei(3*log(sqrt(hBC11*sigd_BC1_i^2 + 1)/bBC11)) + ei(3*log(sqrt(hBC11*y^2 + 1)/bBC11)) + vpa(halfL)*aBC11*hBC11/(2*bBC11^3)));
        yBC11 = vpa(fminsearch(solBC11,1e-6,options));    
        sigdBC1_i_ibs = double(sqrt(yBC11^2 - sigd_BC1_i^2));
    
        sigd_BC1_w = sqrt(sigd_BC1_i^2 + sigdBC1_i_ibs^2);     % Uncorrelated energy spread at the waist of BC1
    
    % IBS-induced fractional rms energy spread with Bane's formula for the second half of BC1
        dispH1_w = beta1_w*etap1^2+2*alpha1_w*eta1*etap1+gamma1_w*eta1^2;

        aBC12 = vpa(2*C1*D/(gammaf1^(3/2)*sqrt(enx*beta1_w)));
        bBC12 = vpa(qmax_BC1_w*enx/(2*sqrt(2)*re));
        hBC12 = vpa(dispH1_w*gammaf1/enx);
        solBC12 = @(y) vpa(abs( -ei(3*log(sqrt(hBC12*sigd_BC1_w^2 + 1)/bBC12)) + ei(3*log(sqrt(hBC12*y^2 + 1)/bBC12)) + vpa(halfL)*aBC12*hBC12/(2*bBC12^3)));
        yBC12 = vpa(fminsearch(solBC12,1e-5,options));    
        sigdBC1_w_ibs = double(sqrt(yBC12^2 - sigd_BC1_w^2));
    
    elseif theta1 == 0       
        sigdBC1_i_ibs = sqrt(2*halfL*D*log(qmax_BC1_i*enx/(re*2*sqrt(2)))./(gammaf1.^2*sqrt(enx*beta1_in./gammaf1)));
        sigd_BC1_w = sqrt(sigd_BC1_i^2 + sigdBC1_i_ibs^2);
        sigdBC1_w_ibs = sqrt(2*halfL*C1*D*log(qmax_BC1_w*enx/(re*2*sqrt(2)))./(gammaf1.^2*sqrt(enx*beta1_w./gammaf1)));
    end
elseif switch_bane1 == 0
    sigdBC1_i_ibs = 0;
    sigdBC1_w_ibs = 0;
    sigd_BC1_w = sqrt(sigd_BC1_i^2 + sigdBC1_i_ibs^2);
end

% IBS-induced energy spread in BC1, in [MeV] 
sigEBC1_i = sigdBC1_i_ibs*Ef1;
sigEBC1_w = sigdBC1_w_ibs*Ef1;


%%%%%%%%%%%%%%
%% BC1: LSC %%
%%%%%%%%%%%%%%

rbBC1_i = 0.8735*(sqrt(enx*beta1_in/gammaf1)+sqrt(eny*beta1_in/gammaf1));
rbBC1_w = 0.8735*(sqrt(enx*beta1_w/gammaf1)+sqrt(eny*beta1_w/gammaf1));
aBC1_i =  @(lambda) k(lambda)*rbBC1_i/gammaf1;
aBC1_w =  @(lambda) k(lambda)*rbBC1_w/gammaf1;
I1_BC1_i = @(lambda) besseli(1,aBC1_i(lambda));
I1_BC1_w = @(lambda) besseli(1,aBC1_w(lambda));
K0_BC1_i = @(lambda) besselk(0,aBC1_i(lambda));
K0_BC1_w = @(lambda) besselk(0,aBC1_w(lambda));
K1_BC1_i = @(lambda) besselk(1,aBC1_i(lambda));
K1_BC1_w = @(lambda) besselk(1,aBC1_w(lambda));
awBC1 =  @(lambda) k(lambda)*rw0/gammaf1;
I0_w_BC1 = @(lambda) besseli(0,awBC1(lambda));
K0_w_BC1 = @(lambda) besselk(0,awBC1(lambda));

% 3-D LSC impedance averaged over transverse dimension (_rr stays for
% chamber shielding)
if switch_shield == 0
    Z_LSC_BC1_i = @(lambda) 1i*(Z0./(pi*k(lambda)*rbBC1_i^2)).*(1-2*I1_BC1_i(lambda).*K1_BC1_i(lambda));
    Z_LSC_BC1_w = @(lambda) 1i*(Z0./(pi*k(lambda)*rbBC1_w^2)).*(1-2*I1_BC1_w(lambda).*K1_BC1_w(lambda));
elseif switch_shield ~= 0
    Z_LSC_rr_BC1_i = @(lambda) 1i*(Z0./(pi*k(lambda)*rbBC1_i^2))*(1-2*(I1_BC1_i(lambda)/I0_w_BC1(lambda))*(I0_w_BC1(lambda)*K1_BC1_i(lambda)+I1_BC1_i(lambda)*K0_w_BC1(lambda)));
    Z_LSC_rr_BC1_w = @(lambda) 1i*(Z0./(pi*k(lambda)*rbBC1_w^2))*(1-2*(I1_BC1_w(lambda)/I0_w_BC1(lambda))*(I0_w_BC1(lambda)*K1_BC1_w(lambda)+I1_BC1_w(lambda)*K0_w_BC1(lambda)));
    Z_LSC_BC1_i = @(lambda) Z_LSC_rr_BC1_i(lambda);
    Z_LSC_BC1_w = @(lambda) Z_LSC_rr_BC1_w(lambda);
end

% LSC-induced energy modulation amplitude in me unit
DgammaBC1_i = @(lambda) -(4*pi/Z0)*(I0/IA)*b0(lambda).*Z_LSC_BC1_i(lambda)*halfL;
DgammaBC1_w = @(lambda) -(4*pi/Z0)*(I0/IA)*b0(lambda).*Z_LSC_BC1_w(lambda)*halfL;
DgammaBC1 = @(lambda) DgammaBC1_i(lambda) + DgammaBC1_w(lambda);
DgBC1_rel = @(lambda) DgammaBC1(lambda)/gammaf1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BC1: Laser Heater effect %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the hypergeometric function in its integral form. The integral
% runs over the ratio of laser/e-beam radius: assume here from 1 to 10) for
% the definite integral.
A1 = @(lambda) abs((2*pi*C1./lambda)*R561*DeltaE_lh/Ef1);
J01_LH = @(r,lambda) besselj(0,A1(lambda).*exp(-r.^2/(4*Blh^2)));
S1_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*J01_LH(r,lambda),0,Inf,'ArrayValued',true);

if sigdE_lh ~= 0 && switch_lh == 1
    S1_LH = @(lambda) S1_LH(lambda);
elseif sigdE_lh == 0 || switch_lh == 0
    S1_LH = @(lambda) 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CSR GAIN in BC1 vs. uncompressed wavelength %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if switch_csr ~= 0 && theta1 ~= 0
%sigT11 = @(lambda) k(lambda)*sqrt(ex1*H11);
sigT11 = @(lambda) k(lambda)*theta1*sqrt(ex1*beta1_in);
erf_T11 = @(lambda) erf(sigT11(lambda));
%sigT12 = @(lambda) k(lambda)*sqrt(ex1*H12);
%erf_T12 = @(lambda) erf(sigT12(lambda));
I1n_csr = @(lambda) (I1*R561*k(lambda).^(4/3)*theta1^(2/3)*Lb1^(1/3)/(IA*gammaf1));
phi1 = 2*DL1/beta1_in;

% Ht1 (Ht2) takes care of transverse Landau damping through Dipole-1 (Dipole-2+Dipole3) of BC1-DS1
Ht11 = @(lambda,t) exp(-sigT11(lambda).^2*((1-2.*t+alpha1_in*phi1.*t).^2+(phi1.*t).^2)/(1+h1*R561*t).^2-(C1*R561*k(lambda).*sigd01).^2*(t.^2+(C1*(1-t)).^2)/(2*(1+h1*R561*t).^2));
%Ht12 = @(lambda,t) exp(-sigT12(lambda).^2*((1-2.*t+alpha1_in*phi1.*t).^2+(phi1.*t).^2)/(1+h1*R561*t).^2-(C1*R561*k(lambda)*sigd01).^2*(t.^2+(C1*(1-t)).^2)/(2*(1+h1*R561*t).^2));

F01 = @(lambda) (exp(-sigT11(lambda).^2)+sigT11(lambda).*sqrt(pi).*erf_T11(lambda)-1)/(2*sigT11(lambda).^2);
F11 = @(lambda) 2*integral(@(t)(1-t).*Ht11(lambda,t)./(1+h1*R561*t).^(4/3),0,1,'ArrayValued',true);
F21 = @(lambda) 2*integral(@(t)(1-t).*t.*Ht11(lambda,t)./(C1*(1+h1*R561*t).^(7/3)),0,1,'ArrayValued',true);

% this is the CSR gain driven by initial bunching, at first order in beam current
b11_csr = @(lambda) -b0(lambda).*A_csr.*I1n_csr(lambda).*(sqrt(pi).*erf_T11(lambda).*exp(-0.5*(C1*R561*k(lambda)*sigd01).^2).*S1_LH(lambda)./(2*sigT11(lambda))+F11(lambda));

% this is the CSR gain driven by initial bunching, at second order in beam current
b12_csr = @(lambda)  b0(lambda).*A_csr^2.*I1n_csr(lambda).^2.*F01(lambda).*F21(lambda);

% this is the CSR gain driven by incoming energy modulation (LSC) in D0+L1, at first order in beam current
b13_csr = @(lambda) -1i*C1*R561*k(lambda).*((Dg1_rel(lambda)*Ef1 + Dg0_rel(lambda)*Ef0)/Ef1).*A_csr.*I1n_csr(lambda).*F21(lambda);

elseif switch_csr == 0 || theta1 == 0
    b11_csr = @(lambda) 0;
    b12_csr = @(lambda) 0;
    b13_csr = @(lambda) 0;
end

Gain1_csr = @(lambda) abs((b11_csr(lambda)+b12_csr(lambda)+b13_csr(lambda))/b0(lambda));

% Delta = Lb1*(Q/e)*re^2/(2*gammaf1^(3/2)*enx^(3/2)*lb*sqrt(beta1_in));
% sigma = 4*sqrt(ex1*beta1_in)*theta1/(R561^2*C1^2);
% klim = sigma/Delta
% 2*pi/klim


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LSC + CSR + IBS GAIN after D0+L1+BC1 vs. uncompressed wavelength %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is the initial bunching with longitudinal damping in BC1, no
% collective effects
if switch_lowgain ~= 0
    b10 = @(lambda) b0(lambda).*exp(-0.5.*C1^2*k(lambda).^2*R561^2*sigd01^2).*S1_LH(lambda);
else 
    b10 = @(lambda) 0;
end

% this is the bunching induced by the energy modulation cumulated in L0+L1, passing
% through DS1, and with damping in DS1
b11 = @(lambda) -1i*C1*k(lambda)*R561.*((Dg1_rel(lambda)*Ef1 + Dg0_rel(lambda)*Ef0)/Ef1).*...
    exp(-0.5*C1^2*k(lambda).^2*R561^2*sigd01^2).*S1_LH(lambda);

b1 = @(lambda) b10(lambda)+b11(lambda)+b11_csr(lambda)+b12_csr(lambda)+b13_csr(lambda);

Gain1 = @(lambda) abs(b1(lambda)./b0(lambda));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Uncorrelated Energy Spread at BC1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncorrelated energy spread RMS induced by LSC through L0+L1 in [keV]
sigdE1_LSC = 1e3*me*sqrt((2*e*c/I0)*integral(@(lambda)abs((Dgamma1(lambda)+Dgamma0(lambda)).^2./(b0(lambda).^2.*lambda.^2)),0.1e-6,200e-6,'ArrayValued',true));

% Total uncorrelated energy spread RMS in [keV] at the EXIT of BC1
sigdE1_up = C1*sqrt(1e6*(sigdE^2+switch_lh*sigdE_lh^2)+sigdE1_LSC^2+...
            (1e3*sigEL0_ibs)^2+(1e3*sigEL1_ibs)^2+(1e3*sigEBC1_i)^2);
sigdE1_down = 1e3*sigEBC1_w;
sigdE1_tot = sqrt(sigdE1_up^2+sigdE1_down^2);


%%%%%%%%%%%%%%%%%
%% LINAC2: IBS %%
%%%%%%%%%%%%%%%%%

% Derived
Ei2 = Ef1;                                     % initial mean energy in [MeV]
G2 = (Ef2-Ei2)/L2;                             % L1 average accelerating gradient in [MeV/m]
gamma2 = @(s2) (Ei2+G2*s2)/me;                 % Lorentz factor for mean energy ALONG L2
gammaf2 = Ef2/me;                              % Lorentz factor for mean energy at the END of L2
gammam2 = (gammaf2 - Ei2/me)/2;                % Mean Lorentz factor for L2
k1 = @(lambda) k(lambda)*C1;                   % wave number compressed by C1 in DS1 in [1/m]

sigd_L2 = sqrt(C1^2*sigd_BC1_w^2 + sigdBC1_w_ibs^2);

% Max angle
qmax_L2 = sqrt(C1*L2*(Q/e)*re^2/(2*enx^(3/2)*lb*sqrt(betax2)));

% Coulomb logarithm
ln2 = log(qmax_L2*enx/(re*2*sqrt(2)));

% IBS-induced rms absolute energy spread in [MeV] with Bane's approximation
% of B-M expression, cumulated through L1
if switch_bane1 == 1
        % Di Mitri-Perosa closed form for energy-dependent Clog 
        sigdL2_DMP1 = sqrt((2*me*C1*D/(3*gammaf2^2*G2*sqrt(enx*betax2)))*(gammaf2^(3/2)-gammaf1^(3/2)+2*ln2*(gammaf2^(3/2)-gammaf1^(3/2))+gammaf1^(3/2)*log(gammaf1^(3/4))-gammaf2^(3/2)*log(gammaf2^(3/4))));
        % Di Mitri-Perosa closed form for energy-INdependent Clog 
        sigdL2_DMP2 = sqrt((4/3)*(re^2*(Q/e)*ln2*me*C1/(8*enx^(3/2)*sqrt(betax2)*lb*G2))*((gammaf2^(3/2)-gammaf1^(3/2))/gammaf2^2));

     sigdL2_ibs = sigdL2_DMP1;
     sigEL2_ibs = sigdL2_ibs*Ef2;
     
elseif switch_bane1 == 0
    sigdL2_ibs = 0;
    sigEL2_ibs = 0;
end


%%%%%%%%%%%%%%%%%
%% LINAC2: LSC %%
%%%%%%%%%%%%%%%%%

% Average beam radius in [m] and related quantities
rb2 = @(s2) 0.8735*(sqrt(enx*betax2./gamma2(s2))+sqrt(eny*betay2./gamma2(s2)));
a2 = @(s2,lambda) k1(lambda).*rb2(s2)./gamma2(s2);
I21 = @(s2,lambda) besseli(1,a2(s2,lambda));
K20 = @(s2,lambda) besselk(0,a2(s2,lambda));
K21 = @(s2,lambda) besselk(1,a2(s2,lambda));
aw2 =  @(s2,lambda) k1(lambda)*rw2/gamma2(s2);
I20_w = @(s2,lambda) besseli(0,aw2(s2,lambda));
K20_w = @(s2,lambda) besselk(0,aw2(s2,lambda));

% 3-D LSC impedance averaged over transverse dimension (_rr stays for
% chamber shielding)
if switch_shield == 0
    Z2_LSC = @(s2,lambda) 1i*(Z0./(pi*k1(lambda).*rb2(s2).^2)).*(1-2*I21(s2,lambda).*K21(s2,lambda));
elseif switch_shield ~= 0
    Z2_LSC_rr = @(s2,lambda) 1i*(Z0./(pi*k1(lambda)*rb2(s2).^2)).*(1-2*(I21(s2,lambda)./I20_w(s2,lambda)).*(I20_w(s2,lambda)*K21(s2,lambda)+I21(s2,lambda)*K20_w(s2,lambda)));
    Z2_LSC = @(s2,lambda) Z2_LSC_rr(s2,lambda);
end
% Approximated expressions
%%Z2_LSC = @(s2,lambda) 1i*(Z0./(pi*k1(lambda)*rb2(s2).^2))*(1-a2(s2,lambda)*K21(s2,lambda));
%%Z2_LSC = @(s2,lambda) 1i*(Z0./(2*pi*gamma2(s2).*rb2(s2)))*sqrt(1+a2(s2,lambda).^2*(K20(s2,lambda).^2-K21(s2,lambda).^2));

Z2_LSC_real = @(s2,lambda) real(Z2_LSC(s2,lambda));
Z2_LSC_imag = @(s2,lambda) imag(Z2_LSC(s2,lambda));
Dgamma2_real = @(lambda) -(4*pi/Z0)*(I0/IA)*b1(lambda).*integral(@(s2)Z2_LSC_real(s2,lambda),0,L2,'ArrayValued',true);  %%%%
Dgamma2_imag = @(lambda) -(4*pi/Z0)*(I0/IA)*b1(lambda).*integral(@(s2)Z2_LSC_imag(s2,lambda),0,L2,'ArrayValued',true);  %%%%
Dgamma2 = @(lambda) Dgamma2_real(lambda)+1i*Dgamma2_imag(lambda);
Dg2_rel = @(lambda) -(4*pi/Z0)*(I0/IA)*b1(lambda).*integral(@(s2)Z2_LSC(s2,lambda)./gamma2(s2),0,L2,'ArrayValued',true); %%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BC2: Compression factor %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigd02 = sqrt((sigd01*Ef1)^2+(sigdBC1_i_ibs^2+sigdBC1_w_ibs^2/C1^2)*Ef1^2+(sigEL2_ibs/C1)^2)/Ef2;     % fractional energy spread rms, normalized at the DS2 energy, uncompressed
sigd12 = C1*sigd02;                                     % fractional energy spread rms, normalized at the DS2 energy, compressed by C1
I2 = C1*C2*I0;                                          % final peak current in [A]
ex2 = enx/gammaf2;                                      % geometric horizontal emittance at chicane in [m rad]
lambda_co_long_2 = 2*pi*abs(R562)*C2*C1.*sqrt(sigdE^2+switch_lh*sigdE_lh^2+sigEL0_ibs^2+sigEL1_ibs^2+(sigdBC1_i_ibs^2+sigdBC1_w_ibs^2/C1^2)*Ef1^2+(sigEL2_ibs/C1)^2)/Ef2;            % Energy Landau damping cutoff wavelength in [m]
lambda_co_tran_2 = 2*pi*sqrt(ex2*H21);                  % Energy Landau damping cutoff wavelength in [m]


%%%%%%%%%%%%%%
%% BC2: IBS %%
%%%%%%%%%%%%%%

sigd_BC2_i = sqrt((sigd_L2*Ef1)^2+(sigEL2_ibs)^2)/Ef2;     % Uncorrelated energy spread at the entrance of BC2
qmax_BC2_i = sqrt(C1*halfL2*(Q/e)*re^2/(2*gammaf2^(3/2)*enx^(3/2)*lb*sqrt(beta2_in)));
qmax_BC2_w = sqrt(C2*C1*halfL2*(Q/e)*re^2/(2*gammaf2^(3/2)*enx^(3/2)*lb*sqrt(beta2_w)));

% optics functions in the first half (_i) and second half (_w) of BC1
etap2 = theta2;
eta2  = theta2*(DL2+Lb2);
gamma2_in = (1+alpha2_in^2)/beta2_in;
gamma2_w  = (1+alpha2_w^2)/beta2_w;

if switch_bane1 == 1
    if theta2 ~= 0
        % IBS-induced fractional rms energy spread with Bane's formula for the first half of BC2
        dispH2_i = beta2_in*etap2^2+2*alpha2_in*eta2*etap2+gamma2_in*eta2^2;
        
        aBC21 = vpa(2*C1*D/(gammaf2^(3/2)*sqrt(enx*beta2_in)));
        bBC21 = vpa(qmax_BC2_i*enx/(2*sqrt(2)*re));
        hBC21 = vpa(dispH2_i*gammaf2/enx);
        solBC21 = @(y) vpa(abs( -ei(3*log(sqrt(hBC21*sigd_BC2_i^2 + 1)/bBC21)) + ei(3*log(sqrt(hBC21*y^2 + 1)/bBC21)) + vpa(halfL2)*aBC21*hBC21/(2*bBC21^3)));
        yBC21 = vpa(fminsearch(solBC21,1e-6,options));
        sigdBC2_i_ibs = double(sqrt(yBC21^2 - sigd_BC2_i^2));
        
        sigd_BC2_w = sqrt(sigd_BC2_i^2 + sigdBC2_i_ibs^2);                      % Uncorrelated energy spread at the waist of BC2
    
        % IBS-induced fractional rms energy spread with Bane's formula for the second half of BC2
        dispH2_w = beta2_w*etap2^2+2*alpha2_w*eta2*etap2+gamma2_w*eta2^2;
    
        aBC22 = vpa(2*C1*C2*D/(gammaf2^(3/2)*sqrt(enx*beta2_w)));
        bBC22 = vpa(qmax_BC2_w*enx/(2*sqrt(2)*re));
        hBC22 = vpa(dispH1_w*gammaf2/enx);
        solBC22 = @(y) vpa(abs( -ei(3*log(sqrt(hBC22*sigd_BC2_w^2 + 1)/bBC22)) + ei(3*log(sqrt(hBC22*y^2 + 1)/bBC22)) + vpa(halfL2)*aBC22*hBC22/(2*bBC22^3)));
        yBC22 = vpa(fminsearch(solBC22,1e-6,options));
        sigdBC2_w_ibs = double(sqrt(yBC22^2 - sigd_BC2_w^2));
        
    elseif theta2 == 0        
        sigdBC2_i_ibs = sqrt(2*halfL2*C1*D*log(qmax_BC2_i*enx/(re*2*sqrt(2)))./(gammaf2.^2*sqrt(enx*beta2_in./gammaf2)));
        sigd_BC2_w = sqrt(sigd_BC2_i^2 + sigdBC2_i_ibs^2);
        sigdBC2_w_ibs = sqrt(2*halfL2*C1*C2*D*log(qmax_BC2_w*enx/(re*2*sqrt(2)))./(gammaf2.^2*sqrt(enx*beta2_w./gammaf2)));
    end
elseif switch_bane1 == 0
    sigdBC2_i_ibs = 0;
    sigdBC2_w_ibs = 0;
    sigd_BC2_w = sqrt(sigd_BC2_i^2 + sigdBC2_i_ibs^2);
end

% IBS-induced energy spread in BC2, in [MeV] 
sigEBC2_i = sigdBC2_i_ibs*Ef2;
sigEBC2_w = sigdBC2_w_ibs*Ef2;


%%%%%%%%%%%%%%
%% BC2: LSC %%
%%%%%%%%%%%%%%

rbBC2_i = 0.8735*(sqrt(enx*beta2_in/gammaf2)+sqrt(eny*beta2_in/gammaf2));
rbBC2_w = 0.8735*(sqrt(enx*beta2_w/gammaf2)+sqrt(eny*beta2_w/gammaf2));
aBC2_i =  @(lambda) k1(lambda)*rbBC2_i/gammaf2;
aBC2_w =  @(lambda) k1(lambda)*rbBC2_w/gammaf2;
I1_BC2_i = @(lambda) besseli(1,aBC2_i(lambda));
I1_BC2_w = @(lambda) besseli(1,aBC2_w(lambda));
K0_BC2_i = @(lambda) besselk(0,aBC2_i(lambda));
K0_BC2_w = @(lambda) besselk(0,aBC2_w(lambda));
K1_BC2_i = @(lambda) besselk(1,aBC2_i(lambda));
K1_BC2_w = @(lambda) besselk(1,aBC2_w(lambda));
awBC2 =  @(lambda) k1(lambda)*rw0/gammaf2;
I0_w_BC2 = @(lambda) besseli(0,awBC2(lambda));
K0_w_BC2 = @(lambda) besselk(0,awBC2(lambda));

% 3-D LSC impedance averaged over transverse dimension (_rr stays for
% chamber shielding)
if switch_shield == 0
    Z_LSC_BC2_i = @(lambda) 1i*(Z0./(pi*k1(lambda)*rbBC2_i^2)).*(1-2*I1_BC2_i(lambda).*K1_BC2_i(lambda));
    Z_LSC_BC2_w = @(lambda) 1i*(Z0./(pi*k1(lambda)*rbBC2_w^2)).*(1-2*I1_BC2_w(lambda).*K1_BC2_w(lambda));
elseif switch_shield ~= 0
    Z_LSC_rr_BC2_i = @(lambda) 1i*(Z0./(pi*k1(lambda)*rbBC2_i^2))*(1-2*(I1_BC2_i(lambda)/I0_w_BC2(lambda))*(I0_w_BC2(lambda)*K1_BC2_i(lambda)+I1_BC2_i(lambda)*K0_w_BC2(lambda)));
    Z_LSC_rr_BC2_w = @(lambda) 1i*(Z0./(pi*k1(lambda)*rbBC2_w^2))*(1-2*(I1_BC2_w(lambda)/I0_w_BC2(lambda))*(I0_w_BC2(lambda)*K1_BC2_w(lambda)+I1_BC2_w(lambda)*K0_w_BC2(lambda)));
    Z_LSC_BC2_i = @(lambda) Z_LSC_rr_BC2_i(lambda);
    Z_LSC_BC2_w = @(lambda) Z_LSC_rr_BC2_w(lambda);
end

% LSC-induced energy modulation amplitude in me unit
DgammaBC2_i = @(lambda) -(4*pi/Z0)*(I0/IA)*b1(lambda).*Z_LSC_BC2_i(lambda)*halfL2;  %%%%
DgammaBC2_w = @(lambda) -(4*pi/Z0)*(I0/IA)*b1(lambda).*Z_LSC_BC2_w(lambda)*halfL2;  %%%%
DgammaBC2 = @(lambda) DgammaBC2_i(lambda) + DgammaBC2_w(lambda);
DgBC2_rel = @(lambda) DgammaBC2(lambda)/gammaf2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BC2: Laser Heater effect %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the hypergeometric function in its integral form. The integral
% runs over the ratio of laser/e-beam radius: assume here from 1 to 10) for
% the definite integral.

A2_0 = @(lambda) abs((2*pi*C1*C2./lambda)*R562*C1*DeltaE_lh/Ef2);
J02_LH_0 = @(r,lambda) besselj(0,A2_0(lambda).*exp(-r.^2/(4*Blh^2)));
S2_LH_0 = @(lambda) integral(@(r)r.*exp(-r.^2/2).*J02_LH_0(r,lambda),0,Inf,'ArrayValued',true);

A2 = @(lambda) abs((2*pi*C1*C2./lambda)*R562*C1*DeltaE_lh/Ef2);
J02_LH = @(r,lambda) besselj(0,A2(lambda).*exp(-r.^2/(4*Blh^2)));
S2_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*J02_LH(r,lambda),0,Inf,'ArrayValued',true);

if switch_lh == 1 && theta2 ~= 0 && sigdE_lh ~= 0
    S2_LH_0 = @(lambda) S2_LH_0(lambda);
    S2_LH = @(lambda) S2_LH(lambda);
elseif switch_lh == 0 || theta2 == 0 || sigdE_lh == 0
    S2_LH_0 = @(lambda) 1;
    S2_LH = @(lambda) 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CSR GAIN in BC2 vs. uncompressed wavelength %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if switch_csr ~= 0 && theta2 ~= 0
%sigT21 = @(lambda) k1(lambda)*sqrt(ex2*H21); %come sopra
sigT21 = @(lambda) k1(lambda)*theta2*sqrt(ex2*beta2_in);
erf_T21 = @(lambda) erf(sigT21(lambda));
%sigT22 = @(lambda) k1(lambda)*sqrt(ex2*H22);
%erf_T22 = @(lambda) erf(sigT22(lambda));
I2n_csr = @(lambda) (I2*R562*k1(lambda).^(4/3)*theta2^(2/3)*Lb2^(1/3)/(IA*gammaf2));
phi2 = 2*DL2/beta2_in;

% Ht21 (Ht22) takes care of transverse Landau damping through Dipole-1 (Dipole-2+Dipole3) of BC2-DS2
Ht21 = @(lambda,t) exp(-sigT21(lambda).^2*((1-2.*t+alpha2_in*phi2.*t).^2+(phi2.*t).^2)./(1+h2*R562*t).^2-(C2*R562*k1(lambda)*sigd12).^2*(t.^2+(C2*(1-t)).^2)/(2*(1+h2*R562*t).^2));
%Ht22 = @(lambda,t) exp(-sigT22(lambda).^2*((1-2.*t+alpha2_in*phi2.*t).^2+(phi2.*t).^2)/(1+h2*R562*t).^2-(C2*R562*k1(lambda)*sigd12).^2*(t.^2+(C2*(1-t)).^2)/(2*(1+h2*R562*t).^2));

F02 = @(lambda) (exp(-sigT21(lambda).^2)+sigT21(lambda)*sqrt(pi).*erf_T21(lambda)-1)./(2*sigT21(lambda).^2);
F12 = @(lambda) 2*integral(@(t)(1-t).*Ht21(lambda,t)./(1+h2*R562*t).^(4/3),0,1,'ArrayValued',true);
F22 = @(lambda) 2*integral(@(t)(1-t).*t.*Ht21(lambda,t)./(C2*(1+h2*R562*t).^(7/3)),0,1,'ArrayValued',true);

% this is the CSR gain driven by initial bunching, at first order in beam current
b21_csr = @(lambda) -b1(lambda).*A_csr.*I2n_csr(lambda).*(sqrt(pi).*erf_T21(lambda).*exp(-0.5.*(C2*R562*k1(lambda)*sigd12).^2).*S2_LH(lambda)./(2*sigT21(lambda))+F12(lambda));

% this is the CSR gain driven by initial bunching, at second order in beam current
b22_csr = @(lambda)  b1(lambda).*A_csr^2.*I2n_csr(lambda).^2.*F02(lambda).*F22(lambda);

% this is the CSR gain driven by incoming energ modulation (LSC) in D0+L1, at first order in beam current
b23_csr = @(lambda) -1i*C2*R562*k1(lambda).*(Dg2_rel(lambda)).*A_csr*I2n_csr(lambda).*F22(lambda);

elseif switch_csr == 0 || theta2 == 0
    b21_csr = @(lambda) 0;
    b22_csr = @(lambda) 0;
    b23_csr = @(lambda) 0;
end

Gain2_csr = @(lambda) abs((b21_csr(lambda)+b22_csr(lambda)+b23_csr(lambda))/b0(lambda));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GAIN after D0+L1+BC1+L2+BC2 vs. uncompressed wavelength %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expression is for any gain regime

if switch_lowgain ~= 0
 
% bunching due to energy modulation cumulated through L0+L1 surviving up to BC2
b20 = @(lambda) -1i*C2*k1(lambda).*R562.*((DgBC1_rel(lambda)*Ef1 + Dg1_rel(lambda)*Ef1 + Dg0_rel(lambda)*Ef0)/Ef2).*...
    exp(-0.5*C2^2*k1(lambda).^2*R562^2*sigd12^2).*S2_LH_0(lambda);
else 
    b20 = @(lambda) 0;
end

% bunching due to the initial bunching arriving at BC2
b21 = @(lambda) b1(lambda).*exp(-0.5*C2^2*k1(lambda).^2*R562^2*sigd12^2).*S2_LH_0(lambda);

% bunching induced by the energy modulation cumulated through L2
b22 = @(lambda) -1i*C2*k1(lambda).*R562.*(Dg2_rel(lambda)).*...
    exp(-0.5*C2^2*k1(lambda).^2*R562^2*sigd12^2).*S2_LH(lambda);

b2 = @(lambda) b20(lambda)+b21(lambda)+b22(lambda)+b21_csr(lambda)+b22_csr(lambda)+b23_csr(lambda);

Gain2 = @(lambda) abs(b2(lambda)./b0(lambda));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Uncorrelated Energy Spread at BC2 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncorrelated energy spread RMS in induced by LSC through L2+L3 in [keV]
if R561 ~= 0
    sigdE2_LSC = 1e3*me*sqrt((2*e*c/I0)*integral(@(lambda)abs(Gain1(lambda).^2.*(Dgamma2(lambda)).^2./(b1(lambda).^2.*lambda.^2)),5e-6,200e-6,'ArrayValued',true));
elseif R561 == 0
    sigdE2_LSC = 1e3*me*sqrt((2*e*c/I0)*integral(@(lambda)abs((Dgamma0(lambda)+Dgamma1(lambda)+Dgamma2(lambda)+DgammaBC1(lambda)).^2./(b0(lambda).^2.*lambda.^2)),1e-6,200e-6,'ArrayValued',true));
end

% Total uncorrelated energy spread RMS in [keV], at the EXIT of BC2
sigdE2_up = C2*sqrt(C1^2*1e6*(sigdE^2+switch_lh*sigdE_lh^2+sigEL0_ibs^2+sigEL1_ibs^2+sigEBC1_i^2)+1e6*sigEBC1_w^2+...
            sigdE2_LSC^2+(1e3*sigEL2_ibs)^2+(1e3*sigEBC2_i)^2);
sigdE2_down = 1e3*sigEBC2_w;
sigdE2_tot = sqrt(sigdE2_up^2+sigdE2_down^2);


%%%%%%%%%%%%%%%%%
%% LINAC3: IBS %%
%%%%%%%%%%%%%%%%%

% Derived
Ei3 = Ef2;                        % initial mean energy in [MeV]
G3 = (Ef3-Ei3)/L3;                % L3 average accelerating gradient in [MeV/m]
gamma3 = @(s3) (Ei3+G3*s3)/me;    % Lorentz factor for mean energy ALONG L3
gammaf3 = Ef3/me;                 % Lorentz factor for mean energy at the END of L3
gammam3 = (gammaf3 - Ei3/me)/2;   % Mean Lorentz factor for L3
k2 = @(lambda) k1(lambda)*C2;  % wave number compressed by C1+C2 in DS1+DS2 in [1/m]

sigd_L3 = sqrt(C2^2*sigd_BC2_w^2 + sigdBC2_w_ibs^2);      % Uncorrelated energy spread at L3

% Max angle
qmax_L3 = sqrt(C2*C1*L3*(Q/e)*re^2/(2*enx^(3/2)*lb*sqrt(betax3)));

% Coulomb logarithm
ln3 = log(qmax_L3*enx/(re*2*sqrt(2)));

% IBS-induced energy spread with Bane's formula
if switch_bane1 == 1     
        % Di Mitri-Perosa closed form for energy-dependent Clog 
        sigdL3_DMP1 = sqrt((2*me*D*C1*C2/(3*gammaf3^2*G3*sqrt(enx*betax3)))*(gammaf3^(3/2)-gammaf2^(3/2)+2*ln3*(gammaf3^(3/2)-gammaf2^(3/2))+gammaf2^(3/2)*log(gammaf2^(3/4))-gammaf3^(3/2)*log(gammaf3^(3/4))));
        % Di Mitri-Perosa closed form for energy-INdependent Clog 
        sigdL3_DMP2 = sqrt((4/3)*(re^2*(Q/e)*ln3*me*C1*C2/(8*enx^(3/2)*sqrt(betax3)*lb*G3))*((gammaf3^(3/2)-gammaf2^(3/2))/gammaf3^2));

     sigdL3_ibs = sigdL3_DMP1;
     sigEL3_ibs = sigdL3_ibs*Ef3;

elseif switch_bane1 == 0
    sigdL3_ibs = 0;
    sigEL3_ibs = 0;
end


%%%%%%%%%%%%%%%%%
%% LINAC3: LSC %%
%%%%%%%%%%%%%%%%%

% Average beam radius in [m] and related quantities
rb3 = @(s3) 0.8735*(sqrt(enx*betax3./gamma3(s3))+sqrt(eny*betay3./gamma3(s3)));
a3 = @(s3,lambda) k2(lambda).*rb3(s3)./gamma3(s3);
I31 = @(s3,lambda) besseli(1,a3(s3,lambda));
K30 = @(s3,lambda) besselk(0,a3(s3,lambda));
K31 = @(s3,lambda) besselk(1,a3(s3,lambda));
aw3 =  @(s3,lambda) k2(lambda)*rw3./gamma3(s3);
I30_w = @(s3,lambda) besseli(0,aw3(s3,lambda));
K30_w = @(s3,lambda) besselk(0,aw3(s3,lambda));

% 3-D LSC impedance averaged over transverse dimension (_rr stays for
% chamber shielding)
if switch_shield == 0
    Z3_LSC = @(s3,lambda) 1i*(Z0./(pi*k2(lambda).*rb3(s3).^2)).*(1-2*I31(s3,lambda).*K31(s3,lambda));
elseif switch_shield ~= 0
    Z3_LSC_rr = @(s3,lambda) 1i*(Z0./(pi*k2(lambda)*rb3(s3).^2)).*(1-2*(I31(s3,lambda)./I30_w(s3,lambda)).*(I30_w(s3,lambda)*K31(s3,lambda)+I31(s3,lambda)*K30_w(s3,lambda)));
    Z3_LSC = @(s3,lambda) Z3_LSC_rr(s3,lambda);
end
% Approximated expressions
%Z3_LSC = @(s3,lambda) 1i*(Z0./(pi*k1(lambda)*rb3(s3).^2))*(1-a3(s3,lambda)*K31(s3,lambda));
%Z3_LSC = @(s3,lambda) 1i*(Z0./(2*pi*gamma3(s3).*rb3(s3)))*sqrt(1+a3(s3,lambda).^2*(K30(s3,lambda).^2-K31(s3,lambda).^2));

% LSC-induced relative energy modulation amplitude
Z3_LSC_real = @(s3,lambda) real(Z3_LSC(s3,lambda));
Z3_LSC_imag = @(s3,lambda) imag(Z3_LSC(s3,lambda));

Dgamma3_real = @(lambda) -(4*pi/Z0)*(I0/IA)*b2(lambda).*integral(@(s3)Z3_LSC_real(s3,lambda),0,L3,'ArrayValued',true); 
Dgamma3_imag = @(lambda) -(4*pi/Z0)*(I0/IA)*b2(lambda).*integral(@(s3)Z3_LSC_imag(s3,lambda),0,L3,'ArrayValued',true);  
Dgamma3 = @(lambda) Dgamma3_real(lambda)+1i*Dgamma3_imag(lambda);
Dg3_rel = @(lambda) -(4*pi/Z0)*(I0/IA)*b2(lambda).*integral(@(s3)Z3_LSC(s3,lambda)./gamma3(s3),0,L3,'ArrayValued',true); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Uncorrelated Energy Spread at Linac End %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncorrelated energy spread RMS in [keV] induced by LSC 
if R562 ~=0
    sigdE3_LSC = 1e3*me*sqrt((2*e*c/I0)*integral(@(lambda)abs(Gain2(lambda).^2.*(Dgamma3(lambda)).^2./(b2(lambda).^2.*lambda.^2)),5e-6,200e-6,'ArrayValued',true));
elseif R561 ~= 0 && R562 == 0
    sigdE3_LSC = 1e3*me*sqrt((2*e*c/I0)*integral(@(lambda)abs(Gain1(lambda).^2.*...
        (Dgamma2(lambda)+DgammaBC2(lambda)+Dgamma3(lambda)).^2./(b1(lambda).^2.*lambda.^2)),1e-6,200e-6,'ArrayValued',true));
elseif R561 == 0 && R562 == 0
    sigdE3_LSC = 1e3*me*sqrt((2*e*c/I0)*integral(@(lambda)abs((Dgamma0(lambda)+Dgamma1(lambda)+Dgamma2(lambda)+Dgamma3(lambda)+...
        DgammaBC1(lambda)+DgammaBC2(lambda)).^2./(b0(lambda).^2.*lambda.^2)),1e-6,200e-6,'ArrayValued',true));
end

% Total uncorrelated energy spread RMS in [keV],at the linac end
sigdE3_tot = sqrt(C2^2*C1^2*1e6*(sigdE^2+switch_lh*sigdE_lh^2+sigEL0_ibs^2+sigEL1_ibs^2+sigEBC1_i^2)+C2^2*1e6*(sigEBC1_w^2+sigEL2_ibs^2+sigEBC2_i^2)+...
                 1e6*sigEBC2_w^2+1e6*sigEL3_ibs^2+sigdE3_LSC^2);
             
     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Energy modulation ampltiudes along the linac %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dgamma1_keV = @(lambda) abs((Dgamma0(lambda)+Dgamma1(lambda))*me*1000);
Dgamma2_keV = @(lambda) abs((Dgamma0(lambda)+Dgamma1(lambda)+DgammaBC1(lambda))*me*1000);
Dgamma3_keV = @(lambda) abs((Dgamma0(lambda)+Dgamma1(lambda)+DgammaBC1(lambda)+Dgamma2(lambda))*me*1000);
Dgamma4_keV = @(lambda) abs((Dgamma0(lambda)+Dgamma1(lambda)+DgammaBC1(lambda)+Dgamma2(lambda)+DgammaBC2(lambda))*me*1000);
Dgamma5_keV = @(lambda) abs((Dgamma0(lambda)+Dgamma1(lambda)+DgammaBC1(lambda)+Dgamma2(lambda)+DgammaBC2(lambda)+Dgamma3(lambda))*me*1000);




%%%%%%%%%%%
%% Plots %%
%%%%%%%%%%%

% Mean Energy in MeV
f1=figure(1);
xlabel('s [m]','FontSize',16)
ylabel('Mean Energy [MeV]','FontSize',16)
set(gca,'FontSize',16)
set(gcf,'DefaultLineLineWidth',3)
energy0 = gamma0*me;
energy1 = @(s1) gamma1(s1)*me;
energy2 = @(s2) gamma2(s2)*me;
energy3 = @(s3) gamma3(s3)*me;
u1=0:0.1:L1;
u2=0:0.1:L2;
u3=0:0.1:L3;
line([0,L0],[energy0,energy0]);
hold on
plot(u1+L0,energy1(u1),'r-')
hold on
plot(u2+L0+L1,energy2(u2),'m-')
hold on
plot(u3+L0+L1+L2,energy3(u3),'b-')
xlabel('s [m]','FontSize',16)
ylabel('Mean Energy [MeV]','FontSize',16)
set(gca,'FontSize',16)
saveas(f1,['./HK_mean_energy_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_','.jpg'])

% Effective beam radius in micron
f2=figure(2);
xlabel('s [m]','FontSize',16)
ylabel('Effective Beam Radius [\mum]','FontSize',16)
set(gca,'FontSize',16)
set(gcf,'DefaultLineLineWidth',3)
rb1_um = @(s1) rb1(s1)*1e6;
rb2_um = @(s2) rb2(s2)*1e6;
rb3_um = @(s3) rb3(s3)*1e6;
v1=0:0.1:L1;
v2=0:0.1:L2;
v3=0:0.1:L3;
line([0,L0],[rb0*1e6,rb0*1e6]);
hold on
plot(v1+L0,rb1_um(v1),'r-')
hold on
plot(v2+L0+L1,rb2_um(v2),'m-')
hold on
plot(v3+L0+L1+L2,rb3_um(v3),'b-')
saveas(f2,['./HK_beam_size_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_','.jpg'])

% LSC-induced energy modulation amplitude in keV vs. uncompressed
% wavelength in micron
% Dg3_keV = @(lambda) abs((Dgamma2(lambda))/b1(lambda));
f3=figure(3);
fplot(Dgamma1_keV,[1e-6,200e-6],'g:','LineWidth',3)
hold on
fplot(Dgamma2_keV,[1e-6,200e-6],'r:','LineWidth',3)
hold on
fplot(Dgamma3_keV,[1e-6,200e-6],'c-','LineWidth',3)
hold on
fplot(Dgamma4_keV,[1e-6,200e-6],'m-.','LineWidth',3)
hold on
fplot(Dgamma5_keV,[1e-6,200e-6],'b-','LineWidth',3)
legend('Entrance of BC1','Exit of BC1','Entrance of BC2','Exit of BC2','Linac End')
xlabel('\lambda_0 [\mum]','FontSize',16)
ylabel('\DeltaE_{MBI} [keV]','FontSize',16)
set(gca,'XTickLabel',{'50','100','150','200'})
set(gca,'FontSize',16)
saveas(f3,['./HK_energy_mod_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_','.jpg'])

% Bunching factor at several stages vs. uncompressed wavelength in micron
abs_b0 = @(lambda) abs(b0(lambda))*100;
abs_b1 = @(lambda) abs(b1(lambda))*100;
abs_b2 = @(lambda) abs(b2(lambda))*100;

f4=figure(4);
fplot(abs_b0,[1e-6,200e-6],'g-','LineWidth',3)
hold on
fplot(abs_b1,[1e-6,200e-6],'r:.','LineWidth',3)
hold on
fplot(abs_b2,[1e-6,200e-6],'b-','LineWidth',3)
legend('Initial','Exit of BC1','Exit of BC2')
xlabel('\lambda_0 [\mum]','FontSize',16)
ylabel('Bunching Factor [%]','FontSize',16)
set(gca,'XTickLabel',{'50','100','150','200'})
set(gca,'FontSize',16)
saveas(4,['./HK_bunching_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_','.jpg'])

% Spectral gain (LSC+CSR)
f5=figure(5);
fplot(Gain1,[1e-6,200e-6],'r:','LineWidth',3)
hold on
fplot(Gain2,[1e-6,200e-6],'b-','LineWidth',3)
legend('Exit of BC1','Exit of BC2')
xlabel('\lambda_0 [\mum]','FontSize',16)
ylabel('Gain','FontSize',16)
set(gca,'XTickLabel',{'50','100','150','200'})
set(gca,'FontSize',16)
saveas(f5,['./HK_gain_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_','.jpg'])    

% Spectral CSR gain (LSC+CSR, but plotting only CSR bunching)
f6=figure(6);
fplot(Gain1_csr,[1e-6,200e-6],'r:','LineWidth',3)
hold on
fplot(Gain2_csr,[1e-6,200e-6],'b-','LineWidth',3)
legend('Exit of BC1','Exit of BC2')
xlabel('\lambda_0 [\mum]','FontSize',16)
ylabel('Gain_{CSR}','FontSize',16)
set(gca,'XTickLabel',{'50','100','150','200'})
set(gca,'FontSize',16)
saveas(f6,['./HK_gain_csr_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_','.jpg'])    



%%%%%%%%%%%
%% Print %%
%%%%%%%%%%%

question='Is the LH turned ON ?';
yes='yes';
no='no';
if switch_lh == 1
    fprintf('%s %s \n', question, yes);
elseif switch_lh == 0
    fprintf('%s %s \n', question, no);
end

question='Is the CSR effect included ?';
yes='yes';
no='no';
if switch_csr == 1
    fprintf('%s %s \n', question, yes);
elseif switch_csr == 0
    fprintf('%s %s \n', question, no);
end

question='Is low gain bunching included ?';
yes='yes';
no='no';
if switch_lowgain == 1
    fprintf('%s %s \n', question, yes);
elseif switch_lowgain == 0
    fprintf('%s %s \n', question, no);
end

question='Is IBS included ?';
yes='yes';
no='no';
if switch_bane1 == 1
    fprintf('%s %s \n', question, yes);
elseif switch_bane1 == 0
    fprintf('%s %s \n', question, no);
end

question='Is ZLSC shielding included ?';
yes='yes';
no='no';
if switch_shield == 1
    fprintf('%s %s \n', question, yes);
elseif switch_shield == 0
    fprintf('%s %s \n', question, no);
end

label='sigE_gun =';
value=1e3*sigdE;
fprintf('%s %d keV \n', label, value);

label='sigE_LH =';
value=1e3*sigdE_lh;
fprintf('%s %d keV \n', label, value);

label='I0 =';
value=I0;
fprintf('%s %d A \n', label, value);

label='I1 =';
value=I1;
fprintf('%s %d A \n', label, value);

label='I2 =';
value=I2;
fprintf('%s %d A \n', label, value);

label='C1 =';
value=C1;
fprintf('%s %d \n', label, value);

label='C2 =';
value=C2;
fprintf('%s %d \n', label, value);

label='R56_BC1 =';
value=1e3*R561;
fprintf('%s %d mm \n', label, value);

label='R56_BC2 =';
value=1e3*R562;
fprintf('%s %d mm \n', label, value);

label='sigDelta_corr(BC1) =';
value=sigd_BC1_cor;
fprintf('%s %d \n', label, value);

label='sigDelta_corr(BC2) =';
value=sigd_BC2_cor;
fprintf('%s %d \n', label, value);

label='lambda_cutoff (delta,BC1) = ';
value=1e6*lambda_co_long_1;
fprintf('%s %d um \n', label, value);

label='lambda_cutoff(delta,BC2) = ';
value=1e6*lambda_co_long_2;
fprintf('%s %d um \n', label, value);

label='lambda_cutoff (ex,BC1) = ';
value=1e6*lambda_co_tran_1;
fprintf('%s %d um \n', label, value);

label='lambda_cutoff (ex,BC2) = ';
value=1e6*lambda_co_tran_2;
fprintf('%s %d um \n', label, value);


label='sigE_IBS (L0+L1) = ';
value=1e3*(sigEL0_ibs+sigEL1_ibs);
fprintf('%s %d keV \n', label, value);

label='sigE_IBS (BC1) = ';
value=1e3*sqrt(sigEBC1_i^2+sigEBC1_w^2);
fprintf('%s %d keV \n', label, value);

label='sigE_IBS (L2+L3) = ';
value=1e3*sigEL2_ibs;
fprintf('%s %d keV \n', label, value);

label='sigE_IBS (BC2) = ';
value=1e3*sqrt(sigEBC2_i^2+sigEBC2_w^2);
fprintf('%s %d keV \n', label, value);

label='sigE_IBS (L4) = ';
value=1e3*sigEL3_ibs;
fprintf('%s %d keV \n', label, value);


if R561 == 0 && R562 == 0
    label='sigE_COLL (L0+L1+BC1+L2+L3+BC2+L4) = ';
    value=sigdE3_LSC;
    fprintf('%s %d keV \n', label, value);
    
    label='sigE_tot (@L4) = ';
    value=sigdE3_tot;
    fprintf('%s %d keV \n', label, value);
end

if R561 ~= 0 && R562 == 0  
    label='sigE_COLL (L0+L1) = ';
    value=sigdE1_LSC;
    fprintf('%s %d keV \n', label, value);
    
    label='sigE_tot (@BC1) = ';
    value=sigdE1_tot;
    fprintf('%s %d keV \n', label, value);
    
    label='sigE_COLL (L2+L3+BC2+L4) = ';
    value=sigdE3_LSC;
    fprintf('%s %d keV \n', label, value);
    
    label='sigE_tot (@L4) = ';
    value=sigdE3_tot;
    fprintf('%s %d keV \n', label, value);
end

if R561 == 0 && R562 ~= 0
    label='sigE_COLL (L0+L1+BC1+L2) = ';
    value=sigdE2_LSC;
    fprintf('%s %d keV \n', label, value);
    
    label='sigE_tot (@BC2) = ';
    value=sigdE2_tot;
    fprintf('%s %d keV \n', label, value);
    
    label='sigE_COLL (L4) = ';
    value=sigdE3_LSC;
    fprintf('%s %d keV \n', label, value);
    
    label='sigE_tot (@L4) = ';
    value=sigdE3_tot;
    fprintf('%s %d keV \n', label, value);
end

if R561 ~= 0 && R562 ~= 0
    label='sigE_COLL (L0+L1) = ';
    value=sigdE1_LSC;
    fprintf('%s %d keV \n', label, value);
    
    label='sigE_tot (@BC1) = ';
    value=sigdE1_tot;
    fprintf('%s %d keV \n', label, value);
    
    label='sigE_COLL (L2+L3) = ';
    value=sigdE2_LSC;
    fprintf('%s %d keV \n', label, value)
    
    label='sigE_tot (@BC2) = ';
    value=sigdE2_tot;
    fprintf('%s %d keV \n', label, value)
    
    label='sigE_COLL (L4) = ';
    value=sigdE3_LSC;
    fprintf('%s %d keV \n', label, value)
    
    label='sigE_tot (@L4) = ';
    value=sigdE3_tot;
    fprintf('%s %d keV \n', label, value)
end

