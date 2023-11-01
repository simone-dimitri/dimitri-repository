%% MODELING: based on the Bosch-Kleman formalism for the microbunching
%% Beam size and energy vary with length. Normalized emittances are constant. Betatron functions are average values along linac sections.
%% 3-D LSC averaged (note: sign of bunching is consistent with R56<0 in a chicane). Staged LSC-microbunching.
%% Z_LSC is averaged over the trasverse distribution, whose average area is constant and = pi*r^2, in free space [M.Venturini+J.Qiang]
%% Effective beam radius is 1.747(sigx+sigy) [M.Venturini]
%% Z_LSC_rr is like Z_LSC but for a round beam in a round vacuum chamber (only important at long wavelengths) [K.Wang,Y.Li]
%% Laser Heater effect is for a beam initial Gaussian energy and transverse distribution -> hypergeometric function [Z.Huang et al.]
%% CSR and transverse Landau damping is for a 4-dipoles chicane [Z.Huang,K.Kim]
%% IBS in all straight sections (constant and varying energy) and chicanes (half+half) [Bane's approx. of B-M formalism]
%% Blog includes a correction factor of 2 w.r.t. Bane, and matches Z. Huang. Raubenheimer's cut with damping time is replaced by time of traveling.
%% Low gain terms are included by default.

clc
clear variables
close all

%% Layout:
% Drift(D0) -> Linac(L1) -> Disp.Section(BC1) -> Linac(L2+L3) -> Disp.Section(BC2) -> Linac(L4) -> Spreader -> Modulator(MOD1) -> Drift(D1) -> Disp.Section(DL) -> Modulator(MOD2)


%% Constant
A_csr = 1.63-1i*0.9454;   % CSR free vacuum impedance
Z0 = 120*pi;              % vacuum impedance n [Ohm]
e = 1.6e-19;              % electron charge in [C]
me = 0.511e6;               % electron rest mass in [eV]
re = 2.818e-15;           % electron classical radius in [m]
IA = 17045;               % Alfven current in [A]
c = 2.998e8;              % light speed in vacuum in [m/s]

%%%%%%%%%%%%%%
%% SWITCHES %%
%%%%%%%%%%%%%%

switch_lh = 0;            % 1 means heating depends on transverse size, 0 means heating is purely Gaussian
switch_R561_spec = 1;     % 1 means user-specified R561, 0 means R561 of a 4-dipole chicane (theta1)
switch_R562_spec = 1;     % 1 means user-specified R561, 0 means R562 of a 4-dipole chicane (theta2)
switch_csr = 0;           % 1 means CSR ON, 0 means CSR OFF
switch_cer = 0;           % 1 means CER ON, 0 means CER OFF
switch_bane1 = 0;         % 1 means IBS effects ON, 0 means IBS effects OFF
switch_shield = 0;        % 1 means LSC impedance with shielding, 0 means LSC impedance in vacuum

switch_spreader = 0;      % 1 means spreader ON, 0 means spreader OFF (identity matrix)
switch_FEL = 1;           % 1 means Spreader-FEL1, 2 means Spreader-FEL2 
switch_Hx = 0;            % 1 means transverse damping is ON in Spreader, 0 means transverse damping is OFF in Spreader
switch_EEHG = 0;          % 1 means EEHG modulators and DL included, 0 means Spreader only

flagLH = {'noLH','LH'};
flagIBS = {'noIBS','IBS'};
flagSpreader = {'woS','wS','wS','wS'};


%%%%%%%%%%%%%%%%%%%%%%
%% Input Parameters %%
%%%%%%%%%%%%%%%%%%%%%%

%% General
Q = 50e-12;          % total bunch charge in [C]
I0 = 7;              % initial peak current in [A]
E0 = 100e6;            % initial mean energy in [eV]
enx = 0.5e-6;           % normalized horizontal emittance rms in [m rad]
eny = 0.5e-6;           % normalized vertical emittance rms in [m rad]
sigdE = 2e3;          % natural initial uncorrelated energy spread rms in [eV]

%% Laser Heater (Hypergeometric Function)
Blh = 1.2;         % ratio of laser over electron beam radius (>=1) in the LH
sigdE_lh = 6e3;    % RMS energy spread induced by the LH in [eV]

%% LH Area
L0 = 5;            % drift path length in [m]
Ei0 = E0;           % initial mean energy in [eV]
Ef0 = E0;           % final mean energy in [MeV]
betax0 = 8;         % average horizontal betatron function in [m]
betay0 = 8;         % average vertical betatron function in [m]
rw0 = 10e-3;        % average inner radius of round vacuum chamber in [m]

%% LINAC1
L1 = 95;            % Linac1 path length in [m]
Ei1 = Ef0;          % initial mean energy in [eV]
Ef1 = 580e6;          % final mean energy in [eV]
betax1 = 10;        % average horizontal betatron function in [m]
betay1 = 10;        % average vertical betatron function in [m]
rw1 = 10e-3;        % average inner radius of round vacuum chamber in [m]

%% BC1
C1 = 1;            % linear compression factor
theta1 = 0.;     % dipole bending angle in [rad]
Lb1 = 0.67;        % dipole length in [m]
DL1 = 1;         % drift length between outer and inner dipole of chicane, in [m]
D112 = 1.0;         % drift length after dipole-2
D113 = DL1;
D114 = D112;
alpha1_in = 1;      % horizontal alpha-function at entrance of chicane
alpha1_w = 0;
beta1_in = 5;      % horizontal betatron function at entrance of chicane, in [m]
beta1_w = 5;        % horizontal betatron function at waist in the second half of chicane, in [m]
R561_spec = -0;      % R561 user-defined in [m], replaces R561 calculated for a chicane from theta1

%% LINAC2+LINAC3
L2 = 130;            % Linac(2+3) path length in [m]
Ef2 = 1000e6;        % final mean energy in [eV]
betax2 = 7;        % average horizontal betatron function in [m]
betay2 = 7;        % average vertical betatron function in [m]
rw2 = 10e-3;        % average inner radius of round vacuum chamber in [m]

%% BC2
C2 = 100;             % linear compression factor
theta2 = 0.19;         % dipole bending angle in [rad]
Lb2 = 0.67;        % dipole length in [m]
DL2 = 2;         % drift length between outer and inner dipole of chicane, in [m]
D212 = 1.0;         % drift length after dipole-2
D213 = DL2;
D214 = D212;
alpha2_in = 1;      % horizontal alpha-function at entrance of chicane
alpha2_w = 0;
beta2_in = 3;      % horizontal betatron function at entrance of chicane, in [m]
beta2_w = 3;        % horizontal betatron function at waist in the second half of chicane, in [m]
R562_spec = -0.19;      % R562 user-defined in [m], replaces R561 calculated for a chicane from theta2

%% LINAC4
L3 = 30;            % Linac4 path length in [m]
Ef3 = 1000e6;         % final mean energy in [MeV]
betax3 = 15;        % average horizontal betatron function in [m]
betay3 = 15;        % average vertical betatron function in [m]
rw3 = 10e-3;        % average inner radius of round vacuum chamber in [m]

%% Spreader
Lbs = 0.4;                      % Spreader dipole magnets' length in [m]
thetas = 6*pi/180;              % Spreader dipole magnets' bending angle in [rad]
betaxdr0 = 25;                  % horizontal betatron function along the Spreader drift sections, in [m]
betaydr0 = 25;                  % vertical betatron function along the Spreader drift sections, in [m]
betaxs1 = 50;                   % horizontal betatron function at the Spreader dipoles' magnets, in [m] 
betays1 = 100;                  % vertical betatron function at the Spreader dipoles' magnets, in [m]
betaxdr1 = 25;                  % horizontal betatron function along the Spreader drift sections, in [m]
betaydr1 = 25;                  % vertical betatron function along the Spreader drift sections, in [m]
betaxs2 = 5;                    % horizontal betatron function at the Spreader dipoles' magnets, in [m] 
betays2 = 10;                   % vertical betatron function at the Spreader dipoles' magnets, in [m]
betaxdr2 = 25;                  % horizontal betatron function along the Spreader drift sections, in [m]
betaydr2 = 25;                  % vertical betatron function along the Spreader drift sections, in [m]
betaxs3 = 10;                   % horizontal betatron function at the Spreader dipoles' magnets, in [m] 
betays3 = 10;                   % vertical betatron function at the Spreader dipoles' magnets, in [m]
betaxdr3 = 25;                  % horizontal betatron function along the Spreader drift sections, in [m]
betaydr3 = 25;                  % vertical betatron function along the Spreader drift sections, in [m]
betaxs4 = 5;                    % horizontal betatron function at the Spreader dipoles' magnets, in [m] 
betays4 = 10;                   % vertical betatron function at the Spreader dipoles' magnets, in [m]
betaxdr4 = 25;                  % horizontal betatron function along the Spreader drift sections, in [m]
betaydr4 = 25;                  % vertical betatron function along the Spreader drift sections, in [m]

R55_S1 = 1;
R55_S2 = 1;
R55_S3 = 1;
R55_S4 = 1;

R56_S1 = Lbs*(1-sin(thetas)/thetas);    % R56 of B_SCL.01 in [m]
R56_S2 = -1.66e-3;                      %-Lbs*(1-sin(thetas)/thetas);  %1.34e-3;    % R56 of B_SCL.02 in [m]
R56_S3 = -1.66e-3;                      %-Lbs*(1-sin(thetas)/thetas);  %1.34e-3;    % R56 of B_SFEL02.01 or B_SFEL01.01 in [m]
R56_S4 = Lbs*(1-sin(thetas)/thetas);    % R56 of B_SFEL02.02 or B_SFEL01.02 in [m]

LS0 = 5;                                % drift length before spreader
if switch_FEL == 1
    LS1 = 4;                            % drift length after B_SCL.01, in [m]
    LS2 = 13.5;                         % drift length after B_SCL.02, in [m]
    LS3 = 4;                            % drift length after B_SFEL01.01, in [m]
    LS4 = 9;                            % drift length from B_SFEL01.02 to undulator, in [m]
elseif switch_FEL == 2
    LS1 = 4;                            % drift length after B_SCL.01, in [m]
    LS2 = 4;                            % drift length after B_SCL.02, in [m]
    LS3 = 4;                            % drift length after B_SFEL02.01, in [m]
    LS4 = 8;                            % drift length from B_SFEL02.02 to undulator, in [m]
end

%% Modulator 1
Lmod1 = 3.2;                            % Modulator length in [m]
Kmod1 = 4.8542;
betaxm1 = 15;
betaym1 = 15;
DeltaE_sl1 = 0;                         % Energy Modulations induced by Seed Laser 1 [eV]

%% Undulators
Lu = 16;
betaxu = 15;
betayu = 15;

%% Delay Line
betax4 = 10;
betay4 = 10;
theta3 = 0.056;
Lb3 = 0.18;         % dipole length in [m]
DL3 = 0.24;         % drift length between outer and inner dipole of chicane, in [m]
D312 = 0.1;         % drift length after dipole-2
D313 = DL3;
D314 = D312;
alpha3_in = 3;      % horizontal alpha-function at entrance of chicane
alpha3_w = 0;
beta3_in = 10;      % horizontal betatron function at entrance of chicane, in [m]
beta3_w = 2;        % horizontal betatron function at waist in the second half of chicane, in [m]

%% Modulator 2
Lmod2 = 1.47;
Kmod2 = 4.5538;
betaxm2 = 20;
betaym2 = 20;

%% ------------------------------------------------------------------------- %%

%%%%%%%%%%%%%
%% Derived %%
%%%%%%%%%%%%%

sigd0 = sigdE/E0;                                            % initial fractional uncorrelated energy spread rms
DeltaE_lh = 2*sigdE_lh;                                      % half-Amplitude of the LH-induced energy modulation, in [MeV]
lb = c*Q/(5.5*I0);                                           % RMS bunch length
halfL = DL1+Lb1;                                             % half length of the compressor, BC1=BC2
halfL_DL = DL3+Lb3;                                          % half length of the Delay Line chicane
k = @(lambda) 2*pi./lambda;                                  % uncompressed wave number in [1/m]
bk0 = @(lambda) sqrt(2*e*c./(I0*lambda));                    % initial bunching factor for density modulation, from shot noise
bm0 = @(lambda) 0;                                           % initial bunching factor for energy modulation, from shot noise
bf0 = @(lambda) [bk0(lambda) 0; bm0(lambda) 0];
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

R56_Spreader = R56_S1+R56_S2+R56_S3+R56_S4;                  % R56 of spreader transfer line in [m]
R56_DL = -2*theta3^2*((2/3)*Lb3+DL3);                        % R56 of Delay Line chicane in [m]

if R561 < 0
    h1 = abs((1-1/C1)/R561);                                 % linear energy chirp at BC1
elseif R561 == 0
    h1 = 0;
end
if R562 < 0
    h2 = abs((1-1/C2)/R562);                                 % linear energy chirp at BC2
elseif R562 == 0
    h2 = 0;
end
R1 = Lb1/theta1;                                             % BC1 dipoles bending radius in [m]
R2 = Lb2/theta2;                                             % BC2 dipoles bending radius in [m]
RS = Lbs/thetas;                                             % Spreader dipoles bending radius in [m]
R3 = Lb3/theta3;                                             % Delay Line chicane bending radius in [m]

H11 = (theta1^2*beta1_in);                                   % horizontal H-function at the exit of first dipole of BC1-DS1, in [m]
H12 = (theta1^2*beta1_w+2*eta1_max^2/beta1_w);               % horizontal H-function at the exit of second+third dipole of BC1-DS1, in [m]
H21 = (theta2^2*beta2_in);                                   % horizontal H-function at the exit of first dipole of BC2-DS2, in [m]
H22 = (theta2^2*beta2_w+2*eta2_max^2/beta2_w);               % horizontal H-function at the exit of second+third dipole of BC2-DS2, in [m]

D = re^2*(Q/e)/(8*lb*enx);                                   % constant factor in Bane's approximation                                 

options = optimset('MaxFunEvals',1000,'TolFun',1e-15,'TolX',1e-15);


%%%%%%%%
%% D0 %%
%%%%%%%%

%-------------------------------------------------------------------------%
%-------------------------------% IBS %-----------------------------------%
%-------------------------------------------------------------------------%

% Lorentz factor for mean energy along D0
gamma0 = E0/me;  

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

%-------------------------------------------------------------------------%
%-------------------------------% LSC %-----------------------------------%
%-------------------------------------------------------------------------%

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
Z0_LSC_i = @(lambda) Z0_LSC(lambda)*L0;
Dgamma0 = @(lambda) -(4*pi/Z0)*(I0/IA)*bk0(lambda).*Z0_LSC(lambda)*L0;

% S-matrix
S_L0 = @(lambda) [1 0;-Z0_LSC_i(lambda)*I0/Ef0 1];


%%%%%%%%%%%%
%% LINAC1 %%
%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%-------------------------------% IBS %-----------------------------------%
%-------------------------------------------------------------------------%

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

%-------------------------------------------------------------------------%
%-------------------------------% LSC %-----------------------------------%
%-------------------------------------------------------------------------%

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

Z1_LSC_real = @(s1,lambda) real(Z1_LSC(s1,lambda));
Z1_LSC_imag = @(s1,lambda) imag(Z1_LSC(s1,lambda));

Dgamma1_real = @(lambda) -(4*pi/Z0)*(I0/IA)*bk0(lambda)*integral(@(s1)Z1_LSC_real(s1,lambda),0,L1);
Dgamma1_imag = @(lambda) -(4*pi/Z0)*(I0/IA)*bk0(lambda)*integral(@(s1)Z1_LSC_imag(s1,lambda),0,L1);
Dgamma1 = @(lambda) Dgamma1_real(lambda)+1i*Dgamma1_imag(lambda);

% Z_LSC impedance integrated through L1
Z1_LSC_i = @(lambda) integral(@(s1) Z1_LSC(s1,lambda),0,L1,'ArrayValued',true);

% S-matrix of Linac1
S_L1 = @(lambda) [1 0; -Z1_LSC_i(lambda)*I0/Ef1 Ef0/Ef1];


%%%%%%%%%
%% BC1 %%
%%%%%%%%%

%-------------------------------------------------------------------------%
%--------------------------------% LH %-----------------------------------%
%-------------------------------------------------------------------------%

% Calculate the hypergeometric function in its integral form. The integral
% runs over the ratio of laser/e-beam radius: assume here from 1 to 10) for
% the definite integral.
A1 = @(lambda) abs((2*pi*C1./lambda)*R561*DeltaE_lh/Ef1);
J01_LH = @(r,lambda) besselj(0,A1(lambda).*exp(-r.^2/(4*Blh^2)));
J11_LH = @(r,lambda) besselj(1,A1(lambda).*exp(-r.^2/(4*Blh^2)));
S01_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*J01_LH(r,lambda),0,Inf);
S11_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*exp(-r.^2/(4*Blh^2)).*J11_LH(r,lambda),0,Inf,'ArrayValued',true);

if sigdE_lh ~= 0 && switch_lh == 1
    S01_LH = @(lambda) S01_LH(lambda);
    S11_LH = @(lambda) S11_LH(lambda);
elseif sigdE_lh == 0 || switch_lh == 0
    S01_LH = @(lambda) 1;
    S11_LH = @(lambda) 1;
end

sigd01 = sqrt((sigd0*Ef0)^2+sigEL0_ibs^2+sigEL1_ibs^2)/Ef1;                                     % fractional energy spread rms, normalized at the BC1 energy
I1 = I0*C1;                                                                                     % peak current compressed by C1 in BC1, in [A]
ex1 = enx/gammaf1;                                                                              % geometric horizontal emittance at chicane in [m rad]
lambda_co_long_1 = 2*pi*abs(R561)*C1.*sqrt(sigdE^2+sigdE_lh^2+sigEL0_ibs^2+sigEL1_ibs^2)/Ef1;   % energy Landau damping cutoff wavelength in [m]
lambda_co_tran_1 = 2*pi*sqrt(ex1*H11);                                                          % energy Landau damping cutoff wavelength in [m]
k1 = @(lambda) k(lambda)*C1;                                                                    % compressed modulation wave number

%-------------------------------------------------------------------------%
%-------------------------------% IBS %-----------------------------------%
%-------------------------------------------------------------------------%

% Derived
sigd_BC1_i = sqrt(sigd_L1^2*Ef0^2 + sigEL1_ibs^2)/Ef1;      % uncorrelated energy spread at the entrance of BC1
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
    sigd_BC1_w = sqrt(sigd_BC1_i^2);
    
end

% IBS-induced energy spread in BC1, in [MeV] 
sigEBC1_i = sigdBC1_i_ibs*Ef1;
sigEBC1_w = sigdBC1_w_ibs*Ef1;

%-------------------------------------------------------------------------%
%-------------------------------% LSC %-----------------------------------%
%-------------------------------------------------------------------------%

% Average beam radius in [m] and related quantities
rbBC1_i = 0.8735*(sqrt(enx*beta1_in/gammaf1) + sqrt(eny*beta1_in/gammaf1));
sBC1_i = sqrt(enx*beta1_in/gammaf1) + sqrt(eny*beta1_in/gammaf1);
spBC1_i = sqrt(enx/(gammaf1*beta1_in)) + sqrt(eny/(gammaf1*beta1_in));
rbBC1_w = 0.8735*(sqrt(enx*beta1_w/gammaf1) + sqrt(eny*beta1_w/gammaf1));
sBC1_w = sqrt(enx*beta1_w/gammaf1) + sqrt(eny*beta1_w/gammaf1);
spBC1_w = sqrt(enx/(gammaf1*beta1_w)) + sqrt(eny/(gammaf1*beta1_w));
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
DgammaBC1_i = @(lambda) -(4*pi/Z0)*(I0/IA)*bk0(lambda).*Z_LSC_BC1_i(lambda)*halfL;
DgammaBC1_w = @(lambda) -(4*pi/Z0)*(I0/IA)*bk0(lambda).*Z_LSC_BC1_w(lambda)*halfL;
DgammaBC1 = @(lambda) DgammaBC1_i(lambda) + DgammaBC1_w(lambda);


if switch_csr ~= 0 && theta1 ~= 0
   
    % CSR impedance of BC1 first half 
    
    % CSR in dipole-1
    ZCSR_BC1_1 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lb1*k(lambda).^(1/3)/R1^(2/3)).*exp(-(k(lambda)*theta1*sBC1_i).^2);     
    % CSR in dipole-2
    ZCSR_BC1_2_i = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(C1/(1+C1))*(Lb1*(((2*C1)/(1+C1))*k(lambda)).^(1/3)/R1^(2/3)).*exp(-(((2*C1)/(1+C1))*k(lambda)*DL1*theta1*spBC1_i).^2); 
    % CSR from dipole-2
    ZCSR_BC1_2_e = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(1/(1+C1))*(Lb1*((2/(1+C1))*(k(lambda))).^(1/3)/R1^(2/3)).*exp(-((2/(1+C1))*k(lambda)*D112*theta1*spBC1_i).^2);          

    % CER impedance of BC1 first half 
    
    % CER from dipole-1
    ZCER_BC1_1 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(DL1,lambda*gammaf1^2/(2*pi))/(R1^(2/3)*lambda.^(1/3)))).*exp(-(k(lambda)*theta1*sBC1_i).^2);      
    % CER from dipole-2
    ZCER_BC1_2_i = @(lambda) switch_cer*(Z0/(2*pi))*(C1/(1+C1))*log((min(D112,((1+C1)/(2*C1))*(lambda*gammaf1^2)/(2*pi))/(R1^(2/3)*(lambda*(1+C1)/(2*C1)).^(1/3)))).*exp(-(((2*C1)/(1+C1))*k(lambda)*DL1*theta1*spBC1_i).^2);   
    % CER from exit of dipole-2
    ZCER_BC1_2_e = @(lambda) switch_cer*(Z0/(2*pi))*(1/(1+C1))*log((min(D112,((1+C1)/2).*(lambda*gammaf1^2)/(2*pi))/((R1^(2/3).*(lambda*(1+C1)/2).^(1/3))))).*exp(-((2/(1+C1))*k(lambda)*D112*theta1*spBC1_w).^2);             

    % LSC impedance of BC1 first half
    ZBC1_LSC_fh = @(lambda) Z_LSC_BC1_i(lambda)*halfL;

    % Total impedance of BC1 first half
    ZBC1_fh = @(lambda) ZCSR_BC1_1(lambda) + ZCSR_BC1_2_i(lambda) + ZCSR_BC1_2_e(lambda) + ZCER_BC1_1(lambda) + ZCER_BC1_2_i(lambda) + ZCER_BC1_2_e(lambda) + ZBC1_LSC_fh(lambda);
    
   
    % CSR impedance of BC1 second half
    
    % CSR in dipole-3
    ZCSR_BC1_3 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lb1*k1(lambda).^(1/3)/R1^(2/3)).*exp(-(k1(lambda)*theta1*sBC1_w).^2);  
    % CSR in dipole-4
    ZCSR_BC1_4 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lb1*k1(lambda).^(1/3)/R1^(2/3));        
    
    % CER impedance of BC1 second half 
    
    % CER from dipole-3
    ZCER_BC1_3 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(D113,lambda*gammaf1^2/(2*pi))/(R1^(2/3)*lambda.^(1/3)))).*exp(-(k1(lambda)*theta1*sBC1_w).^2);     
    % CER from dipole-4
    ZCER_BC1_4 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(D114,lambda*gammaf1^2/(2*pi))/(R1^(2/3)*lambda.^(1/3)))); 
    
    % LSC impedance of BC1 second half
    ZBC1_LSC_sh = @(lambda) Z_LSC_BC1_w(lambda)*halfL;

    % Total impedance of BC1 second half
    ZBC1_sh = @(lambda) ZCSR_BC1_3(lambda) + ZCSR_BC1_4(lambda) + ZCER_BC1_3(lambda) + ZCER_BC1_4(lambda) + ZBC1_LSC_sh(lambda);
    
    
elseif switch_csr == 0 || theta1 == 0
    
    % LSC impedance of BC1 first half
    ZBC1_fh = @(lambda) Z_LSC_BC1_i(lambda)*halfL;
    
    % LSC impedance of BC1 second half
    ZBC1_sh = @(lambda) Z_LSC_BC1_w(lambda)*halfL;
       
end


if theta1 ~= 0

 % Amplification & damping functions

    F1_ld = @(lambda) exp(-0.5*(C1*R561*k(lambda)*sigd_BC1_w).^2).*S01_LH(lambda);
    G1_ld = @(lambda) exp(-0.5*(C1*R561*k(lambda)*sigd_BC1_w).^2).*S11_LH(lambda)*Ef1*(C1*R561*k(lambda)*sigd_BC1_w^2);
    
    % Compression matrix
    S_BC1_fh = @(lambda) [1 0; -ZBC1_fh(lambda)*I0/Ef1 1];
    S_BC1_sh = @(lambda) [1 0; -ZBC1_sh(lambda)*I0/Ef1 1];
    S_BC1_c = @(lambda) [F1_ld(lambda) 1i*F1_ld(lambda)*C1*R561.*k(lambda);...
                          1i*G1_ld(lambda)*C1/Ef1 (C1*F1_ld(lambda)-C1^2*G1_ld(lambda)*R561.*k(lambda)/Ef1)]; 
    
elseif theta1 == 0

    % Compression matrix
    S_BC1_fh = @(lambda) [1 0; -ZBC1_fh(lambda)*I0/Ef1 1];
    S_BC1_sh = @(lambda) [1 0; -ZBC1_sh(lambda)*I0/Ef1 1];
    S_BC1_c = @(lambda) [1 0; 0 1];
    
end

% S-matrix of BC1
S_BC1 = @(lambda) S_BC1_c(lambda)*S_BC1_sh(lambda)*S_BC1_fh(lambda);


%%%%%%%%%%%%
%% LINAC2 %%
%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%-------------------------------% IBS %-----------------------------------%
%-------------------------------------------------------------------------%

% Derived
Ei2 = Ef1;                                     % initial mean energy in [MeV]
G2 = (Ef2-Ei2)/L2;                             % L1 average accelerating gradient in [MeV/m]
gamma2 = @(s2) (Ei2+G2*s2)/me;                 % Lorentz factor for mean energy ALONG L2
gammaf2 = Ef2/me;                              % Lorentz factor for mean energy at the END of L2
gammam2 = (gammaf2 - Ei2/me)/2;                % Mean Lorentz factor for L2

sigd_L2 = sqrt(C1^2*sigd_BC1_w^2 + sigdBC1_w_ibs^2);    % uncorrelated relative enrgy spread at the entrance of Linac2

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

%-------------------------------------------------------------------------%
%-------------------------------% LSC %-----------------------------------%
%-------------------------------------------------------------------------%

% Average beam radius in [m] and related quantities
rb2 = @(s2) 0.8735*(sqrt(enx*betax2./gamma2(s2))+sqrt(eny*betay2./gamma2(s2)));
a2 = @(s2,lambda) k1(lambda).*rb2(s2)./gamma2(s2);
I21 = @(s2,lambda) besseli(1,a2(s2,lambda));
K20 = @(s2,lambda) besselk(0,a2(s2,lambda));
K21 = @(s2,lambda) besselk(1,a2(s2,lambda));
aw2 =  @(s2,lambda) k1(lambda)*rw2./gamma2(s2);
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

Dgamma2_real = @(lambda) -(4*pi/Z0)*(I0/IA)*b1(lambda)*integral(@(s2)Z2_LSC_real(s2,lambda),0,L2);
Dgamma2_imag = @(lambda) -(4*pi/Z0)*(I0/IA)*b1(lambda)*integral(@(s2)Z2_LSC_imag(s2,lambda),0,L2);
Dgamma2 = @(lambda) Dgamma2_real(lambda)+1i*Dgamma2_imag(lambda);

% Z_LSC impedance integrated through Linac2
Z2_LSC_i = @(lambda) integral(@(s2)Z2_LSC(s2,lambda),0,L2,'ArrayValued',true);

% S-matrix of Linac2
S_L2 = @(lambda) [1 0; -Z2_LSC_i(lambda)*I1/Ef2 Ef1/Ef2];   %%%%%% test


%%%%%%%%%
%% BC2 %%
%%%%%%%%%

%-------------------------------------------------------------------------%
%--------------------------------% LH %-----------------------------------%
%-------------------------------------------------------------------------%

% Calculate the hypergeometric function in its integral form. The integral
% runs over the ratio of laser/e-beam radius: assume here from 1 to 10) for
% the definite integral.
A2 = @(lambda) abs((2*pi*C1*C2./lambda)*R562*C1*DeltaE_lh/Ef2);
J02_LH = @(r,lambda) besselj(0,A2(lambda)*exp(-r.^2/(4*Blh^2)));
J12_LH = @(r,lambda) besselj(1,A2(lambda)*exp(-r.^2/(4*Blh^2)));
S02_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*J02_LH(r,lambda),0,Inf,'ArrayValued',true);
S12_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*exp(-r.^2/(4*Blh^2)).*J12_LH(r,lambda),0,Inf,'ArrayValued',true);

if switch_lh == 1 && theta2 ~= 0 && sigdE_lh ~= 0
    S02_LH = @(lambda) S02_LH(lambda);
    S12_LH = @(lambda) S12_LH(lambda);
elseif switch_lh == 0 || theta2 == 0 || sigdE_lh == 0
    S02_LH = @(lambda) 1;
    S12_LH = @(lambda) 1;
end

sigd02 = sqrt((sigd01*Ef1)^2+(sigdBC1_i_ibs^2+sigdBC1_w_ibs^2/C1^2)*Ef1^2+(sigEL2_ibs/C1)^2)/Ef2;       % fractional energy spread rms, normalized at the BC2 energy, uncompressed
sigd12 = C1*sigd02;                                                                                     % fractional energy spread rms, normalized at the DS2 energy, compressed by C1
I2 = C1*C2*I0;                                                                                          % final peak current in [A]
ex2 = enx/gammaf2;                                                                                      % geometric horizontal emittance at chicane in [m rad]
k2 = @(lambda) k1(lambda)*C2;                                                                           % wave number compressed by C1+C2 in BC1+BC2, in [1/m]

 % Energy Landau damping cutoff wavelength in [m]
lambda_co_long_2 = 2*pi*abs(R562)*C2.*2*pi*abs(R562)*C2*C1.*sqrt(sigdE^2+sigdE_lh^2+sigEL0_ibs^2+sigEL1_ibs^2+(sigdBC1_i_ibs^2+sigdBC1_w_ibs^2/C1^2)*Ef1^2+(sigEL2_ibs/C1)^2)/Ef2;  
% Energy Landau damping cutoff wavelength in [m]
lambda_co_tran_2 = 2*pi*sqrt(ex2*H21);   

%-------------------------------------------------------------------------%
%-------------------------------% IBS %-----------------------------------%
%-------------------------------------------------------------------------%

sigd_BC2_i = sqrt((sigd_L2*Ef1)^2+(sigEL2_ibs)^2)/Ef2;     % Uncorrelated energy spread at the entrance of BC2
qmax_BC2_i = sqrt(C1*halfL*(Q/e)*re^2/(2*gammaf2^(3/2)*enx^(3/2)*lb*sqrt(beta2_in)));
qmax_BC2_w = sqrt(C2*C1*halfL*(Q/e)*re^2/(2*gammaf2^(3/2)*enx^(3/2)*lb*sqrt(beta2_w)));

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
        solBC21 = @(y) vpa(abs( -ei(3*log(sqrt(hBC21*sigd_BC2_i^2 + 1)/bBC21)) + ei(3*log(sqrt(hBC21*y^2 + 1)/bBC21)) + vpa(halfL)*aBC21*hBC21/(2*bBC21^3)));
        yBC21 = vpa(fminsearch(solBC21,1e-6,options));
        sigdBC2_i_ibs = double(sqrt(yBC21^2 - sigd_BC2_i^2));

    % Uncorrelated relative energy spread at the waist of BC2
    sigd_BC2_w = sqrt(sigd_BC2_i^2 + sigdBC2_i_ibs^2);                      
    
    % IBS-induced fractional rms energy spread with Bane's formula for the second half of BC2
    dispH2_w = beta2_w*etap2^2+2*alpha2_w*eta2*etap2+gamma2_w*eta2^2;
        aBC22 = vpa(2*C1*C2*D/(gammaf2^(3/2)*sqrt(enx*beta2_w)));
        bBC22 = vpa(qmax_BC2_w*enx/(2*sqrt(2)*re));
        hBC22 = vpa(dispH1_w*gammaf2/enx);
        solBC22 = @(y) vpa(abs( -ei(3*log(sqrt(hBC22*sigd_BC2_w^2 + 1)/bBC22)) + ei(3*log(sqrt(hBC22*y^2 + 1)/bBC22)) + vpa(halfL)*aBC22*hBC22/(2*bBC22^3)));
        yBC22 = vpa(fminsearch(solBC22,1e-6,options));
        sigdBC2_w_ibs = double(sqrt(yBC22^2 - sigd_BC2_w^2));

    elseif theta2 == 0        
        
        sigdBC2_i_ibs = sqrt(2*halfL*C1*D*log(qmax_BC2_i*enx/(re*2*sqrt(2)))./(gammaf2.^2*sqrt(enx*beta2_in./gammaf2)));
        sigd_BC2_w = sqrt(sigd_BC2_i^2 + sigdBC2_i_ibs^2);
        sigdBC2_w_ibs = sqrt(2*halfL*C1*C2*D*log(qmax_BC2_w*enx/(re*2*sqrt(2)))./(gammaf2.^2*sqrt(enx*beta2_w./gammaf2)));
    
    end
    
elseif switch_bane1 == 0
    
    sigdBC2_i_ibs = 0;
    sigdBC2_w_ibs = 0;
    sigd_BC2_w = sqrt(sigd_BC2_i^2);

end

% IBS-induced energy spread in BC2, in [MeV] 
sigEBC2_i = sigdBC2_i_ibs*Ef2;
sigEBC2_w = sigdBC2_w_ibs*Ef2;

%-------------------------------------------------------------------------%
%-------------------------------% LSC %-----------------------------------%
%-------------------------------------------------------------------------%

% Average beam radius in [m] and related quantities
rbBC2_i = 0.8735*(sqrt(enx*beta2_in/gammaf2) + sqrt(eny*beta2_in/gammaf2));
sBC2_i = sqrt(enx*beta2_in/gammaf2) + sqrt(eny*beta2_in/gammaf2);
spBC2_i = sqrt(enx/(gammaf2*beta2_in)) + sqrt(eny/(gammaf2*beta2_in));
rbBC2_w = 0.8735*(sqrt(enx*beta2_w/gammaf2) + sqrt(eny*beta2_w/gammaf2));
sBC2_w = sqrt(enx*beta2_w/gammaf2) + sqrt(eny*beta2_w/gammaf2);
spBC2_w = sqrt(enx/(gammaf2*beta2_w)) + sqrt(eny/(gammaf2*beta2_w));
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
DgammaBC2_i = @(lambda) -(4*pi/Z0)*(I0/IA)*bk0(lambda).*Z_LSC_BC2_i(lambda)*halfL;
DgammaBC2_w = @(lambda) -(4*pi/Z0)*(I0/IA)*bk0(lambda).*Z_LSC_BC2_w(lambda)*halfL;
DgammaBC2 = @(lambda) DgammaBC2_i(lambda) + DgammaBC2_w(lambda);

if switch_csr ~= 0 && theta2 ~= 0

    % CSR impedance of BC2 first half 

    % CSR in dipole-1
    ZCSR_BC2_1 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lb2*k1(lambda).^(1/3)/R2^(2/3)).*exp(-(k1(lambda)*theta2*sBC2_i).^2);
    % CSR in dipole-2
    ZCSR_BC2_2 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lb2*k1(lambda).^(1/3)/R2^(2/3));

    % CSR in dipole-2
    ZCSR_BC2_2_i = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(C2/(1+C2))*(Lb2*(((2*C2)/(1+C2))*k1(lambda)).^(1/3)/R2^(2/3)).*exp(-(((2*C2)/(1+C2))*k1(lambda)*DL2*theta2*spBC2_i).^2);  
    % CSR from dipole-2
    ZCSR_BC2_2_e = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(1/(1+C2))*(Lb2*((2/(1+C2))*(k1(lambda))).^(1/3)/R2^(2/3)).*exp(-((2/(1+C2))*k1(lambda)*D212*theta2*spBC2_i).^2);         

    % CER impedance of BC2 first half 

    % CER from dipole-1
    ZCER_BC2_1 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(DL2,lambda*gammaf2^2/(2*pi))/(R2^(2/3)*lambda.^(1/3)))).*exp(-(k1(lambda)*theta2*sBC2_i).^2);      
    % CER from dipole-2
    ZCER_BC2_2_i = @(lambda) switch_cer*(Z0/(2*pi))*(C2/(1+C2))*log((min(D212,((1+C2)/(2*C2))*(lambda*gammaf2^2)/(2*pi))/(R2^(2/3)*(lambda*(1+C2)/(2*C2)).^(1/3)))).*exp(-(((2*C2)/(1+C2))*k1(lambda)*DL2*theta2*spBC2_i).^2);  
    % CER from exit of dipole-2
    ZCER_BC2_2_e = @(lambda) switch_cer*(Z0/(2*pi))*(1/(1+C2))*log((min(D212,((1+C2)/2).*(lambda*gammaf2^2)/(2*pi))/((R2^(2/3).*(lambda*(1+C2)/2).^(1/3))))).*exp(-((2/(1+C2))*k1(lambda)*D212*theta2*spBC2_w).^2);             

    % LSC impedance of BC2 first half
    ZBC2_LSC_fh = @(lambda) Z_LSC_BC2_i(lambda)*halfL;

    % Total impedance of BC2 first half
    ZBC2_fh = @(lambda) ZCSR_BC2_1(lambda) + ZCSR_BC2_2_i(lambda) + ZCSR_BC2_2_e(lambda) + ZCER_BC2_1(lambda) + ZCER_BC2_2_i(lambda) + ZCER_BC2_2_e(lambda) + ZBC2_LSC_fh(lambda);
    
   
    % CSR impedance of BC2 second half 

    % CSR in dipole-3
    ZCSR_BC2_3 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lb1*k2(lambda).^(1/3)/R2^(2/3)).*exp(-(k2(lambda)*theta2*sBC2_w).^2);      
    % CSR in dipole-4
    ZCSR_BC2_4 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lb1*k2(lambda).^(1/3)/R2^(2/3));        

    % CER impedance of BC2 second half 

    % CER from dipole-3
    ZCER_BC2_3 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(D213,lambda*gammaf2^2/(2*pi))/(R2^(2/3)*lambda.^(1/3)))).*exp(-(k2(lambda)*theta2*sBC2_w).^2);     
    % CER from dipole-4
    ZCER_BC2_4 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(D214,lambda*gammaf2^2/(2*pi))/(R2^(2/3)*lambda.^(1/3)))); 

    % LSC impedance of BC2 second half
    ZBC2_LSC_sh = @(lambda) Z_LSC_BC2_w(lambda)*halfL;

    % Total impedance of BC2 second half
    ZBC2_sh = @(lambda) ZCSR_BC2_3(lambda) + ZCSR_BC2_4(lambda) + ZCER_BC2_3(lambda) + ZCER_BC2_4(lambda) + ZBC2_LSC_sh(lambda);
    
elseif switch_csr == 0 || theta2 == 0

    % LSC impedance of BC2 first half
    ZBC2_fh = @(lambda) Z_LSC_BC2_i(lambda)*halfL;
    
    % LSC impedance of BC2 second half
    ZBC2_sh = @(lambda) Z_LSC_BC2_w(lambda)*halfL;
    
end


if theta2 ~= 0

    % Amplification & damping functions
    F2_ld = @(lambda) exp(-0.5*(C1*C2*R562*k(lambda)*sigd_BC2_w).^2).*S02_LH(lambda);
    G2_ld = @(lambda) exp(-0.5*(C1*C2*R562*k(lambda)*sigd_BC2_w).^2).*S12_LH(lambda)*Ef2.*(C2*R562*k(lambda)*sigd_BC2_w^2);
    F3_bc2_ld = @(lambda) exp(-0.5*(k(lambda)*sigd_BC2_w*(R561*Ef2/Ef1+C1*C2*R562)).^2).*S02_LH(lambda);
    G3_bc2_ld = @(lambda) exp(-0.5*(k(lambda)*sigd_BC2_w*(R561*Ef2/Ef1+C1*C2*R562)).^2).*S12_LH(lambda).*(k(lambda)*sigd_BC2_w^2*(R561*Ef2^2/(C1*Ef1)+C2*R562*Ef2));

    % Compression matrix
    S_BC2_fh = @(lambda) [1 0; -ZBC2_fh(lambda)*I1/Ef2 1];    %%%% test
    S_BC2_sh = @(lambda) [1 0; -ZBC2_sh(lambda)*I1/Ef2 1];    %%%% test

    if theta1 ~= 0
        % Higher orders in C1,C2 to take into account upstream microbunch
        % rotation in BC1
        S_BC2_c = @(lambda) [F3_bc2_ld(lambda) 1i*F3_bc2_ld(lambda)*R561.*k(lambda)+1i*F3_bc2_ld(lambda)*C1*C2*R562.*k(lambda)*Ef1/Ef2;...
            1i*G3_bc2_ld(lambda)*C1*C2/Ef2 F3_bc2_ld(lambda)*C2*Ef1/Ef2-G3_bc2_ld(lambda).*k(lambda)*C1*C2*R561/Ef2-G3_bc2_ld(lambda).*k(lambda)*C1^2*C2^2*R562*Ef1/Ef2^2];

    elseif theta1 == 0
        % Single-stage matrix
        S_BC2_c = @(lambda) [F2_ld(lambda) 1i*F2_ld(lambda)*C2*R562.*k1(lambda);...
            1i*G2_ld(lambda)*C2/Ef2 (C2*F2_ld(lambda)-C2^2*G2_ld(lambda)*R562.*k1(lambda)/Ef2)];
    end

elseif theta2 == 0

    S_BC2_fh = @(lambda) [1 0; -ZBC2_fh(lambda)*I1/Ef2 1];    %%%% test
    S_BC2_sh = @(lambda) [1 0; -ZBC2_sh(lambda)*I1/Ef2 1];    %%%% test
    S_BC2_c = @(lambda) [1 0; 0 1];

end

% S-matrix of BC2
S_BC2 = @(lambda) S_BC2_c(lambda)*S_BC2_sh(lambda)*S_BC2_fh(lambda);


%%%%%%%%%%%%
%% LINAC3 %%
%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%-------------------------------% IBS %-----------------------------------%
%-------------------------------------------------------------------------%

% Derived
Ei3 = Ef2;                        % initial mean energy in [MeV]
G3 = (Ef3-Ei3)/L3;                % L3 average accelerating gradient in [MeV/m]
gamma3 = @(s3) (Ei3+G3*s3)/me;    % Lorentz factor for mean energy ALONG L3
gammaf3 = Ef3/me;                 % Lorentz factor for mean energy at the END of L3
gammam3 = (gammaf3 - Ei3/me)/2;   % Mean Lorentz factor for L3

sigd_L3 = sqrt(C2^2*sigd_BC2_w^2 + sigdBC2_w_ibs^2);      % Uncorrelated relative energy spread at the entrance of L3

% Max angle
qmax_L3 = sqrt(C2*C1*L2*(Q/e)*re^2/(2*enx^(3/2)*lb*sqrt(betax3)));

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


%-------------------------------% LSC %-----------------------------------%

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

Dgamma3_real = @(lambda) -(4*pi/Z0)*(I0/IA)*b2(lambda)*integral(@(s3)Z3_LSC_real(s3,lambda),0,L3);
Dgamma3_imag = @(lambda) -(4*pi/Z0)*(I0/IA)*b2(lambda)*integral(@(s3)Z3_LSC_imag(s3,lambda),0,L3);
Dgamma3 = @(lambda) Dgamma3_real(lambda)+1i*Dgamma3_imag(lambda);

% Z_LSC impedance integrated through Linac3
Z3_LSC_i = @(lambda) integral(@(s3)Z3_LSC(s3,lambda),0,L3,'ArrayValued',true);

% S-matric of Linac3
%S_L3 = @(lambda) [1 0; -Z3_LSC_i(lambda)*I2/Ef3 Ef2/Ef3];   %%%% test
S_L3 = @(lambda) [1 0; -Z3_LSC_i(lambda)*I2/(C1*C2*Ef3) Ef2/Ef3];   %%%% test


%%%%%%%%%%%%%%
%% Spreader %%
%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%--------------------------------% LH %-----------------------------------%
%-------------------------------------------------------------------------%

% The LH hpergeometric function is defined for, and applied to, each dipoole
% magnet individually
A31 = @(lambda) abs((2*pi*C1*C2./lambda)*R56_S1*DeltaE_lh/Ef3);
J031_LH = @(r,lambda) besselj(0,A31(lambda).*exp(-r.^2/(4*Blh^2)));
J131_LH = @(r,lambda) besselj(1,A31(lambda).*exp(-r.^2/(4*Blh^2)));
S031_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*J031_LH(r,lambda),0,Inf);
S131_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*exp(-r.^2/(4*Blh^2)).*J131_LH(r,lambda),0,Inf,'ArrayValued',true);

A32 = @(lambda) abs((2*pi*C1*C2./lambda)*R56_S2*DeltaE_lh/Ef3);
J032_LH = @(r,lambda) besselj(0,A32(lambda).*exp(-r.^2/(4*Blh^2)));
J132_LH = @(r,lambda) besselj(1,A32(lambda).*exp(-r.^2/(4*Blh^2)));
S032_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*J032_LH(r,lambda),0,Inf);
S132_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*exp(-r.^2/(4*Blh^2)).*J132_LH(r,lambda),0,Inf,'ArrayValued',true);

A33 = @(lambda) abs((2*pi*C1*C2./lambda)*R56_S3*DeltaE_lh/Ef3);
J033_LH = @(r,lambda) besselj(0,A33(lambda).*exp(-r.^2/(4*Blh^2)));
J133_LH = @(r,lambda) besselj(1,A33(lambda).*exp(-r.^2/(4*Blh^2)));
S033_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*J033_LH(r,lambda),0,Inf);
S133_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*exp(-r.^2/(4*Blh^2)).*J133_LH(r,lambda),0,Inf,'ArrayValued',true);

A34 = @(lambda) abs((2*pi*C1*C2./lambda)*R56_S4*DeltaE_lh/Ef3);
J034_LH = @(r,lambda) besselj(0,A34(lambda).*exp(-r.^2/(4*Blh^2)));
J134_LH = @(r,lambda) besselj(1,A34(lambda).*exp(-r.^2/(4*Blh^2)));
S034_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*J034_LH(r,lambda),0,Inf);
S134_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*exp(-r.^2/(4*Blh^2)).*J134_LH(r,lambda),0,Inf,'ArrayValued',true);

if sigdE_lh ~= 0 && switch_lh == 1
    S031_LH = @(lambda) S031_LH(lambda);
    S032_LH = @(lambda) S032_LH(lambda);
    S033_LH = @(lambda) S033_LH(lambda);
    S034_LH = @(lambda) S034_LH(lambda);
    S131_LH = @(lambda) S131_LH(lambda);
    S132_LH = @(lambda) S132_LH(lambda);
    S133_LH = @(lambda) S133_LH(lambda);
    S134_LH = @(lambda) S134_LH(lambda);
elseif sigdE_lh == 0 || switch_lh == 0
    S031_LH = @(lambda) 1;
    S032_LH = @(lambda) 1;
    S033_LH = @(lambda) 1;
    S034_LH = @(lambda) 1;
    S131_LH = @(lambda) 1;
    S132_LH = @(lambda) 1;
    S133_LH = @(lambda) 1;
    S134_LH = @(lambda) 1;
end

%-------------------------------------------------------------------------%
%-------------------------------% IBS %-----------------------------------%
%-------------------------------------------------------------------------%

% Relative uncorrelated energy spread at the Spreader entrance
sigd_S = sqrt((sigd_L3^2*Ef2^2 + sigEL3_ibs^2))/Ef3;   

% Max angle
qmax_S0 = sqrt(LS0*C1*C2*(Q/e)*re^2/(2*gammaf3^(3/2)*enx^(3/2)*lb*sqrt(betaxdr0)));

% IBS-induced energy spread with Bane's formula
if switch_bane1 == 1
     sigdS0_ibs = sqrt(2*C1*C2*LS0*D*log(qmax_S0*enx/(re*2*sqrt(2)))/(gammaf3^2*sqrt(enx*betaxdr0./gammaf3)));
     sigES0_ibs = sigdS0_ibs*Ef3;
elseif switch_bane1 == 0
    sigdS0_ibs = 0;
end

% Relative uncorrelated after the first dipole
sigd_Sd1 = sqrt(sigd_S^2 + sigdS0_ibs^2);

% Max angle
qmax_S1 = sqrt(LS1*C1*C2*(Q/e)*re^2/(2*gammaf3^(3/2)*enx^(3/2)*lb*sqrt(betaxdr1)));

% IBS-induced energy spread with Bane's formula
if switch_bane1 == 1
     sigdS1_ibs = sqrt(2*C1*C2*LS1*D*log(qmax_S1*enx/(re*2*sqrt(2)))/(gammaf3^2*sqrt(enx*betaxdr1./gammaf3)));
     sigES1_ibs = sigdS1_ibs*Ef3;
elseif switch_bane1 == 0
    sigdS1_ibs = 0;
end

% Relative uncorrelated after the second dipole
sigd_Sd2 = sqrt(sigd_Sd1^2 + sigdS1_ibs^2);

% Max angle
qmax_S2 = sqrt(LS2*C1*C2*(Q/e)*re^2/(2*gammaf3^(3/2)*enx^(3/2)*lb*sqrt(betaxdr2)));

% IBS-induced energy spread with Bane's formula
if switch_bane1 == 1
     sigdS2_ibs = sqrt(2*C1*C2*LS2*D*log(qmax_S2*enx/(re*2*sqrt(2)))/(gammaf3^2*sqrt(enx*betaxdr2./gammaf3)));
     sigES2_ibs = sigdS2_ibs*Ef3;
elseif switch_bane1 == 0
    sigdS2_ibs = 0;
end

% Relative uncorrelated after the third dipole
sigd_Sd3 = sqrt(sigd_Sd2^2 + sigdS2_ibs^2);

% Max angle
qmax_S3 = sqrt(LS3*C1*C2*(Q/e)*re^2/(2*gammaf3^(3/2)*enx^(3/2)*lb*sqrt(betaxdr3)));

% IBS-induced energy spread with Bane's formula
if switch_bane1 == 1
     sigdS3_ibs = sqrt(2*C1*C2*LS3*D*log(qmax_S3*enx/(re*2*sqrt(2)))/(gammaf3^2*sqrt(enx*betaxdr3./gammaf3)));
     sigES3_ibs = sigdS3_ibs*Ef3;
elseif switch_bane1 == 0
    sigdS3_ibs = 0;
end

% Relative uncorrelated after the fourth dipole
sigd_Sd4 = sqrt(sigd_Sd3^2 + sigdS3_ibs^2);

% Max angle
qmax_S4 = sqrt(LS4*C1*C2*(Q/e)*re^2/(2*gammaf3^(3/2)*enx^(3/2)*lb*sqrt(betaxdr4)));

% IBS-induced energy spread with Bane's formula
if switch_bane1 == 1
     sigdS4_ibs = sqrt(2*C1*C2*LS4*D*log(qmax_S4*enx/(re*2*sqrt(2)))/(gammaf3^2*sqrt(enx*betaxdr4./gammaf3)));
     sigES34_ibs = sigdS4_ibs*Ef3;
elseif switch_bane1 == 0
    sigdS4_ibs = 0;
end

%-------------------------------------------------------------------------%
%-------------------------------% LSC %-----------------------------------%
%-------------------------------------------------------------------------%

% Average beam radius in [m] and related quantities
rbs0 = 0.8735*(sqrt(betaxdr0*enx/gammaf3) + sqrt(betaydr0*eny/gammaf3));    % effective radius along drifts of the Spreader, in [m]
rbs1 = 0.8735*(sqrt(betaxdr1*enx/gammaf3) + sqrt(betaydr1*eny/gammaf3));
rbs2 = 0.8735*(sqrt(betaxdr2*enx/gammaf3) + sqrt(betaydr2*eny/gammaf3));
rbs3 = 0.8735*(sqrt(betaxdr3*enx/gammaf3) + sqrt(betaydr3*eny/gammaf3));
rbs4 = 0.8735*(sqrt(betaxdr4*enx/gammaf3) + sqrt(betaydr4*eny/gammaf3));

sigS1 = sqrt(enx*betaxs1/gammaf3) + sqrt(eny*betays1/gammaf3);              % rms beam size at the Spreader dipoles, in [m] 
sigS2 = sqrt(enx*betaxs2/gammaf3) + sqrt(eny*betays2/gammaf3);
sigS3 = sqrt(enx*betaxs3/gammaf3) + sqrt(eny*betays3/gammaf3);
sigS4 = sqrt(enx*betaxs4/gammaf3) + sqrt(eny*betays4/gammaf3);

sigpS1 = sqrt(enx/(betaxs1*gammaf3)) + sqrt(eny/(betays1*gammaf3));         % rms beam angular divergence  at the Spreader dipoles, in [m] 
sigpS2 = sqrt(enx/(betaxs2*gammaf3)) + sqrt(eny/(betays2*gammaf3));
sigpS3 = sqrt(enx/(betaxs3*gammaf3)) + sqrt(eny/(betays3*gammaf3));
sigpS4 = sqrt(enx/(betaxs4*gammaf3)) + sqrt(eny/(betays4*gammaf3));

a0_sp = @(lambda) k2(lambda).*rbs0/gammaf3;
a1_sp = @(lambda) k2(lambda).*rbs1/gammaf3;
a2_sp = @(lambda) k2(lambda).*rbs2/gammaf3;
a3_sp = @(lambda) k2(lambda).*rbs3/gammaf3;
a4_sp = @(lambda) k2(lambda).*rbs4/gammaf3;

K01_SP = @(lambda) besselk(1,a0_sp(lambda));
I01_SP = @(lambda) besseli(1,a0_sp(lambda));
K11_SP = @(lambda) besselk(1,a1_sp(lambda));
I11_SP = @(lambda) besseli(1,a1_sp(lambda));
K21_SP = @(lambda) besselk(1,a2_sp(lambda));
I21_SP = @(lambda) besseli(1,a2_sp(lambda));
K31_SP = @(lambda) besselk(1,a3_sp(lambda));
I31_SP = @(lambda) besseli(1,a3_sp(lambda));
K41_SP = @(lambda) besselk(1,a4_sp(lambda));
I41_SP = @(lambda) besseli(1,a4_sp(lambda));

% CSR impedance of Spreader dipole magnets
ZCSR_S1 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lbs*k2(lambda).^(1/3)/RS^(2/3)).*exp(-switch_Hx*(k2(lambda)*thetas*sigS1).^2);              % CSR in dipole-1
ZCSR_S2 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lbs*k2(lambda).^(1/3)/RS^(2/3)).*exp(-switch_Hx*(k2(lambda)*Lbs*thetas*sigpS2).^2);         % CSR from dipole-2
ZCSR_S3 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lbs*k2(lambda).^(1/3)/RS^(2/3)).*exp(-switch_Hx*(k2(lambda)*thetas*sigS3).^2);              % CSR in dipole-3
ZCSR_S4 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lbs*k2(lambda).^(1/3)/RS^(2/3)).*exp(-switch_Hx*(k2(lambda)*Lbs*thetas*sigpS4).^2);         % CSR in dipole-4

% CER impedance of Spreader dipole magnets
ZCER_S1 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(LS1,lambda*gammaf3^2/(2*pi))/(RS^(2/3)*lambda.^(1/3)))).*exp(-switch_Hx*(k2(lambda)*thetas*sigS1).^2);          % CER from dipole-1
ZCER_S2 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(LS2,lambda*gammaf3^2/(2*pi))/(RS^(2/3)*lambda.^(1/3)))).*exp(-switch_Hx*(k2(lambda)*Lbs*thetas*sigpS2).^2);     % CER from dipole-2
ZCER_S3 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(LS3,lambda*gammaf3^2/(2*pi))/(RS^(2/3)*lambda.^(1/3)))).*exp(-switch_Hx*(k2(lambda)*thetas*sigS3).^2);          % CER from dipole-3
ZCER_S4 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(LS4,lambda*gammaf3^2/(2*pi))/(RS^(2/3)*lambda.^(1/3)))).*exp(-switch_Hx*(k2(lambda)*Lbs*thetas*sigpS4).^2);     % CER from dipole-4

% LSC impedances of the Spreader drift sections
ZLSC_S1 = @(lambda) (1i*Z0*LS0./(pi*k2(lambda)*rbs0^2)).*(1-2*I01_SP(lambda).*K01_SP(lambda)).*exp(-switch_Hx*(k2(lambda)*thetas*sigS1).^2);
ZLSC_S2 = @(lambda) (1i*Z0*LS1./(pi*k2(lambda)*rbs1^2)).*(1-2*I11_SP(lambda).*K11_SP(lambda)).*exp(-switch_Hx*(k2(lambda)*thetas*sigS2).^2);             % LSC in drift of Spreader, after dipole1
ZLSC_S3 = @(lambda) (1i*Z0*LS2./(pi*k2(lambda)*rbs2^2)).*(1-2*I21_SP(lambda).*K21_SP(lambda)).*exp(-switch_Hx*(k2(lambda)*thetas*sigS3).^2);             % LSC in drift of Spreader, after dipole2
ZLSC_S4 = @(lambda) (1i*Z0*LS3./(pi*k2(lambda)*rbs3^2)).*(1-2*I31_SP(lambda).*K31_SP(lambda)).*exp(-switch_Hx*(k2(lambda)*thetas*sigS4).^2);             % LSC in drift of Spreader, after dipole3
ZLSC_last = @(lambda) (1i*Z0*LS4./(pi*k2(lambda)*rbs4^2)).*(1-2*I41_SP(lambda).*K41_SP(lambda));                                                           % LSC in drift of Spreader, after dipole4

% Matrices of CSR impedance from spreader dipoles
S_Sdip1_csr = @(lambda) [1 0; -ZCSR_S1(lambda)*I2/Ef3 1];
S_Sdip2_csr = @(lambda) [1 0; -ZCSR_S2(lambda)*I2/Ef3 1];
S_Sdip3_csr = @(lambda) [1 0; -ZCSR_S3(lambda)*I2/Ef3 1];
S_Sdip4_csr = @(lambda) [1 0; -ZCSR_S4(lambda)*I2/Ef3 1];

% Matrices of CER impedance from spreader dipoles
S_Sdip1_cer = @(lambda) [1 0; -ZCER_S1(lambda)*I2/Ef3 1];
S_Sdip2_cer = @(lambda) [1 0; -ZCER_S2(lambda)*I2/Ef3 1];
S_Sdip3_cer = @(lambda) [1 0; -ZCER_S3(lambda)*I2/Ef3 1];
S_Sdip4_cer = @(lambda) [1 0; -ZCER_S4(lambda)*I2/Ef3 1];

% Matrices of LSC impedance for spreader drifts
S_Sd0 = @(lambda) [1 0; -ZLSC_S1(lambda)*I2/Ef3 1];
S_Sd1 = @(lambda) [1 0; -ZLSC_S2(lambda)*I2/Ef3 1];
S_Sd2 = @(lambda) [1 0; -ZLSC_S3(lambda)*I2/Ef3 1];
S_Sd3 = @(lambda) [1 0; -ZLSC_S4(lambda)*I2/Ef3 1];
S_Sd4 = @(lambda) [1 0; -ZLSC_last(lambda)*I2/Ef3 1];

% Growth suppression for a gaussian distribution
FS1 = @(lambda) exp(-0.5*(R55_S1^(-1)*R56_S1*k2(lambda)*sigd_S).^2).*S031_LH(lambda);
GS1 = @(lambda) exp(-0.5*(R55_S1^(-1)*R56_S1*k2(lambda)*sigd_S).^2).*S131_LH(lambda)*Ef3*(R55_S1^(-1)*R56_S1*k2(lambda)*sigd_S^2);
FS2 = @(lambda) exp(-0.5*(R55_S2^(-1)*R56_S2*k2(lambda)*sigd_Sd1).^2).*S032_LH(lambda);
GS2 = @(lambda) exp(-0.5*(R55_S2^(-1)*R56_S2*k2(lambda)*sigd_Sd1).^2).*S132_LH(lambda)*Ef3*(R55_S2^(-1)*R56_S2*k2(lambda)*sigd_Sd1^2);
FS3 = @(lambda) exp(-0.5*(R55_S3^(-1)*R56_S3*k2(lambda)*sigd_Sd2).^2).*S033_LH(lambda);
GS3 = @(lambda) exp(-0.5*(R55_S3^(-1)*R56_S3*k2(lambda)*sigd_Sd2).^2).*S133_LH(lambda)*Ef3*(R55_S3^(-1)*R56_S3*k2(lambda)*sigd_Sd2^2);
FS4 = @(lambda) exp(-0.5*(R55_S4^(-1)*R56_S4*k2(lambda)*sigd_Sd3).^2).*S034_LH(lambda);
GS4 = @(lambda) exp(-0.5*(R55_S4^(-1)*R56_S4*k2(lambda)*sigd_Sd3).^2).*S134_LH(lambda)*Ef3*(R55_S4^(-1)*R56_S4*k2(lambda)*sigd_Sd3^2);

S_Sdip1_c = @(lambda) [FS1(lambda) 1i*FS1(lambda)*R55_S1^(-1)*R56_S1.*k2(lambda);...
                          1i*GS1(lambda)*R55_S1^(-1)/Ef3 (R55_S1^(-1)*FS1(lambda)-R55_S1^(-2)*GS1(lambda)*R56_S1.*k2(lambda)/Ef3)];
S_Sdip2_c = @(lambda) [FS2(lambda) 1i*FS2(lambda)*R55_S2^(-1)*R56_S2.*k2(lambda);...
                          1i*GS2(lambda)*R55_S2^(-1)/Ef3 (R55_S2^(-1)*FS2(lambda)-R55_S2^(-2)*GS2(lambda)*R56_S2.*k2(lambda)/Ef3)];
S_Sdip3_c = @(lambda) [FS3(lambda) 1i*FS3(lambda)*R55_S3^(-1)*R56_S3.*k2(lambda);...
                          1i*GS3(lambda)*R55_S3^(-1)/Ef3 (R55_S3^(-1)*FS3(lambda)-R55_S3^(-2)*GS3(lambda)*R56_S3.*k2(lambda)/Ef3)];
S_Sdip4_c = @(lambda) [FS4(lambda) 1i*FS4(lambda)*R55_S4^(-1)*R56_S4.*k2(lambda);...
                          1i*GS4(lambda)*R55_S4^(-1)/Ef3 (R55_S4^(-1)*FS4(lambda)-R55_S4^(-2)*GS4(lambda)*R56_S4.*k2(lambda)/Ef3)];
                      
if switch_spreader == 0
    S_Sdip1 = @(lambda) [1 0;0 1];
    S_Sdip2 = @(lambda) [1 0;0 1];
    S_Sdip3 = @(lambda) [1 0;0 1];
    S_Sdip4 = @(lambda) [1 0;0 1];
elseif switch_spreader == 1 
    S_Sdip1 = @(lambda) S_Sdip1_c(lambda)*S_Sdip1_csr(lambda)*S_Sdip1_cer(lambda);
    S_Sdip2 = @(lambda) S_Sdip2_c(lambda)*S_Sdip2_csr(lambda)*S_Sdip2_cer(lambda);
    S_Sdip3 = @(lambda) S_Sdip3_c(lambda)*S_Sdip3_csr(lambda)*S_Sdip3_cer(lambda);
    S_Sdip4 = @(lambda) S_Sdip4_c(lambda)*S_Sdip4_csr(lambda)*S_Sdip4_cer(lambda);
end

% S-matrix of the whole Spreader line
S_spreader = @(lambda) S_Sd4(lambda)*S_Sdip4(lambda)*S_Sd3(lambda)*S_Sdip3(lambda)*S_Sd2(lambda)*S_Sdip2(lambda)*S_Sd1(lambda)*S_Sdip1(lambda)*S_Sd0(lambda);


%%%%%%%%%%%%%%%%%
%% Modulator 1 %%
%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%-------------------------------% IBS %-----------------------------------%
%-------------------------------------------------------------------------%

% Q max
qmax_Lmod1 = sqrt(Lmod1*C1*C2*(Q/e)*re^2/(2*gammaf3^(3/2)*enx^(3/2)*lb*sqrt(betaxm1)));

% IBS-induced rms relative energy spread with Bane's approx. of the B-M expression
if switch_bane1 == 1
    sigdMOD1_ibs = sqrt(2*Lmod1*C1*C2*D*log(qmax_Lmod1*enx/(re*2*sqrt(2)))/(gammaf3^2*sqrt(enx*betaxm1./gammaf3)));
elseif switch_bane1 == 0
    sigdMOD1_ibs = 0;
end

% IBS-induced absolute energy spread in [MeV] cumulated through L0
sigEMOD1_ibs = sigdMOD1_ibs*Ef3;

%-------------------------------------------------------------------------%
%-------------------------------% LSC %-----------------------------------%
%-------------------------------------------------------------------------%

% Average beam radius in [m] and related quantities
rbm1 = 0.8735*(sqrt(enx*betaxm1/gammaf3)+sqrt(eny*betaym1/gammaf3));
am1 =  @(lambda) k2(lambda)*rbm1*sqrt(1 + Kmod1^2/2)/gammaf3;  
I1_m1 = @(lambda) besseli(1,am1(lambda));
K0_m1 = @(lambda) besselk(0,am1(lambda));
K1_m1 = @(lambda) besselk(1,am1(lambda));
awm1 =  @(lambda) k2(lambda)*rwm1*sqrt(1 + Kmod1^2/2)/gammaf3;
I0_m1_w = @(lambda) besseli(0,awm1(lambda));
K0_m1_w = @(lambda) besselk(0,awm1(lambda));

% 3-D LSC impedance averaged over transverse dimension (_rr stays for
% chamber shielding)
Z_MOD1_LSC = @(lambda) 1i*(Z0./(pi*k2(lambda)*rbm1^2)).*(1-2*I1_m1(lambda).*K1_m1(lambda));
Z_MOD1_LSC_rr = @(lambda) 1i*(Z0./(pi*k2(lambda)*rbm1^2))*(1-2*(I1_m1(lambda)/I0_m1_w(lambda))*(I0_m1_w(lambda)*K1_m1(lambda)+I1_m1(lambda)*K0_m1_w(lambda)));

% LSC-induced energy modulation amplitude in me unit
Z_MOD1_LSC_i = @(lambda) Z_MOD1_LSC(lambda)*Lmod1;
Dgamma_MOD1 = @(lambda) -(4*pi/Z0)*(I0/IA)*b3(lambda).*Z_MOD1_LSC(lambda)*Lmod1;

% S-matrix
S_MOD1 = @(lambda) [1 0;-Z_MOD1_LSC_i(lambda)*I2/Ef3 1];


%%%%%%%%%%%%%%%%%
%% Undulators  %%
%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%-------------------------------% IBS %-----------------------------------%
%-------------------------------------------------------------------------%

% Q max
qmax_U = sqrt(Lu*C1*C2*(Q/e)*re^2/(2*gammaf3^(3/2)*enx^(3/2)*lb*sqrt(betaxu)));

% IBS-induced rms relative energy spread with Bane's approx. of the B-M expression
if switch_bane1 == 1
    sigdU_ibs = sqrt(2*Lu*C1*C2*D*log(qmax_U*enx/(re*2*sqrt(2)))/(gammaf3^2*sqrt(enx*betaxu./gammaf3)));
elseif switch_bane1 == 0
    sigdU_ibs = 0;
end

% IBS-induced absolute energy spread in [MeV] cumulated through L0
sigEU_ibs = sigdU_ibs*Ef3;

%-------------------------------------------------------------------------%
%-------------------------------% LSC %-----------------------------------%
%-------------------------------------------------------------------------%

% Average beam radius in [m] and related quantities
rbu = 0.8735*(sqrt(enx*betaxu/gammaf3)+sqrt(eny*betayu/gammaf3));
au =  @(lambda) k2(lambda)*rbu/gammaf3;  
I1_u = @(lambda) besseli(1,au(lambda));
K0_u = @(lambda) besselk(0,au(lambda));
K1_u = @(lambda) besselk(1,au(lambda));
awu =  @(lambda) k2(lambda)*rwu/gammaf3;
I0_u_w = @(lambda) besseli(0,awu(lambda));
K0_u_w = @(lambda) besselk(0,awu(lambda));

% 3-D LSC impedance averaged over transverse dimension (_rr stays for
% chamber shielding)
Z_u_LSC = @(lambda) 1i*(Z0./(pi*k2(lambda)*rbu^2)).*(1-2*I1_u(lambda).*K1_u(lambda));
Z_u_LSC_rr = @(lambda) 1i*(Z0./(pi*k2(lambda)*rbu^2))*(1-2*(I1_u(lambda)/I0_u_w(lambda))*(I0_u_w(lambda)*K1_u(lambda)+I1_u(lambda)*K0_u_w(lambda)));

% LSC-induced energy modulation amplitude in me unit
Z_u_LSC_i = @(lambda) Z_u_LSC(lambda)*Lu;
Dgamma_u = @(lambda) -(4*pi/Z0)*(I0/IA)*b3(lambda).*Z_u_LSC(lambda)*Lu;

% S-matrix
S_U = @(lambda) [1 0;-Z_u_LSC_i(lambda)*I2/Ef3 1];


%%%%%%%%%%%%%%%%
%% Delay Line %%
%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%--------------------------------% LH %-----------------------------------%
%-------------------------------------------------------------------------%

A4 = @(lambda) abs((2*pi*C1*C2./lambda)*R56_DL*DeltaE_lh/Ef3);
J04_LH = @(r,lambda) besselj(0,A4(lambda).*exp(-r.^2/(4*Blh^2)));
J14_LH = @(r,lambda) besselj(1,A4(lambda).*exp(-r.^2/(4*Blh^2)));
S04_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*J04_LH(r,lambda),0,Inf);
S14_LH = @(lambda) integral(@(r)r.*exp(-r.^2/2).*exp(-r.^2/(4*Blh^2)).*J14_LH(r,lambda),0,Inf,'ArrayValued',true);

if sigdE_lh ~= 0 && switch_lh == 1
    S04_LH = @(lambda) S04_LH(lambda);
    S14_LH = @(lambda) S14_LH(lambda);
elseif sigdE_lh == 0 || switch_lh == 0
    S04_LH = @(lambda) 1;
    S14_LH = @(lambda) 1;
end

%-------------------------------------------------------------------------%
%----------------------------% Seed Laser 1 %-----------------------------%
%-------------------------------------------------------------------------%

B4 = @(lambda) R56_DL*k2(lambda)*DeltaE_sl1/Ef3;
S04_SL = @(r,lambda) besselj(0,B4(lambda));
S14_SL = @(r,lambda) besselj(1,B4(lambda));

if DeltaE_sl1 ~= 0
    S04_SL = @(lambda) S04_SL(lambda);
    S14_SL = @(lambda) S14_SL(lambda);
elseif DeltaE_sl1 == 0
    S04_SL = @(lambda) 1;
    S14_SL = @(lambda) 1;
end

%-------------------------------------------------------------------------%
%-------------------------------% IBS %-----------------------------------%
%-------------------------------------------------------------------------%

sigd_DL = sqrt(sigd_Sd4^2 + sigdS4_ibs^2 + sigdMOD1_ibs^2 + sigdU_ibs^2);   % relative uncorrelated energy spread at the entrance of Delay line

qmax_DL_i = sqrt(C2*C1*halfL_DL*(Q/e)*re^2/(2*gammaf3^(3/2)*enx^(3/2)*lb*sqrt(beta3_in)));
qmax_DL_w = sqrt(C2*C1*halfL_DL*(Q/e)*re^2/(2*gammaf3^(3/2)*enx^(3/2)*lb*sqrt(beta3_w)));

% optics functions in the first half (_i) and second half (_w) of EEHG chicane
etap3 = theta3;
eta3  = theta3*(DL3+Lb3);
gamma3_in = (1+alpha3_in^2)/beta3_in;
gamma3_w  = (1+alpha3_w^2)/beta3_w;

if switch_bane1 == 1
    
    if theta3 ~= 0
        
    % IBS-induced fractional rms energy spread with Bane's formula for the first half of EEHG chicane
    dispH3_i = beta3_in*etap3^2+2*alpha3_in*eta3*etap3+gamma3_in*eta3^2;  
    
        aHG1 = vpa(2*C1*C2*D/(gammaf3^(3/2)*sqrt(enx*beta3_in)));
        bHG1 = vpa(qmax_DL_i*enx/(2*sqrt(2)*re));
        hHG1 = vpa(dispH3_i*gammaf3/enx);
        solHG1 = @(y) vpa(abs( -ei(3*log(sqrt(hHG1*sigd_DL^2 + 1)/bHG1)) + ei(3*log(sqrt(hHG1*y^2 + 1)/bHG1)) + vpa(halfL_DL)*aHG1*hHG1/(2*bHG1^3)));
        yHG1 = vpa(fminsearch(solHG1,1e-6,options));
        sigdDL_i_ibs = double(sqrt(yHG1^2 - vpa(sigd_DL)^2));
    
    % Uncorrelated energy spread at the waist of EEHG
    sigd_DL_w = sqrt(sigd_DL^2 + sigdDL_i_ibs^2);       
    
     % IBS-induced fractional rms energy spread with Bane's formula for the second half of EEHG chicane
     dispH3_w = beta3_w*etap3^2+2*alpha3_w*eta3*etap3+gamma3_w*eta3^2;
     
        aHG2 = vpa(2*C1*C2*D/(gammaf3^(3/2)*sqrt(enx*beta3_w)));
        bHG2 = vpa(qmax_DL_w*enx/(2*sqrt(2)*re));
        hHG2 = vpa(dispH3_w*gammaf3/enx);
        solHG2 = @(y) vpa(abs( -ei(3*log(sqrt(hHG2*sigd_DL_w^2 + 1)/bHG2)) + ei(3*log(sqrt(hHG2*y^2 + 1)/bHG2)) + vpa(halfL_DL)*aHG2*hHG2/(2*bHG2^3)));
        yHG2 = vpa(fminsearch(solHG2,1e-6,options));
        sigdDL_w_ibs = double(sqrt(yHG2^2 - vpa(sigd_DL_w)^2));
    
    elseif theta3 == 0   
        
        sigdDL_i_ibs = sqrt(2*halfL_DL*C2*C1*D*log(qmax_DL_i*enx/(re*2*sqrt(2)))./(gammaf3.^2*sqrt(enx*beta3_in./gammaf3)));
        sigd_DL_w = sqrt(sigd_DL^2 + sigdDL_i_ibs^2);
        sigdDL_w_ibs = sqrt(2*halfL_DL*C1*C2*D*log(qmax_DL_w*enx/(re*2*sqrt(2)))./(gammaf3.^2*sqrt(enx*beta3_w./gammaf3)));
        
    end
    
elseif switch_bane1 == 0
    
    sigdDL_i_ibs = 0;
    sigdDL_w_ibs = 0;
    sigd_DL_w = sqrt(sigd_DL^2);
    
end

%-------------------------------------------------------------------------%
%-------------------------------% LSC %-----------------------------------%
%-------------------------------------------------------------------------%

sDL_i = sqrt(enx*beta3_in/gammaf3) + sqrt(eny*beta3_in/gammaf3);                                                 
spDL_i = sqrt(enx/(beta3_in*gammaf3)) + sqrt(eny/(beta3_in*gammaf3));
sDL_w = sqrt(enx*beta3_w/gammaf3) + sqrt(eny*beta3_w/gammaf3);
spDL_w = sqrt(enx/(beta3_w*gammaf3)) + sqrt(eny/(beta3_w*gammaf3));

% CSR impedance of EEHG first half 
ZCSR_DL_1 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lb3*k2(lambda).^(1/3)/R3^(2/3)).*exp(-(k2(lambda)*theta3*sDL_i).^2);                     % CSR in dipole-1
ZCSR_DL_2_i = @(lambda) switch_csr*(Z0/(8*pi))*A_csr*(Lb3*(k2(lambda)).^(1/3)/R3^(2/3)).*exp(-(k2(lambda)*DL3*theta3*spDL_i).^2);            % CSR in dipole-2
ZCSR_DL_2_e = @(lambda) switch_csr*(Z0/(8*pi))*A_csr*(Lb3*(k2(lambda)).^(1/3)/R3^(2/3)).*exp(-(k2(lambda)*D312*theta3*spDL_i).^2);           % CSR from dipole-2

% CER impedance of EEHG first half 
ZCER_DL_1 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(DL3,lambda*gammaf3^2/(2*pi))/(R3^(2/3)*lambda.^(1/3)))).*exp(-(k2(lambda)*theta3*sDL_i).^2);                  % CER from dipole-1
ZCER_DL_2_i = @(lambda) switch_cer*(Z0/(4*pi))*log((min(D312,(lambda*gammaf3^2)/(2*pi))/(R3^(2/3)*(lambda).^(1/3)))).*exp(-(k2(lambda)*DL3*theta3*spDL_i).^2);      % CER from dipole-2
ZCER_DL_2_e = @(lambda) switch_cer*(Z0/(4*pi))*log((min(D312,(lambda*gammaf3^2)/(2*pi))/((R3^(2/3).*(lambda).^(1/3))))).*exp(-(k2(lambda)*D312*theta3*spDL_w).^2);  % CER from exit of dipole-2

% Total impedance of EEHG first half
ZDL_fh = @(lambda) ZCSR_DL_1(lambda) + ZCSR_DL_2_i(lambda) + ZCSR_DL_2_e(lambda) + ZCER_DL_1(lambda) + ZCER_DL_2_i(lambda) + ZCER_DL_2_e(lambda);
S_DL_fh = @(lambda) [1 0; -ZDL_fh(lambda)*I2/Ef3 1];

% Amplification & damping functions
F3_ld = @(lambda) exp(-0.5*(R56_DL*k2(lambda)*sigd_DL_w).^2).*S04_LH(lambda).*S04_SL(lambda);
G3_ld = @(lambda) exp(-0.5*(R56_DL*k2(lambda)*sigd_DL_w).^2).*S14_LH(lambda).*S14_SL(lambda)*Ef3*(R56_DL*k2(lambda)*sigd_DL_w^2);

% Compression matrix
S_DL_c = @(lambda) [F3_ld(lambda) 1i*F3_ld(lambda)*R56_DL.*k2(lambda);...
                      1i*G3_ld(lambda)/Ef3 (F3_ld(lambda)-G3_ld(lambda)*R56_DL.*k2(lambda)/Ef3)]; 

% CSR impedance of EEHG second half 
ZCSR_DL_3 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lb3*k2(lambda).^(1/3)/R3^(2/3)).*exp(-(k2(lambda)*theta3*sDL_w).^2);     % CSR in dipole-3
ZCSR_DL_4 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lb3*k2(lambda).^(1/3)/R3^(2/3));                                         % CSR in dipole-4

% CER impedance of EEHG second half 
ZCER_DL_3 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(D313,lambda*gammaf3^2/(2*pi))/(R3^(2/3)*lambda.^(1/3)))).*exp(-(k2(lambda)*theta3*sDL_w).^2);     % CER from dipole-3
ZCER_DL_4 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(D314,lambda*gammaf3^2/(2*pi))/(R3^(2/3)*lambda.^(1/3))));                                         % CER from dipole-4

% Total impedance of EEHG second half
ZDL_sh = @(lambda) ZCSR_DL_3(lambda) + ZCSR_DL_4(lambda) + ZCER_DL_3(lambda) + ZCER_DL_4(lambda);
S_DL_sh = @(lambda) [1 0; -ZDL_sh(lambda)*I2/Ef3 1];

% S-matrix of the EEHG line (incuding big chicane)
S_DL = @(lambda) S_DL_sh(lambda)*S_DL_c(lambda)*S_DL_fh(lambda);


%%%%%%%%%%%%%%%%%
%% Modulator 2 %%
%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%-------------------------------% IBS %-----------------------------------%
%-------------------------------------------------------------------------%

% Q max
qmax_Lmod2 = sqrt(Lmod2*C1*C2*(Q/e)*re^2/(2*gammaf3^(3/2)*enx^(3/2)*lb*sqrt(betaxm2)));

% IBS-induced rms relative energy spread with Bane's approx. of the B-M expression
if switch_bane1 == 1
    sigdMOD2_ibs = sqrt(2*Lmod2*C1*C2*D*log(qmax_Lmod2*enx/(re*2*sqrt(2)))/(gammaf3^2*sqrt(enx*betaxm2./gammaf3)));
elseif switch_bane1 == 0
    sigdMOD2_ibs = 0;
end

% IBS-induced absolute energy spread in [MeV] cumulated through L0
sigEMOD2_ibs = sigdMOD2_ibs*Ef3;

%-------------------------------------------------------------------------%
%-------------------------------% LSC %-----------------------------------%
%-------------------------------------------------------------------------%

% Average beam radius in [m] and related quantities
rbm2 = 0.8735*(sqrt(enx*betaxm2/gammaf3)+sqrt(eny*betaym2/gammaf3));
am2 =  @(lambda) k2(lambda)*rbm2*sqrt(1 + Kmod2^2/2)/gammaf3;  
I1_m2 = @(lambda) besseli(1,am2(lambda));
K0_m2 = @(lambda) besselk(0,am2(lambda));
K1_m2 = @(lambda) besselk(1,am2(lambda));
awm2 =  @(lambda) k2(lambda)*rwm2*sqrt(1 + Kmod2^2/2)/gammaf3;
I0_m2_w = @(lambda) besseli(0,awm2(lambda));
K0_m2_w = @(lambda) besselk(0,awm2(lambda));

% 3-D LSC impedance averaged over transverse dimension (_rr stays for
% chamber shielding)
Z_MOD2_LSC = @(lambda) 1i*(Z0./(pi*k2(lambda)*rbm2^2)).*(1-2*I1_m2(lambda).*K1_m2(lambda));
Z_MOD2_LSC_rr = @(lambda) 1i*(Z0./(pi*k2(lambda)*rbm2^2))*(1-2*(I1_m2(lambda)/I0_m2_w(lambda))*(I0_m2_w(lambda)*K1_m2(lambda)+I1_m2(lambda)*K0_m2_w(lambda)));

% LSC-induced energy modulation amplitude in me unit
Z_MOD2_LSC_i = @(lambda) Z_MOD2_LSC(lambda)*Lmod2;
Dgamma_MOD2 = @(lambda) -(4*pi/Z0)*(I0/IA)*b4(lambda).*Z_MOD2_LSC(lambda)*Lmod2;

% S-matrix
S_MOD2 = @(lambda) [1 0;-Z_MOD2_LSC_i(lambda)*I2/Ef3 1];


%%----------------------------------------------------------------------------------------------------------------------------------------------%%


%%%%%%%%%%%%%%%%
%% MBI Matrix %%
%%%%%%%%%%%%%%%%

% Transfer matrix from injection to BC1 (included)
M_BC1 = @(lambda) S_BC1(lambda)*S_L1(lambda)*S_L0(lambda);

% Transfer matrix from injection to BC2 (included)
M_BC2 = @(lambda) S_BC2(lambda)*S_L2(lambda)*M_BC1(lambda);

% Transfer matrix from injection to Linac end
M_linac = @(lambda) S_L3(lambda)*M_BC2(lambda);

% Transfer matrix from injection to spreader (included)
M_spreader = @(lambda) S_spreader(lambda)*M_linac(lambda);

% Transfer matrix from injection to first modulator (included)
M_MOD1 = @(lambda) S_MOD1(lambda)*M_spreader(lambda);

% Transfer matrix from injection to Delay Line (included)
M_DL = @(lambda) S_DL(lambda)*S_U(lambda)*M_MOD1(lambda);

% Transfer matrix from injection to second Modulator (included)
M_MOD2 = @(lambda) S_MOD2(lambda)*M_DL(lambda);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MBI matrices entries %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

select = @(M,row,col) M(row,col);

% Gain function after D0+L1+BC1
Gain1 = @(lambda) abs(select(M_BC1(lambda),1,1));

% Energy modulation induced after D0+L1+BC1
Emod1 = @(lambda) abs(select(M_BC1(lambda),2,1));

% Gain function after D0+L1+BC1+L2+BC2
Gain2 = @(lambda) abs(select(M_BC2(lambda),1,1));

% Energy modulation induced after D0+L1+BC1+L2+BC2
Emod2 = @(lambda) abs(select(M_BC2(lambda),2,1));

% Energy modulation induced after D0+L1+BC1+L2+BC2+L3
Emod_linac = @(lambda) abs(select(M_linac(lambda),2,1));

% Gain function after D0+L1+BC1+L2+BC2+L3+Spreader
Gain3 = @(lambda) abs(select(M_spreader(lambda),1,1));

Gain_spreader = @(lambda) abs(select(S_spreader(lambda),1,1));

% Energy modulation induced after D0+L1+BC1+L2+BC2+L3+Spreader
Emod_spreader = @(lambda) abs(select(M_spreader(lambda),2,1));

% Gain function after D0+L1+BC1+L2+BC2+L3+Spreader+MOD1
Gain4 = @(lambda) abs(select(M_MOD1(lambda),1,1));

% Energy modulation induced after D0+L1+BC1+L2+BC2+L3+Spreader+MOD1
Emod_mod1 = @(lambda) abs(select(M_MOD1(lambda),2,1));

% Gain function after D0+L1+BC1+L2+BC2+L3+Spreader+MOD1+Undulators+DL
Gain5 = @(lambda) abs(select(M_DL(lambda),1,1));

% Energy modulation induced after D0+L1+BC1+L2+BC2+L3+Spreader+MOD1+Undulators+DL
Emod_dl = @(lambda) abs(select(M_DL(lambda),2,1));

% Gain function after D0+L1+BC1+L2+BC2+L3+Spreader+MOD1+Undulators+DL+MOD2
Gain6 = @(lambda) abs(select(M_MOD2(lambda),1,1));

% Energy modulation induced after D0+L1+BC1+L2+BC2+L3+Spreader+MOD1+Undulators+DL+MOD2
Emod_mod2 = @(lambda) abs(select(M_MOD2(lambda),2,1));


%%%%%%%%%%%%%%%%
%% MBI Curves %%
%%%%%%%%%%%%%%%%

% Bunching factors (LSC + CSR)
abs_b0 = @(lambda) abs(bk0(lambda))*100;
abs_b1 = @(lambda) abs(select(M_BC1(lambda)*bf0(lambda),1,1))*100;
abs_b2 = @(lambda) abs(select(M_BC2(lambda)*bf0(lambda),1,1))*100;
abs_b3 = @(lambda) abs(select(M_linac(lambda)*bf0(lambda),1,1))*100;
abs_b4 = @(lambda) abs(select(M_spreader(lambda)*bf0(lambda),1,1))*100;
abs_b5 = @(lambda) abs(select(M_MOD1(lambda)*bf0(lambda),1,1))*100;
abs_b6 = @(lambda) abs(select(M_DL(lambda)*bf0(lambda),1,1))*100;
abs_b7 = @(lambda) abs(select(M_MOD2(lambda)*bf0(lambda),1,1))*100;

% Energy Modulation 
Emod1_keV = @(lambda) abs(select(M_BC1(lambda)*bf0(lambda),2,1))*Ef1*1e-3;
Emod2_keV = @(lambda) abs(select(M_BC2(lambda)*bf0(lambda),2,1))*Ef2*1e-3;
Emod_linac_keV = @(lambda) abs(select(M_linac(lambda)*bf0(lambda),2,1))*Ef3*1e-3;
Emod_spreader_keV = @(lambda) abs(select(M_spreader(lambda)*bf0(lambda),2,1))*Ef3*1e-3;
Emod_mod1_keV = @(lambda) abs(select(M_MOD1(lambda)*bf0(lambda),2,1))*Ef3*1e-3;
Emod_dl_keV = @(lambda) abs(select(M_DL(lambda)*bf0(lambda),2,1))*Ef3*1e-3;
Emod_mod2_keV = @(lambda) abs(select(M_MOD2(lambda)*bf0(lambda),2,1))*Ef3*1e-3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Uncorrelated Energy spread %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncorrelated energy spread RMS in induced by LSC in [keV] after BC1
sigE1_LSC = 1e-3*Ef1*sqrt((2*e*c/I0)*integral(@(lambda)abs(Emod1(lambda).^2./(lambda.^2)),1e-6,200e-6,'ArrayValued',true));

% Total energy spread after BC1
sigE1_tot = sqrt((C1^2*sigd_BC1_w^2 + sigdBC1_w_ibs^2)*Ef1^2*1e-6 + sigE1_LSC^2 + switch_lh*C1^2*sigdE_lh^2*1e-6);

% Uncorrelated energy spread RMS in induced by LSC in [keV] after BC2
sigE2_LSC = 1e-3*Ef2*sqrt((2*e*c/I0)*integral(@(lambda)abs(Emod2(lambda).^2./(lambda.^2)),1e-6,200e-6,'ArrayValued',true));

% Total energy spread after BC2
sigE2_tot = sqrt((C2^2*sigd_BC2_w^2 + sigdBC2_w_ibs^2)*Ef2^2*1e-6 + sigE2_LSC^2 + switch_lh*C1^2*C2^2*sigdE_lh^2*1e-6);

% Uncorrelated energy spread RMS in induced by LSC in [keV] at linac end
sigE3_LSC = 1e-3*Ef3*sqrt((2*e*c/I0)*integral(@(lambda)abs(Emod_linac(lambda).^2./(lambda.^2)),1e-6,200e-6,'ArrayValued',true));

% Total energy spread at linac end
sigE3_tot = sqrt((sigd_L3*Ef2)^2*1e-6 + (sigdL3_ibs*Ef3)^2*1e-6 + sigE3_LSC^2 + switch_lh*C1^2*C2^2*sigdE_lh^2*1e-6);

% Uncorrelated energy spread RMS in induced by LSC in [keV] at spreader end
sigES_LSC = 1e-3*Ef3*sqrt((2*e*c/I0)*integral(@(lambda)abs(Emod_spreader(lambda).^2./(lambda.^2)),1e-6,200e-6,'ArrayValued',true));

% Total energy spread at spreader end
sigES_tot = sqrt((sigd_Sd4^2 + sigdS4_ibs^2)*Ef3^2*1e-6 + sigES_LSC^2 + switch_lh*C1^2*C2^2*sigdE_lh^2*1e-6);

% Uncorrelated energy spread RMS in induced by LSC in [keV] at first modulator
sigEMOD1_LSC = 1e-3*Ef3*sqrt(2)*sqrt((e*c/I0)*integral(@(lambda)abs(Emod_mod1(lambda).^2./(lambda.^2)),1e-6,200e-6,'ArrayValued',true));

% Total energy spread at first modulator end
sigEMOD1_tot = sqrt((sigd_Sd4^2 + sigdS4_ibs^2 + sigdMOD1_ibs^2)*Ef3^2*1e-6 + sigEMOD1_LSC^2 + switch_lh*C1^2*C2^2*sigdE_lh^2*1e-6);

% Uncorrelated energy spread RMS in induced by LSC in [keV] at Delay line end
sigEDL_LSC = 1e-3*Ef3*sqrt(2)*sqrt((e*c/I0)*integral(@(lambda)abs(Emod_dl(lambda).^2./(lambda.^2)),1e-6,200e-6,'ArrayValued',true));

% Total energy spread at Delay line end
sigEDL_tot = sqrt((sigd_DL_w^2 + sigdDL_w_ibs^2)*Ef3^2*1e-6 + sigEDL_LSC^2 + switch_lh*C1^2*C2^2*sigdE_lh^2*1e-6);

% Uncorrelated energy spread RMS in induced by LSC in [keV] at second modulator
sigEMOD2_LSC = 1e-3*Ef3*sqrt(2)*sqrt((e*c/I0)*integral(@(lambda)abs(Emod_mod2(lambda).^2./(lambda.^2)),1e-6,200e-6,'ArrayValued',true));

% Total energy spread at second modulator end
sigEF_tot = sqrt((sigd_DL_w^2 + sigdDL_w_ibs^2 + sigdMOD2_ibs^2)*Ef3^2*1e-6 + sigEMOD2_LSC^2 + switch_lh*C1^2*C2^2*sigdE_lh^2*1e-6);



%%%%%%%%%%
%% Plot %%
%%%%%%%%%%


% Mean Energy in MeV
f10=figure(10);
xlabel('s [m]','FontSize',16)
ylabel('Mean Energy [MeV]','FontSize',16)
set(gca,'FontSize',16)
set(gcf,'DefaultLineLineWidth',4)
me=0.511;
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
plot(u2+L0+L1,energy2(u2),'m-')
plot(u3+L0+L1+L2,energy3(u3),'b-')
xlabel('s [m]','FontSize',16)
ylabel('Mean Energy [MeV]','FontSize',16)
set(gca,'FontSize',16)
saveas(f10,['./BK_mean_energy_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])

% Effective beam radius in micron
f20=figure(20);
xlabel('s [m]','FontSize',16)
ylabel('Effective Beam Radius [\mum]','FontSize',16)
set(gca,'FontSize',16)
set(gcf,'DefaultLineLineWidth',4)
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
saveas(f20,['./BK_beam_size_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])

if switch_spreader == 0 && switch_EEHG == 0
    
    % LSC-induced energy modulation amplitude in keV vs. uncompressed
    % wavelength in micron
    f30=figure(30);
    fplot(Emod1_keV,[1e-6,200e-6],'r:','LineWidth',3)
    hold on
    fplot(Emod2_keV,[1e-6,200e-6],'m-.','LineWidth',3)
    fplot(Emod_linac_keV,[1e-6,200e-6],'b-','LineWidth',3)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('\DeltaE_{MBI} [keV]','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    legend('Exit of BC1','Exit of BC2','Linac End')
    set(gca,'FontSize',16)
    saveas(f30,['./BK_energy_mod_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])

    % Bunching factor
    f40=figure(40);
    fplot(abs_b0,[1e-6,200e-6],'g-','LineWidth',3)
    hold on
    fplot(abs_b1,[1e-6,200e-6],'r:','LineWidth',3)
    fplot(abs_b2,[1e-6,200e-6],'m-..','LineWidth',3)
    fplot(abs_b3,[1e-6,200e-6],'b-','LineWidth',3)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('Bunching factor [%]','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    legend('Initial','Exit of BC1','Exit of BC2','Linac End')
    set(gca,'FontSize',16)
    saveas(40,['./BK_bunching_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])

    % Spectral gain
    f50=figure(50);
    fplot(Gain1,[1e-6,200e-6],'r:','LineWidth',3)
    hold on
    fplot(Gain2,[1e-6,200e-6],'m-.','LineWidth',3)
    hold on
    fplot(Gain2,[1e-6,200e-6],'b-','LineWidth',3)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('Gain','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    legend('Exit of BC1','Exit of BC2','Linac End')
    set(gca,'FontSize',16)
    saveas(f50,['./BK_gain_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])    

elseif switch_spreader ~= 0 && switch_EEHG == 0

    % LSC-induced energy modulation amplitude in keV vs. uncompressed
    % wavelength in micron
    f31=figure(31);
    fplot(Emod1_keV,[1e-6,200e-6],'r:','LineWidth',3)
    hold on
    fplot(Emod2_keV,[1e-6,200e-6],'m-.','LineWidth',3)
    fplot(Emod_linac_keV,[1e-6,200e-6],'b-','LineWidth',3)
    fplot(Emod_spreader_keV,[1e-6,200e-6],'k','LineWidth',3)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('\DeltaE_{MBI} [keV]','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    legend('Exit of BC1','Exit of BC2','Linac End','Spreader End')
    set(gca,'FontSize',16)
    saveas(f31,['./BK_energy_mod_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])

    % Bunching factor
    f41=figure(41);
    fplot(abs_b0,[1e-6,200e-6],'g-','LineWidth',3)
    hold on
    fplot(abs_b1,[1e-6,200e-6],'r:','LineWidth',3)
    fplot(abs_b2,[1e-6,200e-6],'m-.','LineWidth',3)
    fplot(abs_b3,[1e-6,200e-6],'b-','LineWidth',3)
    fplot(abs_b4,[1e-6,200e-6],'k-','LineWidth',3)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('Bunching factor [%]','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    legend('Initial','Exit of BC1','Exit of BC2','Linac End','Spreader End')
    set(gca,'FontSize',16)
    saveas(f41,['./BK_bunching_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])
    
    % Spectral gain
    f51=figure(51);
    fplot(Gain1,[1e-6,200e-6],'r:','LineWidth',3)
    hold on
    fplot(Gain2,[1e-6,200e-6],'m-.','LineWidth',3)
    fplot(Gain3,[1e-6,200e-6],'k-','LineWidth',3)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('Gain','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})   
    legend('Exit of BC1','Exit of BC2','Spreader End')
    set(gca,'FontSize',16)
    saveas(f51,['./BK_gain_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])
    
    % Spectral gain of Spreader only
    f61=figure(61);
    fplot(Gain_spreader,[1e-6,200e-6],'LineWidth',3)
    xlabel('\lambda [\mum]','FontSize',16)
    ylabel('Gain','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    set(gca,'FontSize',16)
    saveas(f61,['./BK_gain_spreader_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])
    
    
elseif switch_spreader ~= 0 && switch_EEHG ~= 0

    % LSC-induced energy modulation amplitude in keV vs. uncompressed
    % wavelength in micron
    f32=figure(32);
    fplot(Emod1_keV,[1e-6,200e-6],'r:','LineWidth',3)
    hold on
    fplot(Emod2_keV,[1e-6,200e-6],'m-.','LineWidth',3)
    fplot(Emod_linac_keV,[1e-6,200e-6],'b-','LineWidth',3)
    fplot(Emod_spreader_keV,[1e-6,200e-6],'k-','LineWidth',3)
    fplot(Emod_mod1_keV,[1e-6,200e-6],'y--','LineWidth',3)
    fplot(Emod_dl_keV,[1e-6,200e-6],'g--','LineWidth',3)
    fplot(Emod_mod2_keV,[1e-6,200e-6],'c--','LineWidth',3)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('\DeltaE_{MBI} [keV]','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    legend('Exit of BC1','Exit of BC2','Linac End','Spreader End','Exit of MOD1','Exit of DL','Exit of MOD2')
    set(gca,'FontSize',16)
    saveas(f32,['./BK_energy_mod_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])

    % Bunching factor
    f42=figure(42);
    fplot(abs_b0,[1e-6,200e-6],'g-','LineWidth',3)
    hold on
    fplot(abs_b1,[1e-6,200e-6],'r:','LineWidth',3)
    fplot(abs_b2,[1e-6,200e-6],'m-.','LineWidth',3)
    fplot(abs_b3,[1e-6,200e-6],'b-','LineWidth',3)
    fplot(abs_b4,[1e-6,200e-6],'k-','LineWidth',3)
    fplot(abs_b5,[1e-6,200e-6],'y--','LineWidth',3)
    fplot(abs_b6,[1e-6,200e-6],'g--','LineWidth',3)
    fplot(abs_b7,[1e-6,200e-6],'c--','LineWidth',3)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('Bunching factor [%]','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    legend('Initial','Exit of BC1','Exit of BC2','Linac End','Spreader End','Exit of MOD1','Exit of DL','Exit of MOD2')
    set(gca,'FontSize',16)
    saveas(f42,['./BK_bunching_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])
    
    % Spectral gain 
    f52=figure(52);
    fplot(Gain1,[1e-6,200e-6],'r:','LineWidth',3)
    hold on
    fplot(Gain2,[1e-6,200e-6],'m-.','LineWidth',3)
    fplot(Gain3,[1e-6,200e-6],'b-','LineWidth',3)
    fplot(Gain4,[1e-6,200e-6],'y--','LineWidth',3)
    fplot(Gain5,[1e-6,200e-6],'g--','LineWidth',3)
    fplot(Gain6,[1e-6,200e-6],'c--','LineWidth',3)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('Gain','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})   
    legend('Exit of BC1','Exit of BC2','Spreader End','Exit of MOD1','Exit of DL','Exit of MOD2')
    set(gca,'FontSize',16)
    saveas(f52,['./BK_gain_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])
    
    % Spectral gain of Spreader only
    f62=figure(62);
    fplot(Gain_spreader,[1e-6,200e-6],'LineWidth',3)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('Gain','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    set(gca,'FontSize',16)
    saveas(f62,['./BK_gain_spreader_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])
    
end
    
    % f70 = figure(70);
    % l = (1:10:500)*1e-6;
    % s = (1:10:500)*1e-6;
    % [X,Y] = meshgrid(l,s);
    % e = exp(-switch_Hx*(k2(X)*thetas.*Y).^2);
    % p = pcolor(X,Y,e);
    % p.FaceColor = 'interp';
    % fplot(e,[1e-6,200e-6],'LineWidth',3)
    % set(gca,'FontSize',16)
    % xlabel('\lambda_0 [\mum]','FontSize',16)
    % ylabel('\sigma_x','FontSize',16)
    % set(gca,'XTickLabel',{'50','100','150','200'})
    


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
fprintf('%s %s \n', question, yes);

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
value=1e-3*sigdE;
fprintf('%s %d keV \n', label, value);

label='sigE_LH =';
value=1e-3*sigdE_lh;
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

if switch_spreader == 0 && switch_EEHG == 0
elseif switch_spreader == 0 && switch_EEHG ~= 0
    label='R56_DL =';
    value=1e3*R56_DL;
    fprintf('%s %d mm \n', label, value);
elseif switch_spreader ~= 0 && switch_EEHG == 0
    label='R56_Spreader =';
    value=1e3*R56_Spreader;
    fprintf('%s %d mm \n', label, value);
elseif switch_spreader ~= 0 && switch_EEHG ~= 0    
    label='R56_Spreader =';
    value=1e3*R56_Spreader;
    fprintf('%s %d mm \n', label, value);
    label='R56_DL =';
    value=1e3*R56_DL;
    fprintf('%s %d mm \n', label, value);
end

label='sigE_IBS(L0+L1) = ';
value=1e-3*(sigEL0_ibs+sigEL1_ibs);
fprintf('%s %d keV \n', label, value);

label='sigE_IBS(BC1) = ';
value=1e-3*sqrt(sigEBC1_i^2+sigEBC1_w^2);
fprintf('%s %d keV \n', label, value);

label='sigE_IBS(L2+L3) = ';
value=1e-3*sigEL2_ibs;
fprintf('%s %d keV \n', label, value);

label='sigE_IBS(BC2) = ';
value=1e-3*sqrt(sigEBC2_i^2+sigEBC2_w^2);
fprintf('%s %d keV \n', label, value);

label='sigE_IBS(L4) = ';
value=1e-3*sigEL3_ibs;
fprintf('%s %d keV \n', label, value);

if switch_spreader ~= 0 
    label='sigE_IBS(SPRD) = ';
    value=1e-3*Ef3*sigdS4_ibs;
    fprintf('%s %d keV \n', label, value);
end

if switch_EEHG ~= 0
    label='sigE_IBS(MOD1) = ';
    value=1e-3*sigEMOD1_ibs;
    fprintf('%s %d keV \n', label, value);

    label='sigE_IBS(DL) = ';
    value=1e-3*Ef3*sqrt(sigdDL_i_ibs^2 + sigdDL_w_ibs^2);
    fprintf('%s %d keV \n', label, value);

    label='sigE_IBS(MOD2) = ';
    value=1e-3*sigEMOD2_ibs;
    fprintf('%s %d keV \n', label, value);
end


if switch_spreader == 0 && switch_EEHG == 0

    if R561 == 0 && R562 == 0

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

    elseif R561 ~= 0 && R562 == 0

        label='sigE_coll at BC1 end = ';
        value=sigE1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at BC1 end = ';
        value=sigE1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

    elseif R561 == 0 && R562 ~= 0

        label='sigE_coll after L0+L1+BC1+L2) = ';
        value=sigE2_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot after BC2 = ';
        value=sigE2_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

    elseif R561 ~= 0 && R562 ~= 0

        label='sigE_coll (L0+L1) = ';
        value=sigE1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at BC1 end = ';
        value=sigE1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll (L2+L3) = ';
        value=sigE2_LSC;
        fprintf('%s %d keV \n', label, value)

        label='sigE_tot at BC2 end = ';
        value=sigE2_tot;
        fprintf('%s %d keV \n', label, value)

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);
    end

elseif switch_spreader ~= 0 && switch_EEHG == 0

    if R561 == 0 && R562 == 0

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at spreader end = ';
        value=sigES_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at spreader end = ';
        value=sigES_tot;
        fprintf('%s %d keV \n', label, value);

    elseif R561 ~= 0 && R562 == 0

        label='sigE_coll at BC1 end = ';
        value=sigE1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at BC1 end = ';
        value=sigE1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at spreader end = ';
        value=sigES_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at spreader end = ';
        value=sigES_tot;
        fprintf('%s %d keV \n', label, value);

    elseif R561 == 0 && R562 ~= 0

        label='sigE_coll (L0+L1+BC1+L2)  = ';
        value=sigE2_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at BC2 end = ';
        value=sigE2_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at spreader end = ';
        value=sigES_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at spreader end = ';
        value=sigES_tot;
        fprintf('%s %d keV \n', label, value);

    elseif R561 ~= 0 && R562 ~= 0

        label='sigE_coll at BC1 end = ';
        value=sigE1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at BC1 end = ';
        value=sigE1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll (L2+L3) = ';
        value=sigE2_LSC;
        fprintf('%s %d keV \n', label, value)

        label='sigE_tot at BC2 end = ';
        value=sigE2_tot;
        fprintf('%s %d keV \n', label, value)

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at spreader end = ';
        value=sigES_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at spreader end = ';
        value=sigES_tot;
        fprintf('%s %d keV \n', label, value);

    end

elseif switch_spreader ~= 0 && switch_EEHG ~= 0

    if R561 == 0 && R562 == 0

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at spreader end = ';
        value=sigES_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at spreader end = ';
        value=sigES_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at first modulator end = ';
        value=sigEMOD1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at first modulator end = ';
        value=sigEMOD1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at Delay line end = ';
        value=sigEDL_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at Delay line end = ';
        value=sigEDL_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at second modulator end = ';
        value=sigEMOD2_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at second modulator end = ';
        value=sigEF_tot;
        fprintf('%s %d keV \n', label, value);

    elseif R561 ~= 0 && R562 == 0

        label='sigE_coll at BC1 end = ';
        value=sigE1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at BC1 end = ';
        value=sigE1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at spreader end = ';
        value=sigES_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at spreader end = ';
        value=sigES_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at first modulator end = ';
        value=sigEMOD1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at first modulator end = ';
        value=sigEMOD1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at Delay line end = ';
        value=sigEDL_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at Delay line end = ';
        value=sigEDL_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at second modulator end = ';
        value=sigEMOD2_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at second modulator end = ';
        value=sigEF_tot;
        fprintf('%s %d keV \n', label, value);

    elseif R561 == 0 && R562 ~= 0

        label='sigE_coll (L0+L1+BC1+L2)  = ';
        value=sigE2_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at BC2 end = ';
        value=sigE2_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at spreader end = ';
        value=sigES_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at spreader end = ';
        value=sigES_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at first modulator end = ';
        value=sigEMOD1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at first modulator end = ';
        value=sigEMOD1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at Delay line end = ';
        value=sigEDL_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at Delay line end = ';
        value=sigEDL_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at second modulator end = ';
        value=sigEMOD2_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at second modulator end = ';
        value=sigEF_tot;
        fprintf('%s %d keV \n', label, value);

    elseif R561 ~= 0 && R562 ~= 0

        label='sigE_coll at BC1 end = ';
        value=sigE1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at BC1 end = ';
        value=sigE1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll (L2+L3) = ';
        value=sigE2_LSC;
        fprintf('%s %d keV \n', label, value)

        label='sigE_tot at BC2 end = ';
        value=sigE2_tot;
        fprintf('%s %d keV \n', label, value)

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at spreader end = ';
        value=sigES_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at spreader end = ';
        value=sigES_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at first modulator end = ';
        value=sigEMOD1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at first modulator end = ';
        value=sigEMOD1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at Delay line end = ';
        value=sigEDL_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at Delay line end = ';
        value=sigEDL_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at second modulator end = ';
        value=sigEMOD2_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at second modulator end = ';
        value=sigEF_tot;
        fprintf('%s %d keV \n', label, value);
    end

elseif switch_spreader == 0 && switch_EEHG ~= 0

    if R561 == 0 && R562 == 0

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at spreader end = ';
        value=sigES_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at spreader end = ';
        value=sigES_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at first modulator end = ';
        value=sigEMOD1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at first modulator end = ';
        value=sigEMOD1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at Delay line end = ';
        value=sigEDL_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at Delay line end = ';
        value=sigEDL_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at second modulator end = ';
        value=sigEMOD2_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at second modulator end = ';
        value=sigEF_tot;
        fprintf('%s %d keV \n', label, value);

    elseif R561 ~= 0 && R562 == 0

        label='sigE_coll at BC1 end = ';
        value=sigE1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at BC1 end = ';
        value=sigE1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at spreader end = ';
        value=sigES_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at spreader end = ';
        value=sigES_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at first modulator end = ';
        value=sigEMOD1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at first modulator end = ';
        value=sigEMOD1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at Delay line end = ';
        value=sigEDL_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at Delay line end = ';
        value=sigEDL_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at second modulator end = ';
        value=sigEMOD2_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at second modulator end = ';
        value=sigEF_tot;
        fprintf('%s %d keV \n', label, value);

    elseif R561 == 0 && R562 ~= 0

        label='sigE_coll (L0+L1+BC1+L2)  = ';
        value=sigE2_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at BC2 end = ';
        value=sigE2_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at spreader end = ';
        value=sigES_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at spreader end = ';
        value=sigES_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at first modulator end = ';
        value=sigEMOD1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at first modulator end = ';
        value=sigEMOD1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at Delay line end = ';
        value=sigEDL_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at Delay line end = ';
        value=sigEDL_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at second modulator end = ';
        value=sigEMOD2_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at second modulator end = ';
        value=sigEF_tot;
        fprintf('%s %d keV \n', label, value);

    elseif R561 ~= 0 && R562 ~= 0

        label='sigE_coll at BC1 end = ';
        value=sigE1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at BC1 end = ';
        value=sigE1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll (L2+L3) = ';
        value=sigE2_LSC;
        fprintf('%s %d keV \n', label, value)

        label='sigE_tot at BC2 end = ';
        value=sigE2_tot;
        fprintf('%s %d keV \n', label, value)

        label='sigE_coll at linac end = ';
        value=sigE3_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at linac end = ';
        value=sigE3_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at spreader end = ';
        value=sigES_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at spreader end = ';
        value=sigES_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at first modulator end = ';
        value=sigEMOD1_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at first modulator end = ';
        value=sigEMOD1_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at Delay line end = ';
        value=sigEDL_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at Delay line end = ';
        value=sigEDL_tot;
        fprintf('%s %d keV \n', label, value);

        label='sigE_coll at second modulator end = ';
        value=sigEMOD2_LSC;
        fprintf('%s %d keV \n', label, value);

        label='sigE_tot at second modulator end = ';
        value=sigEF_tot;
        fprintf('%s %d keV \n', label, value);
    end

end


