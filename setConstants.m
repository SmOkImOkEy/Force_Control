%% NO FRICTION
% This Script creates User Defined Constant and add to workspace

global g mu0 % constants
% constants
% parameters (GUI): ------------------------------------
    
    %  field constant
        mu0=4*pi*1e-7; % [Hy/m]
        g=9.8;
    % geometric constants:
        xCoil.Radius=.5; % mean coil Radius
        xCoil.Turns=600; % number of turns
        
        yCoil.Radius=.5; % mean coil Radius
        yCoil.Turns=600; % number of turns

        Rs=0; %[Ohm] source resistance
        AWG=10;
        switch AWG
            case 13
                copperRes=6.57e-3;  %[Ohm/m] copper  wire AWG13 (r= 0.914mm) (I_fuse=198A)
                preece=198; % ~10 sec fusing current
                yCoil.Width=0.08; % coil width (R2-R1)
                yCoil.Height=0.03; % coil width (w)
                xCoil.Width=0.08; % coil width (R2-R1)   
                xCoil.Height=0.03; % coil height (w)
                
            case 10
                copperRes=3.28e-3; %[Ohm/m] copper  wire AWG10 (r= 1.2941) (I_fuse=333A)
                preece=333; % ~10 sec fusing current
                yCoil.Width=0.12; % coil width (R2-R1)
                yCoil.Height=0.045; % coil width (w)
                xCoil.Width=0.12; % coil width (R2-R1)   
                xCoil.Height=0.045; % coil height (w)
        end
        
    maxCurr=inf%0.3*preece; % third of preece (~10sec)current   
% Medium Parameters
    Medium.mu=0;%0.3;
    Medium.alpha=0; % liquid friction coefficient for moment    
    Medium.viscosity=0;%0.2;
    gain=1; % conpensator is regular PID controller on this version
    alpha=Medium.alpha;
% Capsule Parameters    
    capLen=30e-3; % capsule length (30mm)
    capDia=11e-3; % capsule diameter (11mm) (height on 2d space)
    Dipole.moment=0.394; % capsule magnetic moment magnitude%
    Dipole.mass=4e-3;    % capsule mass (4gr)
    Dipole.inertia=pi/4*(capLen/2)*(capDia/2)^3; % moment of inertia
% -----------------------------------------------------------------

 % some more parameters
    Misc.maxCurr=maxCurr;    
    Misc.Bmag=1e-4; 
    BmW=Misc.Bmag;
 
 % coil Model Parameters
    xCoil.R=copperRes*2*pi*xCoil.Radius*xCoil.Turns;
    xCoil.L=1e-6*0.8*0.0254*xCoil.Radius^2*xCoil.Turns^2/(6*xCoil.Radius+9*xCoil.Height+10*xCoil.Width); % wheelers approximation
    yCoil.R=copperRes*2*pi*yCoil.Radius*yCoil.Turns;
    yCoil.L=1e-6*0.8*0.0254*yCoil.Radius^2*yCoil.Turns^2/(6*yCoil.Radius+9*yCoil.Height+10*yCoil.Width); % wheelers approximation

% loading previously calculated field magnitude and gradient
% for each coil maxwell and helmholtz
%     load('Bmat_file.mat')
Rx=xCoil.R+Rs;
Ry=yCoil.R+Rs;
Lx=xCoil.L;
Ly=yCoil.L;
M=Dipole.moment;
I=Dipole.inertia;

Nx=xCoil.Turns;
Ny=yCoil.Turns;
m=Dipole.mass;