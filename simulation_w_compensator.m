% SET RIGHT POLES

clear all
% close all
clc
global h kappa 

disp('setting constants...')
 setConstants; % Creates User Defined Constant and add to workspace

    h=4.0693e-06;% ratio between helmholtz current to field magnitude
    % it is the field magnitude by a 1A current runs in x coil on r=(0,0)
    
    kappa=-1.1700e-05;% ratio between maxwell current to field gradient
    % it is the field d/dx made by a 1A current runs in x coil on r=(0,0)

% --------------------------------------
   
    en=11;
% sample rate:
    Sample.dt=en/1001;
    dt=Sample.dt;
% controler options
    ClosedLoop=1; 
    % select action (closed loop or open loop comparison)

disp('predefined path')
    sigma=en/2;
    f=0;
    Route.start_time=0;
    Route.end_time=en;      
    Route.xfun=@(t) .1*cos(5*pi/en*t) ; 
    Route.yfun=@(t) .1*sin(5*pi/en*t);
    
    tt=(Route.start_time:dt:Route.end_time);  
    x=Route.xfun(tt)-Route.xfun(0); % route x samples
    y=Route.yfun(tt)-Route.yfun(0); % route y samples 
    figure(40)
plot(x,y)
%%
% ------------------- build spline ----------------------
    %x(t)=an(1,t)*t^2+an(2,t)*t+an(3,t)
    %y(t)=bn(1,t)*t^2+bn(2,t)*t+bn(3,t)
disp('spline building...')
    [an, bn]=splineBuild(x,y,tt);
% --------------------------------------------------------    

% -------------- Calculating F(t) ------------------------
% force on time t depends on acceleration -> 2nd derivative of route
    F=2*Dipole.mass*[an(1,:);bn(1,:)]+Medium.viscosity*[2*an(1,:).*tt+an(2,:);
                                                        2*bn(1,:).*tt+bn(2,:)];
%         F(2,:)=zeros(1,numel(tt));
%         F(1,:)=[10*ones(1,1) zeros(1,499) zeros(1,numel(tt)-500)];
%         
%         v=cumtrapz(tt.',F.'/m);
%         v(1,:)=v(1,:)+0;
%         v(2,:)=v(2,:)+0;
%         
%         r=cumtrapz(tt.',v);
%         x=r(:,1).';
%         y=r(:,2).'+0;

% ---------------------------------------------------------

% -------------- B Field Direction Calculating ------------
% field direction depends on direction of velocity -> 1st derivative of route
    phiB=(atan2(2*bn(1,:).*tt+bn(2,:),2*an(1,:).*tt+an(2,:)));
% ---------------------------------------------------------
    Bmag=BmW*ones(1,numel(tt));


disp('setting simulation...')
    simIn=timeseries([F;phiB],tt,'name','Predefined Route');
 % Initial Values
    initialValues.R0=[x(1);y(1)];
    initialValues.V0=[diff(x(1:2))/dt;diff(y(1:2))/dt];
    initialValues.phiB=phiB(1);
 % -----------
    save('par','simIn','-v7.3')
    
 % Create Compensator:

 % linear model:
disp('Calculating Linear Model..')

disp('loading TF')
% loading general transfer function
    Gstruct=load('linMat','A','B','C');
    Amatrix=Gstruct.A;
    Bmatrix=Gstruct.B;
    Cmatrix=Gstruct.C;

% substituting constants:
    Amatrix=subs(Amatrix,{'Rx','Lx','Ry','Ly','h'},{Rx,Lx,Rx,Lx,h}); % substituting coil parameters
    Amatrix=subs(Amatrix,{'M','I','alpha'},{Dipole.moment,Dipole.inertia,alpha});
    Bmatrix=subs(Bmatrix,{'M','Rx','Lx','Ry','Ly','h','kappa'},{Dipole.moment,Rx,Lx,Rx,Lx,h,kappa}); % substituting coil parameters
    Cmatrix=subs(Cmatrix,{'Rx','Lx','Ry','Ly','h','kappa'},{Rx,Lx,Rx,Lx,h,kappa}); % substituting coil parameters
% ----------------------
        disp('initializing')
        phiB_ss=mean(phiB); 
        F0=F(:,1);
        Bm0=BmW;
    % initial Currents:    
        Ihx0=Bm0/h*cos(phiB_ss);
        Ihy0=Bm0/h*sin(phiB_ss);
        Imx0=2/(3*M*kappa)*((2*F0(1))/cos(phiB_ss) +(F0(2))/sin(phiB_ss));
        Imy0=2/(3*M*kappa)*((F0(1))/cos(phiB_ss) +(2*F0(2))/sin(phiB_ss));
        
    % substituting initial values    
    disp('substituting part parameters...')
    Amatrix=double(subs(Amatrix,{'Ihx','Ihy','Imx','Imy'},{Ihx0,Ihy0,Imx0,Imy0}));
    Bmatrix=subs(Bmatrix,{'Fx_hat','Fy_hat'},{F0(1),F0(2)}); % substituting F0
    Bmatrix=double(subs(Bmatrix,{'Bm_hat','phiB_hat'},{Bm0,phiB_ss})); % substituting initial B
    Cmatrix=subs(Cmatrix,'M',Dipole.moment,kappa);
    Cmatrix=double(subs(Cmatrix,{'Ihx','Ihy','Imx','Imy'},{Ihx0,Ihy0,Imx0,Imy0}));
    % converting to tf
    disp('Creating open loop system')
    StateSpace=(ss(Amatrix,Bmatrix,Cmatrix,0));
    G_ol=minreal(tf(StateSpace));
% -------------- Open loop system ready ---------------
if ClosedLoop
     loopstr='closed';

        % act 1: calculating left factors (non-coprime)
        disp('Calculate Left Fraction')
        [LN,LD]=left_poly_fractions(G_ol);
    % act 2: Creating coefficient matrix Ngal,Dgal (not co-prime)
        disp('Coefficient Matrixes')
        [Ngal, Dgal]=Coefficient_Matrixes(LN,LD);
    % act 3: Co-prime Factorization (pole zero cancelation)
        disp('Co-prime Factorization...')
        [Ngal_cp,Dgal_cp]=coprime_Factorization(Ngal,Dgal,1e-8);
        [Ngal_cp,Dgal_cp]=TFSimplify(Ngal_cp,Dgal_cp,1e-10);
        [Ntf, Dtf]=create_tf(Ngal_cp,Dgal_cp); % needed only for debugging
    % act 4: create compensator:
        disp('Compensator Design')
        psi=[pi pi-pi/40 pi+pi/40];
        poles=1e4*exp(1j*psi);
        poles=1e4*[-1 -1-.001j -1+.001j];
    %     poles=[-.01 -.01-.1j -.01+.1j];
        [Agal, Bgal, Ftf]=Compensator_design(Ngal_cp,Dgal_cp,poles);
        [Bgal,Agal]=TFSimplify(Bgal,Agal,1e-10);
        Agal=real(Agal);
        Bgal=real(Bgal);
        [Atf, Btf]=create_tf(Agal,Bgal); % needed only for debugging

    % ----------------------Linear Simulation:--------------------
        subMat_len=size(G_ol,1);
        Ctf=(Btf*Atf^-1); % Compensator
        disp('Compensator ready')
        G_cl=minreal(Ntf*(Atf*Dtf+Btf*Ntf)^-1*Btf); % closed loop - without amp
%         G_cl=minreal((eye(3)+Ntf*Dtf^-1*Ctf)^-1*Ntf*Dtf^-1*Ctf);
        [ystp_cl, ~]=step(G_cl,tt);
        amp_cl=diag(diag(1./squeeze(real(ystp_cl(end,:,:)))));
    % ----------------

        % multiplying systems by current AMP
        G_cl=amp_cl*G_cl; % closed loo 
            % ---------
    disp('Close loop system ready')
     % running simulation with input:
        % set input vector:
           u=[F; mod(phiB+pi,2*pi)-pi].';
        % simulating closed and open loop:
            % closed loop:
         [y_sim_cl,~]=lsim(G_cl,u,tt);
          y_sim_cl=real(y_sim_cl);
          y_sim_cl(:,3)=mod(y_sim_cl(:,3)+pi,2*pi)-pi; % setting limits of phi fo plotting

%      % comparison between Open and Closed loop:
%          for n=1:subMat_len       
%              switch n
%                 case 1 
%                     strTit='F_x'; 
%                 case 2 
%                     strTit='F_y';
%                 case 3 
%                     strTit='\phi_B';
%              end   
%              figure(5)
%              % plot CL results
%                subplot(subMat_len, 1,n)
%                plot(tt,u(:,n),'--b',tt,y_sim_cl(:,n),'-k')
%                legend('Input','Output(CL)')
%                title (['Linear system: Closed Loop ',strTit])
%          end      
else % not closed loop
    loopstr='Open';
    Ctf=eye(3);
    amp_cl=eye(3);
    G_cl=zeros(3);
end

% Run nonlinear simulation of closed loop with compensator:
 disp('starting nonlinear simulation... this may take a while')
    openModels = find_system('SearchDepth', 0);
    tic
    if ~any(strcmp(openModels,'ForceControl_closed_loop'))
        open 'ForceControl_closed_loop.slx'
    end
    set_param('ForceControl_closed_loop', 'StopTime', num2str(tt(end)))
    set_param('ForceControl_closed_loop', 'StartTime', num2str(tt(1)))
    %%
    sim('ForceControl_closed_loop')
    toc
 disp('plotting results, actuall model')
simsamp=floor(linspace(1,numel(simForce.time),size(simForce.data,1)/10));

figure(10*ClosedLoop+1)
    subplot(2,1,1)
        fig=plot(tt,F(1,:),simForce.time(simsamp),squeeze(simForce.data(simsamp,1)));
        set(fig(1),'lineWidth',3)
        legend('F wanted','F got')
        % axis([0 en -7*abs(max(F(1,:))) 7*abs(max(F(1,:)))])
        title (['Fx comparison: ',loopstr])

        subplot(2,1,2)
        fig=plot(tt,F(2,:),simForce.time(simsamp),squeeze(simForce.data(simsamp,2)));
        set(fig(1),'lineWidth',3)
        legend('F wanted','F got')
        title (['Fy comparison: ',loopstr])

figure(10*ClosedLoop+2)
    subplot(2,2,1)
        plot(simCurrent.time(simsamp),squeeze(simCurrent.data(simsamp,1)))
        title (['Ih_x: ',loopstr])
        xlabel 't [sec]'
        ylabel '[A]'
    subplot(2,2,2)
        plot(simCurrent.time(simsamp),squeeze(simCurrent.data(simsamp,2)))
        title (['Ih_y: ',loopstr])
        xlabel 't [sec]'
        ylabel '[A]'
    subplot(2,2,3)
        plot(simCurrent.time(simsamp),squeeze(simCurrent.data(simsamp,3)))
        title (['Im_x: ',loopstr])
        xlabel 't [sec]'
        ylabel '[A]'
    subplot(2,2,4)
        plot(simCurrent.time(simsamp),squeeze(simCurrent.data(simsamp,4)))
        title (['Im_y: ',loopstr])
        xlabel 't [sec]'
        ylabel '[A]'
figure(10*ClosedLoop+3)
    fig=plot(x,y,simDipoleLoc.data(simsamp,1),simDipoleLoc.data(simsamp,2));
    legend('Desired','Simulated')
    title (['Route Comparison: ',loopstr])
    set(fig,'LineWidth',1.5)
    grid

 
 
% axis([-1 1 -1 1])
T=simDipoleLoc.time(simsamp);
xSim=simDipoleLoc.data(simsamp,1);
ySim=simDipoleLoc.data(simsamp,2);
phiSim=simPhiB.data(simsamp);
% plot(tt,phiB,'b*',simphiB.time(simsamp),simphiB.data(simsamp),'r')
figure(10*ClosedLoop+4)
    subplot(3,1,1)
        fig=plot(tt,x,'k',T,xSim,'b--');
        set(fig(1),'lineWidth',1.7)
        legend ('Desired','Simulated')
        title (['x(t) full model - ',loopstr,' loop'])
    subplot(3,1,2)
        fig=plot(tt,y,'k',T,ySim,'b--');
        set(fig(1),'lineWidth',1.7)
        title (['y(t) full model - ',loopstr,' loop'])
        legend ('Desired','Simulated')
    subplot(3,1,3)
        fig=plot(tt,mod(phiB,2*pi),'k',simPhiM.time,mod(simPhiM.data,2*pi),'b:');
        title (['phiM(t) full model - ',loopstr,' loop'])
        legend ('Desired','Simulated')
        set(fig(1),'lineWidth',1.7)
%
% forceControMatl_linearModel  
