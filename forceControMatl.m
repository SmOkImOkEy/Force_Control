% main Function
% clear all;
% global B dBx dBy X Y  %#ok<NUSED>
clear all
    close all
    clc
global h kappa 
% h = waitbar(0,'Please wait...');

% End_Time = 100;

% myoptions = simset('SrcWorkspace','current','DstWorkspace','current','FixedStep', 0.1);
% sim('model',End_Time,myoptions)

% close(h)
disp('setting constants...')
 setConstants; % Creates User Defined Constant and add to workspace
% B shape: B(X,Y,coil,Direction)
% coil=1->x axis helmholtz
% coil=2->y axis helmholtz
% coil=3->x axis maxwell
% coil=4->y axis maxwell

% calculating linearization parameters of field control:
%     [Bi, ~, ~]=B_Field_Interp(0,0,[1;0;0;0]);
    h=4.0693e-06;% ratio between helmholtz current to field magnitude
    % it is the field magnitude by a 1A current runs in x coil on r=(0,0)
    
%     [~, dBxi, ~]=B_Field_Interp(0,0,[0;0;1;0]);
    kappa=-1.1700e-05;% ratio between maxwell current to field gradient
    % it is the field d/dx made by a 1A current runs in x coil on r=(0,0)
% Open or closed loop? 0-open; 1-closed
    state=0;
    if state==0
        loopstr='open';
    else
        loopstr='closed';
    end
    disp(loopstr)
% --------------------------------------
   
    en=10;
% sample rate:
    Sample.dt=en/1001;
    dt=Sample.dt;
 
disp('predefined path')
    Route.start_time=0;
    Route.end_time=en;      
    Route.xfun=@(t) 1/(2*m)*t.^2;   
    Route.yfun=@(t) 2/(2*m)*t.^2; 
    tt=(Route.start_time:dt:Route.end_time);  
    x=Route.xfun(tt)-Route.xfun(Route.start_time); % route x samples
    y=Route.yfun(tt)-Route.yfun(Route.start_time); % route y samples 
    
  mov=0;
    plov % plot route (if mov==1 than plot movie to examine velocity)
 
    if state==0
        feedbackGain=0;
        Control.Pcon=1;
        Control.Icon=0;
        Control.Dcon=0;
    elseif state==1
 % Control Parameters
        feedbackGain=1;
        Control.Pcon=.5;
        Control.Icon=1000;
        Control.Dcon=0;
    end
        
% ------------------------------------------------------
tlen=numel(tt);
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
%     F=[x;y];
                                                    %      F=[0.03*tt;(tt>=0.3)];
% ---------------------------------------------------------

% -------------- B Field Direction Calculating ------------
% field direction depends on direction of velocity -> 1st derivative of route
    phiB=(atan2(2*bn(1,:).*tt+bn(2,:),2*an(1,:).*tt+an(2,:)));
%     phiB=pi/4*ones(size(phiB));
% ---------------------------------------------------------
    Bmag=BmW*ones(1,numel(tt));


disp('setting simulation...')
    simIn=timeseries([F;phiB],tt,'name','Predefined Route');
    % Initial Values
    initialValues.R0=[x(1);y(1)];
    initialValues.V0=[diff(x(1:2))/dt;diff(y(1:2))/dt];
    initialValues.phiB=phiB(1);

    save('par','simIn','-v7.3')

%     run specific_Compensator
% f=step(G0);
% gain=1/f(end);
 
%% Simulation:
 disp('starting simulation... this may take a while')

    tic
    set_param('ForceControl_closed_loop', 'StopTime', num2str(tt(end)))
    set_param('ForceControl_closed_loop', 'StartTime', num2str(tt(1)))
    sim('ForceControl_closed_loop')
    toc
    save('result')
 %% Plot Results
%  clear all
%  close all
%  clc
%  load('result')
 disp('plotting results, actuall model')
simsamp=floor(linspace(1,numel(simForce.time),size(simForce.data,1)/10));

figure(51)
    subplot(2,1,1)
        fig=plot(tt,F(1,:),simForce.time(simsamp),squeeze(simForce.data(simsamp,1)));
        set(fig(1),'lineWidth',3)
        legend('F wanted','F got')
        % axis([0 en -7*abs(max(F(1,:))) 7*abs(max(F(1,:)))])
        subplot(2,1,2)
        title 'Fx comparison'
        fig=plot(tt,F(2,:),simForce.time(simsamp),squeeze(simForce.data(simsamp,2)));
        set(fig(1),'lineWidth',3)
        legend('F wanted','F got')
        title 'Fy comparison'
figure(52)
    subplot(2,2,1)
        plot(simCurrent.time(simsamp),squeeze(simCurrent.data(simsamp,1)))
        title 'Ih_x'
        xlabel 't [sec]'
        ylabel '[A]'
    subplot(2,2,2)
        plot(simCurrent.time(simsamp),squeeze(simCurrent.data(simsamp,2)))
        title 'Ih_y'
        xlabel 't [sec]'
        ylabel '[A]'
    subplot(2,2,3)
        plot(simCurrent.time(simsamp),squeeze(simCurrent.data(simsamp,3)))
        title 'Im_x'
        xlabel 't [sec]'
        ylabel '[A]'
    subplot(2,2,4)
        plot(simCurrent.time(simsamp),squeeze(simCurrent.data(simsamp,4)))
        title 'Im_y'
        xlabel 't [sec]'
        ylabel '[A]'
figure(53)
    fig=plot(x,y,simDipoleLoc.data(simsamp,1),simDipoleLoc.data(simsamp,2));
    legend('Desired','Simulated')
    title 'Route Comparison'
    set(fig,'LineWidth',1.5)
    grid

 
 
% axis([-1 1 -1 1])
T=simDipoleLoc.time(simsamp);
xSim=simDipoleLoc.data(simsamp,1);
ySim=simDipoleLoc.data(simsamp,2);
phiSim=simPhiB.data(simsamp);
% plot(tt,phiB,'b*',simphiB.time(simsamp),simphiB.data(simsamp),'r')
figure(54)
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