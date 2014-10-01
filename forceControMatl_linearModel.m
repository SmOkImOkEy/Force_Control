% linear model:
disp('Calculating Linear Model..')

% reseting matrixes
parts=1;
padnum=parts-mod(numel(phiB),parts); % amount of zero padding to divide to parts
phiBl=[phiB phiB(end)*ones(1,padnum)]; % padded phi
Bmagl=[Bmag phiB(end)*ones(1,padnum)]; % padded Bm
Fl=[F zeros(2,padnum)]; % padded F
ttl=[tt tt(end)+dt*(1:padnum)]; % padded t

phiBmat=reshape(phiBl',[],parts)'; % each part on row phi
Bmmat=reshape(Bmagl',[],parts)'; % each part on row phi

Fmat=reshape(Fl,2,[],parts); % each part on row F
Fsp=zeros([size(Fmat,2),3,parts]); 
Tsp=zeros([size(Fmat,2),parts]);
ttmat=reshape(ttl,size(Fmat,2),[]);



disp('loading TF')
% loading general transfer function
    Gstruct=load('linMat','A','B','C');
    Amatrix=Gstruct.A;
    Bmatrix=Gstruct.B;
    Cmatrix=Gstruct.C;
disp('Starting loop:')    
for part=1:parts
    disp(['PART ',num2str(part)])
    disp('initializing')
        phiB_ss=mean(phiBmat(part,:)); 
        F0=Fmat(:,1,part);
        Bm0=Bmmat(part,1);
    % initial Currents:    
        Ihx0=Bm0/h*cos(phiB_ss);
        Ihy0=Bm0/h*sin(phiB_ss);
        Imx0=2/(3*M*kappa)*((2*F0(1))/cos(phiB_ss) +(F0(2))/sin(phiB_ss));
        Imy0=2/(3*M*kappa)*((F0(1))/cos(phiB_ss) +(2*F0(2))/sin(phiB_ss));
        
  % REGULAR PARAMETERS GIVES LOW CONDITION NUMBER OF Smt(needed to solve for compensator)
 % setting easy values to improve condition number of Sm
             Rx=1;Lx=1.3;
             Dipole.moment=1;
             Dipole.inertia=0.002;
             Ihx0=1;
             Ihy0=2;
             Imx0=-1;
             Imy0=-3;
             F0=[1.5;1];
             phiB_ss=pi/4;
             h=1;kappa=1;
             alpha=1;Bm0=2;
    % -------------------------------------------DEBUGING ONLY
    % substituting initial values    
    disp('substituting part parameters...')
    Amatrix=subs(Amatrix,{'Rx','Lx','Ry','Ly','h'},{Rx,Lx,Rx,Lx,h}); % substituting coil parameters
    Amatrix=subs(Amatrix,{'M','I','alpha'},{Dipole.moment,Dipole.inertia,alpha});
    Amatrix=double(subs(Amatrix,{'Ihx','Ihy','Imx','Imy'},{Ihx0,Ihy0,Imx0,Imy0}));

    Bmatrix=subs(Bmatrix,{'Rx','Lx','Ry','Ly','h','kappa'},{Rx,Lx,Rx,Lx,h,kappa}); % substituting coil parameters
    Bmatrix=subs(Bmatrix,{'Fx_hat','Fy_hat'},{F0(1),F0(2)}); % substituting F0
    Bmatrix=subs(Bmatrix,{'Bm_hat','phiB_hat'},{Bm0,phiB_ss}); % substituting initial B
    Bmatrix=double(subs(Bmatrix,'M',Dipole.moment));

    Cmatrix=subs(Cmatrix,'M',Dipole.moment,kappa);
    Cmatrix=subs(Cmatrix,{'Ihx','Ihy','Imx','Imy'},{Ihx0,Ihy0,Imx0,Imy0});
    Cmatrix=double(subs(Cmatrix,{'Rx','Lx','Ry','Ly','h','kappa'},{Rx,Lx,Rx,Lx,h,kappa})); % substituting coil parameters
    % converting to tf
    disp('Creating open loop system')
    StateSpace=(ss(Amatrix,Bmatrix,Cmatrix,0));
    Gtf=minreal(tf(StateSpace));
%     save('tfMat','Gtf');
% act 1: calculating left factors (non-coprime)
    disp('Calculate Left Fraction')
    [LN,LD]=left_poly_fractions(Gtf);
% act 2: Creating coefficient matrix Ngal,Dgal (not co-prime)
    disp('Coefficient Matrixes')
    [Ngal, Dgal]=Coefficient_Matrixes(LN,LD);
   % set tolerance function    
     tol=@(A,B) min(max([size(A) size(B)])*eps([norm3Dmat(A) norm3Dmat(B)]));

% act 3: Co-prime Factorization (pole zero cancelation)
    disp('Co-prime Factorization...')
    [Ngal_cp,Dgal_cp]=coprime_Factorization(Ngal,Dgal);
    [Ngal_cp,Dgal_cp]=TFSimplify(Ngal_cp,Dgal_cp,1e-10);
    [Ntf, Dtf]=create_tf(Ngal_cp,Dgal_cp); % needed only for debugging
% act 4: create compensator:
    disp('Compensator Design')
    psi=[pi pi-pi/4 pi+pi/4];
    
%     poles=[-1 -2-1j -2+1j];
    poles=0.01*exp(1j*psi);
    [Agal, Bgal, Ftf]=Compensator_design(Ngal_cp,Dgal_cp,poles);
    [Bgal,Agal]=TFSimplify(Bgal,Agal,1e-10);
    Agal=real(Agal);
    Bgal=real(Bgal);
    [Atf, Btf]=create_tf(Agal,Bgal); % needed only for debugging
% act 5: testing
    disp('Testing Compensator')
    [~,~,kf]=zpkdata(Ftf);
 [~, ~, k]=zpkdata(minreal(Atf*Dtf+Btf*Ntf-Ftf));
   if norm(k)<tol(kf,kf)
       disp('VERY GOOD')
   elseif norm(k)<tol(Ngal_cp,Dgal_cp)
       disp('GOOD')
   elseif norm(k)<1e-8;
       disp('Fine')
   end
    
% ----------------------Simulation:--------------------
disp('Simulating...')
subMat_len=size(Gtf,1);
    Ctf=minreal(Atf^-1*Btf); % Compensator
    G_ol_cp=Gtf;%Ntf_cp*Dtf_cp^-1; % Open loop system
    G_cl_cp=(eye(3)+G_ol_cp*Ctf)*G_ol_cp*Ctf; % Closed Loop System
    
    Den=(Atf*Dtf+Btf*Ntf)^-1; % dirty denominator
   % cleaning denominator:
      [z, p, k]=zpkdata(Den); 
      k(abs(k)<sqrt(subMat_len*eps(norm(k))))=0;
      Den=zpk(z,p,k);
   % ----
    
    Gcl1=minreal(Ntf*Den*Btf);
    Gclf=minreal(Ntf*Ftf^-1*Btf);
%      [y_imp, t_imp]=impulse(Gcl1) ; 
     [ystp_t,~]=step(Gcl1,ttl);
      ystp_t=real(ystp_t);
      amp=diag(diag(1./squeeze(ystp_t(end,:,:))));
      [ystp,tstp]=step(amp*Gcl1,ttl);
      ystp=real(ystp);

      clear ystp_t
   
     [ystpf,tstpf]=step(amp*Gclf,ttl);
      ystpf=real(ystpf);

    u=[Fmat(:,:,part); mod(phiBl+pi,2*pi)-pi].';
     [y_sim,t_sim]=lsim(amp*Gcl1,u,ttl);
     y_sim=real(y_sim);
     y_sim(:,3)=mod(y_sim(:,3)+pi,2*pi)-pi;
     


    % check if compensator gives F poles:
    nvec=(reshape(1:subMat_len^2,subMat_len,subMat_len)');
    n=0;
     for in=1:subMat_len
         for out=1:subMat_len
             n=n+1;
             figure(10)
             subplot(subMat_len,subMat_len,nvec(n))
             fig=plot(tstp,ystp(:,n),'-b',tstpf,ystpf(:,n),':r');
             grid on
             legend('Simulated','Desired')
             title (['Step input:',num2str(in),' output:',num2str(out)])
             xlabel 't [sec]'
             ylabel 'step response'
             set(fig,'linewidth',2)
            stepIn=stepinfo(ystp,tstp);
            stepIn.RiseTime
         end
     end
disp(['F error norm: ',num2str(norm(ystp(:)-ystpf(:))/norm(ystpf(:))*100),'%'])
     
     for n=1:subMat_len  
         figure(part)
           subplot(subMat_len, 1,n)
           plot(ttl,u(:,n),'--b',ttl,(y_sim(:,n)),'-k')
           legend('Input','Output(CL)')
        switch n
            case 1 
                title ('F_x')
                
            case 2 
                title ('F_y')
            case 3 
                title ('\phi_B')
        end
     end
end
% ----------------------------------------------------------