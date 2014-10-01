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

% substituting constants:
    Amatrix=subs(Amatrix,{'Rx','Lx','Ry','Ly','h'},{Rx,Lx,Rx,Lx,h}); % substituting coil parameters
    Amatrix=subs(Amatrix,{'M','I','alpha'},{Dipole.moment,Dipole.inertia,alpha});
    Bmatrix=subs(Bmatrix,{'M','Rx','Lx','Ry','Ly','h','kappa'},{Dipole.moment,Rx,Lx,Rx,Lx,h,kappa}); % substituting coil parameters
    Cmatrix=subs(Cmatrix,{'Rx','Lx','Ry','Ly','h','kappa'},{Rx,Lx,Rx,Lx,h,kappa}); % substituting coil parameters
% ----------------------

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
%              Rx=1;Lx=1.3;
%              Dipole.moment=1;
%              Dipole.inertia=0.002;
%              Ihx0=1;
%              Ihy0=2;
%              Imx0=-1;
%              Imy0=-3;
%              F0=[1.5;1];
%              phiB_ss=pi/4;
%              h=1;kappa=1;
%              alpha=1;Bm0=2;
    % -------------------------------------------DEBUGING ONLY
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

 % act 1: calculating left factors (non-coprime)
    disp('Calculate Left Fraction')
    [LN,LD]=left_poly_fractions(G_ol);
% act 2: Creating coefficient matrix Ngal,Dgal (not co-prime)
    disp('Coefficient Matrixes')
    [Ngal, Dgal]=Coefficient_Matrixes(LN,LD);
% act 3: Co-prime Factorization (pole zero cancelation)
    disp('Co-prime Factorization...')
    [Ngal_cp,Dgal_cp]=coprime_Factorization(Ngal,Dgal);
    [Ngal_cp,Dgal_cp]=TFSimplify(Ngal_cp,Dgal_cp,1e-10);
    [Ntf, Dtf]=create_tf(Ngal_cp,Dgal_cp); % needed only for debugging
% act 4: create compensator:
    disp('Compensator Design')
    psi=[pi pi-pi/4 pi+pi/4];
    poles=0.01*exp(1j*psi);
    [Agal, Bgal, Ftf]=Compensator_design(Ngal_cp,Dgal_cp,poles);
    [Bgal,Agal]=TFSimplify(Bgal,Agal,1e-10);
    Agal=real(Agal);
    Bgal=real(Bgal);
    [Atf, Btf]=create_tf(Agal,Bgal); % needed only for debugging
    
% ----------------------Simulation:--------------------
disp('Simulating...')
subMat_len=size(G_ol,1);
    Ctf=minreal(Atf^-1*Btf); % Compensator
    
    Den=(Atf*Dtf+Btf*Ntf)^-1; % dirty denominator
   % cleaning denominator:
      [z, p, k]=zpkdata(Den); 
      k(abs(k)<sqrt(subMat_len*eps(norm(k))))=0;
      Den=zpk(z,p,k);
   % ----
    
    G_cl=minreal(Ntf*Den*Btf); % closed loop - without amp
    G_clf=minreal(Ntf*Ftf^-1*Btf); % imaged cl with poles we want
    % simulating Closed Loop to find AMP that get unity response
    [ystp_cl, ~]=step(G_cl,ttl);
    amp_cl=diag(diag(1./squeeze(real(ystp_cl(end,:,:)))));
    % ----------------------
% system with F poles only: (debugginh only)
    % simulating Imaged Closed Loop to find AMP that get unity response
    [ystpf_cl, ~]=step(G_clf,ttl);
    amp_clf=diag(diag(1./squeeze(real(ystpf_cl(end,:,:)))));
    % ----------------------
% ----------------
    % simulating Open Loop to find AMP that get unity response
      [ystp_ol, ~]=step(G_ol,ttl);
      amp_ol=diag(diag(1./squeeze(real(ystp_ol(end,:,:)))));
    % ----------------------
    % multiplying systems by current AMP
      G_cl=amp_cl*G_cl; % closed loop
      G_clf=amp_clf*G_clf; % imaged closed loop
      G_ol=amp_ol*G_ol; % open loop
    % ---------
 % running simulation with input:
    % set input vector:
       u=[Fmat(:,:,part); mod(phiBl+pi,2*pi)-pi].';
    % simulating closed and open loop:
        % closed loop:
     [y_sim_cl,~]=lsim(G_cl,u,ttl);
      y_sim_cl=real(y_sim_cl);
      y_sim_cl(:,3)=mod(y_sim_cl(:,3)+pi,2*pi)-pi; % setting limits of phi fo plotting
        % open loop:
      [y_sim_ol,~]=lsim(G_ol,u,ttl);
      y_sim_ol=real(y_sim_ol);
      y_sim_ol(:,3)=mod(y_sim_ol(:,3)+pi,2*pi)-pi; % setting limits of phi fo plotting
   % --------------------------------------------


    % check step response poled of G is poles of F:
        nvec=(reshape(1:subMat_len^2,subMat_len,subMat_len)');
        n=0;
%          for in=1:subMat_len
%              for out=1:subMat_len
%                  n=n+1;
%                  figure(10)
%                  subplot(subMat_len,subMat_len,nvec(n))
%                  fig=plot(ttl,ystp_cl(:,n),'-b',ttl,ystpf_cl(:,n),':r');
%                  grid on
%                  legend('Simulated','Desired')
%                  title (['Step input:',num2str(in),' output:',num2str(out)])
%                  xlabel 't [sec]'
%                  ylabel 'step response'
%                  set(fig,'linewidth',2)
% %                 stepIn=stepinfo(ystp,tstp);
% %                 stepIn.RiseTime
%              end
%          end 
%          disp(['F error norm: ',num2str(norm(ystp_cl(:)-ystpf_cl(:))/norm(ystpf_cl(:))*100),'%'])
    % -------------------------
   % comparison between Open and Closed loop:
     for n=1:subMat_len       
         switch n
            case 1 
                strTit='F_x'; 
            case 2 
                strTit='F_y';
            case 3 
                strTit='\phi_B';
         end   
         figure(part)
         % plot CL results
           subplot(subMat_len, 2,2*n-1)
           plot(ttl,u(:,n),'--b',ttl,y_sim_cl(:,n),'-k')
           legend('Input','Output(CL)')
           title (['Closed Loop ',strTit])
           
         % plot OL results
           subplot(subMat_len, 2,2*n)
           plot(ttl,u(:,n),'--b',ttl,(y_sim_ol(:,n)),'-r')
           legend('Input','Output(OL)')
           title (['Open Loop ',strTit])  
     end
end
% ----------------------------------------------------------