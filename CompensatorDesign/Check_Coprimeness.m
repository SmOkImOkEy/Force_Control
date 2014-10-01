% monitor Coprime error:
    [y1, t1]=step(Ntf*Dtf^-1,10) ; 
     [y2,t2]=step(LD^-1*LN,10);
     [y3,t3]=step(Gtf,10);
     for n=1:subMat_len^2     
           subplot(subMat_len,subMat_len,n)
           fig=plot(t1,y1(:,n),'--b',t2,y2(:,n),'-k',t3,y3(:,n),':r');
           legend('cp','ncp','original')
           set(fig(1),'lineWidth',2)
     end