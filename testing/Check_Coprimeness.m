% monitor Coprime error:
    [y1, t1]=step(Ntf*Dtf^-1,tt) ; 
%      [y2,t2]=step(LD^-1*LN,10);
     [y3,t3]=step(G_ol,tt);
     for n=1:subMat_len^2     
           subplot(subMat_len,subMat_len,n)
           fig=plot(t1,y1(:,n),'--b',t3,y3(:,n),':r');
           legend('cp','original')
           set(fig(1),'lineWidth',2)
           err(n)=norm(y1(:,n)-y3(:,n));
     end
     
     [znd,pnd,knd]=zpkdata(Ntf*Dtf^-1);
     [zgo,pgo,kgo]=zpkdata(G_ol);
