function [XK,YK,ZK,nframes]=omlsa6(noise_ori,noise_est1,noise_est2)

nwin=512;		
inc=0.25*nwin;
win=hamming(nwin);
M21=nwin/2+1;
  
xx=enframe(noise_ori,nwin,inc);
yy=enframe(noise_est1,nwin,inc);
zz=enframe(noise_est2,nwin,inc);
nframes=size(xx,1);
XK=zeros(1,nframes);
YK=zeros(1,nframes);
ZK=zeros(1,nframes);
for i=1:nframes
    x=(xx(i,:))';
    y=(yy(i,:))';
    z=(zz(i,:))';
    X=fft(win.*x);
    Y=fft(win.*y);
    Z=fft(win.*z);
    X2=abs(X(1:M21)).^2;
    Y2=abs(Y(1:M21)).^2;
    Z2=abs(Z(1:M21)).^2;
    XK(i)=X2(32);
    YK(i)=Y2(32);
    ZK(i)=Z2(32);
end

    
        