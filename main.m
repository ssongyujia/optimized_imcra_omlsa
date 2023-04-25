clear;close all;clc;
noisy_file='noisy_snr0_马路';
enhanced1_file='noisy_snr0_马路_enhanced1';
enhanced2_file='noisy_snr0_马路_enhanced2';
tstart1=tic;
omlsa1(noisy_file,enhanced1_file);
tend1=toc(tstart1)
tstart2=tic;
omlsa7(noisy_file,enhanced2_file);
tend2=toc(tstart2)
 
clean_path='clean.wav';
noisy_path=[noisy_file,'.wav'];
enhanced1_path=[enhanced1_file,'.wav'];
enhanced2_path=[enhanced2_file,'.wav'];
[clean,~]=audioread(clean_path);
[noisy,~]=audioread(noisy_path);
[enhanced1,~]=audioread(enhanced1_path);
[enhanced2,fs]=audioread(enhanced2_path);
N1=length(clean);
N2=length(noisy);
N3=length(enhanced1);
N4=length(enhanced2);
N=min([N1,N2,N3,N4]);
clean=clean(1:N);
noisy=noisy(1:N);
enhanced1=enhanced1(1:N);
enhanced2=enhanced2(1:N);
 
snr1=SNR_singlech(clean,noisy);
snr2=SNR_singlech(clean,enhanced1);
snr3=SNR_singlech(clean,enhanced2);
ssnr1=segsnr(clean,noisy,fs);
ssnr2=segsnr(clean,enhanced1,fs);
ssnr3=segsnr(clean,enhanced2,fs);
llr1=comp_llr(clean_path,noisy_path);
llr2=comp_llr(clean_path,enhanced1_path);
llr3=comp_llr(clean_path,enhanced2_path);
lsd1=LogSpectralDistance(clean,noisy,fs);
lsd2=LogSpectralDistance(clean,enhanced1,fs);
lsd3=LogSpectralDistance(clean,enhanced2,fs);
pesq1=pesq(clean_path,noisy_path);
pesq2=pesq(clean_path,enhanced1_path);
pesq3=pesq(clean_path,enhanced2_path);
wss1=comp_wss(clean_path,noisy_path);
wss2=comp_wss(clean_path,enhanced1_path);
wss3=comp_wss(clean_path,enhanced2_path);
 
excel{2,1}='——';excel{3,1}='omlsa';excel{4,1}='omlsa1';
excel{1,2}='snr';excel{2,2}=snr1;excel{3,2}=snr2;excel{4,2}=snr3;
excel{1,3}='ssnr';excel{2,3}=ssnr1;excel{3,3}=ssnr2;excel{4,3}=ssnr3;
excel{1,4}='llr';excel{2,4}=llr1;excel{3,4}=llr2;excel{4,4}=llr3;
excel{1,5}='lsd';excel{2,5}=lsd1;excel{3,5}=lsd2;excel{4,5}=lsd3;
excel{1,6}='pesq';excel{2,6}=pesq1;excel{3,6}=pesq2;excel{4,6}=pesq3;
excel{1,7}='wss';excel{2,7}=wss1;excel{3,7}=wss2;excel{4,7}=wss3;
 
x=0:1:5;
x1={'snr','ssnr','llr','lsd','pesq','wss'};
plot(x,cell2mat(excel(3,2:7))-cell2mat(excel(2,2:7)),'bx');hold on;plot(x,cell2mat(excel(3,2:7))-cell2mat(excel(2,2:7)),'b-');hold on;
plot(x,cell2mat(excel(4,2:7))-cell2mat(excel(2,2:7)),'rx');hold on;plot(x,cell2mat(excel(4,2:7))-cell2mat(excel(2,2:7)),'r-');
set(gca,'xtick',x);
set(gca,'xticklabel',x1);
legend(' ','omlsa',' ','omlsa1');
