% 一次降噪
clc;clear;close all;
noise_name='马路';
fs=16000;
excel{1,1}='';excel{1,2}='snr';excel{1,3}='llr';excel{1,4}='lsd';excel{1,5}='pesq';
 
%% 
clean_path='clean.wav';
noisy_name=['noisy_snr0_',noise_name];
noisy_path=[noisy_name,'.wav'];
 
enhanced1_path=[noisy_name,'_enhanced1.wav'];
enhanced2_path=[noisy_name,'_enhanced2.wav'];
enhanced1_name=[noisy_name,'_enhanced1'];
enhanced2_name=[noisy_name,'_enhanced2'];
omlsa1(noisy_name,enhanced1_name);
omlsa7(noisy_name,enhanced2_name);
[clean,~]=audioread(clean_path);
[noisy,~]=audioread(noisy_path);
[enhanced1,~]=audioread(enhanced1_path);
[enhanced2,~]=audioread(enhanced2_path);
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
llr1=comp_llr(clean_path,noisy_path);
llr2=comp_llr(clean_path,enhanced1_path);
llr3=comp_llr(clean_path,enhanced2_path);
lsd1=LogSpectralDistance(clean,noisy,fs);
lsd2=LogSpectralDistance(clean,enhanced1,fs);
lsd3=LogSpectralDistance(clean,enhanced2,fs);
pesq1=pesq(clean_path,noisy_path);
pesq2=pesq(clean_path,enhanced1_path);
pesq3=pesq(clean_path,enhanced2_path);
 
i=2;
j=1;
excel{i,j}='——';j=j+1;
excel{i,j}=snr2-snr1;j=j+1;
excel{i,j}=llr1-llr2;j=j+1;
excel{i,j}=lsd1-lsd2;j=j+1;
excel{i,j}=pesq2-pesq1;
i=i+1;
j=1;
excel{i,j}='——';j=j+1;
excel{i,j}=snr3-snr1;j=j+1;
excel{i,j}=llr1-llr3;j=j+1;
excel{i,j}=lsd1-lsd3;j=j+1;
excel{i,j}=pesq3-pesq1;
 
%% 
clean_path='clean.wav';
noisy_name=['noisy_snr5_',noise_name];
noisy_path=[noisy_name,'.wav'];
 
enhanced1_path=[noisy_name,'_enhanced1.wav'];
enhanced2_path=[noisy_name,'_enhanced2.wav'];
enhanced1_name=[noisy_name,'_enhanced1'];
enhanced2_name=[noisy_name,'_enhanced2'];
omlsa1(noisy_name,enhanced1_name);
omlsa7(noisy_name,enhanced2_name);
[clean,~]=audioread(clean_path);
[noisy,~]=audioread(noisy_path);
[enhanced1,~]=audioread(enhanced1_path);
[enhanced2,~]=audioread(enhanced2_path);
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
llr1=comp_llr(clean_path,noisy_path);
llr2=comp_llr(clean_path,enhanced1_path);
llr3=comp_llr(clean_path,enhanced2_path);
lsd1=LogSpectralDistance(clean,noisy,fs);
lsd2=LogSpectralDistance(clean,enhanced1,fs);
lsd3=LogSpectralDistance(clean,enhanced2,fs);
pesq1=pesq(clean_path,noisy_path);
pesq2=pesq(clean_path,enhanced1_path);
pesq3=pesq(clean_path,enhanced2_path);
 
i=i+1;
j=1;
excel{i,j}='——';j=j+1;
excel{i,j}=snr2-snr1;j=j+1;
excel{i,j}=llr1-llr2;j=j+1;
excel{i,j}=lsd1-lsd2;j=j+1;
excel{i,j}=pesq2-pesq1;
i=i+1;
j=1;
excel{i,j}='——';j=j+1;
excel{i,j}=snr3-snr1;j=j+1;
excel{i,j}=llr1-llr3;j=j+1;
excel{i,j}=lsd1-lsd3;j=j+1;
excel{i,j}=pesq3-pesq1;
 
%% 
clean_path='clean.wav';
noisy_name=['noisy_snr10_',noise_name];
noisy_path=[noisy_name,'.wav'];
 
enhanced1_path=[noisy_name,'_enhanced1.wav'];
enhanced2_path=[noisy_name,'_enhanced2.wav'];
enhanced1_name=[noisy_name,'_enhanced1'];
enhanced2_name=[noisy_name,'_enhanced2'];
omlsa1(noisy_name,enhanced1_name);
omlsa7(noisy_name,enhanced2_name);
[clean,~]=audioread(clean_path);
[noisy,~]=audioread(noisy_path);
[enhanced1,~]=audioread(enhanced1_path);
[enhanced2,~]=audioread(enhanced2_path);
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
llr1=comp_llr(clean_path,noisy_path);
llr2=comp_llr(clean_path,enhanced1_path);
llr3=comp_llr(clean_path,enhanced2_path);
lsd1=LogSpectralDistance(clean,noisy,fs);
lsd2=LogSpectralDistance(clean,enhanced1,fs);
lsd3=LogSpectralDistance(clean,enhanced2,fs);
pesq1=pesq(clean_path,noisy_path);
pesq2=pesq(clean_path,enhanced1_path);
pesq3=pesq(clean_path,enhanced2_path);
 
i=i+1;
j=1;
excel{i,j}='——';j=j+1;
excel{i,j}=snr2-snr1;j=j+1;
excel{i,j}=llr1-llr2;j=j+1;
excel{i,j}=lsd1-lsd2;j=j+1;
excel{i,j}=pesq2-pesq1;
i=i+1;
j=1;
excel{i,j}='——';j=j+1;
excel{i,j}=snr3-snr1;j=j+1;
excel{i,j}=llr1-llr3;j=j+1;
excel{i,j}=lsd1-lsd3;j=j+1;
excel{i,j}=pesq3-pesq1;
 
%% 
clean_path='clean.wav';
noisy_name=['noisy_snr15_',noise_name];
noisy_path=[noisy_name,'.wav'];
 
enhanced1_path=[noisy_name,'_enhanced1.wav'];
enhanced2_path=[noisy_name,'_enhanced2.wav'];
enhanced1_name=[noisy_name,'_enhanced1'];
enhanced2_name=[noisy_name,'_enhanced2'];
omlsa1(noisy_name,enhanced1_name);
omlsa7(noisy_name,enhanced2_name);
[clean,~]=audioread(clean_path);
[noisy,~]=audioread(noisy_path);
[enhanced1,~]=audioread(enhanced1_path);
[enhanced2,~]=audioread(enhanced2_path);
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
llr1=comp_llr(clean_path,noisy_path);
llr2=comp_llr(clean_path,enhanced1_path);
llr3=comp_llr(clean_path,enhanced2_path);
lsd1=LogSpectralDistance(clean,noisy,fs);
lsd2=LogSpectralDistance(clean,enhanced1,fs);
lsd3=LogSpectralDistance(clean,enhanced2,fs);
pesq1=pesq(clean_path,noisy_path);
pesq2=pesq(clean_path,enhanced1_path);
pesq3=pesq(clean_path,enhanced2_path);
 
i=i+1;
j=1;
excel{i,j}='——';j=j+1;
excel{i,j}=snr2-snr1;j=j+1;
excel{i,j}=llr1-llr2;j=j+1;
excel{i,j}=lsd1-lsd2;j=j+1;
excel{i,j}=pesq2-pesq1;
i=i+1;
j=1;
excel{i,j}='——';j=j+1;
excel{i,j}=snr3-snr1;j=j+1;
excel{i,j}=llr1-llr3;j=j+1;
excel{i,j}=lsd1-lsd3;j=j+1;
excel{i,j}=pesq3-pesq1;
 
x=[0,5,10,15];
subplot(2,2,1);
plot(x,cell2mat(excel(2:2:8,2)),'k-');hold on;
plot(x,cell2mat(excel(3:2:9,2)),'r-');hold on;
plot(x,cell2mat(excel(2:2:8,2)),'kx');hold on;
plot(x,cell2mat(excel(3:2:9,2)),'rx');
xlabel('合成信噪比/dB');
ylabel('△SNR');
title('△SNR');
axis([0 15 min([cell2mat(excel(2:2:8,2));cell2mat(excel(3:2:9,2));0]) max([cell2mat(excel(2:2:8,2));cell2mat(excel(3:2:9,2))])])
 
subplot(2,2,2);
plot(x,cell2mat(excel(2:2:8,3)),'k-');hold on;
plot(x,cell2mat(excel(3:2:9,3)),'r-');hold on;
plot(x,cell2mat(excel(2:2:8,3)),'kx');hold on;
plot(x,cell2mat(excel(3:2:9,3)),'rx');
legend('IMCRA-OMLSA','CI\_IMCRA-OMLSA');
xlabel('合成信噪比/dB');
ylabel('△LLR');
title('△LLR');
axis([0 15 min([cell2mat(excel(2:2:8,3));cell2mat(excel(3:2:9,3));0]) max([cell2mat(excel(2:2:8,3));cell2mat(excel(3:2:9,3))])])
 
subplot(2,2,3);
plot(x,cell2mat(excel(2:2:8,4)),'k-');hold on;
plot(x,cell2mat(excel(3:2:9,4)),'r-');hold on;
plot(x,cell2mat(excel(2:2:8,4)),'kx');hold on;
plot(x,cell2mat(excel(3:2:9,4)),'rx');
xlabel('合成信噪比/dB');
ylabel('△LSD');
title('△LSD');
axis([0 15 min([cell2mat(excel(2:2:8,4));cell2mat(excel(3:2:9,4));0]) max([cell2mat(excel(2:2:8,4));cell2mat(excel(3:2:9,4))])])
 
subplot(2,2,4);
plot(x,cell2mat(excel(2:2:8,5)),'k-');hold on;
plot(x,cell2mat(excel(3:2:9,5)),'r-');hold on;
plot(x,cell2mat(excel(2:2:8,5)),'kx');hold on;
plot(x,cell2mat(excel(3:2:9,5)),'rx');
xlabel('合成信噪比/dB');
ylabel('△PESQ');
title('△PESQ');
axis([0 15 min([cell2mat(excel(2:2:8,5));cell2mat(excel(3:2:9,5));0]) max([cell2mat(excel(2:2:8,5));cell2mat(excel(3:2:9,5))])])

