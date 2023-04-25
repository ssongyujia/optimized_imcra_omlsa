clc;clear;close all;

snr = 15;

clean_name = 'clean';
noise_name = '油烟机';
file.clean = [clean_name '.wav']; 
file.noise = [noise_name '.wav'];
noisy_name = ['noisy_snr' num2str(snr) '_' noise_name];
file.noisy = [noisy_name '.wav'];      
addnoise_asl(file.clean, file.noise, file.noisy , snr);    
% 
% clean_name = 'clean';
% noise_name = '排气扇';
% file.clean = [clean_name '.wav']; 
% file.noise = [noise_name '.wav'];
% noisy_name = ['noisy_snr' num2str(snr) '_' noise_name];
% file.noisy = [noisy_name '.wav'];      
% addnoise_asl(file.clean, file.noise, file.noisy , snr);
% 
% clean_name = 'clean';
% noise_name = '工厂';
% file.clean = [clean_name '.wav']; 
% file.noise = [noise_name '.wav'];
% noisy_name = ['noisy_snr' num2str(snr) '_' noise_name];
% file.noisy = [noisy_name '.wav'];      
% addnoise_asl(file.clean, file.noise, file.noisy , snr);
% 
% clean_name = 'clean';
% noise_name = '马路';
% file.clean = [clean_name '.wav']; 
% file.noise = [noise_name '.wav'];
% noisy_name = ['noisy_snr' num2str(snr) '_' noise_name];
% file.noisy = [noisy_name '.wav'];      
% addnoise_asl(file.clean, file.noise, file.noisy , snr);

% clean_name = 'clean';
% noise_name = '键盘';
% file.clean = [clean_name '.wav']; 
% file.noise = [noise_name '.wav'];
% noisy_name = ['noisy_snr' num2str(snr) '_' noise_name];
% file.noisy = [noisy_name '.wav'];      
% addnoise_asl(file.clean, file.noise, file.noisy , snr);
% 
% clean_name = 'clean';
% noise_name = '风扇';
% file.clean = [clean_name '.wav']; 
% file.noise = [noise_name '.wav'];
% noisy_name = ['noisy_snr' num2str(snr) '_' noise_name];
% file.noisy = [noisy_name '.wav'];      
% addnoise_asl(file.clean, file.noise, file.noisy , snr);