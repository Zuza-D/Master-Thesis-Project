a%% Reading audio file
clear all; close all; clc;
% path = 'C:\Users\zuzanna.Dlugosz\OneDrive - Aalto University\Thesis\scripts\sample to render\'; % define the working folder
path = 'C:\Users\zuzanna.Dlugosz\OneDrive - Aalto University\Thesis\Old things\sample to render\';
pattern =  fullfile(path, '*.wav'); % create a patter for .wav files
audiofiles = dir(pattern); % all .wav files in folder
for k = 1 : length(audiofiles)
  basename = audiofiles(k).name;
  fullname = fullfile(path, basename);
  [S{k},Fs] = audioread(fullname); % store the files as struct
end
figure(1)
plot(S{1})
grid on
%% Equalize the loudness before
target = -23;
for k = 1 : length(audiofiles)
    [ac_loudness{k}, specificLoudness{k}] = acousticLoudness(S{k},Fs);
    ac_loudness{k} = round(ac_loudness{k},0);
    [loudness{k}] = integratedLoudness(S{k},Fs);    
    gaindB(k) = target - double(loudness{k});
    gain(k) = 10.^(gaindB(k)/20);
    S{k} = S{k}.*gain(k);
    ac_loudness_new{k} = acousticLoudness(S{k},Fs);
    ac_loudness_new{k} = round(ac_loudness_new{k},0);
    loudness_new{k} = integratedLoudness(S{k},Fs);
end
% 
% figure(1)
% plot(S{1})
% grid on
figure(2)
plot(S{1})
grid on


%% Equalize the RMS before
% RMS = [];
% for k = 1 : length(audiofiles)
%     RMS = [RMS rms(S{k})];
% end
% 
% RMS_ONE = min(RMS);
% RMS_EQ = [];
% % change loudness of the files
% for k = 1 : length(audiofiles)
%     coeff(k) = RMS_ONE ./ RMS(k);
%     S{k} = S{k}.* coeff(k);
%     RMS_EQ = [RMS_EQ rms(S{k})];
% end    
%% Auralization
Folder = 'C:\Users\zuzanna.Dlugosz\OneDrive - Aalto University\Thesis\SRIR\';
names = dir([Folder '*.mat']);
SDM = cell(1,length(names)/4);
NORMALIZATION = 0;
MINUS = 1; %(1 - default array ; 0 - changed array(with -))
for audio = 1 : length(S)
    for j = 1 : length(SDM)
%         if NORMALIZATION == 1
%             SDM{j} = load([Folder names(j + 3).name]);
%         else
%             SDM{j} = load([Folder names(j + 3).name]);
%         end
        if NORMALIZATION == 1
          if MINUS == 1
              SDM{j} = load([Folder names(j + 3).name]);
          else
              SDM{j} = load([Folder names(j + 9).name]);
          end
        else
            if MINUS == 1
                SDM{j} = load([Folder names(j).name]);
            else
                SDM{j} = load([Folder names(j + 6).name]);
            end
        end
    
%         SDM{j} = load([Folder names(j).name]);
        S{audio} = [zeros(length(SDM{j}.H),1); S{audio} ; zeros(length(SDM{j}.H),1)]; % zero-padding the signal
    
        numOfLsp = size(SDM{j}.H,2); % number of loudspeakers 
        Y = zeros(size(S{audio},1),numOfLsp); % prepare array for each audio file
    

        for lsp = 1:numOfLsp
            % Convolution with Matlab's overlap-add
            Y(:,lsp) = Y(:,lsp) +  fftfilt(SDM{j}.H(:,lsp),S{audio});
        end
    
        % Y contains the auralization of the spatial IRs with S
        %% Auralization 

        disp('Started Auralization');

        % define the names of the files and the destination folder
        if NORMALIZATION == 1
            if MINUS == 1
                savefolder = 'C:\Users\zuzanna.Dlugosz\OneDrive - Aalto University\Thesis\scripts\sample to render\rendered\NORM\-GRAVI25\';
                savename = [savefolder names(j + 3).name(16 : 17) '_' audiofiles(audio).name];
            else
                savefolder = 'C:\Users\zuzanna.Dlugosz\OneDrive - Aalto University\Thesis\scripts\sample to render\rendered\NORM\GRAVI25\';
                savename = [savefolder names(j + 9).name(10 : 11) '_' audiofiles(audio).name];
            end
            
        else
            if MINUS == 1
%                 savefolder = 'C:\Users\zuzanna.Dlugosz\OneDrive - Aalto University\Thesis\scripts\sample to render\rendered\-GRAVI25\';
                  savefolder = 'C:\Users\zuzanna.Dlugosz\OneDrive - Aalto University\Thesis\scripts\MyTests\test_samples\test\-GRAVI25\';
                  savename = [savefolder names(j).name(11 : 12) '_' audiofiles(audio).name];
            else
%                 savefolder = 'C:\Users\zuzanna.Dlugosz\OneDrive - Aalto University\Thesis\scripts\sample to render\rendered\GRAVI25\';   
                  savefolder = 'C:\Users\zuzanna.Dlugosz\OneDrive - Aalto University\Thesis\scripts\MyTests\test_samples\test\';
                  savename = [savefolder names(j + 6).name(5 : 6) '_' audiofiles(audio).name];
            end
        end
        
        

        % Sound normalization (if clicks)
        Y = Y/max(abs(Y(:)))*.9;
        disp('Sound normalized, since otherwise would have clipped')

        % Saving aurialized sample as .wav
        audiowrite(savename,Y,Fs)
        info = audioinfo(savename);
    end
end
%% Checking the loudness
% clear all
% savefolder = 'C:\Users\zuzanna.Dlugosz\OneDrive - Aalto University\Thesis\scripts\sample to render\rendered\-GRAVI25\';
% pattern =  fullfile(savefolder, '*.wav'); % create a patter for .wav files
% audiofiles = dir(pattern); % all .wav files in folder
% for k = 1 : length(audiofiles) 
%   basename = audiofiles(k).name;
%   fullname = fullfile(savefolder, basename);
%   [A{k},Fs] = audioread(fullname); % store the files as struct
%     
%     for i = 1 : length(A{k})
%         mono{k}(i) = sum(A{k}(i,:));
%     end
%     ac_loudness_after{k} = acousticLoudness(mono{k}',Fs);
%     ac_loudness_after{k} = round(ac_loudness_after{k},0);
%     loudness_after{k} = integratedLoudness(mono{k}',Fs);
% end