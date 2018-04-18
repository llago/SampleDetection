[original, fs_o] = audioread('D:\Faculdade\Mono\SD_Repository\Base Teste\PFNMF_o4_puro.wav');
[sample, fs_s] = audioread('D:\Faculdade\Mono\SD_Repository\Base Teste\PFNMF_o4_s32.wav');

% window = 4096; % Play around with this
% hop = 1024;

original = mean(original, 2);
sample = mean(sample, 2);

hop = 256;
sizeFrame = 512;
sizeSTFT = 2048;

Vst = imq.STFT(original, @imq.windows.HammingWindow, sizeFrame, hop, sizeSTFT);
V = abs(Vst);
Xst = imq.STFT(sample, @imq.windows.HammingWindow, sizeFrame, hop, sizeSTFT);
X = abs(Xst);

% 
% Xo = spectrogram(original, window, window-hop);
% Xs = spectrogram(sample, window, window-hop);

% [M, N] = size(V);
% min = min(size(Xo));

k = 5;

[Wo, Ho] = nnmf(V, k);
[~, Ho_hypo, ~, ~, err] = PfNmf(X, Wo, [], [], [], 0, 0);

for i = 1:k
    Ho(i,:) = Ho(i,:)/norm(Ho(i,:),1);
end

for i = 1:k
    Ho_hypo(i,:) = Ho_hypo(i,:)/norm(Ho_hypo(i,:),1);
end

Ho_gm = geomean(Ho,1);
Ho_hypo_gm = geomean(Ho_hypo,1);

% Ho_gm = Ho_gm/norm(Ho_gm,2);
% Ho_hypo_gm = Ho_hypo_gm/norm(Ho_hypo_gm,2);

cor = xcorr(Ho_gm, Ho_hypo_gm);
figure;
plot(cor)

% save('D:\Faculdade\Mono\SD_Repository\Resultados Teste\PFNMF_o4_s32_spec_1.mat', 'Wo', 'Ho', 'Ho_hypo', 'Ho_gm', 'Ho_hypo_gm')