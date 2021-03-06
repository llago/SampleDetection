[original, fs_o] = audioread('D:\Faculdade\Mono\SD_Repository\Base Teste\PFNMF_boss_curto.wav');
[copy, fs_s] = audioread('D:\Faculdade\Mono\Gururani\Dataset\Copied\47.mp3');
[s, r] = size(original);

dur = fs_s*120;

original = mean(original, 2);
copy = mean(copy, 2);
copy = copy(1: dur);

[orig_size, ~] = size(original);
[copy_size, ~] = size(copy);

hop = 128;
sizeFrame = 512;
sizeSTFT = 4096;

Vst = imq.STFT(original, @imq.windows.HammingWindow, sizeFrame, hop, sizeSTFT);
V = abs(Vst);
Xst = imq.STFT(copy, @imq.windows.HammingWindow, sizeFrame, hop, sizeSTFT);
X = abs(Xst);

% 
% Xo = spectrogram(original, window, window-hop);
% Xs = spectrogram(sample, window, window-hop);

% [M, N] = size(V);
% min = min(size(Xo));Ho_hypo_gm

k = 5;
L = 10;

[Wo, Ho] = nnmf(V, k);
[~, Ho_hypo, ~, ~, err] = PfNmf(X, Wo, [], [], [], L, 0);

for i = 1:k
    Ho(i,:) = Ho(i,:)./max(Ho(i,:),1);
end

for i = 1:k
    Ho_hypo(i,:) = Ho_hypo(i,:)./max(Ho_hypo(i,:),1);
end

Ho_gm = geomean(Ho,1);
Ho_hypo_gm = geomean(Ho_hypo,1);

% Ho_gm = Ho_gm./norm(Ho_gm,2);
% Ho_hypo_gm = Ho_hypo_gm./norm(Ho_hypo_gm,2);

cor = xcorr(Ho_gm, Ho_hypo_gm);
figure;
plot(cor)

% save('D:\Faculdade\Mono\SD_Repository\Resultados Teste\PFNMF_boss_orig_curto.mat', 'Wo', 'Ho', 'Ho_hypo', 'Ho_gm', 'Ho_hypo_gm', 'cor', 'copy_size', 'orig_size')