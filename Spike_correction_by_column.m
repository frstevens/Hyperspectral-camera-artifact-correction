
function Spike_correction_by_column(image_path)
% INPUT
% image_path: image name without the extension .mat.
% EXPLANATION
% The function treats the problems of artifacts, also called spikes, that
% may appear hyperspectral images measured with line scan type (pushbroom)
% cameras. The function searches for abnormal values at single wavelengths.
% It exploits the fact that these problems generally show vertical
% patterns: when a pixel is impacted, usually all the pixels of the
% column of the image, or a large part of it, are also impacted. 

%% Load image
im = load([image_path '.mat']);
data = im.spc.data; % The name of the fields (spc, data, ...) may vary for different images. Adapt the script if necessary.

%% Computation of the LPI (local peak index) 
% The matrix of LPI has the same dimension and size as the hypercube. For
% each pixel (i, j) and each wavelength k, LPI is equal to the double of 
% this wavelength minus the sum of the adjacent spectral neighbours :
% LPI(i, j, k) = 2*im(i, j, k) - im(i, j, k-1) - im(i, j, k+1)
nBand = size(data, 3);
LPI = zeros(size(data));
diffWl = diff(data, 1, 3); % First order difference along dimension 3 (spectral dimension)
LPI(:, :, 2:nBand-1) = LPI(:, :, 2:nBand-1) + diffWl(:, :, 1:end-1); % LPI computation (efficient)
LPI(:, :, 2:nBand-1) = LPI(:, :, 2:nBand-1) - diffWl(:, :, 2:end); % LPI computation (efficient)

%% Compuation of LPI statistics
% mean and median by column, of size nRow x nBand
medLPI = squeeze(median(LPI, 1));
meanLPI = squeeze(mean(LPI, 1));
clear diffWl LPI

%% Computation of statistics based on the values in "im"
% For each column and each wavelength, we compute
% 
% 1/ Two times the mean of the column, minus the sum of
% the mean of the left column and the mean of the right column 
meanIm = squeeze(mean(data, 1));
diffIm = diff(meanIm, 1, 1);
diffLR = zeros(size(meanIm));
diffLR(2:end, :) = diffLR(2:end, :) + diffIm;
diffLR(1:end-1, :) = diffLR(1:end-1, :) - diffIm;
clear meanIm diffIm

% 2/ Two times the standard deviation of the column, minus the sum of
% the standard deviation of the left column and the standard deviation of the right column 
stdIm = squeeze(std(data, [], 1));
diffImStd = diff(stdIm, 1, 1);
diffLRStd = zeros(size(stdIm));
diffLRStd(2:end, :) = diffLRStd(2:end, :) + diffImStd;
diffLRStd(1:end-1, :) = diffLRStd(1:end-1, :) - diffImStd;
clear stdIm diffImStd

%% Selection of indicator and computation + thresholding of Mahananobis distance
X = [meanLPI(:) medLPI(:) diffLR(:) diffLRStd(:)]; % we combine the four indicators of potential anomaly. Each row of X corresponds to a combination column + wavelength.
D = sqrt(mahal(X, X)); % the Mahalanobis distance accounts for the difference of variance of the four indicators and the covariance between them
[~, order] = sort(D, 'descend'); % sorting of the cases from the most abnormal to the most normal
nProb = sum(D > mean(D)+3*std(D));
[colError, bandError] = ind2sub(size(data, [2 3]), order); % retrienving of the indices of abnormal combinations column + wavelength 
colError = colError(1:nProb);
bandError = bandError(1:nProb);

%% Correction by spline interpolation of the spectrum
dataCorr = data;
for col = unique(colError)'
    NOKband = bandError(colError==col);
    OKband = setdiff(1:nBand, NOKband);
    dataCorr(:, col, NOKband) = interp1(OKband, squeeze(data(:, col, OKband))' , NOKband, 'linear')';
end

%% Export correction
im.spc.data = dataCorr;
image_path_corr = [image_path '_CORR'];
spc = im.spc;
save([image_path_corr '.mat'], 'spc')

%% plot
% figure
% tiledlayout(3, 2)
% hold on
% nexttile
% plot(D(order))
% xline(nProb)
% nexttile
% hist(D(:), 100)
% xline(D(order(nProb)))
% nexttile
% hist(meanLPI(:), 100)
% nexttile
% hist(medLPI(:), 100)
% nexttile
% hist(diffLR(:), 100)
% nexttile
% hist(diffLRStd(:), 100)

end

