function Summary = fla_calcium6d2(directory, fps)

% Flagellar GCaMP analysis (190327)
% Updated on 12/29/2021
% Updated on 5/17/2022 changed calculation methods


%% Parameters for kymograph analysis

% px_size = 0.104667; % pixel size in (um/px)
% fps = 19.33; %frame per second
% threshold = 0.9; % between 0 and 1, larger is severe?

%% Read kymograph files

directory01 = dir(directory);
% directory01 = dir;
val1 = length(directory01);

for m01 = 1:val1
    if regexp(directory01(m01).name, 'kymograph\d.tif') > 0
        filename01 = directory01(m01).name;
        image_ca01 = im2double(imread(filename01));
        %imsho(image_km01), title(filename01);
    elseif regexp(directory01(m01).name, 'kymograph_\d.labels.png') > 0
        image_label01 = imread(directory01(m01).name);
    end
end

% [~,name1,~] = fileparts(directory01(1).folder);
% figure, imshow(imrotate(image_ca01, 90), []), title('GCaMP');

clear val1 m01 filename01

%% Background correction

siz_im01 = size(image_ca01);
image_label02 = imbinarize(image_label01);
flaRegion1 = sum(image_label02,2);
flaRegion2 = flaRegion1;
flaRegion2(flaRegion1==0) = NaN;

image_caA4 = imgaussfilt(image_ca01, 1, 'FilterSize', 3);
% image_caA5 = image_caA4 .* image_label02;
% image_caA5(image_caA5==0) = NaN;

% SE_line = strel('line', 11,0);
% image_labelB1 = imdilate(image_label02, SE_line);
% image_labelB2 = imcomplement(image_labelB1);
% image_labelC1 = false(siz_im01);
% image_labelC1(:, round(siz_im01(2)/2):siz_im01(2)-6) = 1;
% image_labelD1 = image_labelB2 & image_labelC1; % use the background of the tip side
% % image_labelC2 = false(siz_im01);
% % image_labelC2(:, 7:round(siz_im01(2)/2)) = 1;
% % image_labelD2 = image_labelB2 & image_labelC2;
% 
% back1 = sum(image_labelD1(:), 'omitnan');
% 
% if back1 == 0
%     image_labelC2 = false(siz_im01);
%     image_labelC2(:, round(siz_im01(2)/2):siz_im01(2)-1) = 1;
%     image_labelD2 = image_labelB2 & image_labelC2;
%     image_caB1 = image_caA4 .* image_labelD2;
% else
%     image_caB1 = image_caA4 .* image_labelD1;
% end
% 
% image_caB1(image_caB1==0) = NaN;
% gcampBackground1 = mean(image_caB1(:), 'omitnan');

min_imcaA4 = min(image_caA4(:));
gcampBackground1 = min_imcaA4;

median_length1 = median(sum(image_label02, 2)); % the length of flagella which was used for analysis
time1 = siz_im01(1);

mean_imcaA4 = mean(image_caA4(:), 'omitnan');

if mean_imcaA4 < gcampBackground1
    min_imcaA4 = min(image_caA4(:));
    image_ca02 = image_caA4 - min_imcaA4;
else
    image_ca02 = image_caA4 - gcampBackground1; % background corrected image
end

image_ca03 = image_ca02 .* image_label02; %exclude non-flagella position
image_ca03(image_ca03==0) = NaN;

% min_imca02 = min(image_ca02(:));
% mean_imca03 = mean(image_ca03(:), 'omitnan');
% figure, imshow(imrotate(image_ca03,90),[]), title('GCaMP analysis');



%% Baseline correction

gcampIntensity2 = mean(image_ca03, 2, 'omitnan');
gcampIntensity3a = gcampIntensity2(1:100); % first 100 pixels
gcampIntensity3b = gcampIntensity2(end-99:end); % last 100 pixels

mean_gcampIntensity2 = mean(gcampIntensity2);
mean_gcampIntensity3a = mean(gcampIntensity3a);
std_gcampIntensity3a = std(gcampIntensity3a);
mean_gcampIntensity3b = mean(gcampIntensity3b);
std_gcampIntensity3b = std(gcampIntensity3b);

[pks1a, ~] = findpeaks(-gcampIntensity3a, 'MinPeakHeight', -mean_gcampIntensity3a-std_gcampIntensity3a*0.5); %, 'MinPeakProminence', std_gcampIntensity3a*0.1);
[pks1b, ~] = findpeaks(-gcampIntensity3b, 'MinPeakHeight', -mean_gcampIntensity3b-std_gcampIntensity3b*0.5); %, 'MinPeakProminence', std_gcampIntensity3b*0.1);
[pks1c, loc1c] = findpeaks(-gcampIntensity2, 'MinPeakHeight', -mean_gcampIntensity2*1.2);


if length(pks1a) <= 2 || length(pks1b) <= 2
    base1 = 0;
    dataA3b = gcampIntensity2;
else
    baseFirst1 = mean(pks1a);
    baseLast1 = mean(pks1b);
    baseMean1 = mean(pks1c);
    base1 = (baseFirst1 - baseLast1) / baseMean1;
    
    [p1,s1,mu1] = polyfit(loc1c,abs(pks1c),2);
    f_y = polyval(p1,(1:numel(gcampIntensity2))',[],mu1);
    f_y2 = f_y - abs(baseFirst1);
    dataA3b = gcampIntensity2 - f_y2;
end

if base1 > 0    
    image_base01 = ones(siz_im01);
    image_base01 = image_base01 .* (base1 + 1);
    base2 = (1:base1/300:301)';
    base2 = base2(1:301,:) * ones(1,siz_im01(2));
    image_base01(1:301,:) = base2;
elseif baseFirst1 < 0 || baseLast1 < 0 || baseMean1 < 0
    image_base01 = ones(siz_im01);
    
elseif base1 < -1
    image_base01 = ones(siz_im01);
    image_base01 = image_base01 .* (base1 + 1);
    base2 = (-1:base1/300:-301)';
    base3 = ones(301,1) .* 2;
    base2 = (base3 + base2(1:301,:)) * ones(1,siz_im01(2));
    image_base01(1:301,:) = base2;
else
    image_base01 = ones(siz_im01);
end

% check whether the kymograph has peak or not
mean_dataA3b = mean(dataA3b);
[~,locs3] = findpeaks(dataA3b, 'MinPeakHeight', mean_dataA3b, 'MinPeakProminence', mean_dataA3b*0.5);

if isempty(locs3) % peaks not deteced
    f_y3 = f_y ./ abs(baseFirst1);
    image_ca04 = image_ca03 ./ f_y3;
    peaks1 = 0;
elseif numel(locs3) == 1 % peaks deteced
    image_ca04 = image_ca03 .* image_base01; % backgrond corrected image
    peaks1 = 1;
else
    image_ca04 = image_ca03 .* image_base01; % backgrond corrected image
    peaks1 = 2;
end


%figure, imshow(imrotate(image_ca04,90),[]), title('Background corrected');



%% GCaMP intensity

intensity1 = sum(image_ca04, 2, 'omitnan');

dataA1 = smooth(intensity1, 5, 'sgolay', 1); 
min_dataA1 = min(dataA1);
dataA4 = dataA1 - min_dataA1;
mean_dataA1 = mean(dataA1);

if peaks1 == 2
    threshPeakHeight1 = 0.9 * mean_dataA1;
    threshProminence1 = 0.2 * mean_dataA1;
elseif peaks1 == 1
    threshPeakHeight1 = 1.25 * mean_dataA1;
    threshProminence1 = 0.4 * mean_dataA1;
elseif peaks1 == 0
    threshPeakHeight1 = 1.25 * mean_dataA1;
    threshProminence1 = 0.5 * mean_dataA1;
end
    
    
[pks1,locs1,w1,~] = findpeaks(dataA1, 'WidthReference','halfheight', 'MinPeakHeight', threshPeakHeight1, 'MinPeakProminence', threshProminence1);
PeakData1 = horzcat(locs1, pks1, w1);
PeakArea1 = w1 .* pks1; %estimate area of peak
PeakData1 = horzcat(PeakData1, PeakArea1);

if isempty(PeakData1)
    PeakData1(1,1:6) = NaN;
else
    for n06 = 1:length(locs1)
        if n06 > 1
            PeakData1(n06,5) = locs1(n06) - locs1(n06-1);
            PeakData1(n06,6) = PeakData1(n06,5) / fps;
        else
            PeakData1(n06,5) = NaN;
            PeakData1(n06,6) = NaN;
        end
    end
end

mean_PeakData1 = mean(PeakData1, 1, 'omitnan');
%std_PeakData1 = std(PeakData1, 1, 'omitnan');
Total_calcium = sum(PeakArea1) / (siz_im01(1,1) / fps);

if isnan(PeakData1(1,1))
    numPeaks = 0;
else
    numPeaks = size(PeakData1,1);
end

% figure, findpeaks(dataA4, 'WidthReference','halfheight', 'MinPeakHeight', threshPeakHeight, 'MinPeakProminence', threshProminence);

intensity2 = intensity1 ./ flaRegion2;
unitCal1 = sum(intensity2(:)) / (siz_im01(1,1) / fps);
totalCal1 = sum(dataA1) / (siz_im01(1,1) / fps);
frequency1 = numPeaks / (siz_im01(1,1) / fps);




%% Extract data

Summary = [numPeaks, mean_PeakData1(1,2), mean_PeakData1(1,3), mean_PeakData1(1,4), mean_PeakData1(1,5), mean_PeakData1(1,6), Total_calcium, median_length1, time1/fps, sum(dataA4), totalCal1, unitCal1, frequency1]; 
% number of peaks, mean peak height, mean peak wide, mean peak area, mean
% peak interval (px), mean peak interval (s), total peak calcium intensity,
% calculated length (px), total time (s), total calcium amount, total
% calcium/s,  unit intensity(/s), frequency
%imwrite(image_bw03,'image_bw03.tif');
csvwrite('PeakData1d.csv', PeakData1);
csvwrite('Summary1d.csv', Summary);
