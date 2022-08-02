function [Summary_IFT, Summary_Cal] = fla_calciumD5f(directory, px_size, fps)

% 3/18/2021
% 5/14/2021 Updated
% 5/20/2021 Updated minor bug fix
% 8/24/2021 Updated for new version D2, not calculate calcium concentration
% 9/1/2021 Updated for peak analysis (changed condition for peak
% detection) and background correction
% 11/17/2021 Updated for peak analysis (changed condition for peak
% detection) and background correction
% 12/7/2021 Updated: changed peak detection and calculation method 
% 12/14/2021 Updated: changed calculation of high calcium position
% 2/3/2022 Updated: 
% 3/2/2022 Updated: calculate calcium & IFT
% 5/27/2022 Updated: change calculation of GCaMP intensity and IFT data
% etc.

% Flagellar GCaMP analysis for CalciumIFT.m 
% Flagellum position was manually selected as file "kymograph1_labels.png"

% Need the following functions: formatKymoButlerC2.m, remove_traj_overlap4.m, findPeaks2.m, (VisualizeTrajectories1.m), 

%% Parameters

% fps = 9.66; % frames per second
% px_size = 0.104667; % pixel size in (um/px)
gcampThresh1 = 18; % A threshold for judging have peaks or not
gcampThresh2 = 3; % A threshold for judging have peaks or not
gcampThresh3 = 0.14; % 
gcampThresh4 = 0.22;
duration_speed = 5; % (px)
interval1s = 5; % analysis time (s) after peak
speedTresh1 = 0.5; % (um/s)

%% Read Kymograph images

%directory01 = dir;
directory01 = dir(directory);
val1 = length(directory01);

for m01 = 1:val1
    if regexp(directory01(m01).name, 'kymograph\dI.tif') > 0
        image_km01 = im2double(imread(directory01(m01).name));
    elseif regexp(directory01(m01).name, 'kymograph\dC.tif') > 0
        image_ca01 = im2double(imread(directory01(m01).name));
    elseif regexp(directory01(m01).name, 'kymograph\d_labels.png') > 0
        image_label01 = imread(directory01(m01).name);
    elseif regexp(directory01(m01).name, 'coordinatesA.txt') > 0
        coordinatesA1 = importdata(directory01(m01).name);
    elseif regexp(directory01(m01).name, 'coordinatesR.txt') > 0
        coordinatesR1 = importdata(directory01(m01).name);
    elseif regexp(directory01(m01).name, 'trajectories_Antero.tif') > 0
        image_Antero1 = imread(directory01(m01).name);
    elseif regexp(directory01(m01).name, 'trajectories_Retro.tif') > 0
        image_Retro1 = imread(directory01(m01).name);
    end
end

% figure, imshow(imrotate(image_ca01, 90), []), title('GCaMP');
% figure, imshow(imrotate(image_km01, 90), []), title('IFT');
% figure, imshow(imrotate(image_label01, 90), []), title('label');

clear val1 m01

%% Subtract IFT trajectories from GCaMP

%image_kmA3 = vertcat(image_km01(1,:), image_km01(1:end-1,:));
image_kmA1 = imtranslate(image_km01, [0,0.5]); % slightly shift the IFT image
image_kmA1(1,:) = image_km01(1,:);
%image_kmA3 = imgaussfilt(image_kmA1, 1, 'FilterSize', 3);
image_caA3 = image_ca01 - (image_kmA1 * 0.2); % IFT subtracted GCaMP image
%figure, imshow(imrotate(image_caA3, 90), []), title('GCaMP subtracted');

% figure, imshowpair(imrotate(image_caA3,90), imrotate(image_km01,90));


%% Baseline correction

siz_im01 = size(image_km01);
image_label02 = imbinarize(image_label01);
flaRegion1 = sum(image_label02,2);
flaRegion2 = flaRegion1;
flaRegion2(flaRegion1==0) = NaN;

image_caA4 = imgaussfilt(image_caA3, 1, 'FilterSize', 3);
image_caA5 = image_caA4 .* image_label02;
image_caA5(image_caA5==0) = NaN;
gcampIntensity1 = sum(image_caA5, 2, 'omitnan');

gcampIntensity2 = gcampIntensity1 ./ flaRegion2;
gcampIntensity3a = gcampIntensity2(1:100); % first 100 pixels
gcampIntensity3b = gcampIntensity2(end-99:end); % last 100 pixels

mean_gcampIntensity2 = mean(gcampIntensity2);
mean_gcampIntensity3a = mean(gcampIntensity3a);
std_gcampIntensity3a = std(gcampIntensity3a);
mean_gcampIntensity3b = mean(gcampIntensity3b);
std_gcampIntensity3b = std(gcampIntensity3b);

[pks1a, ~] = findpeaks(-gcampIntensity3a, 'MinPeakHeight', -mean_gcampIntensity3a-std_gcampIntensity3a*0.5); %, 'MinPeakProminence', std_gcampIntensity3a*0.1);
[pks1b, ~] = findpeaks(-gcampIntensity3b, 'MinPeakHeight', -mean_gcampIntensity3b-std_gcampIntensity3b*0.5); %, 'MinPeakProminence', std_gcampIntensity3b*0.1);
[pks1c, ~] = findpeaks(-gcampIntensity2, 'MinPeakHeight', -mean_gcampIntensity2*1.2);


if length(pks1a) <= 2 || length(pks1b) <= 2
    base1 = 0;
else
    baseFirst1 = mean(pks1a);
    baseLast1 = mean(pks1b);
    baseMean1 = mean(pks1c);
    base1 = (baseFirst1 - baseLast1) / baseMean1;
end

if base1 > 1    
    image_base01 = ones(siz_im01);
    image_base01 = image_base01 .* (base1 + 1);
    base2 = (1:base1/300:301)';
    base2 = base2(1:301,:) * ones(1,siz_im01(2));
    image_base01(1:301,:) = base2;
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

image_caA6 = image_caA4 .* image_base01; % backgrond corrected image
%figure, imshow(imrotate(image_caA6,90),[]), title('Background corrected');

%% Background correction

SE_line = strel('line', 9,0);
image_labelB1 = imdilate(image_label02, SE_line);
image_labelB2 = imcomplement(image_labelB1);
image_labelC1 = false(siz_im01);
image_labelC1(:, round(siz_im01(2)/2):siz_im01(2)-6) = 1;
image_labelD1 = image_labelB2 & image_labelC1;
% image_labelC2 = false(siz_im01);
% image_labelC2(:, 7:round(siz_im01(2)/2)) = 1;
% image_labelD2 = image_labelB2 & image_labelC2;

%figure, imshow(imrotate(image_labelC2, 90), [])

image_caB1 = image_caA6 .* image_labelD1;
image_caB1(image_caB1==0) = NaN;
% mean_imcaB1 = mean(image_caB1(:),'omitnan');
% std_imcaB1 = std(image_caB1(:),'omitnan');
% gcampBackground1 = mean_imcaB1 - std_imcaB1;
min_imcaB1 = min(image_caB1(:));
gcampBackground1 = min_imcaB1;

% image_caB2 = image_caA6 .* image_labelD2;
% image_caB2(image_caB2==0) = NaN;
% mean_imcaB2 = mean(image_caB2(:),'omitnan');

% if mean_imcaB1 > mean_imcaB2
%     gcampBackground1 = mean_imcaB2;
% else
%     gcampBackground1 = mean_imcaB1;
% end

image_ca02 = image_caA6 .* image_label02; %exclude non-flagella position
image_ca02(image_ca02==0) = NaN;
%min_imca02 = min(image_ca02(:));
image_ca03 = image_ca02 - gcampBackground1;
% mean_imca03 = mean(image_ca03(:), 'omitnan');
% figure, imshow(imrotate(image_ca03,90),[]), title('GCaMP analysis');

%% Calcium peak intensity

intensity1 = sum(image_ca03, 2, 'omitnan');

dataA1 = smooth(intensity1, 5, 'sgolay', 1); 
minDataA1 = min(dataA1);
maxDataA1 = max(dataA1);
meanDataA1 = mean(dataA1, 'omitnan');
stdDataA1 = std(dataA1, 'omitnan');
intensityRatio1 = (meanDataA1-minDataA1) / (maxDataA1-minDataA1);

intensity2 = intensity1 ./ flaRegion2; % intensity/pixel
dataA2 = smooth(intensity2, 5, 'sgolay', 1); 
maxDataA2 = max(dataA2);
minDataA2 = min(dataA2);
meanDataA2 = mean(dataA2, 'omitnan');
stdDataA2 = std(dataA2, 'omitnan');
% intensityRatio1 = meanDataA2 / maxDataA2;

image_labelA1 = false(siz_im01);
if (maxDataA2-minDataA2) < gcampThresh1 && stdDataA2 < gcampThresh2 % no peaks
    PeakData2 = double.empty; 
elseif intensityRatio1 < gcampThresh3 % small peaks
    [Ypk1,Xpk1,wxPk1,Wpk1,Ppk1] = findPeaks2(dataA1, 'MinPeakHeight', meanDataA1+stdDataA1, 'MinPeakProminence', stdDataA1*0.5);
    PeakData2 = horzcat(Ypk1,Xpk1,wxPk1,Wpk1,Ppk1);
    for n01a = 1:size(PeakData2,1)
        PeakData2(n01a,7) = PeakData2(n01a,1) * PeakData2(n01a,5);
        image_labelA1(round(PeakData2(n01a,3)):round(PeakData2(n01a,4)),:) = 1;
    end
elseif intensityRatio1 < gcampThresh4 % mediam peaks
    [Ypk1,Xpk1,wxPk1,Wpk1,Ppk1] = findPeaks2(dataA1, 'MinPeakHeight', meanDataA1, 'MinPeakProminence', stdDataA1*0.1);
    PeakData2 = horzcat(Ypk1,Xpk1,wxPk1,Wpk1,Ppk1);
    for n01a = 1:size(PeakData2,1)
        PeakData2(n01a,7) = PeakData2(n01a,1) * PeakData2(n01a,5);
        image_labelA1(round(PeakData2(n01a,3)):round(PeakData2(n01a,4)),:) = 1;
    end
else % high peaks
    [Ypk1,Xpk1,wxPk1,Wpk1,Ppk1] = findPeaks2(dataA1, 'MinPeakHeight', meanDataA1-stdDataA1, 'MinPeakProminence', stdDataA1*0.1);
    PeakData2 = horzcat(Ypk1,Xpk1,wxPk1,Wpk1,Ppk1);
    for n01a = 1:size(PeakData2,1)
        PeakData2(n01a,7) = PeakData2(n01a,1) * PeakData2(n01a,5);
        image_labelA1(round(PeakData2(n01a,3)):round(PeakData2(n01a,4)),:) = 1;
    end
end

%% Visualize Calcium peaks on GCaMP image

max_imcaA3 = max(image_caA3(:));
min_imcaA3 = min(image_caA3(:));
image_ca06 = (image_caA3 - min_imcaA3) ./ (max_imcaA3 - min_imcaA3);

% image_B = labeloverlay(image_ca06,image_labelA1, 'Colormap', [1,0,0],'Transparency', 0.8);
% figure, imshow(imrotate(image_B,90))

%% Identify Highly Calcium Influx Time 

image_labelA2 = false(siz_im01); 

if maxDataA2 > gcampThresh1 || stdDataA2 > gcampThresh2 % with peak
    %meanPeakHeight2 = mean(PeakData2(:,1));
    meanCal1 = smooth(dataA2, 9 ,'moving');
    
    if intensityRatio1 < gcampThresh3 % less peaks
        temp1 = meanCal1(:,1) > meanDataA2+stdDataA2*0.5;
        image_labelA2(temp1,:) = 1;
    elseif intensityRatio1 > gcampThresh4 && size(PeakData2,1) > 30 % many peaks
        [~,~,wxPk2,~,~] = findPeaks2(-dataA2, 'MinPeakHeight', -meanDataA2);
        for n02a = 1:size(wxPk2,1)
            image_labelA2(round(wxPk2(n02a,1)):round(wxPk2(n02a,2)),:) = 1;
        end
        image_labelA2 = not(image_labelA2);
        [~,~,wxPk3,~,~] = findPeaks2(dataA2, 'MinPeakHeight', meanDataA2-stdDataA2);
        for n02b = 1:size(wxPk3,1)
            image_labelA2(round(wxPk3(n02b,1)):round(wxPk3(n02b,2)),:) = 1;
        end
    else % 
        temp1 = find(meanCal1(:,1)>meanDataA2-stdDataA2*0.2);
        image_labelA2(temp1,:) = 1;         
    end
end

%image_C = labeloverlay(image_ca06,image_labelA2, 'Colormap', [1,0,0],'Transparency', 0.8);
%figure, imshow(imrotate(image_C,90))

%% Format KymoButler Coordinates

if exist('image_Antero1', 'var')
    max_Antero1 = max(image_Antero1(:));
    trajectory_listA2 = cell(max_Antero1,2);
    for n02a = 1:max_Antero1
        trajectory_listA2{n02a,2} = find(image_Antero1==n02a);
        [x_val1,y_val1] = ind2sub(siz_im01, trajectory_listA2{n02a,2});
        trajectory_listA2{n02a,1} = [x_val1,y_val1];
        trajectory_listA2{n02a,1} = sortrows(trajectory_listA2{n02a,1});
    end

else
    trajectory_listA1 = formatKymoButlerC2(siz_im01, coordinatesA1);

    % Exclude trajectories outside of flagella region
    trajectory_listA2 = cell(size(trajectory_listA1,1),2);
    for n03a = 1:size(trajectory_listA1,1)
        positive1 = image_label02(trajectory_listA1{n03a,2})>0;
        trajectory_listA2{n03a,2} = trajectory_listA1{n03a,2}(positive1);
        [x_val1,y_val1] = ind2sub(siz_im01, trajectory_listA2{n03a,2});
        trajectory_listA2{n03a,1} = [x_val1,y_val1];
        trajectory_listA2{n03a,1} = sortrows(trajectory_listA2{n03a,1});
    end
end

temp_Idx2a = cellfun('length',trajectory_listA2(:,1));
temp_Idx2b = find(temp_Idx2a < 10);
trajectory_listA2(temp_Idx2b,:) = [];

% fix flipped trajectories of retrograde IFT

if exist('image_Antero1', 'var')
    max_Retro1 = max(image_Retro1(:));
    trajectory_listR2 = cell(max_Retro1,2);
    for n02b = 1:max_Retro1
        trajectory_listR2{n02b,2} = find(image_Retro1==n02b);
        [x_val2,y_val2] = ind2sub(siz_im01, trajectory_listR2{n02b,2});
        trajectory_listR2{n02b,1} = [x_val2,y_val2];
        trajectory_listR2{n02b,1} = sortrows(trajectory_listR2{n02b,1},[1,2],{'ascend' 'descend'});
    end

else
    trajectory_listR1 = formatKymoButlerC2(siz_im01, coordinatesR1);

    trajectory_listR2 = cell(size(trajectory_listR1,1),1);
    image_label02b = fliplr(image_label02);

    for n03b = 1:size(trajectory_listR1,1)
        positive2 = image_label02b(trajectory_listR1{n03b,2})>0;
        temp_Idx1 = trajectory_listR1{n03b,2}(positive2);
        [x_val2,y_val2] = ind2sub(siz_im01, temp_Idx1);
        temp_Sub1 = [x_val2,y_val2];
        y_val3 = siz_im01(1,2) - temp_Sub1(:,2) +1;
        x_val3 = temp_Sub1(:,1);
        trajectory_listR2{n03b,1} = [x_val3, y_val3];
        trajectory_listR2{n03b,1} = sortrows(trajectory_listR2{n03b,1},[1,2],{'ascend' 'descend'});
        trajectory_listR2{n03b,2} = sub2ind(siz_im01, x_val3, y_val3);
    end
end

temp_Idx3a = cellfun('length',trajectory_listR2(:,1));
temp_Idx3b = find(temp_Idx3a < 10);
trajectory_listR2(temp_Idx3b,:) = [];

%VisualizeTrajectories1(trajectory_listA2,siz_im01);
%VisualizeTrajectories1(trajectory_listR2,siz_im01);

%% IFT velocity


% Anterograde
sizeTrajectoryA1 = size(trajectory_listA2);
for n05a1 = 1:sizeTrajectoryA1(1,1)
    duration1 = size(trajectory_listA2{n05a1,1},1);
    speedA1 = zeros(duration1, 1);
    for n05a2 = 1:duration1-duration_speed
            startPosition1 = trajectory_listA2{n05a1,1}(n05a2,2);
            endPosition1 = trajectory_listA2{n05a1,1}(find(trajectory_listA2{n05a1,1}(:,1) == trajectory_listA2{n05a1,1}(n05a2,1)+duration_speed,1,'last'),2);
            if isempty(endPosition1)
                approxIdx1 = knnsearch(trajectory_listA2{n05a1,1}(:,1),trajectory_listA2{n05a1,1}(n05a2,1)+duration_speed);
                endPosition1 = trajectory_listA2{n05a1,1}(approxIdx1,2);
                speedA1(n05a2,1) = (endPosition1-startPosition1) * px_size / ((trajectory_listA2{n05a1,1}(approxIdx1,1) - trajectory_listA2{n05a1,1}(n05a2,1)) / fps);
            else
                speedA1(n05a2,1) = (endPosition1-startPosition1) * px_size / (duration_speed / fps);
            end
    end
    speedA1(duration1-duration_speed+1:end,1) = (trajectory_listA2{n05a1,1}(duration1,2) - trajectory_listA2{n05a1,1}(duration1-duration_speed,2)) * px_size * fps / (trajectory_listA2{n05a1,1}(duration1,1) - trajectory_listA2{n05a1,1}(duration1-duration_speed,1));
    speedA1(isinf(speedA1)) = NaN;
    trajectory_listA2{n05a1,3} = speedA1;
end
    
% Retrograde
sizeTrajectoryR1 = size(trajectory_listR2);
for n05b1 = 1:sizeTrajectoryR1(1,1)
    duration1 = size(trajectory_listR2{n05b1,1},1);
    speedR1 = zeros(duration1, 1);
    for n05b2 = 1:duration1-duration_speed
            startPosition1 = trajectory_listR2{n05b1,1}(n05b2,2);
            endPosition1 = trajectory_listR2{n05b1,1}(find(trajectory_listR2{n05b1,1}(:,1) == trajectory_listR2{n05b1,1}(n05b2,1)+duration_speed,1,'last'),2);
            if isempty(endPosition1)
                approxIdx1 = knnsearch(trajectory_listR2{n05b1,1}(:,1),trajectory_listR2{n05b1,1}(n05b2,1)+duration_speed);
                endPosition1 = trajectory_listR2{n05b1,1}(approxIdx1,2);
                speedR1(n05b2,1) = (endPosition1-startPosition1) * px_size / ((trajectory_listR2{n05b1,1}(approxIdx1,1) - trajectory_listR2{n05b1,1}(n05b2,1)) / fps);
            else
                speedR1(n05b2,1) = (endPosition1-startPosition1) * px_size / (duration_speed / fps);
            end
    end
    speedR1(duration1-duration_speed+1:end,1) = (trajectory_listR2{n05b1,1}(duration1,2) - trajectory_listR2{n05b1,1}(duration1-duration_speed,2)) * px_size * fps / (trajectory_listR2{n05b1,1}(duration1,1) - trajectory_listR2{n05b1,1}(duration1-duration_speed,1));
    speedR1(isinf(speedR1)) = NaN;
    trajectory_listR2{n05b1,3} = speedR1;
end


%% IFT analysis 


image_km02 = image_km01 .* image_base01; %baseline correction

% Anterograde
trajectoryDataA1 = double.empty(sizeTrajectoryA1(1,1),0);
for n06a = 1:sizeTrajectoryA1(1,1)
    trajectoryDataA1(n06a,1) = n06a; % ID
    trajectoryDataA1(n06a,2) = sum(image_km02(trajectory_listA2{n06a,2}),'omitnan'); % total intensity
    trajectoryDataA1(n06a,3) = mean(image_km02(trajectory_listA2{n06a,2}),'omitnan'); % mean intensity
    trajectoryDataA1(n06a,4) = mean(trajectory_listA2{n06a,3},'omitnan'); % mean velocity (um/s)
    trajectoryDataA1(n06a,5) = median(trajectory_listA2{n06a,3},'omitnan'); % median velocity (um/s)

    if trajectory_listA2{n06a,3}(1,1) < speedTresh1
        startA1 = find(trajectory_listA2{n06a,3}(:,1)>speedTresh1,1);
        if ~isempty(startA1)
            trajectoryDataA1(n06a,6) = trajectory_listA2{n06a,1}(startA1+duration_speed-1,2); % start position (px)
            trajectoryDataA1(n06a,7) = trajectory_listA2{n06a,1}(startA1+duration_speed-1,1); % start time (px)
        else
            trajectoryDataA1(n06a,6) = NaN;
            trajectoryDataA1(n06a,7) = NaN;
        end
    else
        trajectoryDataA1(n06a,6) = trajectory_listA2{n06a,1}(1,2); % start position (px)
        trajectoryDataA1(n06a,7) = trajectory_listA2{n06a,1}(1,1); % start time (px)
        startA1 = 1;
    end
    
    if ~isnan(trajectoryDataA1(n06a,6))
        stopA1 = find(trajectory_listA2{n06a,3}(startA1:end,1)<speedTresh1/2,1);
        if ~isempty(stopA1)
            restartA1 = find(trajectory_listA2{n06a,3}(stopA1+startA1:end,1)>speedTresh1,1);
            if ~isempty(restartA1)
                trajectoryDataA1(n06a,8) = trajectory_listA2{n06a,1}(restartA1+stopA1+startA1+duration_speed-3,2); % restart position
                trajectoryDataA1(n06a,9) = trajectory_listA2{n06a,1}(restartA1+stopA1+startA1+duration_speed-3,1); % restart time (px)
            else
                trajectoryDataA1(n06a,8) = NaN; % restart position (px)
                trajectoryDataA1(n06a,9) = NaN; % restart time (px)
            end
        else
            trajectoryDataA1(n06a,8) = NaN; % restart position (px)
            trajectoryDataA1(n06a,9) = NaN; % restart time (px)
        end
    else
        trajectoryDataA1(n06a,8) = NaN; % restart position (px)
        trajectoryDataA1(n06a,9) = NaN; % restart time (px)
    end
    trajectoryDataA1(n06a,10) = find(image_label02(trajectory_listA2{n06a,1}(1,1),:),1); % flagella position
end

% Retrograde
trajectoryDataR1 = double.empty(sizeTrajectoryR1(1,1),0);
for n06b = 1:sizeTrajectoryR1(1,1)
    trajectoryDataR1(n06b,1) = n06b; % ID
    trajectoryDataR1(n06b,2) = sum(image_km02(trajectory_listR2{n06b,2}),'omitnan'); % total intensity
    trajectoryDataR1(n06b,3) = mean(image_km02(trajectory_listR2{n06b,2}),'omitnan'); % mean intensity
    trajectoryDataR1(n06b,4) = mean(trajectory_listR2{n06b,3},'omitnan'); % mean velocity (um/s)
    trajectoryDataR1(n06b,5) = median(trajectory_listR2{n06b,3},'omitnan'); % median velocity (um/s)
    if trajectory_listR2{n06b,3}(1,1) > -speedTresh1
        startR1 = find(trajectory_listR2{n06b,3}(:,1)<-speedTresh1,1);
        if ~isempty(startR1)
            trajectoryDataR1(n06b,6) = trajectory_listR2{n06b,1}(startR1+duration_speed-1,2); % start position (px)
            trajectoryDataR1(n06b,7) = trajectory_listR2{n06b,1}(startR1+duration_speed-1,1); % start time (px)
        else
            trajectoryDataR1(n06b,6) = NaN; % start position (px)
            trajectoryDataR1(n06b,7) = NaN; % start time (px)
        end
    else
        trajectoryDataR1(n06b,6) = trajectory_listR2{n06b,1}(1,2); % start position (px)
        trajectoryDataR1(n06b,7) = trajectory_listR2{n06b,1}(1,1); % start time (px)
        startR1 = 1;
    end
    if ~isnan(trajectoryDataR1(n06b,6))
        stopR1 = find(trajectory_listR2{n06b,3}(startR1:end,1)>-speedTresh1/2,1);
        if ~isempty(stopR1)
            restartR1 = find(trajectory_listR2{n06b,3}(stopR1+startR1:end,1)<-speedTresh1,1);
            if ~isempty(restartR1)
                trajectoryDataR1(n06b,8) = trajectory_listR2{n06b,1}(restartR1+stopR1+startR1+duration_speed-3,2); % start position (px)
                trajectoryDataR1(n06b,9) = trajectory_listR2{n06b,1}(restartR1+stopR1+startR1+duration_speed-3,1); % strat time (px)
            else
                trajectoryDataR1(n06b,8) = NaN; % restart position (px)
                trajectoryDataR1(n06b,9) = NaN; % restart time (px)
            end
        else
            trajectoryDataR1(n06b,8) = NaN; % restart position (px)
            trajectoryDataR1(n06b,9) = NaN; % restart time (px)
        end
    else
        trajectoryDataR1(n06b,8) = NaN; % restart position (px)
        trajectoryDataR1(n06b,9) = NaN; % restart time (px)
    end
    
    trajectoryDataR1(n06b,10) = find(image_label02(trajectory_listR2{n06b,1}(1,1),:),1,'last'); % flagella position
    
    trajectoryDataR1(n06b,12) = trajectory_listR2{n06b,1}(end,2); % end position (px)
    trajectoryDataR1(n06b,13) = trajectory_listR2{n06b,1}(end,1); % end time (px)
    trajectoryDataR1(n06b,14) = find(image_label02(trajectoryDataR1(n06b,13),:),1);
end


trajectorySummaryA1 = zeros(1,6);
trajectorySummaryA1(1,1) = size(trajectoryDataA1,1); % numbers of IFT trajectories
mean_trajDataA1 = mean(trajectoryDataA1,'omitnan'); 
trajectorySummaryA1(1,2:5) = mean_trajDataA1(1,2:5); % mean IFT intensity, average IFT intensity, speed, median speed

trajectorySummaryR1 = zeros(1,6);
trajectorySummaryR1(1,1) = size(trajectoryDataR1,1);
mean_trajDataR1 = mean(trajectoryDataR1,'omitnan');
trajectorySummaryR1(1,2:5) = mean_trajDataR1(1,2:5);


%% Image for IFT intensity

image_IFTA1 = false(siz_im01);

for m02a = 1:length(trajectory_listA2)
    image_IFTA1(trajectory_listA2{m02a,2}) = 1;
end

% SE_line3 = strel('line',3,90);
SE_disk1 = strel('disk',1);
image_IFTA2 = imdilate(image_IFTA1, SE_disk1);
image_IFTA3 = image_km02 .* image_IFTA2 .* image_label02;


% figure, imshow(imrotate(image_IFTA1,90))
unitIntA1 = sum(image_IFTA3(:)) / (siz_im01(1,1) / fps) / (mean(sum(image_label02,2)) * px_size);

image_IFTR1 = false(siz_im01);

for m02b = 1:length(trajectory_listR2)
    image_IFTR1(trajectory_listR2{m02b,2}) = 1;
end
image_IFTR2 = imdilate(image_IFTR1, SE_disk1);
image_IFTR3 = image_km02 .* image_IFTR2 .* image_label02;
unitIntR1 = sum(image_IFTR3(:)) / (siz_im01(1,1) / fps) / (mean(sum(image_label02,2)) * px_size);

% background
image_IFTinv1 = image_km02 .* imcomplement(image_IFTA2) .* imcomplement(image_IFTR2) .* image_label02; % background (exclude IFT trajectories)
image_IFTinv2 = imcomplement(image_IFTA2) .* imcomplement(image_IFTR2) .* image_label02; % background location

%% Detect IFT trajectories which start from the base and back to the base

% select anterograde trajectories which start from the base
for n07a1 = 1:sizeTrajectoryA1(1,1)
    if trajectoryDataA1(n07a1,6) <= trajectoryDataA1(n07a1,10)+5
        trajectoryDataA1(n07a1,11) = 1;
    else
        trajectoryDataA1(n07a1,11) = 0;
    end
end

trajectoryDataA2 = sortrows(trajectoryDataA1, [7,6]);
temp_Idx4a = find(~trajectoryDataA2(:,11));
trajectoryDataA2(temp_Idx4a,:) = [];

for n07a2 = 2:size(trajectoryDataA2,1)
    trajectoryDataA2(n07a2,12) = sum(dataA1(trajectoryDataA2(n07a2-1,7):trajectoryDataA2(n07a2,7))); % amount of calcium injection from previous IFT injection
    trajectoryDataA2(n07a2,13) = trajectoryDataA2(n07a2,7) - trajectoryDataA2(n07a2-1,7); % time from previous IFT injection (px)
end

trajectorySummaryA1(1,6) = size(trajectoryDataA2,1) / (siz_im01(1,1) / fps); % frequency

% select retrograde trajectories which start from the tip
for n07b1 = 1:sizeTrajectoryR1(1,1)
    if trajectoryDataR1(n07b1,6) >= trajectoryDataR1(n07b1,10)-15
        trajectoryDataR1(n07b1,11) = 1;
    else
        trajectoryDataR1(n07b1,11) = 0;
    end
end

trajectoryDataR2 = sortrows(trajectoryDataR1, [7,6]);
temp_Idx4b = find(~trajectoryDataR2(:,11));
trajectoryDataR2(temp_Idx4b,:) = [];

for n07b2 = 2:size(trajectoryDataR2,1)
    trajectoryDataR2(n07b2,12) = sum(dataA1(trajectoryDataR2(n07b2-1,7):trajectoryDataR2(n07b2,7))); % amount of calcium injection from previous IFT injection
    trajectoryDataR2(n07b2,13) = trajectoryDataR2(n07b2,7) - trajectoryDataR2(n07b2-1,7); % time from previous IFT injection (px)
end

trajectorySummaryR1(1,6) = size(trajectoryDataR2,1) / (siz_im01(1,1) / fps);

%figure, histogram(trajectoryDataA2(:,11))
%figure, scatter(trajectoryDataA2(:,11),trajectoryDataA2(:,2)),
%xlabel('calcium'), ylabel('IFT intensity')


% select retrograde trajectories which are back to the base

for n07c1 = 1:sizeTrajectoryR1(1,1)
    if trajectoryDataR1(n07c1,12) <= trajectoryDataR1(n07c1,14)+6
        trajectoryDataR1(n07c1,15) = 1;
    else
        trajectoryDataR1(n07c1,15) = 0;
    end
end

trajectoryDataR3 = sortrows(trajectoryDataR1, [7,6]); % retrograde IFT which returns to the base
temp_Idx4c = find(~trajectoryDataR3(:,15));
trajectoryDataR3(temp_Idx4c,:) = [];



%% Calculate background intensity

backgroundIntensity1 = sum(image_IFTinv1,2);
backgroundRegion1 = sum(image_IFTinv2,2);
backgroundIntensity2 = backgroundIntensity1 ./ backgroundRegion1;
backgroundIntensity3 = smooth(backgroundIntensity2, 5, 'sgolay', 1); 

% figure, plot(backgroundIntensity3);

% difference1 = diff(backgroundIntensity3);
% 
% mean_back3 = mean(backgroundIntensity3);
% figure, plot(difference1)


%% IFT injection graph

iftA1 = zeros(siz_im01(1),1);

for n11a = 1:size(trajectoryDataA2,1)
    iftA1(trajectoryDataA2(n11a,7),1) = trajectoryDataA2(n11a,3);
end

iftA2 = smoothdata(iftA1,'gaussian',10);
%figure, plot(iftA2);

iftR1 = zeros(siz_im01(1),1);
for n11b = 1:size(trajectoryDataR3,1)
    iftR1(trajectoryDataR3(n11b,13),1) = trajectoryDataR3(n11b,3);
end
iftR2 = smoothdata(iftR1,'gaussian',10);
%figure, plot(iftR2);

%% Calcium & IFT analysis 1 (Calcium vs. non-calcium time)

% IFT speed
calPositive1 = find(image_labelA2(:,1));
calNegative1 = find(~image_labelA2(:,1));

%Anterograde
for n09a = 1:size(trajectoryDataA1,1)
    Lia1a = ismember(trajectory_listA2{n09a,1}(:,1),calPositive1);
    Lia2a = ismember(trajectory_listA2{n09a,1}(:,1),calNegative1);
    trajectoryDataA1(n09a,12) = mean(trajectory_listA2{n09a,3}(Lia1a,1));
    trajectoryDataA1(n09a,13) = mean(trajectory_listA2{n09a,3}(Lia2a,1));
end

%Retrograde
for n09b = 1:size(trajectoryDataR1,1)
    Lia1b = ismember(trajectory_listR2{n09b,1}(:,1),calPositive1);
    Lia2b = ismember(trajectory_listR2{n09b,1}(:,1),calNegative1);
    trajectoryDataR1(n09b,16) = mean(trajectory_listR2{n09b,3}(Lia1b,1));
    trajectoryDataR1(n09b,17) = mean(trajectory_listR2{n09b,3}(Lia2b,1));
end

% Check trajectories start with calcium influx or not
for n08a = 1:size(trajectoryDataA2,1)
    if ~isnan(trajectoryDataA2(n08a,7))
        if find(calPositive1==trajectoryDataA2(n08a,7))
            trajectoryDataA2(n08a,14) = 1; % IFT trajectories start with calcium influx
        else
            trajectoryDataA2(n08a,14) = 0;
        end
    else
        trajectoryDataA2(n08a,14) = 0;
    end
end


for n08b = 1:size(trajectoryDataR2,1)
    if ~isnan(trajectoryDataR2(n08b,7))
        if find(calPositive1==trajectoryDataR2(n08b,7))
            trajectoryDataR2(n08b,15) = 1; % IFT trajectories start with calcium influx
        else
            trajectoryDataR2(n08b,15) = 0;
        end
    else
        trajectoryDataR2(n08b,15) = 0;
    end
end

fla_length1 = mean(sum(image_label02(image_labelA2(:,1),:),2)); % mean flagella length with calcium influx region
fla_length2 = mean(sum(image_label02(~image_labelA2(:,1),:),2)); % mean flagella length with non-calcium influx region
CalIFT1 = double.empty(0,13);

%Anterograde
CalIFT1(1,1) = length(find(trajectoryDataA2(:,14))); % numbers of IFT trajectories start with calcium influx
CalIFT1(2,1) = length(find(~trajectoryDataA2(:,14))); % numbers of IFT trajectories start with calcium influx
CalIFT1(1,2) = length(find(image_labelA2(:,1))) / fps; % total time with calcium influx (s)
CalIFT1(2,2) = length(find(~image_labelA2(:,1))) / fps; % total time without calcium influx (s)
CalIFT1(1,3) = CalIFT1(1,1) / CalIFT1(1,2); % IFT frequency (mean number of IFT trajectories) with calcium influx
CalIFT1(2,3) = CalIFT1(2,1) / CalIFT1(2,2); % IFT frequency (mean number of IFT trajectories) witouth calcium influx
CalIFT1(1,4) = sum(trajectoryDataA2(find(trajectoryDataA2(:,14)),3)); % total intensity of IFT trajectories with calcium influx
CalIFT1(2,4) = sum(trajectoryDataA2(find(~trajectoryDataA2(:,14)),3)); % total intensity of IFT trajectories witouth calcium influx
CalIFT1(1,5) = CalIFT1(1,4) / CalIFT1(1,2); % mean IFT intensity with calcium influx
CalIFT1(2,5) = CalIFT1(2,4) / CalIFT1(2,2); % mean IFT intensity without calcium influx

CalIFT1(1,6) = mean(trajectoryDataA1(:,12),'omitnan'); % mean IFT velosity with calcium influx
CalIFT1(2,6) = mean(trajectoryDataA1(:,13),'omitnan'); % mean IFT velosity without calcium influx
image_temp1a = image_IFTA3 .* image_labelA2;
image_temp1b = image_IFTA3 .* ~image_labelA2;
CalIFT1(1,7) = sum(image_temp1a(:)) / CalIFT1(1,2) / fla_length1; % unit intensity of IFT with calcium influx
CalIFT1(2,7) = sum(image_temp1b(:)) / CalIFT1(2,2) / fla_length2; % unit intensity of IFT with calcium influx

%Retrograde
CalIFT1(1,8) = length(find(trajectoryDataR2(:,15))); % numbers of IFT trajectories start with calcium influx
CalIFT1(2,8) = length(find(~trajectoryDataR2(:,15))); % numbers of IFT trajectories start with calcium influx
CalIFT1(1,9) = CalIFT1(1,8) / CalIFT1(1,2); % IFT frequency (mean number of IFT trajectories) with calcium influx
CalIFT1(2,9) = CalIFT1(2,8) / CalIFT1(2,2); % IFT frequency (mean number of IFT trajectories) witouth calcium influx
CalIFT1(1,10) = sum(trajectoryDataR2(find(trajectoryDataR2(:,15)),3)); % total intensity of IFT trajectories with calcium influx
CalIFT1(2,10) = sum(trajectoryDataR2(find(~trajectoryDataR2(:,15)),3)); % total intensity of IFT trajectories witouth calcium influx
CalIFT1(1,11) = CalIFT1(1,10) / CalIFT1(1,2); % mean IFT intensity with calcium influx
CalIFT1(2,11) = CalIFT1(2,10) / CalIFT1(2,2); % mean IFT intensity without calcium influx
CalIFT1(1,12) = mean(trajectoryDataR1(:,16),'omitnan'); % mean IFT velosity with calcium influx
CalIFT1(2,12) = mean(trajectoryDataR1(:,17),'omitnan'); % mean IFT velosity without calcium influx

image_temp1c = image_IFTR3 .* image_labelA2;
image_temp1d = image_IFTR3 .* ~image_labelA2;
CalIFT1(1,13) = sum(image_temp1c(:)) / CalIFT1(1,2) / fla_length1; % unit intensity of IFT with calcium influx
CalIFT1(2,13) = sum(image_temp1d(:)) / CalIFT1(2,2) / fla_length2; % unit intensity of IFT without calcium influx

image_temp1e = image_IFTinv1 .* image_labelA2;
image_temp1f = image_IFTinv2 .* image_labelA2;
image_temp1g = image_IFTinv1 .* ~image_labelA2;
image_temp1h = image_IFTinv2 .* ~image_labelA2;
CalIFT1(1,14) = sum(image_temp1e(:)) / sum(image_temp1f(:)); % background intensity per pixel with calcium influx
CalIFT1(2,14) = sum(image_temp1g(:)) / sum(image_temp1h(:)); % background intensity per pixel without calcium influx

%% Calcium & IFT analysis 2 

% calculates IFT data by calcium influx group

if ~isempty(PeakData2)
    interval1p = round(interval1s * fps);
    interval2s = interval1p / fps;
    sizePeakData2 = size(PeakData2);
    
    group1 = bwlabel(image_labelA2(:,1),4);
    
    for n10a = 1:sizePeakData2(1)
        PeakData2(n10a,8) = group1(PeakData2(n10a,2),1); % group ID
        if n10a == 1
            PeakData2(n10a,9) = NaN;
        else
            PeakData2(n10a,9) = PeakData2(n10a,2) - PeakData2(n10a-1,2); % interval (px) from the previous peak
            PeakData2(n10a,10) = PeakData2(n10a,9) / fps; % interval (s) from the previous peak
        end
    end
    
    CalPosition1 = zeros(PeakData2(end,8),6);
    group3a = unique(PeakData2(:,8));
    group3b = group3a(find(group3a));
    for n10b = 1:length(group3b)
        group2 = find(PeakData2(:,8)==group3b(n10b),1);
        if isempty(group2)
            CalPosition1(group3b(n10b),1:6) = NaN;
        else
            CalPosition1(group3b(n10b),1) = PeakData2(group2,2) - interval1p *2;
            CalPosition1(group3b(n10b),2) = PeakData2(group2,2) - interval1p;
            CalPosition1(group3b(n10b),3) = PeakData2(group2,2);
            CalPosition1(group3b(n10b),4) = PeakData2(group2,2) + interval1p;
            CalPosition1(group3b(n10b),5) = PeakData2(group2,2) + interval1p *2;
            CalPosition1(group3b(n10b),6) = PeakData2(group2,2) + interval1p *3; 
        end
    end
%     outValue1 = CalPosition1<1;
%     outValue2 = find(CalPosition1>=siz_im01(1));
    CalPosition1(CalPosition1<1) = NaN;
    CalPosition1(CalPosition1>=siz_im01(1)) = NaN;
    maxNum1 = max(PeakData2(:,8));
    
    DataCalInt1 = zeros(maxNum1,6);
    DataIftUnitIntA1 = zeros(maxNum1,6);
    DataIftUnitIntR1 = zeros(maxNum1,6);
    DataIftFrequencyA1 = zeros(maxNum1,6);
    DataIftFrequencyR1 = zeros(maxNum1,6);
    DataIftTraIntA1 = zeros(maxNum1,6);
    DataIftTraIntR1 = zeros(maxNum1,6);
    DataIftBackground1 = zeros(maxNum1,6);
    DataIftReturnInt1 = zeros(maxNum1,6);
    DataIftVelosA1 = zeros(maxNum1,6);
    DataIftVelosR1 = zeros(maxNum1,6);
    for n10c = 1:maxNum1
        DataCalInt1(n10c,1) = n10c;
        DataIftUnitIntA1(n10c,1) = n10c;
        DataIftUnitIntR1(n10c,1) = n10c;
        DataIftFrequencyA1(n10c,1) = n10c;
        DataIftFrequencyR1(n10c,1) = n10c;
        DataIftTraIntA1(n10c,1) = n10c;
        DataIftTraIntR1(n10c,1) = n10c;
        DataIftBackground1(n10c,1) = n10c;
        DataIftReturnInt1(n10c,1) = n10c;
        DataIftVelosA1(n10c,1) = n10c;
        DataIftVelosR1(n10c,1) = n10c;
        for n10d = 2:6
            if isnan(CalPosition1(n10c,n10d-1)) || isnan(CalPosition1(n10c,n10d))
                DataCalInt1(n10c,n10d) = NaN;
                DataIftUnitIntA1(n10c,n10d) = NaN;
                DataIftUnitIntR1(n10c,n10d) = NaN;
                DataIftFrequencyA1(n10c,n10d) = NaN;
                DataIftFrequencyR1(n10c,n10d) = NaN;
                DataIftTraIntA1(n10c,n10d) = NaN;
                DataIftTraIntR1(n10c,n10d) = NaN;
                DataIftBackground1(n10c,n10d) = NaN;
                DataIftReturnInt1(n10c,n10d) = NaN;
                DataIftVelosA1(n10c,n10d) = NaN;
                DataIftVelosR1(n10c,n10d) = NaN;
            else 
                DataCalInt1(n10c,n10d) = sum(dataA1(CalPosition1(n10c,n10d-1):CalPosition1(n10c,n10d)-1));
                flaLengthTemp1 = mean(sum(image_label02(CalPosition1(n10c,n10d-1):CalPosition1(n10c,n10d)-1,:),2));
                imageCalTemp1 = false(siz_im01);
                imageCalTemp1(CalPosition1(n10c,n10d-1):CalPosition1(n10c,n10d)-1,:) = 1;
                imageTemp1a = image_IFTA3 .* imageCalTemp1;
                DataIftUnitIntA1(n10c,n10d) = sum(imageTemp1a(:)) / interval2s / flaLengthTemp1; 
                imageTemp1b = image_IFTR3 .* imageCalTemp1;
                DataIftUnitIntR1(n10c,n10d) = sum(imageTemp1b(:)) / interval2s / flaLengthTemp1; 
                trajTempA1 = trajectoryDataA2(trajectoryDataA2(:,7)>CalPosition1(n10c,n10d-1) & trajectoryDataA2(:,7)<CalPosition1(n10c,n10d)-1,:);
                DataIftFrequencyA1(n10c,n10d) = size(trajTempA1,1) / interval2s;
                DataIftTraIntA1(n10c,n10d) = sum(trajTempA1(:,3)) / interval2s;
                trajTempR1 = trajectoryDataR2(trajectoryDataR2(:,7)>CalPosition1(n10c,n10d-1) & trajectoryDataR2(:,7)<CalPosition1(n10c,n10d)-1,:);
                DataIftFrequencyR1(n10c,n10d) = size(trajTempR1,1) / interval2s;
                DataIftTraIntR1(n10c,n10d) = sum(trajTempR1(:,3)) / interval2s;
                trajTempR2 = trajectoryDataR3(trajectoryDataR3(:,13)>CalPosition1(n10c,n10d-1) & trajectoryDataR3(:,13)<CalPosition1(n10c,n10d)-1,:);
                DataIftReturnInt1(n10c,n10d) = sum(trajTempR2(:,3)) / interval2s;
                imageTemp2a = image_IFTinv1 .* imageCalTemp1;
                imageTemp2b = image_IFTinv2 .* imageCalTemp1;
                DataIftBackground1(n10c,n10d) = sum(imageTemp2a(:)) / sum(imageTemp2b(:));
                trajTempA3 = NaN(size(trajectory_listA2,1),1);
                for n10e = 1:size(trajectory_listA2,1)
                    Lia1c = ismember(trajectory_listA2{n10e,1}(:,1),CalPosition1(n10c,n10d-1):CalPosition1(n10c,n10d)-1);
                    trajTempA3(n10e,1) = mean(trajectory_listA2{n10e,3}(Lia1c,1), 'omitnan');
                end
                DataIftVelosA1(n10c,n10d) = mean(trajTempA3(:,1), 'omitnan');
                trajTempR3 = NaN(size(trajectory_listR2,1),1);
                for n10f = 1:size(trajectory_listR2,1)
                    Lia2c = ismember(trajectory_listR2{n10f,1}(:,1),CalPosition1(n10c,n10d-1):CalPosition1(n10c,n10d)-1);
                    trajTempR3(n10f,1) = mean(trajectory_listR2{n10f,3}(Lia2c,1), 'omitnan');
                end
                DataIftVelosR1(n10c,n10d) = mean(trajTempR3(:,1), 'omitnan');
                
            end
        end
    end
end



      

%% Summarize data

if ~isempty(PeakData2)
    TotalPeakCal = sum(PeakData2(:,7)) / (siz_im01(1,1) / fps);
    MeanPeakCal = mean(PeakData2(:,7));
else
    TotalPeakCal = 0;
    MeanPeakCal = 0;
end

Summary_IFT = horzcat(trajectorySummaryA1, unitIntA1, trajectorySummaryR1, unitIntR1);
UnitCal1 = sum(dataA2) / (siz_im01(1,1) / fps); % unit intensity of calcium
TotalCal1 = sum(dataA1) / (siz_im01(1,1) / fps); % total calcium/sec
Summary_Cal = horzcat(UnitCal1, TotalCal1, TotalPeakCal, MeanPeakCal, CalIFT1(1,1:14),CalIFT1(2,1:14)); % summary of calcium & IFT analysis 1




%%

% figure, subplot(3,1,1)
% plot(dataA1)
% title('GCaMP')
% 
% subplot(3,1,2)
% plot(iftA2)
% title('IFT injection')
% 
% subplot(3,1,3)
% plot(iftR2)
% title('IFT return')

%%

% savefig('graphs1.fig')

writematrix(trajectoryDataA1, 'Trajectory_Antero1.xls');
writematrix(trajectoryDataR1, 'Trajectory_Retro1.xls');
writematrix(trajectoryDataA2, 'Trajectory_Antero2.xls'); % anterograde trajectories which start from the base
writematrix(trajectoryDataR2, 'Trajectory_Retro2.xls'); % retrograde trajectories which start from the tip
writematrix(trajectoryDataR3, 'Trajectory_Retro3.xls'); % retrograde trajectories which reach to the base
writematrix(CalIFT1, 'Analysis_Calcium-IFT.xls');
writematrix(dataA1, 'calcium_influx.txt');
writematrix(iftA1, 'IFT_antero.txt');
writematrix(iftR1, 'IFT_retro.txt');




if ~isempty(PeakData2)
    writematrix(PeakData2, 'calcium_peaks.xls');
    writematrix(DataCalInt1, 'DataCalInt1.xls');
    writematrix(DataIftUnitIntA1, 'DataIftUnitIntA1.xls');
    writematrix(DataIftUnitIntR1, 'DataIftUnitIntR1.xls');
    writematrix(DataIftFrequencyA1, 'DataIftFrequencyA1.xls');
    writematrix(DataIftFrequencyR1, 'DataIftFrequencyR1.xls');
    writematrix(DataIftTraIntA1, 'DataIftTraIntA1.xls');
    writematrix(DataIftTraIntR1, 'DataIftTraIntR1.xls');
    writematrix(DataIftBackground1, 'DataIftBackground1.xls');
    writematrix(DataIftReturnInt1, 'DataIftReturnInt1.xls');
    writematrix(DataIftVelosA1, 'DataIftVelosA1.xls');
    writematrix(DataIftVelosR1, 'DataIftVelosR1.xls');
end

if exist('image_Antero1', 'var') == 0
    image_Antero2 = uint16(zeros(siz_im01));
    for n11a = 1:size(trajectory_listA2,1)
        ind1 = sub2ind(siz_im01,trajectory_listA2{n11a,1}(:,1),trajectory_listA2{n11a,1}(:,2));
        image_Antero2(ind1) = n11a;
    end
    imwrite(image_Antero2, 'trajectories_Antero.tif');
end

if exist('image_Retro1', 'var') == 0
    image_Retro2 = uint16(zeros(siz_im01));
    for n11b = 1:size(trajectory_listR2,1)
        ind2 = sub2ind(siz_im01,trajectory_listR2{n11b,1}(:,1),trajectory_listR2{n11b,1}(:,2));
        image_Retro2(ind2) = n11b;
    end
    imwrite(image_Retro2, 'trajectories_Retro.tif');
end

% extract parameters
fileID = fopen('parameter_matlab5.txt','a');
script_name1 = 'fla_calciumD5f.m\n';
fprintf(fileID, script_name1);
dt1 = datestr(now,'mmmm dd, yyyy HH:MM:SS\n');
fprintf(fileID, dt1);
px_size1 = 'pixel size = %6.6f\n';
fprintf(fileID, px_size1, px_size);
fps1 = 'fps = %4.2f\n';
fprintf(fileID, fps1, fps);
gcampThresh1a = 'GCaMP Threshold1 = %4.2f\n';
fprintf(fileID, gcampThresh1a, gcampThresh1);
gcampThresh2a = 'GCaMP Threshold2 = %4.2f\n';
fprintf(fileID, gcampThresh2a, gcampThresh2);
duration_speed1 = 'Duration for speed caclucation (px) = %d\n';
fprintf(fileID, duration_speed1, duration_speed);
interval1 = 'Duration for interval analysis (S) = %d\n';
fprintf(fileID, interval1, interval1s);
speedTresh1a = 'Threshold for speed analysis (um/s) = %4.2f\n\n';
fprintf(fileID, speedTresh1a, speedTresh1);

