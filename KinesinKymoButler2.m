function Summary_IFT = KinesinKymoButler2(directory, px_size, fps)

% 8/15/2021
% Updated on 5/2/2022 remove inaccurate trajectories
% Updated on 5/4/2022 calculate IFT intensity from trajectories

% Analayze KAP-GFP kymographs
% Need KymoButler coordinate data "coordinatesA.txt"
% Flagellum position was manually selected as file "kymograph1_labels.png"
% Need the following scripts, "formatKymoButlerC1", 


%% Parameters

% fps = 19.33; % frames per second
% px_size = 0.104667; % pixel size in (um/px)

duration_speed = 5; % (px)
speedTresh1 = 0.5; % (um/s)
speedTresh2 = 4.5; % (um/s) cut out more than this speed

%% Read Kymograph images

% directory01 = dir;
directory01 = dir(directory);
val1 = length(directory01);

for m01 = 1:val1
    if regexp(directory01(m01).name, 'kymograph\dI.tif') > 0
        image_km01 = im2double(imread(directory01(m01).name));
    elseif regexp(directory01(m01).name, 'kymograph\d_labels.png') > 0
        image_label01 = imread(directory01(m01).name);
    elseif regexp(directory01(m01).name, 'coordinatesA.txt') > 0
        coordinatesA1 = importdata(directory01(m01).name);
%     elseif regexp(directory01(m01).name, 'coordinatesR.txt') > 0
%         coordinatesR1 = importdata(directory01(m01).name);
    elseif regexp(directory01(m01).name, 'trajectories_Antero2.tif') > 0
        image_Antero1 = imread(directory01(m01).name);
    end
end

% figure, imshow(imrotate(image_km01, 90), []), title('Kymo');
% figure, imshow(imrotate(image_label01, 90), []), title('label');

clear val1 m01

siz_im01 = size(image_km01);
image_label02 = imbinarize(image_label01);

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
    trajectory_listA1 = formatKymoButlerC1(siz_im01, coordinatesA1);

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

%VisualizeTrajectories1(trajectory_listA2,siz_im01);

%% Baseline correction (only for trajectory selection)

image_km02 = imgaussfilt(image_km01, 1, 'FilterSize', 3); 
image_km03a = image_label02 .* image_km02;
mean_km03a = mean(image_km03a(:), 'omitnan');
image_km03b = image_km03a;
image_label03a = false(siz_im01);

for n04a = 1:size(trajectory_listA2,1)
    image_label03a(trajectory_listA2{n04a,2}) = 1;
end

SE_suare3 = strel('square',3);
image_label03b = imdilate(image_label03a, SE_suare3);
image_km03b(image_label03b) = mean_km03a;
image_km03b = image_km03b .* image_label02;

mean_km03b = mean(image_km03b, 2, 'omitnan');
mean_km03c = smooth(mean_km03b, 5, 'sgolay', 1);
max_meankm03c = max(mean_km03c);
mean_km03d = mean_km03c ./ max_meankm03c;
%figure, plot(mean_km03c);

image_base01 =ones(siz_im01);
image_base01 = image_base01 .* mean_km03d;

image_km04 = image_km02 ./ image_base01;
%figure, imshow(imrotate(image_km04, 90), [])

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

%% IFT analysis 


% Anterograde
trajectoryDataA1 = double.empty(sizeTrajectoryA1(1,1),0);
for n06a = 1:sizeTrajectoryA1(1,1)
    trajectoryDataA1(n06a,1) = n06a; % ID
    trajectoryDataA1(n06a,2) = sum(image_km01(trajectory_listA2{n06a,2})); % total intensity
    trajectoryDataA1(n06a,3) = mean(image_km01(trajectory_listA2{n06a,2})); % mean intensity
    trajectoryDataA1(n06a,4) = mean(trajectory_listA2{n06a,3},'omitnan'); % mean velocity (um/s)
    trajectoryDataA1(n06a,5) = median(trajectory_listA2{n06a,3},'omitnan'); % median velocity (um/s)
    trajectoryDataA1(n06a,11) = sum(image_km04(trajectory_listA2{n06a,2})); % total intensity
    trajectoryDataA1(n06a,12) = mean(image_km04(trajectory_listA2{n06a,2})); % mean intensity
    slow_speed = abs(trajectory_listA2{n06a,3}) < speedTresh1;
    trajectoryDataA1(n06a,13) = sum(slow_speed) / fps; % duration of stopping time

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



%% Remove inaccurate trajectories

trajectoryDataA2 = trajectoryDataA1;
mean_intensity1 = mean(image_km04(image_label02(:)));
temp_Idx3a = trajectoryDataA1(:,12) < mean_intensity1; % based on mean intensity
trajectoryDataA2(temp_Idx3a,:) = [];
temp_Idx3b = trajectoryDataA2(:,4) > speedTresh2 | trajectoryDataA2(:,5) > speedTresh2; % based on IFT speed
trajectoryDataA2(temp_Idx3b,:) = [];

trajectorySummaryA1 = zeros(1,7);
trajectorySummaryA1(1,1) = size(trajectoryDataA2,1); % numbers of IFT trajectories
mean_trajDataA2 = mean(trajectoryDataA2,'omitnan'); 
trajectorySummaryA1(1,2:5) = mean_trajDataA2(1,2:5); % mean IFT intensity, average IFT intensity, speed, median speed
trajectorySummaryA1(1,6) = mean_trajDataA2(1,13);

%% Image for IFT intensity

image_IFTA1 = false(siz_im01);

for m02a = 1:size(trajectoryDataA2,1)
    image_IFTA1(trajectory_listA2{trajectoryDataA2(m02a,1),2}) = 1;
end

% SE_line3 = strel('line',3,90);
SE_disk1 = strel('disk',1);
image_IFTA2 = imdilate(image_IFTA1, SE_disk1);
image_IFTA3 = image_km01 .* image_IFTA2 .* image_label02;

image_bg1 = logical(imcomplement(image_IFTA2) .* image_label02);
image_background1 = image_km01 .* image_bg1;

% figure, imshow(imrotate(image_IFTA1,90))
unitIntA1 = sum(image_IFTA3(:)) / (siz_im01(1,1) / fps) / (mean(sum(image_label02,2)) * px_size);
background1 = sum(image_background1(:)) / (siz_im01(1,1) / fps) / (mean(sum(image_label02,2)) * px_size); %unit background
background2 = mean(image_background1(image_bg1(:))); % mean background (px)
TotalIntensity1 = sum(trajectoryDataA2(:,3)) / (siz_im01(1,1) / fps);

%% Detect IFT trajectories which start from the base or tip

% select anterograde trajectories which start from the base
for n07a1 = 1:size(trajectoryDataA2,1)
    if trajectoryDataA2(n07a1,6) <= trajectoryDataA2(n07a1,10)+5
        trajectoryDataA2(n07a1,14) = 1;
    else
        trajectoryDataA2(n07a1,14) = 0;
    end
end

trajectoryDataA3 = sortrows(trajectoryDataA2, [7,6]);
temp_Idx4a = ~trajectoryDataA3(:,14);
trajectoryDataA3(temp_Idx4a,:) = [];

trajectorySummaryA1(1,7) = size(trajectoryDataA3,1) / (siz_im01(1,1) / fps); % frequency

Summary_IFT = horzcat(trajectorySummaryA1, unitIntA1, background1, background2, TotalIntensity1);

%% Extract data

writematrix(trajectoryDataA2, 'Trajectory_Antero2.xls');

if exist('image_Antero1', 'var') == 0
    image_Antero2 = uint16(zeros(siz_im01));
    for n11a = 1:size(trajectory_listA2,1)
        ind1 = sub2ind(siz_im01,trajectory_listA2{n11a,1}(:,1),trajectory_listA2{n11a,1}(:,2));
        image_Antero2(ind1) = n11a;
    end
    imwrite(image_Antero2, 'trajectories_Antero2.tif');
end

% extract parameters
fileID = fopen('parameter_matlab.txt','a');
script_name1 = 'KinesinKymoButtler1b.m\n';
fprintf(fileID, script_name1);
dt1 = datestr(now,'mmmm dd, yyyy HH:MM:SS\n');
fprintf(fileID, dt1);
px_size1 = 'pixel size = %6.6f\n';
fprintf(fileID, px_size1, px_size);
fps1 = 'fps = %4.2f\n';
fprintf(fileID, fps1, fps);
duration_speed1 = 'Duration for speed caclucation (px) = %d\n';
fprintf(fileID, duration_speed1, duration_speed);
speedTresh1a = 'Threshold for speed analysis (um/s) = %4.2f\n\n';
fprintf(fileID, speedTresh1a, speedTresh1);
speedTresh2a = 'Threshold for speed analysis (um/s) = %4.2f\n\n';
fprintf(fileID, speedTresh2a, speedTresh2);