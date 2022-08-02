function trajectory_list5 = formatKymoButlerC2(siz_img01, coordinates01)



% 9/7/2020
% 9/6/2021 Updated

% This script formats the coordinates from KymoButler and detects IFT
% trajectories.

% Imput needs coordinates made by KymoButler and size of the kymograph.
% This script needs "remove_traj_overlap4.m"



%% Format KymoButler Coordinates for Matlab

size_coordi01 = size(coordinates01,1);
trajectory_list1 = cell(size_coordi01, 1);

% Format shape of coordinates for Matlab
for n01 = 1:size_coordi01
    wStr1 = strrep(coordinates01{n01,1},"{{","[");
    wStr2 = strrep(wStr1,"}}","]");
    wStr3 = strrep(wStr2,"}, {","; ");
    X = str2num(wStr3);
    trajectory_list1{n01,1} = X;
end

clear wStr1 wStr2 wStr3 X

%% Connect trajectories

trajectory_list2 = cell(size_coordi01, 1);

for n02a = 1:size(trajectory_list1,1)
    image_bw01a = false(siz_img01);
    index01 = sub2ind(siz_img01,trajectory_list1{n02a,1}(:,1),trajectory_list1{n02a,1}(:,2));
    image_bw01a(index01) = true;
    % connect points in a trajectory step by step
    for n02b = 1:size(trajectory_list1{n02a,1},1)-1
        x01 = trajectory_list1{n02a,1}(n02b,:);
        y01 = trajectory_list1{n02a,1}(n02b+1,:);
        dist01 = pdist2(x01,y01);
        if dist01 > 1.5
            numPoints = dist01 * 2;
            xLine = linspace(x01(1,2), y01(1,2), numPoints);
            yLine = linspace(x01(1,1), y01(1,1), numPoints);
            rows = round(yLine);
            columns = round(xLine);
            for k = 1 : length(xLine)
                image_bw01a(rows(k), columns(k)) = true;
            end
        end
    end
    %image_bw01b = bwskel(image_bw01a);
    image_bw01b = bwmorph(image_bw01a,'skel',Inf);
    image_bw01c = bwskel(image_bw01b,'MinBranchLength',1);
    image_bw01c(:,1) = false;
    image_bw01c(:,end) = false;
    trajectory_list2{n02a,2} = find(image_bw01c);
    [trajectory_list2{n02a,1}(:,1), trajectory_list2{n02a,1}(:,2)] = ind2sub(siz_img01, trajectory_list2{n02a,2});
    trajectory_list2{n02a,1} = sortrows(trajectory_list2{n02a,1});
end



%% Remove and connect overlapped trajectories  

% Need "remove_traj_overlap1.m" function
trajectory_list3 = remove_traj_overlap4(trajectory_list2, siz_img01);
trajectory_list4 = remove_traj_overlap4(trajectory_list3, siz_img01);


%% Connect two trajectories whose start and end are next each other

trajectory_list5 = trajectory_list4;
cord_start = zeros(1,2);
cord_end = zeros(1,2);

for n04a = 1:size(trajectory_list5,1)
    cord_start(n04a,:) = trajectory_list5{n04a,1}(1,:);
    cord_end(n04a,:) = trajectory_list5{n04a,1}(end,:);
end


for n04b = 1:size(cord_end,1)
    Idx04b = rangesearch(cord_start,cord_end(n04b,:),1.5);
    Idx04b = cell2mat(Idx04b);
    if ~isempty(Idx04b)
        image_bw04a = false(siz_img01);
        if length(Idx04b) > 1
            k1 = dsearchn(cord_start(Idx04b,:),cord_end(n04b,:));
            image_bw04a(trajectory_list5{n04b,2}) = true;
            image_bw04a(trajectory_list5{Idx04b(k1),2}) = true;
            image_bw04b = bwmorph(image_bw04a, 'skel');
            trajectory_list5{n04b,2} = find(image_bw04b);
            [xB1, yB1] = ind2sub(siz_img01,trajectory_list5{n04b,2});
            trajectory_list5{n04b,1} = [xB1, yB1];
            trajectory_list5{n04b,1} = sortrows(trajectory_list5{n04b,1});
            trajectory_list5{Idx04b(k1),1} = []; 
            trajectory_list5{Idx04b(k1),2} = [];
        else
            image_bw04a(trajectory_list5{n04b,2}) = true;
            image_bw04a(trajectory_list5{Idx04b,2}) = true;
            image_bw04b = bwmorph(image_bw04a, 'skel');
            trajectory_list5{n04b,2} = find(image_bw04b);
            [xB1, yB1] = ind2sub(siz_img01,trajectory_list5{n04b,2});
            trajectory_list5{n04b,1} = [xB1, yB1];
            trajectory_list5{n04b,1} = sortrows(trajectory_list5{n04b,1});
            trajectory_list5{Idx04b,1} = []; 
            trajectory_list5{Idx04b,2} = [];
        end
    end 
end
    
zB1 = cellfun('isempty',trajectory_list5(:,1));
trajectory_list5(zB1,:) = [];
    
