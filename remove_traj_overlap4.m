function trajectory_listA = remove_traj_overlap4(trajectory_list0, siz_img01)

% Remove and connect overlapped trajectories
% 7/13/2020
% Update 7/16/2020 Remove a small trajectory
% Update 8/11/2020 Consider the endpoints for removing overlaps
% Update 8/25/2021 Add overlap patterns
%


%%
trajectory_listA = trajectory_list0;
image_bw06a = zeros(siz_img01);

%SE_SQ3 = strel('square',3);
SE_L3a = strel([1 1 1; 1 0 0; 0 0 0]);
SE_L3b = strel([0 0 0; 0 0 1; 1 1 1]);

%%

% Plot all trajectories
for n06a = 1:size(trajectory_listA,1)
    image_bw06b = false(siz_img01);
    %trajectory_list2{n06a,2} = sub2ind(siz_img01,trajectory_list2{n06a,1}(:,1),trajectory_list2{n06a,1}(:,2));
    image_bw06b(trajectory_listA{n06a,2}) = 1;
    image_bw06a = image_bw06a + image_bw06b;
end

% Detect overlap
image_bw06c = imbinarize(image_bw06a, 1);
CCbw06c = bwconncomp(image_bw06c);
%figure, imshow(imrotate(image_bw06a, 90), [])
check1 = cell(CCbw06c.NumObjects,4);

for n06b = 1:CCbw06c.NumObjects
    % Detect overlap trajectories
    image_bw06A1 = false(siz_img01);
    image_bw06B1 = false(siz_img01);
    overlap_info1 = zeros(1,5); % for overlap trajectories
    n06c = 1;
    for n06d = 1:size(trajectory_listA,1)
        if find(trajectory_listA{n06d,2} == CCbw06c.PixelIdxList{1,n06b}(1,1))
            overlap_info1(n06c,1) = n06d;
            n06c = n06c+1; 
        end
    end
    check1{n06b,1} = CCbw06c.PixelIdxList{1,n06b};

    
    if size(overlap_info1,1) > 1
        image_bw06A1(trajectory_listA{overlap_info1(1,1),2}) = true;
        image_bw06B1(trajectory_listA{overlap_info1(2,1),2}) = true;
        image_bw06AB1 = image_bw06A1 + image_bw06B1; 
        image_bw06AB2 = imbinarize(image_bw06AB1, 1);
        image_bw06AB2b = bwmorph(image_bw06AB2, 'bridge');
        %image_bw06AB3 = imdilate(image_bw06AB2, SE_SQ3);
        image_bw06AB3 = bwmorph(image_bw06AB2b, 'skel');
        image_bw06AB3b = bwareaopen(image_bw06AB3, 4);
        image_bw06AB4 = image_bw06AB3b | image_bw06AB2;
        CCbw06AB4 = bwconncomp(image_bw06AB4);
        check1{n06b,4} = CCbw06AB4.NumObjects;
        
        % detect endpoints of trajectories
        x04a1 = trajectory_listA{overlap_info1(1,1),1}(1,1);
        y04a1 = trajectory_listA{overlap_info1(1,1),1}(1,2);
        point_startA = sub2ind(siz_img01, x04a1, y04a1);
        x04a2 = trajectory_listA{overlap_info1(1,1),1}(end,1);
        y04a2 = trajectory_listA{overlap_info1(1,1),1}(end,2);
        point_endA = sub2ind(siz_img01, x04a2, y04a2);
        x04b1 = trajectory_listA{overlap_info1(2,1),1}(1,1);
        y04b1 = trajectory_listA{overlap_info1(2,1),1}(1,2);
        point_startB = sub2ind(siz_img01, x04b1, y04b1);
        x04b2 = trajectory_listA{overlap_info1(2,1),1}(end,1);
        y04b2 = trajectory_listA{overlap_info1(2,1),1}(end,2);
        point_endB = sub2ind(siz_img01, x04b2, y04b2);
        
        for n06f = 1:CCbw06AB4.NumObjects
            % Detect overlap point
                image_bw06d = false(siz_img01);
                image_bw06e = false(siz_img01);
                [x05, y05] = ind2sub(siz_img01,CCbw06AB4.PixelIdxList{1,n06f});
                [M1,I1] = min(x05,[],'linear');
                [M2,I2a] = max(flipud(x05),[],'linear');
                I2b = length(x05) - I2a +1;
                overlap_first_idx1 = sub2ind(siz_img01, M1,y05(I1));
                overlap_end_idx1 = sub2ind(siz_img01, M2,y05(I2b));

                image_bw06d(overlap_first_idx1) = true; % First overlap point
                image_bw06e(overlap_end_idx1) = true; % Last overlap point
                image_bw06f = imdilate(image_bw06d, SE_L3a); 
                image_bw06g = imdilate(image_bw06e, SE_L3b);
                
                % index of overlap positin in trajectory list
                OverlapIdxA1 = intersect(find(trajectory_listA{overlap_info1(1,1),1}(:,1) == M1), find(trajectory_listA{overlap_info1(1,1),1}(:,2) == y05(I1))); % index of overlap start in trajectory A
                OverlapIdxB1 = intersect(find(trajectory_listA{overlap_info1(2,1),1}(:,1) == M1), find(trajectory_listA{overlap_info1(2,1),1}(:,2) == y05(I1))); % index of overlap start in trajectory B
                OverlapIdxA2 = intersect(find(trajectory_listA{overlap_info1(1,1),1}(:,1) == M2), find(trajectory_listA{overlap_info1(1,1),1}(:,2) == y05(I2b))); % index of overlap end in trajectory A
                OverlapIdxB2 = intersect(find(trajectory_listA{overlap_info1(2,1),1}(:,1) == M2), find(trajectory_listA{overlap_info1(2,1),1}(:,2) == y05(I2b))); % index of overlap end in trajectory B

                % index of merge and branch points for trajectory A & B
                point_mergeA = find(flipud(image_bw06A1&image_bw06f)); % merge point A
                point_mergeB = find(flipud(image_bw06B1&image_bw06f)); % merge point B
                point_branchA = find(flipud(image_bw06A1&image_bw06g)); % branch point A
                point_branchB = find(flipud(image_bw06B1&image_bw06g)); % branch point B

                % information of overlap and end points
                overlap_info1(1,2) = ~isempty(point_mergeA); % merge point A
                overlap_info1(1,3) = ~isempty(point_branchA); % merge point B
                overlap_info1(2,2) = ~isempty(point_mergeB); % branch point A
                overlap_info1(2,3) = ~isempty(point_branchB); % branch point B
                overlap_info1(1,4) = ~isempty(find(CCbw06AB4.PixelIdxList{1,n06f} == point_startA, 1)); % start point A in the overlap
                overlap_info1(1,5) = ~isempty(find(CCbw06AB4.PixelIdxList{1,n06f} == point_endA, 1)); % end point A in the overlap
                overlap_info1(2,4) = ~isempty(find(CCbw06AB4.PixelIdxList{1,n06f} == point_startB, 1)); % start point B in the overlap
                overlap_info1(2,5) = ~isempty(find(CCbw06AB4.PixelIdxList{1,n06f} == point_endB, 1)); % end point B in the overlap
                
                check1{n06b,2} = overlap_info1;

            % 1.Cross (1,1,1,1)
            if overlap_info1(1,2) == 1 && overlap_info1(2,2) == 1 && overlap_info1(1,3) == 1 && overlap_info1(2,3) == 1
                if point_mergeA < point_mergeB % remove overlap from trajectory B
                    x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                    y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                    trajectory_listA{overlap_info1(2,1),2} = [];
                    trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);

                elseif point_mergeA > point_mergeB % remove overlap from trajectory A
                    x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1));
                    y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                    trajectory_listA{overlap_info1(1,1),2} = [];
                    trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                end
                check1{n06b,3}(n06f,1) = 1;
                
            % 2.Merge (1,1,1,0)
            elseif overlap_info1(1,2) == 1 && overlap_info1(2,2) == 1 && overlap_info1(1,3) == 1 && overlap_info1(2,3) == 0
                if overlap_info1(2,5) == 1 % merge
                    if point_mergeA < point_mergeB % remove overlap from trajectory B
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2));
                    elseif point_mergeA > point_mergeB % remove overlap from trajectory A & connect trajectory first half of B and last half of A
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2));
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA1+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA1+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    end
                    check1{n06b,3}(n06f,1) = 2.1;
                else % cross
                    if point_mergeA < point_mergeB % remove overlap from trajectory B
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    elseif point_mergeA > point_mergeB % remove overlap from trajectory A
                        x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                    end
                    check1{n06b,3}(n06f,1) = 2.2;
                end

            % 3.Merge (1,1,0,1)
            elseif overlap_info1(1,2) == 1 && overlap_info1(2,2) == 1 && overlap_info1(1,3) == 0 && overlap_info1(2,3) == 1
                if overlap_info1(1,5) == 1 % merge
                    if point_mergeA < point_mergeB % remove overlap from trajectory B & connect trajectory first half of A and last half of B
                        x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB1+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB1+1:end,2));
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2));
                    elseif point_mergeA > point_mergeB % remove overlap from trajectory A
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2));
                    end
                    check1{n06b,3}(n06f,1) = 3.1;
                else % cross
                    if point_mergeA < point_mergeB % remove overlap from trajectory B
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    elseif point_mergeA > point_mergeB % remove overlap from trajectory A
                        x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                    end
                    check1{n06b,3}(n06f,1) = 3.2;
                end
                

            % 4.Branch (1,0,1,1)
            elseif overlap_info1(1,2) == 1 && overlap_info1(2,2) == 0 && overlap_info1(1,3) == 1 && overlap_info1(2,3) == 1
                if overlap_info1(2,4) == 1 % branch
                    if point_branchA < point_branchB % remove first half and overlap from trajectory A & connect first half of A nad last half of B
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA2,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA2,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    elseif point_branchA > point_branchB % remove overlap from trajectory B
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                    end
                    check1{n06b,3}(n06f,1) = 4.1;
                else % merge
                    if point_branchA > point_branchB % remove overlap from trajectory B
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    elseif point_branchA < point_branchB % remove overlap from trajectory A
                        x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                    end
                    check1{n06b,3}(n06f,1) = 4.2;
                end
                    
            % 5.Union (1,0,1,0)
            elseif overlap_info1(1,2) == 1 && overlap_info1(2,2) == 0 && overlap_info1(1,3) == 1 && overlap_info1(2,3) == 0
                if overlap_info1(2,4) == 1 && overlap_info1(2,5) == 1 % union
                    trajectory_listA{overlap_info1(2,1),1} = 0;
                    trajectory_listA{overlap_info1(2,1),2} = 0;
                    check1{n06b,3}(n06f,1) = 5.1;
                elseif overlap_info1(2,4) == 1 % branch
                    if x04a2 > x04b2 % remove first half and overlap from trajectory A & connect first half of A nad last half of B 
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA2,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA2,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    elseif x04a2 < x04b2 % remove overlap from trajectory B
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                    end                    
                    check1{n06b,3}(n06f,1) = 5.2;
                elseif overlap_info1(2,5) == 1 % merge
                    if x04a1 > x04b1 % remove overlap from trajectory B
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2));
                    elseif x04a1 < x04b1 % remove overlap from trajectory A & connect trajectory first half of B and last half of A
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2));
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA1+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA1+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    end
                    check1{n06b,3}(n06f,1) = 5.3;
                else % cross
                    if x04a1 > x04b1 % remove overlap from trajectory B
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    elseif x04a1 < x04b1 % remove overlap from trajectory A
                        x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                    end
                    check1{n06b,3}(n06f,1) = 5.4;
                end
                
            % 6.Union (1,0,0,1)
            elseif overlap_info1(1,2) == 1 && overlap_info1(2,2) == 0 && overlap_info1(1,3) == 0 && overlap_info1(2,3) == 1
                if overlap_info1(2,4) == 1 && overlap_info1(1,5) == 1 % union
                    x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB1+1:end,1));
                    y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB1+1:end,2));
                    trajectory_listA{overlap_info1(1,1),2} = [];
                    trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                    trajectory_listA{overlap_info1(2,1),1} = 0;
                    trajectory_listA{overlap_info1(2,1),2} = 0;
                    check1{n06b,3}(n06f,1) = 6.1;
                elseif overlap_info1(2,4) == 1 % branch
                    if x04a2 > x04b2 % remove first half and overlap from trajectory A & connect first half of A nad last half of B
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA2,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA2,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    elseif x04a2 < x04b2 % remove overlap from trajectory B
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                    end
                    check1{n06b,3}(n06f,1) = 6.2;
                elseif overlap_info1(1,5) == 1 % merge
                    if x04a1 > x04b1 % remove overlap from trajectory B & connect trajectory first half of A and last half of B
                        x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB1+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB1+1:end,2));
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2));
                    elseif x04a1 < x04b1 % remove overlap from trajectory A
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2));
                    end
                else %cross
                    if x04a1 > x04b1 % remove overlap from trajectory B
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    elseif x04a1 < x04b1 % remove overlap from trajectory A
                        x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                    end
                    check1{n06b,3}(n06f,1) = 6.4;
                end
                
            % 7.Branch (0,1,1,1)
            elseif overlap_info1(1,2) == 0 && overlap_info1(2,2) == 1 && overlap_info1(1,3) == 1 && overlap_info1(2,3) == 1
                if overlap_info1(1,4) == 1 % branch
                    if point_branchA < point_branchB % remove overlap from trajectory A
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                    elseif point_branchA > point_branchB % remove overlap from trajectory B & connect trajectory first half of B and last half of A
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB2,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB2,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                    end
                    check1{n06b,3}(n06f,1) = 7.1;
                else % cross
                    if point_branchA > point_branchB % remove overlap from trajectory B
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    elseif point_branchA < point_branchB % remove overlap from trajectory A
                        x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                    end
                    check1{n06b,3}(n06f,1) = 7.2;
                end
                
            % 8.Union (0,1,1,0)
            elseif overlap_info1(1,2) == 0 && overlap_info1(2,2) == 1 && overlap_info1(1,3) == 1 && overlap_info1(2,3) == 0
                if overlap_info1(1,4) == 1 && overlap_info1(2,5) == 1 % union
                    x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA1+1:end,1));
                    y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA1+1:end,2));
                    trajectory_listA{overlap_info1(2,1),2} = [];
                    trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    trajectory_listA{overlap_info1(1,1),1} = 0;
                    trajectory_listA{overlap_info1(1,1),2} = 0;
                    check1{n06b,3}(n06f,1) = 8.1;
                elseif overlap_info1(1,4) == 1 % merge
                    if x04a1 > x04b1 % remove overlap from trajectory B
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2));
                    elseif x04a1 < x04b1 % remove overlap from trajectory A & connect trajectory first half of B and last half of A
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2));
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA1+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA1+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    end
                    check1{n06b,3}(n06f,1) = 8.2;
                elseif overlap_info1(2,5) == 1 % branch
                    if x04a2 > x04b2 % remove overlap from trajectory A
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                    elseif x04a2 < x04b2 % remove overlap from trajectory B & connect trajectory first half of B and last half of A
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB2,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB2,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                    end
                    check1{n06b,3}(n06f,1) = 8.3;
                else % cross
                    if x04a1 > x04b1 % remove overlap from trajectory B
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    elseif x04a1 < x04b1 % remove overlap from trajectory A
                        x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                    end
                    check1{n06b,3}(n06f,1) = 8.4;
                end
            % 9.Union (0,1,0,1)
            elseif overlap_info1(1,2) == 0 && overlap_info1(2,2) == 1 && overlap_info1(1,3) == 0 && overlap_info1(2,3) == 1
                if overlap_info1(1,4) == 1 && overlap_info1(1,5) == 1 % union
                    trajectory_listA{overlap_info1(1,1),1} = 0;
                    trajectory_listA{overlap_info1(1,1),2} = 0;
                    check1{n06b,3}(n06f,1) = 9.1;
                elseif overlap_info1(1,4) == 1 % branch
                    if x04a2 > x04b2 % remove overlap from trajectory A
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                    elseif x04a2 < x04b2 % remove overlap from trajectory B & connect trajectory first half of B and last half of A
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB2,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB2,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                    end
                    check1{n06b,3}(n06f,1) = 9.2;
                elseif overlap_info1(1,5) == 1 % merge
                    if x04a1 > x04b1 % remove overlap from trajectory B
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2));
                    elseif x04a1 < x04b1 % remove overlap from trajectory A & connect trajectory first half of B and last half of A
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2));
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA1+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA1+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    end
                    check1{n06b,3}(n06f,1) = 9.3;
                else % cross
                    if x04a1 > x04b1 % remove overlap from trajectory B
                        x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    elseif x04a1 < x04b1 % remove overlap from trajectory A
                        x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1));
                        y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                    end
                    check1{n06b,3}(n06f,1) = 9.4;
                end

            % 10. Merge (1,1,0,0)
            elseif overlap_info1(1,2) == 1 && overlap_info1(2,2) == 1 && overlap_info1(1,3) == 0 && overlap_info1(2,3) == 0
                if overlap_info1(1,5) == 1 && overlap_info1(2,5) == 1 % merge
                    if point_mergeA < point_mergeB
                        trajectory_listA{overlap_info1(2,1),2} = [];
                        trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2));
                    elseif point_mergeA > point_mergeB
                        trajectory_listA{overlap_info1(1,1),2} = [];
                        trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2));
                    end
                    check1{n06b,3}(n06f,1) = 10.1;
                else
                    check1{n06b,3}(n06f,1) = 10.2;
                end
            % 11. Union (1,0,0,0)
            elseif overlap_info1(1,2) == 1 && overlap_info1(2,2) == 0 && overlap_info1(1,3) == 0 && overlap_info1(2,3) == 0
                if overlap_info1(1,5) == 1 && overlap_info1(2,4) == 1 && overlap_info1(2,5) == 1 % union
                    trajectory_listA{overlap_info1(2,1),1} = 0;
                    trajectory_listA{overlap_info1(2,1),2} = 0;
                    check1{n06b,3}(n06f,1) = 11.1;
                else
                    check1{n06b,3}(n06f,1) = 11.2;
                end
            % 12. Union (0,1,0,0)
            elseif overlap_info1(1,2) == 0 && overlap_info1(2,2) == 1 && overlap_info1(1,3) == 0 && overlap_info1(2,3) == 0
                if overlap_info1(1,4) == 1 && overlap_info1(1,5) == 1 && overlap_info1(2,5) == 1 % union
                    trajectory_listA{overlap_info1(1,1),1} = 0;
                    trajectory_listA{overlap_info1(1,1),2} = 0;
                    check1{n06b,3}(n06f,1) = 12.1;
                else
                    check1{n06b,3}(n06f,1) = 12.2;
                end  
            % 13. Union (0,0,1,0)
            elseif overlap_info1(1,2) == 0 && overlap_info1(2,2) == 0 && overlap_info1(1,3) == 1 && overlap_info1(2,3) == 0
                if overlap_info1(1,4) == 1 && overlap_info1(2,4) == 1 && overlap_info1(2,5) == 1 % union
                    trajectory_listA{overlap_info1(2,1),1} = 0;
                    trajectory_listA{overlap_info1(2,1),2} = 0;
                    check1{n06b,3}(n06f,1) = 13.1;
                else
                    check1{n06b,3}(n06f,1) = 13.2;
                end  
            % 14. Union (0,0,0,1)
            elseif overlap_info1(1,2) == 0 && overlap_info1(2,2) == 0 && overlap_info1(1,3) == 0 && overlap_info1(2,3) == 1
                if overlap_info1(1,4) == 1 && overlap_info1(1,5) == 1 && overlap_info1(2,4) == 1 % union
                    trajectory_listA{overlap_info1(1,1),1} = 0;
                    trajectory_listA{overlap_info1(1,1),2} = 0;
                    check1{n06b,3}(n06f,1) = 14.1;
                else
                    check1{n06b,3}(n06f,1) = 14.2;
                end 
            % 15. Start? (0,0,1,1)
            elseif overlap_info1(1,2) == 0 && overlap_info1(2,2) == 0 && overlap_info1(1,3) == 1 && overlap_info1(2,3) == 1
                if overlap_info1(1,4) == 1 && overlap_info1(2,4) == 0 % remove overlap from trajectory A 
                    x06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,1), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,1));
                    y06 = vertcat(trajectory_listA{overlap_info1(1,1)}(1:OverlapIdxA1-1,2), trajectory_listA{overlap_info1(1,1)}(OverlapIdxA2+1:end,2));
                    trajectory_listA{overlap_info1(1,1),2} = [];
                    trajectory_listA{overlap_info1(1,1),2} = sub2ind(siz_img01, x06, y06);
                    check1{n06b,3}(n06f,1) = 15.1;
                elseif overlap_info1(1,4) == 0 && overlap_info1(2,4) == 1 % remove overlap from trajectory B
                    x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                    y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                    trajectory_listA{overlap_info1(2,1),2} = [];
                    trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    check1{n06b,3}(n06f,1) = 15.2;
                elseif overlap_info1(1,4) == 1 && overlap_info1(2,4) == 1
                    x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                    y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                    trajectory_listA{overlap_info1(2,1),2} = [];
                    trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    check1{n06b,3}(n06f,1) = 15.3;    
                elseif overlap_info1(1,4) == 0 && overlap_info1(2,4) == 0
                    x06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,1), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,1));
                    y06 = vertcat(trajectory_listA{overlap_info1(2,1)}(1:OverlapIdxB1-1,2), trajectory_listA{overlap_info1(2,1)}(OverlapIdxB2+1:end,2));
                    trajectory_listA{overlap_info1(2,1),2} = [];
                    trajectory_listA{overlap_info1(2,1),2} = sub2ind(siz_img01, x06, y06);
                    check1{n06b,3}(n06f,1) = 15.4;
                else
                    check1{n06b,3}(n06f,1) = 15.5;
                end 
                
            else
                lengthA = length(trajectory_listA{overlap_info1(1,1),2});
                lengthB = length(trajectory_listA{overlap_info1(2,1),2});
                if lengthA > lengthB
                    trajectory_listA{overlap_info1(2,1),2} = [];
                else
                    trajectory_listA{overlap_info1(1,1),2} = [];
                end
                print3 = [num2str(n06b), '	', num2str(overlap_info1(1,1)), '	', num2str(overlap_info1(2,1))];
                disp(print3)
                
                check1{n06b,3}(n06f,1) = 17;
            end


            for n06e = 1:2
                [x07, y07] = ind2sub(siz_img01,trajectory_listA{overlap_info1(n06e),2});
                trajectory_listA{overlap_info1(n06e),1} = [x07, y07];
                trajectory_listA{overlap_info1(n06e),1} = sortrows(trajectory_listA{overlap_info1(n06e),1});
            end  
        end             
    end
end

x3a = cellfun('length',trajectory_listA(:,2));
x3b = x3a < 6;
trajectory_listA(x3b,:) = [];