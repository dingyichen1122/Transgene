% This code serves to calculate the displacement and strain on zebrafish
% embryos in 2D images. The idea of this code originates from the paper:
% JCI Insight, 2019, DOI: 10.1172/jci.insight.125362. Yichen incorporated
% the multi-threshold segmentation, PCA analysis and interactive analysis
% to efficiently extract the principal directions and calculate the strain
% in 2D time-lapse images. This new version incorporates a group of fixed
% PCA axes, allowing for segmentation, division and quantification from
% end-diastole to end-systole, bypassing the registration between two
% conditons.
%
% Disclaimer: This code was initiated by Yichen Ding on 11/15/2019. For any
% information, please reach out to Yichen at ycding@g.ucla.edu.
%
% Copyright (c) 2019, Yichen Ding
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the
%       distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% Initiate the parameters
clear; close all force; clc

% parameters of the data set
row = 512;  % y = 512
col = 512;  % x = 512
ind = 55;   % num of images
cycle = 5;  % num of cycles  
GroupNo = 24;    % num of groups / sections >= 6    to determine the circ and trans direction
thresholdLevel = 5; % set the threshold level of segmentation
time_series = 0;    % 0: dia -> sys; 1: time series

% raw data root folder
RawPath = '.\test';   % data path
SavePath = strcat(RawPath, '\output');  % output folder
if ~exist(SavePath, 'dir')
    mkdir([SavePath '\seg']); % create the output folder
end
% load images
sample = zeros(col, row, ind);
list = dir(strcat(RawPath, '\*.tif'));
for i = 1:ind
    disp(['loading image # ' num2str(i)]);
    sample(:,:,i) = imread(strcat(list(1).folder, '\', list(i).name));
end
clear list;

%% Perform multi-threshold segmentation
% pick the first image as a representative
thresh = multithresh(sample(:,:,1), thresholdLevel); % find the different intensity levels
disp(['threshold = ' num2str(thresh)]);
% choose the level for segmentation
% initate the while loop
move = 'n'; % input 1 to stay in the loop
while move ~= 'y'
    img = sample(:,:,1);
    question_1 = {'Input the level # for segmentation', 'move or not? y/n'};
    answer = inputdlg(question_1);    % chosen one from above; 2 or 3 for current slice
    img = (img >= thresh(str2num(answer{1})));
    imshow(img,[]);
    % prepare for next round
    move = answer{2};
end
hold on;
% segment all the images with the same threshold in sample
bw_sample = (sample >= thresh(str2num(answer{1})));
% Morphological operations
bw_inner = false(size(bw_sample));
for i = 1:ind % ind
    temp = 0;
    buff = zeros(row);
    bw_sample(:,:,i) = bwmorph(bw_sample(:,:,i), 'majority', 'bothat');
    [B,L,N,A] = bwboundaries(bw_sample(:,:,i));  % B for boundary, A for adjacency
    if (nnz(A(:,1)) > 0)
        for k = find(A(:,1))'
            if size(B{k}, 1) > temp
                boundary = B{k};
                temp = size(B{k}, 1);
            end
        end
    end
    index = sub2ind(size(buff), boundary(:,1), boundary(:,2));
    buff(index) = 1;
    bw_inner(:,:,i) = buff;
end
mask_sample = sample .* bw_sample;  % gray scale images in the segmented regions
clear sample bw_sample;
disp('Segmentation DONE');

%% Perform PCA to find the central axis
% extract the coordinates of nonzero pixels
[r,c] = ind2sub(size(bw_inner(:,:,1)), find(bw_inner(:,:,1)));  % find nonzero pixels, and reassign the coordinates
coordinate = cat(2, c, r);  % combine both arrays; c for x; r for y; [c,r] = [x,y]
coord_mean = mean(coordinate);    % get the coordinates for each column
% remove the mean of each column; no need for standardization due to the same unit
deCent = coordinate - repmat(coord_mean, size(coordinate, 1), 1);
[vector, ~, value] = pca(deCent);   % vector is eigenvector, value is eigenvalue
% x = x0+r1*e1x; y = y0+r1*e1y
plot([coord_mean(1), coord_mean(1)+value(1)*vector(1,1)],[coord_mean(2), coord_mean(2)+value(1)*vector(2,1)],...
    'Color','r','LineWidth',3);
plot([coord_mean(1), coord_mean(1)+value(2)*vector(1,2)],[coord_mean(2), coord_mean(2)+value(2)*vector(2,2)],...
    'Color','g','LineWidth',3);
clear r c;
disp('PCA analysis DONE');

%% Generate the boundary for each group
theta = 2*pi/GroupNo;   % rotation angle
% theta has no sign, counter-clockwise rotation, but it appears clockwise
% in Matlab due to the x-y direction and original points location in matrix
rotationA = [cos(theta), -sin(theta); sin(theta), cos(theta)];   % rotation matrix for P
rotationB = [1-cos(theta), sin(theta); -sin(theta), 1-cos(theta)];  % rotation matrix for origin
origin = coord_mean';   % original point
P = zeros(2,GroupNo);   % points on different lines
% initiate the starting point
P(:,1) = [origin(1)+ value(1)* vector(1,1); origin(2)+ value(1)* vector(2,1)];  % 1st point on the main axis
% define the boundary of different sections
for i = 2:GroupNo
    % (x',y') = rotationA*(x,y) + rotationB*(x0,y0)
    P(:,i) = rotationA * P(:,i-1) + rotationB * origin;
    if i == GroupNo
        plot([origin(1), P(1,i)], [origin(2), P(2,i)], 'r--');
    else
        plot([origin(1), P(1,i)], [origin(2), P(2,i)], 'w-');
    end
end
hold off;
% generate directional vectors for cross product
P = [P, P(:,1)]; % move P1 to the end for looping
V = P - repmat(origin, 1, size(P,2));   % directional vectors for each line
V = [V; zeros(1, size(V, 2))];  % cross product requires 3d arrays
disp('Boundary DONE');

%% Assign nonzero pixels into coresponding section
% myocardium
[rz, cz, iz] = ind2sub(size(mask_sample), find(mask_sample));
coord_sample = cat(2, cz, rz, iz); % coordinates of nonzero pixels
coord_sample = coord_sample'; % column vector
% endocardium
[ri, ci, ii] = ind2sub(size(bw_inner), find(bw_inner));
coord_inner = cat(2, ci, ri, ii);
coord_inner = coord_inner';
% generate vectors from origin to target points
V_target = coord_sample(1:2,:) - repmat(origin, 1, size(coord_sample,2));
V_target = [V_target; zeros(1, size(V_target, 2))]; % remove the ind for cross product
V_inner = coord_inner(1:2,:) - repmat(origin, 1, size(coord_inner,2));
V_inner = [V_inner; zeros(1, size(V_inner, 2))];
% calculate the cross product between V_target and V
% initiate a matrix for all segmented sections
section = zeros(col, row, ind);
section_inner = zeros(col, row, ind);
for i = 1:GroupNo
    disp(['relocate pixels to section # = ' num2str(i)]);
    % colone vectors along the boundary
    V_1 = repmat(V(:,i), 1, size(V_target, 2));
    V_2 = repmat(V(:,i+1), 1, size(V_target, 2));
    V_3 = repmat(V(:,i), 1, size(V_inner, 2));
    V_4 = repmat(V(:,i+1), 1, size(V_inner, 2));
    % V_target x V
    cProduct_1 = cross(V_target, V_1, 1);
    cProduct_2 = cross(V_target, V_2, 1);
    cProduct_3 = cross(V_inner, V_3, 1);
    cProduct_4 = cross(V_inner, V_4, 1);
    % determine the direction of cProduct; the original point of a image in Matlab is at top left
    % find the overlap region between two lines
    intersection = find((cProduct_1(3,:) <= 0) .* (cProduct_2(3,:) > 0));   % return the colunm number of overlap region
    seq = coord_sample(:, intersection); % coordinates for the overlapping region
    index = sub2ind(size(section), seq(2,:), seq(1,:), seq(3,:));
    section(index) = i;
    % find the overlap for endocardium
    inter_inner = find((cProduct_3(3,:) <= 0) .* (cProduct_4(3,:) > 0 ));
    seq_inner = coord_inner(:, inter_inner);
    index_inner = sub2ind(size(section_inner), seq_inner(2,:), seq_inner(1,:), seq_inner(3,:));
    section_inner(index_inner) = i;
end
% save the divided section
for i = 1:ind
    filename = [SavePath '\seg\S_' num2str(i) '.tif'];
    imwrite(section_inner(:,:,i), filename, 'tif');
end
% double check the direction
test = section(:,:,1);
test(test ~= 0) = 1;
assert(isequal(test, mask_sample(:,:,1)>0), 'error: Reconstruction WRONG');
figure; imshow(test,[]);  % for testing
title('reconstructed image after division and morphological change');
clear rz cz iz cProduct_1 cProduct_2 cProduct_3 cProduct_4 V_1 V_2 V_3 V_4;
disp('Division DONE');

%% Calculate the centroid of each section in each 2d image
[x,y,~] = meshgrid(1:col, 1:row, 1:ind);
% num of nonzero pixels in each 2D section in a single image
x_cent = zeros(ind, GroupNo);   % x-y coordinates
y_cent = zeros(ind, GroupNo);
% find the centroid for each section in each 2D image
for i = 1:GroupNo
    disp(['calculating the central position for # ' num2str(i)]);
    mask_sample_section = (section == i).* mask_sample;
    mass_x = mask_sample_section.*x;   % mass * position
    mass_y = mask_sample_section.*y;
    for j = 1:ind
        if any(mask_sample_section(:,:,j), 'all')
            total = mask_sample_section(:,:,j);
        else
            total = mask_sample_section(:,:,j-1);
            mask_sample_section(:,:,j) = mask_sample_section(:,:,j-1);
            mass_x(:,:,j) = mass_x(:,:,j-1);
            mass_y(:,:,j) = mass_y(:,:,j-1);
        end
        mass = sum(total(:));
        x_cent(j,i) = sum(sum(mass_x(:,:,j)))/mass;
        y_cent(j,i) = sum(sum(mass_y(:,:,j)))/mass;
    end
end
disp('Center determination DONE');

%% Perform PCA to find the circumferential and transmural direction of each section
vector_section = zeros(2,2,GroupNo);    % 2 coords x 2 vectors x Group
range = zeros(ind, 2, GroupNo); % range along pca1 and 2 axis of different sections
coord_pca = zeros(ind, 2, GroupNo); % pca coordinates, ind x 2 coords x Group
figure;
for i = 1:GroupNo
    disp(['pca analysis of section # ', num2str(i)]);
    % norm of V along the boundary
    norm_1 = V(1:2,i) ./ norm(V(1:2,i));
    norm_2 = V(1:2,i+1) ./ norm(V(1:2,i+1));
    % determine pca axes
    for j = 1:ind
        if j == 1
            % preprocess data before PCA
            [r1, c1] = ind2sub(size((section(:,:,1) == i)), find((section(:,:,1) == i)));
            coord_ref = cat(2, c1, r1);
            coord_ref_mean = mean(coord_ref);
            deCent_ref = coord_ref - repmat(coord_ref_mean, size(coord_ref, 1), 1);
            [eig_vec, ~, eig_val] = pca(deCent_ref); % norm of vectors is 1
        end
        % determine circ and tran directions
        cos11 = abs(dot(eig_vec(:,1), norm_1));
        cos21 = abs(dot(eig_vec(:,2), norm_1));
        cos12 = abs(dot(eig_vec(:,1), norm_2));
        cos22 = abs(dot(eig_vec(:,2), norm_2));
        cosRef = dot(norm_1, norm_2);
        if cos11 >= cosRef && cos12 >= cosRef
            vector_section(:,1,i) = eig_vec(:,2);
            vector_section(:,2,i) = eig_vec(:,1);
            value_section(1) = eig_val(2);
            value_section(2) = eig_val(1);
        else
            vector_section(:,1,i) = eig_vec(:,1);
            vector_section(:,2,i) = eig_vec(:,2);
            value_section(1) = eig_val(1);
            value_section(2) = eig_val(2);
        end
        % find the contour and area of each section on the each plane
        section_cont = (section_inner(:,:,j) == i);
        % extract coordinates of each section on the each plane               
        [ri, ci] = ind2sub(size(section_cont), find(section_cont));
        coord_ind = cat(2, ci, ri);
        coord_ind_mean = mean(coord_ind);
        deCent_ind = coord_ind - repmat(coord_ind_mean, size(coord_ind, 1), 1);
        % new coordinates in the pca system
        coord_ind_new = deCent_ind * vector_section(:,:,i);
        % calculate the length along pca2 axis
        int_1 = round(coord_ind_new);
        new_1 = sortrows(int_1, [1,2], {'ascend', 'ascend'});
        new_2 = sortrows(int_1, [1,2], {'ascend', 'descend'});
        length_1 = new_1 - new_2;
        range(j,2,i) = max(1, max(abs(length_1(:,2)))); % use inner length as the reference
        % calculate the length along pca1 axis
        int_2 = round(coord_ind_new(:, [2,1])); % change the column sequence, and then fix elements to integers
        new_3 = sortrows(int_2, [1,2], {'ascend', 'ascend'});   % The 1st colunm is sorted, and then sort the 2nd 'ascend'
        new_4 = sortrows(int_2, [1,2], {'ascend', 'descend'});  % sort the 2nd 'descend'. The length is along pca1 axis
        length_2 = new_3 - new_4;   % quantify the length perpendicular to pca2
        range(j,1,i) = max(1, max(abs(length_2(:,2))));
    end
    % transfer x_cent and y_cent to the pca coordinates
    coord_cent = [x_cent(:,i), y_cent(:,i)];
    coord_cent_mean = mean(coord_cent);
    deCent_coord_cent = coord_cent - repmat(coord_cent_mean, size(coord_cent, 1), 1);
    coord_pca(:,:,i) = deCent_coord_cent * vector_section(:,:,i); % ind x 2 x Group
    % plot the representative image
    subplot(2, ceil(GroupNo/2), i);
    imshow(section(:,:,1)==i,[]);
    hold on;
    plot([coord_ref_mean(1), coord_ref_mean(1)+value_section(1)*vector_section(1,1,i)],...
        [coord_ref_mean(2), coord_ref_mean(2)+value_section(1)*vector_section(2,1,i)], 'r');
    plot([coord_ref_mean(1), coord_ref_mean(1)+value_section(2)*vector_section(1,2,i)],...
        [coord_ref_mean(2), coord_ref_mean(2)+value_section(2)*vector_section(2,2,i)], 'g');
    title(['section # ' num2str(i)]);
    hold off;
end
clear ri ci;
disp('Circumferential direction DONE');

%% Calculate the velocity and acceleration of each section
% compute the velocity or displacement of each center along pca1 and pca2 axes
v = diff(coord_pca, 1, 1);
v_c = permute(v(:,1,:), [1,3,2]);
v_t = permute(v(:,2,:), [1,3,2]);
% compute the accelaration of each center along pca1 and pca2 axes
a = diff(coord_pca, 2, 1);
a_c = permute(a(:,1,:), [1,3,2]);
a_t = permute(a(:,2,:), [1,3,2]);
% determine the direction to extract the systole and diastole
choice = {num2str(GroupNo)}; % initiate the parameter 'choice'
move = {'n'};
while move{1} ~= 'y'
    % plot the selected section along circ-axis
    figure(4); subplot(2,1,1);
    plot(a_c(:,str2num(choice{1})), 'r');
    hold on; grid on;
    plot(v_c(:,str2num(choice{1})), 'b');
    legend('Acceleration','Velocity','Location','southeast');
    xlabel('Time (frames)');
    ylabel('Change (pixels)');
    title('Pulsation on the circumferential axis');
    hold off;
    % plot the selected section along tran-axis
    subplot(2,1,2);
    plot(a_t(:,str2num(choice{1})), 'r');
    hold on; grid on;
    plot(v_t(:,str2num(choice{1})), 'b');
    legend('Acceleraton','Velocity','Location','southeast');
    xlabel('Time (frames)');
    ylabel('Change (pixels)');
    title('Pulsation on the transmural axis');
    hold off;
    % move or not
    question_2 = {'Move to next step or not? y/n'};
    move = inputdlg(question_2);
    % choose the section
    question_3 = {'Which section to choose? default = last'};
    choice = inputdlg(question_3);   % choice for section #
end
% choose the axis to use
question_4 = {'Which axis is used for phase extraction? c/t'};
choice_2 = inputdlg(question_4);
if choice_2{1} == 'c'
    acc = a_c(:,str2num(choice{1}));
else
    acc = a_t(:,str2num(choice{1}));
end
disp(['Selected section = ' choice{1}]);
disp(['Selected axis = ' choice_2{1}]);

%% Determine the systole and diastole
% find the systole and diastole based on acceleration
period = ind/cycle; % calculate the period for each phase
% find five largest acc
[~, order_1] = sort(acc, 'descend');
% initialize sys
phase_1 = zeros(cycle, 1);
i = 1;  % i for order array
t = 1;  % t for sys array
zone = [];  % forbidden zone for next item in the order array
while phase_1(end) == 0
    if ismember(order_1(i), zone)    % check whether order(i) is in the zone
        i = i+1;
    else
        phase_1(t) = order_1(i);
        zone = cat(2, zone, round(phase_1(t)-0.5*period):round(phase_1(t)+0.5*period));  % connect different zones
        t = t+1;
        i = i+1;
    end
end
phase_1 = sort((phase_1 + 1), 'ascend'); % original # in the raw stack
% find five smallest acc
[~, order_2] = sort(acc, 'ascend');
% initialize sys and dia
phase_2 = zeros(cycle, 1);
i = 1;
t = 1;
zone = [];
while phase_2(end) == 0
    if ismember(order_2(i), zone)
        i = i+1;
    else
        phase_2(t) = order_2(i);
        zone = cat(2, zone, round(phase_2(t)-0.5*period):round(phase_2(t)+0.5*period));
        t = t+1;
        i = i+1;
    end
end
phase_2 = sort((phase_2 + 1), 'ascend'); % original # in the raw stack
% find systole and diastole
% assess the position of last 2 numbers
if phase_1(end) > phase_2(end)  % last image is systole
    dia = phase_2;
    sys = phase_1;
else
    dia = phase_1;
    sys = phase_2;
end

%% Calculate the displacement and strain along two PCA axes
disp('displacement along pca axes');
% circumferential displacement from diastole to systole
circ_dis = abs(permute(coord_pca(dia,1,:) - coord_pca(sys,1,:), [1,3,2])); % define displacement = Dia_cent - Sys_cent
circ_dis_rel = circ_dis ./ repmat(permute(range(1,1,:), [1,3,2]), cycle, 1);    % relative displacement = displacement / inner
circ_dis_rel(isnan(circ_dis_rel(:))|isinf(circ_dis_rel(:))) = 0;    % remove nan and inf
mean_circ_dis = abs(mean(circ_dis_rel));
std_circ_dis = std(circ_dis_rel);
% transmural displacement from diastole to systole
tran_dis = abs(permute(coord_pca(dia,2,:) - coord_pca(sys,2,:), [1,3,2]));
tran_dis_rel = tran_dis ./ repmat(permute(range(1,2,:), [1,3,2]), cycle, 1);
tran_dis_rel(isnan(tran_dis_rel(:))|isinf(tran_dis_rel(:))) = 0;
mean_tran_dis = abs(mean(tran_dis_rel));
std_tran_dis = std(tran_dis_rel);
% circumferential strain from diastole to systole
circ_strain = abs(permute(range(dia,1,:) - range(sys,1,:), [1,3,2]));
circ_strain = circ_strain ./ repmat(permute(range(1,1,:), [1,3,2]), cycle, 1); % defind strain = displacement / sys_range
circ_strain(isnan(circ_strain(:))|isinf(circ_strain(:))) = 0;
mean_circ_strain = abs(mean(circ_strain));   % mean of multiple measurements
std_circ_strain = std(circ_strain); % std of multiple measurements
% radial strain from diastole to systole
tran_strain = abs(permute(range(dia,2,:) - range(sys,2,:), [1,3,2]));
tran_strain = tran_strain ./ repmat(permute(range(1,2,:), [1,3,2]), cycle, 1);
tran_strain(isnan(tran_strain(:))|isinf(tran_strain(:))) = 0;
mean_tran_strain = abs(mean(tran_strain));
std_tran_strain = std(tran_strain);
disp('Strain calculation DONE');

%% Visualization of sections with peudo-color
% assign the color to each section
circ_D = zeros(col, row, GroupNo);  % displacement
tran_D = zeros(col, row, GroupNo);
circ_S = zeros(col, row, GroupNo);  % strain
tran_S = zeros(col, row, GroupNo);
for i = 1:GroupNo
    circ_D(:,:,i) = (section(:,:,1) == i) * mean_circ_dis(i);
    tran_D(:,:,i) = (section(:,:,1) == i) * mean_tran_dis(i);
    circ_S(:,:,i) = (section_inner(:,:,1) == i) * mean_circ_strain(i);
    tran_S(:,:,i) = (section_inner(:,:,1) == i) * mean_tran_strain(i);
end
% print circumferential image
circ_D_scale = sum(circ_D, 3);
circ_S_scale = sum(circ_S, 3);
figure; imshow(circ_D_scale, []);
title('Circumferential displacement');
colormap('hot');
colorbar;
% save the image
filename = [SavePath '\circumferential displacement.tif'];
print(filename, '-dtiff', '-r600');

figure; imshow(circ_S_scale, []);
title('Circumferential strain');
colormap('hot');
colorbar;
% save the image
filename = [SavePath '\circumferential strain.tif'];
print(filename, '-dtiff', '-r600');

% print transmural image
tran_D_scale = sum(tran_D, 3);
tran_S_scale = sum(tran_S, 3);
figure; imshow(tran_D_scale, []);
title('Transmural displacement');
colormap('hot');
colorbar;
% save the image
filename = [SavePath '\transmural displacement.tif'];
print(filename, '-dtiff', '-r600');

figure; imshow(tran_S_scale, []);
title('Transmural strain');
colormap('hot');
colorbar;
% save the image
filename = [SavePath '\transmural strain.tif'];
print(filename, '-dtiff', '-r600');
disp('Analysis DONE');

%% Visualize the displacement and strain changes at different time points
if time_series == 1
    % v_c: displacement along circumferential axis; v_t: displacement along transmural axi
    % assign the sign to different sections
    vc_ts = v_c;
    vt_ts = v_t;
    % strain
    circ_dis_ts = abs(vc_ts ./ permute(range(ind-1,1,:), [1,3,2])); % defind relative displacement = disp / length
    circ_dis_ts(isnan(circ_dis_ts(:))|isinf(circ_dis_ts(:))) = 0;
    
    tran_dis_ts = abs(vt_ts ./ permute(range(ind-1,2,:), [1,3,2]));
    tran_dis_ts(isnan(tran_dis_ts(:))|isinf(tran_dis_ts(:))) = 0;
    
    circ_strain_ts = abs(permute(diff(range(:,1,:), 1, 1) ./ range(ind-1, 1, :), [1,3,2]));    % strain
    if circ_strain_ts(1,:) == 0
        circ_strain_ts(1,:) = mean_circ_strain;
    elseif circ_strain_ts(end,:) == 0
        circ_strain_ts(end,:) = mean_circ_strain;
    else
        for i = 2:ind-2
            if ~any(circ_strain_ts(i,:), 'all')
                circ_strain_ts(i,:) = 0.5*(circ_strain_ts(i-1,:) + circ_strain_ts(i+1,:));
            end
        end
    end
    
    tran_strain_ts = abs(permute(diff(range(:,2,:), 1, 1) ./ range(ind-1, 2, :), [1,3,2]));
    if tran_strain_ts(1,:) == 0
        tran_strain_ts(1,:) = mean_tran_strain;
    elseif tran_strain_ts(end,:) == 0
        tran_strain_ts(end,:) = mean_tran_strain;
    else
        for i = 2:ind-2
            if ~any(tran_strain_ts(i,:), 'all')
                tran_strain_ts(i,:) = 0.5*(tran_strain_ts(i-1,:) + tran_strain_ts(i+1,:));
            end
        end
    end
    
    % maximal and minimal intensity for displacement
    circ_D_max = max(circ_dis_ts(:));
    circ_D_min = min(circ_dis_ts(:));
    circ_D_ts_scale = 255 * (circ_dis_ts / circ_D_max);
    tran_D_max = max(tran_dis_ts(:));
    tran_D_min = min(tran_dis_ts(:));
    tran_D_ts_scale = 255 * (tran_dis_ts / tran_D_max);
    % maximal and minimal intensity for strain
    circ_S_max = max(circ_strain_ts(:));
    circ_S_min = min(circ_strain_ts(:));
    circ_S_ts_scale = 255 * (circ_strain_ts / circ_S_max);
    tran_S_max = max(tran_strain_ts(:));
    tran_S_min = min(tran_strain_ts(:));
    tran_S_ts_scale = 255 * (tran_strain_ts / tran_S_max);
    
    % circ D
    circ_D_ts = zeros(col, row, 1, ind, GroupNo);
    for i = 1:GroupNo
        disp(['printing group # ' num2str(i)]);
        for j = 1:ind-1
            circ_D_ts(:,:,1,j,i) = (section(:,:,j) == i) * circ_D_ts_scale(j,i);
        end
    end
    circ_D_final = sum(circ_D_ts, 5) + 1;
    % visualize circ displacement
    mov_circD = immovie(ceil(circ_D_final), hot);
    h_circD = implay(mov_circD, 10);
    h_circD.Visual.ColorMap.UserRange = 1;
    h_circD.Visual.ColorMap.UserRangeMin = 1;
    h_circD.Visual.ColorMap.UserRangeMax = 255;
    h_circD.Parent.Name = 'Circumferential displacement - time dependent';
    filename = [SavePath '\Circumferential displacement.avi'];
    v1 = VideoWriter(filename);
    v1.FrameRate = 10;
    v1.Quality = 95;
    open(v1);
    writeVideo(v1, mov_circD);
    clear circ_D_ts;
    close(v1);
    
    % tran D
    tran_D_ts = zeros(col, row, 1, ind, GroupNo);
    for i = 1:GroupNo
        disp(['printing group # ' num2str(i)]);
        for j = 1:ind-1
            tran_D_ts(:,:,1,j,i) = (section(:,:,j) == i) * tran_D_ts_scale(j,i);
        end
    end
    tran_D_final = sum(tran_D_ts, 5) + 1;
    % visualize tran displacement
    mov_tranD = immovie(ceil(tran_D_final), hot);
    h_tranD = implay(mov_tranD, 10);
    h_tranD.Visual.ColorMap.UserRange = 1;
    h_tranD.Visual.ColorMap.UserRangeMin = 1;
    h_tranD.Visual.ColorMap.UserRangeMax = 255;
    h_tranD.Parent.Name = 'Transmural displacement - time dependent';
    filename = [SavePath '\Transmural displacement.avi'];
    v2 = VideoWriter(filename);
    v2.FrameRate = 10;
    v2.Quality = 95;
    open(v2);
    writeVideo(v2, mov_tranD);
    clear tran_D_ts;
    close(v2);
    
    % circ S
    circ_S_ts = zeros(col, row, 1, ind, GroupNo);
    for i = 1:GroupNo
        disp(['printing group # ' num2str(i)]);
        for j = 1:ind-1
            circ_S_ts(:,:,1,j,i) = (section_inner(:,:,j) == i) * circ_S_ts_scale(j,i);
        end
    end
    circ_S_final = sum(circ_S_ts, 5) + 1;
    % visualize circ strain
    mov_circS = immovie(ceil(circ_S_final), hot);
    h_circS = implay(mov_circS, 10);
    h_circS.Visual.ColorMap.UserRange = 1;
    h_circS.Visual.ColorMap.UserRangeMin = 1;
    h_circS.Visual.ColorMap.UserRangeMax = 255;
    h_circS.Parent.Name = 'Circumferential strain - time dependent';
    filename = [SavePath '\Circumferential strain.avi'];
    v3 = VideoWriter(filename);
    v3.FrameRate = 10;
    v3.Quality = 95;
    open(v3);
    writeVideo(v3, mov_circS);
    clear circ_S_ts;
    close(v3);
    
    % tran S
    tran_S_ts = zeros(col, row, 1, ind, GroupNo);
    for i = 1:GroupNo
        disp(['printing group # ' num2str(i)]);
        for j = 1:ind-1
            tran_S_ts(:,:,1,j,i) = (section_inner(:,:,j) == i) * tran_S_ts_scale(j,i);
        end
    end
    tran_S_final = sum(tran_S_ts, 5) + 1;
    % visualize tran strain
    mov_tranS = immovie(ceil(tran_S_final), hot);
    h_tranS = implay(mov_tranS, 10);
    h_tranS.Visual.ColorMap.UserRange = 1;
    h_tranS.Visual.ColorMap.UserRangeMin = 1;
    h_tranS.Visual.ColorMap.UserRangeMax = 255;
    h_tranS.Parent.Name = 'Transmural strain - time dependent';
    filename = [SavePath '\Transmural strain.avi'];
    v4 = VideoWriter(filename);
    v4.FrameRate = 10;
    v4.Quality = 95;
    open(v4);
    writeVideo(v4, mov_tranS);
    clear tran_S_ts;
    close(v4);
    
    % plot the time series circ-displacement
    disp('Time lapse processing of circ direction');
    figure('Position', [10 10 1200 800]);
    hold on;
    for i = 1:GroupNo
        if i == str2num(choice{1})
            subplot(GroupNo, 1, i);
            plot(circ_dis_ts(:,i), 'r', 'LineWidth', 3);
        else
            subplot(GroupNo, 1, i)
            plot(circ_dis_ts(:,i), 'k', 'LineWidth', 1.5);
        end
        if i == 1
            title('Relative displacement along the circumferential direction');
        end
        if i == GroupNo
            xlabel('Time (frames)');
        else
            set(gca,'xtick',[]);
        end
    end
    hold off;
    % save the image
    filename = [SavePath '\Time-dependent circumferential displacement.tif'];
    print(filename, '-dtiff', '-r600');
    % plot the time series circ-strain
    figure('Position', [10 10 1200 800]);
    hold on;
    for i = 1:GroupNo
        if i == str2num(choice{1})
            subplot(GroupNo, 1, i);
            plot(circ_strain_ts(:,i), 'r', 'LineWidth', 3);
        else
            subplot(GroupNo, 1, i)
            plot(circ_strain_ts(:,i), 'k', 'LineWidth', 1.5);
        end
        if i == 1
            title('Strain along the circumferential direction');
        end
        if i == GroupNo
            xlabel('Time (frames)');
        else
            set(gca,'xtick',[]);
        end
    end
    hold off;
    % save the image
    filename = [SavePath '\Time-dependent circumferential strain.tif'];
    print(filename, '-dtiff', '-r600');
    
    % plot the time series tran-displacement
    disp('Time lapse processing of tran direction');
    figure('Position', [10 10 1200 800]);
    hold on;
    for i = 1:GroupNo
        if i == str2num(choice{1})
            subplot(GroupNo, 1, i);
            plot(tran_dis_ts(:,i), 'r', 'LineWidth', 3);
        else
            subplot(GroupNo, 1, i)
            plot(tran_dis_ts(:,i), 'k', 'LineWidth', 1.5);
        end
        if i == 1
            title('Relative displacement along the transmural direction');
        end
        if i == GroupNo
            xlabel('Time (frames)');
        else
            set(gca,'xtick',[]);
        end
    end
    hold off;
    % save the image
    filename = [SavePath '\Time-dependent transmural displacement.tif'];
    print(filename, '-dtiff', '-r600');
    % plot the time series tran-strain
    figure('Position', [10 10 1200 800]);
    hold on;
    for i = 1:GroupNo
        if i == str2num(choice{1})
            subplot(GroupNo, 1, i);
            plot(tran_strain_ts(:,i), 'r', 'LineWidth', 3);
        else
            subplot(GroupNo, 1, i)
            plot(tran_strain_ts(:,i), 'k', 'LineWidth', 1.5);
        end
        if i == 1
            title('Strain along the transmural direction');
        end
        if i == GroupNo
            xlabel('Time (frames)');
        else
            set(gca,'xtick',[]);
        end
    end
    hold off;
    % save the image
    filename = [SavePath '\Time-dependent transmural strain.tif'];
    print(filename, '-dtiff', '-r600');
    disp('Time lapse analysis DONE');
end




