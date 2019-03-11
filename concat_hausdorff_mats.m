function [outmat] = concat_hausdorff_mats(exptdir,suffixes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(suffixes)
    finput = char(strcat('hausdorff',num2str(suffixes(i)),'.mat'));
    %fhausdorff = 'hausdorff.csv';
    load(fullfile(exptdir,finput));
    sz = size(hausdorff_mat);
    if sz(1)~=sz(2)
        hausdorff_mat = squareform(hausdorff_mat);
    end
    disp(size(hausdorff_mat));
    if i==1
        outmat = zeros(length(suffixes),size(hausdorff_mat,1),size(hausdorff_mat,2));
    end
   outmat(i,:,:) = hausdorff_mat;
end

