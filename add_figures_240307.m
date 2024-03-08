clear all
folder = 'D:\HPC-SWR project\Processed Data_231113\ripples_mat\ProfilingSheet\R34 (bar only)_SUB\S_';
folder2 = 'D:\HPC-SWR project\Processed Data_231113\ripples_mat\ProfilingSheet\R34 (bar only)_SUB\sum';
files = dir(fullfile(folder, '*.png'));  % 이미지 파일 확장자에 따라 변경
num_files = numel(files);
num_per_page = 3;

for k = 1:num_per_page:num_files
    figure('position',[100 100 1500 2000]);
    for j = 0:num_per_page-1
        if k+j <= num_files
            subplot(num_per_page, 1, j+1);
            img = imread(fullfile(folder, files(k+j).name));
            imshow(img);
            set(gca, 'Position', [0 1-(j+1)/num_per_page 1 1/num_per_page]);  % 간격 없애기
        end
    end
    saveas(gcf, fullfile(folder2, sprintf('Page_%d.jpg', ceil(k/num_per_page))));
    close(gcf);
end