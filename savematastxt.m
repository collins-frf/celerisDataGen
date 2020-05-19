water_level_change = -1;
jj=0;
homeDir = "D:\celerisWorkFlow";
stdOutput='Z:\celerisWorkFlow\output\';
frames_dir=[stdOutput];

for ii=1:99
    cd(homeDir)
    time_mat_name = ['WL_', num2str(water_level_change), '_time_', num2str(jj), '.mat'];
    mat_path = [stdOutput, '\', time_mat_name];
    load(mat_path)
    filepath = [frames_dir, num2str(time_mat_name), '.txt'];
    filepath
    fid = fopen(filepath, 'wt');
    fprintf(fid, '%s\n', timedata{:});
    fclose(fid);
    jj=jj+1

close all
end