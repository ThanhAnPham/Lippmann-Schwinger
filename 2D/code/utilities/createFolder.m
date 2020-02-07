function fpath = createFolder(folder,algo,varargin)

iter = 0;
run_var = strcat('run_',num2str(iter));
while exist(fullfile(folder,algo,run_var),'dir')...
        && ~isempty(dir(fullfile(folder,algo,run_var,'*.mat')))
    iter = iter + 1;
    run_var = strcat('run_',num2str(iter));
end
fpath = fullfile(folder,algo,run_var);
mkdir(fpath);

fprintf('New folder : %s\n',fpath)
end