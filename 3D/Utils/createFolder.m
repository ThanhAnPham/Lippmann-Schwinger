function fpath = createFolder(folder,algo,varargin)

if nargin > 2
    param_fname = varargin{1};
else
    param_fname = 'main_param.m';
end
iter = 0;
run_var = strcat('run_',num2str(iter));
while exist(fullfile(folder,algo,run_var),'dir')...
        && ~isempty(dir(fullfile(folder,algo,run_var,'*.mat')))
    iter = iter + 1;
    run_var = strcat('run_',num2str(iter));
end
fpath = fullfile(folder,algo,run_var);
mkdir(fpath);
if exist(param_fname,'file')
    copyfile(which(param_fname),fpath);
else
    warning('%s not found',param_fname);
end

fprintf('New folder : %s\n',fpath)
end