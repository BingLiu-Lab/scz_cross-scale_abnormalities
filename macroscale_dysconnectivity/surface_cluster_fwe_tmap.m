% load fs5 surf
surf_left_path = '/home/mwang/fsaverage5/surf/lh.pial';
surf_right_path = '/home/mwang/fsaverage5/surf/rh.pial';
surf_mesh = SurfStatReadSurf({surf_left_path, surf_right_path});

% load data Y
y_path = '/home/mwang/sz_nc_array_fwe.mat';
y_data = load(y_path);
y_data = y_data.arr;

% load design matrix
a_path = '/home/mwang/a_design_matrix_fwe.mat';
a_data = load(a_path);
a_data = a_data.arr;

% fitting linear model
Y = y_data';
M = term(a_data);
slm = SurfStatLinMod(Y, M, surf_mesh);

% FWE mask
Y_size = size(Y);
mask = logical(ones(1, Y_size(2)));

% thresh
clusthresh = 0.001;

% t test p
contrast_p = M.a_data8 - M.a_data7;
slm_p = SurfStatT(slm, contrast_p);

[pval_p, peak_p, clus_p, clusid_p] = SurfStatP(slm_p, mask, clusthresh);
save('/home/mwang/fwe_clusid_p', 'clusid_p');

t_map_p = slm_p.t;
save('/home/mwang/fwe_t_map_p', 't_map_p');

% % t test n
% contrast_n = M.a_data7 - M.a_data8;
% slm_n = SurfStatT(slm, contrast_n);
% 
% [pval_n, peak_n, clus_n, clusid_n] = SurfStatP(slm_n, mask, clusthresh);
% save('/home/mwang/fwe_clusid_n', 'clusid_n');
% 
% t_map_n = slm_n.t;
% save('/home/mwang/fwe_t_map_n', 't_map_n');