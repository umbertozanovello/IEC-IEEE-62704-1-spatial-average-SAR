% prepare IEC-IEEE_62704-1 SAR star data for testing spatialAverageSAR.py averaging algorithm and display some figures
% N De Zanche 2024-11-22

% do if running in Octave
if exist('OCTAVE_VERSION', 'builtin')
  pkg load io
  pkg load image
endif

load_data = 1;
load_type = 0; % 0 = text file; 1 = HDF5 file; 2 = MAT file
remove_duplicates = 1;    %remove duplicate average SAR entries for same position (seen at boundaries)
save_data = 0;      % save the original table in .mat format to speed up reading data
save_nrrd = 0;      % for visualization; NRRD is easier to import into 3D Slicer
uniform_grid=1;   % required to display correctly; set manually according to the file type chosen below

datadir = 'D:\Library\Reference\Standards\IEC\SAR_Star_Reference_Results_160318\';
% NB: spatialAverageSAR.py works only with data on a uniform grid and does not accept non-uniform (graded) data
%filename = 'sar_star_graded_ref_10g_2016_03_18.txt'
%filename = 'sar_star_graded_ref_01g_2016_03_18.txt';
%filename = 'sar_star_uniform_ref_01g_2016_03_18_V3.txt';
filename = 'sar_star_uniform_ref_10g_2016_03_18_V3.txt';

% data from IEC-IEEE_62704:
% position 1 = background (air)
% position 2 = core
% position 3 = outer
densities = [0 2000 1100]; % mass density [kg/m^3]
permittivities = [1 12 41];  % relative permittivity
conductivities = [0 0.15 0.8]; % conductivity in S/m

if load_data
% load SAR star data
  disp('loading data...')

  if load_type == 0

    header_lines = 25;

    file = fopen([datadir filename], "r");
    header = cell(header_lines,1);
    for n= 1:header_lines
      header(n) = fgetl(file);
    end
    fclose(file);

    data = dlmread([datadir filename],'\t');
    data = data((header_lines+1):end,:);

    % grid subdivision vectors
    tempx = textscan(char(header{4}),"%s"){1};
    x = str2num(char(tempx(2:end)))';
    tempy = textscan(char(header{5}),"%s"){1};
    y = str2num(char(tempy(2:end)))';
    tempz = textscan(char(header{6}),"%s"){1};
    z = str2num(char(tempz(2:end)))';

    if save_data
      disp('saving data...')
      if exist([datadir filename '.mat'],'file') == 2
        error([datadir filename '.mat already exists!']); % prevent overwriting
      else
        save('-v7', [datadir filename '.mat'], 'x', 'y', 'z', 'data', 'header')
      endif
    endif

  elseif load_type == 1
    load([datadir filename '.hdf5']);
  elseif load_type == 2
    load([datadir filename '.mat']);
  else
    error("unspecified file type");
  endif
endif

% offsets to line up star correctly with axes
x_offset = x + 0.0005;
y_offset = y + 0.0005;
z_offset = z + 0.0005;
[X Y Z] = ndgrid(x_offset, y_offset, z_offset);

array_size = [size(x,2) size(y,2) size(z,2)];
indices = data(:,1:3);

%% find lines that have the same indices (duplicated)
% see https://www.mathworks.com/matlabcentral/answers/336500-finding-the-indices-of-duplicate-values-in-one-array
[temp,idx,n] = unique(indices,"rows");  % "stable" not implemented but not necessary because fortunately indices are ordered
counts = accumarray(n, 1);%, [], @sum);
duplicates = find(counts~=1); % indices into the first duplicated values
% index into first occurence of duplicate k is idx(duplicates(k))

if remove_duplicates
  for k = 1:size(duplicates,1)
    no_duplicates = counts(duplicates(k));    %number of duplicates of current position
    tmp = 1:no_duplicates;
    % find duplicate with max SAR
    [m im] = max(data(idx(duplicates(k)):idx(duplicates(k))+counts(duplicates(k))-1,9));
    for k1 = 1:no_duplicates
      if k1 ~= im
        indices(idx(duplicates(k))+k1-1,:) = [NaN NaN NaN];  % tag unwanted duplicates
      endif
    end
  end
  data(isnan(indices(:,1)),:) = [];    % remove duplicate rows from data and indices arrays
  indices(isnan(indices(:,1)),:) = [];
endif

%% create array to classify voxels belonging to the background (0) star core (1) and star shell (2)

% dimensions according to Fig. 10 of IEC/IEEE 62704-1:2017
embd_cube_half = 0.007; %[m]
cube_hollow_half = 0.012; %[m]
cube_side_half = 0.04; %[m]
inner_peg_length = 0.045; %[m]
inner_peg_radius = 0.01; %[m]

star = zeros(array_size);

% start by setting all star voxels to shell (2)
star(sub2ind(array_size,indices(:,1)'+1,indices(:,2)'+1,indices(:,3)'+1)) = 2;
% +1 needed to fix the offset in addition to grid offsets above

% set central cubes (x,y,z <= 7 and 12 <= x,y,z <= 40 mm) to core material (1)
star(find(abs(X) <= embd_cube_half & abs(Y) <= embd_cube_half & abs(Z) <= embd_cube_half)) = 1;
star(find((abs(X) <= cube_side_half & abs(Y) <= cube_side_half & abs(Z) <= cube_side_half) & !(abs(X) <= cube_hollow_half & abs(Y) <= cube_hollow_half & abs(Z) <= cube_hollow_half))) = 1;

% set 6 inner pegs to core material (1)
star(find(abs(Z) >= cube_side_half & (abs(Z) <= cube_side_half + inner_peg_length) & (X.^2 + Y.^2 <= inner_peg_radius^2))) = 1;
star(find(abs(Y) >= cube_side_half & (abs(Y) <= cube_side_half + inner_peg_length) & (X.^2 + Z.^2 <= inner_peg_radius^2))) = 1;
star(find(abs(X) >= cube_side_half & (abs(X) <= cube_side_half + inner_peg_length) & (Z.^2 + Y.^2 <= inner_peg_radius^2))) = 1;

%% extract local SAR from imported data
local_SAR = zeros(array_size);
%local_SAR(sub2ind(array_size,indices(:,1)',indices(:,2)',indices(:,3)')) = data.data(:,8);
local_SAR(sub2ind(array_size,indices(:,1)'+1,indices(:,2)'+1,indices(:,3)'+1)) = data(:,8);
% +1 needed to fix the offset as above

% use this as a check to verify that support of local_SAR and star are identical
% comparison = xor(local_SAR > 0, star > 0);

%% extract averaged SAR from imported data
average_SAR = zeros(array_size);
average_SAR(sub2ind(array_size,indices(:,1)'+1,indices(:,2)'+1,indices(:,3)'+1)) = data(:,9);
% +1 needed to fix the offset as above

%comparison2 = xor(average_SAR > 0, star > 0);

%% extract voxel status flag from imported data
status = zeros(array_size);
status(sub2ind(array_size,indices(:,1)'+1,indices(:,2)'+1,indices(:,3)'+1)) = data(:,4);
% +1 needed to fix the offset as above

%% extract volume, mass and orientation columns from imported data
mass = zeros(array_size);
volume = zeros(array_size);
orientation = zeros(array_size);
mass(sub2ind(array_size,indices(:,1)'+1,indices(:,2)'+1,indices(:,3)'+1)) = data(:,5);
volume(sub2ind(array_size,indices(:,1)'+1,indices(:,2)'+1,indices(:,3)'+1)) = data(:,6);
orientation(sub2ind(array_size,indices(:,1)'+1,indices(:,2)'+1,indices(:,3)'+1)) = data(:,7);


%% create figures to display results

if uniform_grid

  % create .nrrd files to view in 3D slicer
  if save_nrrd
    disp('saving NRRD file...');
    img.pixelData = star;
    img.ijkToLpsTransform = [ 1 0 0 x(1); 0 1 0 y(1); 0 0 1 z(1); 0 0 0 1];
    if exist([datadir filename 'star.nrrd'],'file') == 2
      error([datadir filename 'star.nrrd already exists!'])
    else
      nrrdwrite([datadir filename 'star.nrrd'], img);
    endif
    img.pixelData = local_SAR;
    if exist([datadir filename 'SAR.nrrd'],'file') == 2
      error([datadir filename 'SAR.nrrd already exists!'])
    else
      nrrdwrite([datadir filename 'SAR.nrrd'], img);
    endif
  endif

  disp('displaying on uniform grid...')

  figure(1)
  set(gcf,'Name',['SAR MIPs: ' filename])
  subplot(1,3,1)
  imagesc([y(1) y(end)], [z(1) z(end)], squeeze(max(local_SAR,[],1)));
  %uimagesc(y, z, squeeze(max(local_SAR,[],1))');
  axis equal
  title('YZ plane');
  subplot(1,3,2)
  imagesc([x(1) x(end)], [z(1) z(end)], squeeze(max(local_SAR,[],2)));
  %uimagesc(x, z, squeeze(max(local_SAR,[],2))');
  axis equal
  title('XZ plane');
  subplot(1,3,3)
  imagesc([x(1) x(end)], [y(1) y(end)], squeeze(max(local_SAR,[],3)));
  %uimagesc([x(1) x(end)], [y(1) y(end)], squeeze(max(local_SAR,[],3))');
  axis equal
  title('XY plane');

  figure(2)
  set(gcf,'Name',['SAR projections: ' filename])
  subplot(1,3,1)
  imagesc([y(1) y(end)], [z(1) z(end)], squeeze(sum(local_SAR,1)));
  %uimagesc(y, z, squeeze(sum(local_SAR,1))');
  axis equal
  title('YZ plane');
  subplot(1,3,2)
  imagesc([x(1) x(end)], [z(1) z(end)], squeeze(sum(local_SAR,2)));
  %uimagesc(x, z, squeeze(sum(local_SAR,2))');
  axis equal
  title('XZ plane');
  subplot(1,3,3)
  imagesc([x(1) x(end)], [y(1) y(end)], squeeze(sum(local_SAR,3)));
  %uimagesc(x, y, squeeze(sum(local_SAR,3))');
  axis equal
  title('XY plane');

  figure(3)
  set(gcf,'Name',['SAR slices: ' filename])
%  colormap gray
  subplot(1,3,1)
  imagesc([y(1) y(end)], [z(1) z(end)], squeeze(local_SAR(ceil(array_size(1)/2),:,:)));
  axis equal
  title('YZ plane');
  subplot(1,3,2)
  imagesc([x(1) x(end)], [z(1) z(end)], squeeze(local_SAR(:,ceil(array_size(2)/2),:)));
  axis equal
  title('XZ plane');
  subplot(1,3,3)
  imagesc([x(1) x(end)], [y(1) y(end)], squeeze(local_SAR(:,:,ceil(array_size(3)/2))));
  axis equal
  title('XY plane');

  figure(4)
  set(gcf,'Name',['star slices: ' filename])
  colormap gray
  subplot(2,3,1)
  imagesc([y_offset(1) y_offset(end)], [z_offset(1) z_offset(end)], squeeze(star(ceil(array_size(1)/2),:,:)));
  axis equal
  title('YZ plane');
  ylabel('star material (background, core, outer)');
  subplot(2,3,2)
  imagesc([x_offset(1) x_offset(end)], [z_offset(1) z_offset(end)], squeeze(star(:,ceil(array_size(2)/2),:)));
  axis equal
  title('XZ plane');
  subplot(2,3,3)
  imagesc([x_offset(1) x_offset(end)], [y_offset(1) y_offset(end)], squeeze(star(:,:,ceil(array_size(3)/2))));
  axis equal
  title('XY plane');
  subplot(2,3,4)
  imagesc([y_offset(1) y_offset(end)], [z_offset(1) z_offset(end)], squeeze(status(ceil(array_size(1)/2),:,:)));
  axis equal
  ylabel('status flag (INVALID, UNUSED, USED, VALID)');
  subplot(2,3,5)
  imagesc([x_offset(1) x_offset(end)], [z_offset(1) z_offset(end)], squeeze(status(:,ceil(array_size(2)/2),:)));
  axis equal
  subplot(2,3,6)
  imagesc([x_offset(1) x_offset(end)], [y_offset(1) y_offset(end)], squeeze(status(:,:,ceil(array_size(3)/2))));
  axis equal

  figure(5)
  set(gcf,'Name',['star and SAR overlay: ' filename])
  colormap gray
  subplot(1,3,1)
  imshowpair(squeeze(star(ceil(array_size(1)/2),:,:)),squeeze(local_SAR(ceil(array_size(1)/2),:,:)),'blend');
  axis equal
  title('YZ plane');
  subplot(1,3,2)
  imshowpair(squeeze(star(:,ceil(array_size(1)/2),:)),squeeze(local_SAR(:,ceil(array_size(1)/2),:)),'blend');
  axis equal
  title('XZ plane');
  subplot(1,3,3)
  imshowpair(squeeze(star(:,:,ceil(array_size(1)/2))),squeeze(local_SAR(:,:,ceil(array_size(1)/2))),'blend');
  axis equal
  title('XY plane');

  figure(6)
  set(gcf,'Name',['SAR comparison: ' filename])
  colormap gray
  subplot(2,3,1)
  imagesc([y_offset(1) y_offset(end)], [z_offset(1) z_offset(end)], squeeze(local_SAR(ceil(array_size(1)/2),:,:)));
  axis equal
  title('YZ plane');
  ylabel('local SAR');
  subplot(2,3,2)
  imagesc([x_offset(1) x_offset(end)], [z_offset(1) z_offset(end)], squeeze(local_SAR(:,ceil(array_size(2)/2),:)));
  axis equal
  title('XZ plane');
  subplot(2,3,3)
  imagesc([x_offset(1) x_offset(end)], [y_offset(1) y_offset(end)], squeeze(local_SAR(:,:,ceil(array_size(3)/2))));
  axis equal
  title('XY plane');
  subplot(2,3,4)
  imagesc([y_offset(1) y_offset(end)], [z_offset(1) z_offset(end)], squeeze(average_SAR(ceil(array_size(1)/2),:,:)));
  axis equal
  ylabel('averaged SAR');
  subplot(2,3,5)
  imagesc([x_offset(1) x_offset(end)], [z_offset(1) z_offset(end)], squeeze(average_SAR(:,ceil(array_size(2)/2),:)));
  axis equal
  subplot(2,3,6)
  imagesc([x_offset(1) x_offset(end)], [y_offset(1) y_offset(end)], squeeze(average_SAR(:,:,ceil(array_size(3)/2))));
  axis equal


  elseif

  disp('displaying on non-uniform grid...')
  %NB: the middle of the star is not exactly in the middle of the array
  % find indices of slices through the axes
  x_center = find(x == 0);
  y_center = find(y == 0);
  z_center = find(z == 0);

% pcolor without "shading flat" makes a lot of black lines
  figure(1)
  set(gcf,'Name',['SAR MIPs: ' filename])
  subplot(1,3,1)
  pcolor(y, z, squeeze(max(local_SAR,[],1))');
  shading flat
  axis equal
  title('YZ plane');
  subplot(1,3,2)
  pcolor(x, z, squeeze(max(local_SAR,[],2))');
  shading flat
  axis equal
  title('XZ plane');
  subplot(1,3,3)
  pcolor(x, y, squeeze(max(local_SAR,[],3))');
  shading flat
  axis equal
  title('XY plane');

  figure(2)
  set(gcf,'Name',['SAR projections: ' filename])
  subplot(1,3,1)
  pcolor(y, z, squeeze(sum(local_SAR,1))');
  shading flat
  axis equal
  title('YZ plane');
  subplot(1,3,2)
  pcolor(x, z, squeeze(sum(local_SAR,2))');
  shading flat
  axis equal
  title('XZ plane');
  subplot(1,3,3)
  pcolor(x, y, squeeze(sum(local_SAR,3))');
  shading flat
  axis equal
  title('XY plane');

  figure(3)
  set(gcf,'Name',['SAR slices: ' filename])
  subplot(1,3,1)
  pcolor(y, z, squeeze(local_SAR(x_center,:,:))');
  shading flat
  axis equal
  title('YZ plane');
  subplot(1,3,2)
  pcolor(x, z, squeeze(local_SAR(:,y_center,:))');
  shading flat
  axis equal
  title('XZ plane');
  subplot(1,3,3)
  pcolor(x, y, squeeze(local_SAR(:,:,z_center))');
  shading flat
  axis equal
  title('XY plane');

  figure(4)
  set(gcf,'Name',['star slices: ' filename])
  colormap gray
  subplot(2,3,1)
  pcolor(y, z, squeeze(star(x_center,:,:))');
  shading flat
  axis equal
  title('YZ plane');
  ylabel('star material (background, core, outer)');
  subplot(2,3,2)
  pcolor(x, z, squeeze(star(:,y_center,:))');
  shading flat
  axis equal
  title('XZ plane');
  subplot(2,3,3)
  pcolor(x, y, squeeze(star(:,:,z_center))');
  shading flat
  axis equal
  title('XY plane');
  subplot(2,3,4)
  pcolor(y, z, squeeze(status(x_center,:,:))');
  shading flat
  axis equal
  ylabel('status flag (INVALID, UNUSED, USED, VALID)');
  subplot(2,3,5)
  pcolor(x, z, squeeze(status(:,y_center,:))');
  shading flat
  axis equal
  subplot(2,3,6)
  pcolor(x, y, squeeze(status(:,:,z_center))');
  shading flat
  axis equal

  figure(6)
  set(gcf,'Name',['SAR comparison: ' filename])
  colormap gray
  subplot(2,3,1)
  pcolor(y, z, squeeze(local_SAR(x_center,:,:))');
  shading flat
  axis equal
  title('YZ plane');
  ylabel('local SAR');
  subplot(2,3,2)
  pcolor(x, z, squeeze(local_SAR(:,y_center,:))');
  shading flat
  axis equal
  title('XZ plane');
  subplot(2,3,3)
  pcolor(x, y, squeeze(local_SAR(:,:,z_center))');
  shading flat
  axis equal
  title('XY plane');
  subplot(2,3,4)
  pcolor(y, z, squeeze(average_SAR(x_center,:,:))');
  shading flat
  axis equal
  ylabel('averaged SAR');
  subplot(2,3,5)
  pcolor(x, z, squeeze(average_SAR(:,y_center,:))');
  shading flat
  axis equal
  subplot(2,3,6)
  pcolor(x, y, squeeze(average_SAR(:,:,z_center))');
  shading flat
  axis equal

endif


%% saving all offset (centred) data arrays as input for spatialAverageSAR.py
if exist([datadir filename '_arrays.mat'],'file') == 2
  error([datadir filename '_arrays.mat already exists!']); % prevent overwriting
else
  save('-v7', [datadir filename '_arrays.mat'], 'x_offset', 'y_offset', 'z_offset', 'local_SAR', 'average_SAR', 'star', 'status', 'densities', 'mass', 'volume', 'orientation');
endif

