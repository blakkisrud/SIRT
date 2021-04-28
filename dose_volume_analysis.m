% Script to perform analysis on SIRT-images

% Define constants - IMPORTANT, verify these - Johan april 2021

TISSUE_DENSITY = 1.03; % Density in g/ml?
time_between_admin_and_scan = 24; % Time in hrs
Y_90_half_life = 64.2; % Half life in hrs
Y_90_decay_const = log(2)/Y_90_half_life; % Decay-constant 90-Y
Y_90_energy_constant = 1.4958e-13 ; % in J per desintegration, calculate
                                    % assuming 0.9336 of beta-energy. This
                                    % have to be replaced by a verified
                                    % value at some point


% Load the PET-data

path_to_pet_data = ""

[pet_ref, pet_matrix, pet_info] = read_dicom_dir(path_to_pet_data);

% Load the mask and pet-matrix

load('data\struct_results.mat')
liver_struct = struct_results.ttumor_auto;
liver_mask = liver_struct.mask;

% Have to have correct units - start with mask being in kBq/ml?

% Voxel volume in ml:
% This is strange - the dZ-value from the read_dicom-function says 0.4, but
% the data indicates that the slice-thickness is 2.79, is it not uniform?
% For now hard-code the slice-thickness, but it should be replaced by the 
% thickness retreived from somewhere else
% Ask LTGM about the slice-thickness

slice_thickness = pet_info.SliceThickness; % Units of mm

voxel_volume = (pet_ref.PixelExtentInWorldX/10)*(pet_ref.PixelExtentInWorldY/10)*(slice_thickness/10);%(pet_ref.PixelExtentInWorldZ/10);
voxel_mass = (voxel_volume/1e3)*TISSUE_DENSITY;

pet_matrix_Bq = pet_matrix*voxel_volume; % Now have matrix in Bq
% Convert the matrix in Bq at imaging to activity at administration
pet_matrix_Bq_time_zero = pet_matrix_Bq/(exp(Y_90_decay_const*-1*time_between_admin_and_scan));
pet_matrix_Bq_s = (pet_matrix_Bq_time_zero/Y_90_decay_const)*3600; % First Bq*hrs, then Bq*s
pet_matrix_energy = pet_matrix_Bq_s*Y_90_energy_constant; % Energy in joule released in each voxel
pet_matrix_absorbed_dose = pet_matrix_energy/voxel_mass; % Absorbed dose in Gy

% Check the sums to see if we get OK results

disp('Total activity at AQ - MBq')
disp(sum(pet_matrix_Bq(:))/1e6)

disp('Total activity at Admin - MBq')
disp(sum(pet_matrix_Bq_time_zero(:))/1e6)

disp('Maximum absorbed dose in voxel (Gy)')
disp(max(pet_matrix_absorbed_dose(:)))

% Check the absorbed dose in the tumour

mean(pet_matrix_absorbed_dose(liver_mask==2))

% Lets have a look at the mask and pet-image

pet_img = pet_matrix(:,:,50);
mask_img = liver_mask(:,:,50);

imshowpair(pet_img, mask_img)

% Now if we assume everything works - we should be able to plot historgrams
% of the distributions
%%
target_voxels = pet_matrix_absorbed_dose(liver_mask==2);
max_dose = max(target_voxels);
D_vec = 1:max_dose;
fractions = zeros(length(D_vec),1);
num_target_voxels = length(target_voxels);

for i = 1:length(D_vec)
    less_than_i = target_voxels>D_vec(i);
    fractions(i) = sum(less_than_i)./num_target_voxels;
    
end

plot(D_vec, fractions)

