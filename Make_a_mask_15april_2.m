% Pixel dump to make mask
% MLS 
% time 16:51
clear all 
close all 

%% Load the PET 

Path_PET='C:\Users\Maria\Documents\MATLAB\NUK\SIRT_all\SIRT_2021\Data_fra_Database_SIRINO_14_04_2021\YM68_PET\20190710';
[pet_ref, pet_matrix, pet_info] = read_dicom_dir(Path_PET);

PETdim=size(pet_matrix);
makeMask=zeros(PETdim);
%% This is the path where the filedump is in, and where results are saved
Path='C:\Users\Maria\Documents\MATLAB\NUK\SIRT_all\SIRT_2021\ym_pixel_dump_pet_only_enh\';
%'C:\Users\blakk\Documents\BibtexLibs\SIRT_voxeldosimetry\data\';
cd(Path)

%Import Excel
global_list = dir('*.xlsx');
global_N_File=numel(global_list);
excelfile_name=global_list(1).name;
disp(['Working on: ' global_list(1).name])

T2 = readtable(excelfile_name);
SizeTable= size(T2);
rows=SizeTable(1);
T2Variables=T2.Properties.VariableNames;

% Find unique VOI names/types 
Column1=T2(:,1);
[C, ic]=unique(Column1);
number_of_VOis=length(C.Properties.DimensionNames);
VOInames=C.Variables;

%% Mask Creation Here 
% For each voi name
for kk=1:number_of_VOis
    currentVOI=VOInames(kk);
    disp(['Working for ' currentVOI])
    currentvoi_str=string(currentVOI);
%% Find indices
    CurrentVOI_ind =[];
    for jj=1:rows
        if strcmp(T2.VoiName_Region__string_(jj), currentVOI)==1
            CurrentVOI_ind=[CurrentVOI_ind, jj];
        end
    end
    clear jj

    %% Create Mask
    makeMask_current=zeros(PETdim);
    takeXYZ=[];
    for jj=1:length(CurrentVOI_ind)
        current_row=CurrentVOI_ind(jj);
        takeXYZ_temp=[T2.Y_pixel_(current_row) T2.X_pixel_(current_row) 71-T2.Z_pixel_(current_row)];
%         T2.Y_pixel_(current_row);
%         T2.X_pixel_(current_row);
%         T2.Z_pixel_(current_row);
        
        X_temp=takeXYZ_temp(1);
        Y_temp=takeXYZ_temp(2);
        Z_temp=takeXYZ_temp(3);
        
        %% Ask the coordinates so that they belong to PET
        pet_value=pet_matrix(Y_temp,X_temp, Z_temp);
        %T2.Y_pixel_(current_row), T2.X_pixel_(current_row), 71-T2.Z_pixel_(current_row))
        T2.Value_kBq_cc_(current_row);
        %%
        % Give a value to the mask 
        if kk==1; 
            mask_value=1;
            elseif kk==2, mask_value=2;
        end
        
        makeMask_current(X_temp, Y_temp, Z_temp)=mask_value;
        makeMask(X_temp, Y_temp, Z_temp)=mask_value;
        
        takeXYZ=[takeXYZ; takeXYZ_temp];
        clear current_row X_temp Y_temp Z_temp
    end
    clear jj 
    
    %% Visualize mask
    figure('Name',currentvoi_str);
    implay(mat2gray(makeMask_current, [0, 2]))
    
    %% Save some results in the original Path with the pixel dump 
    cd(Path)
    save makeMask makeMask
    
%% Make a structure to save 
    struct_results.(currentvoi_str).mask=makeMask_current;
    struct_results.(currentvoi_str).takeXYZ=takeXYZ;
    struct_results.(currentvoi_str).CurrentVOI_indices=CurrentVOI_ind;
    struct_results.T2=T2;
    struct_results.Mask=makeMask;
    save struct_results struct_results
    % exit first loop
    clear currentVOI makeMask_current currentVOI_ind
    
end

save VOInames VOInames
%% Play unified Mask 
figure(100); implay(mat2gray(makeMask))
