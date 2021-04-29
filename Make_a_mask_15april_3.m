% Pixel dump to make mask
% MLS 
% time 17:19
% 29.04.2021, fixed : sorted the images out of PMOD 

clear all 
close all 

%% Load the PET 
% YM 
Path_PET='C:\Users\Maria\Documents\MATLAB\NUK\SIRT_all\SIRT_2021\Data_fra_Database_SIRINO_14_04_2021\YM68_PET\20190710';
% NK
%Path_PET='C:\Users\Maria\Documents\MATLAB\NUK\SIRT_all\SIRT_2021\Data_fra_Database_SIRINO_14_04_2021\NK72_PET\20190710'
[pet_ref, pet_matrix, pet_info] = read_dicom_dir_290421maria(Path_PET);

figure; imshow(pet_matrix(:,:,50), []);
figure; implay(mat2gray(pet_matrix));

% Get out the imageindex
for jj=1:71, 
    sequence(jj,:)=pet_info(jj).ImageIndex(:);
end

% reorganize a new PET matrix
%TODO

%% Run now for PET images
% path= 'C:\Users\Maria\Documents\MATLAB\NUK\EARL_2019_images\PET AC';
% [PETimref, PET_IM, PET_Info] = bildlesing_earl_func(path)

PETdim=size(pet_matrix);
makeMask=zeros(PETdim);
%% This is the path where the filedump is in, and where results are saved
% YM 
Path='C:\Users\Maria\Documents\MATLAB\NUK\SIRT_all\SIRT_2021\ym_pixel_dump_pet_only_enh\';
%'C:\Users\blakk\Documents\BibtexLibs\SIRT_voxeldosimetry\data\';
% NK 
%Path='C:\Users\Maria\Documents\MATLAB\NUK\SIRT_all\SIRT_2021\nk_pet_only_140421\'
cd(Path)

%Du må gå og åpner xl file selv manuelt og lager en xlsx file. 



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
        % From PET matrix 
        pet_value=pet_matrix(Y_temp,X_temp, Z_temp);
        %T2.Y_pixel_(current_row), T2.X_pixel_(current_row), 71-T2.Z_pixel_(current_row))
        pet_sjekk_value=str2mat(T2.Value_kBq_cc_(current_row));
       % disp(['At line: ' num2str(jj) 'PET matrix:' num2str(pet_value) ' t2 value: ' num2str(pet_sjekk_value)])
        wrong_jjs=[];
        if strcmp(pet_value, pet_sjekk_value)==0
            wrong_jjs=[wrong_jjs, jj];
        end
        % temporary end to check what is happening
    %end
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
