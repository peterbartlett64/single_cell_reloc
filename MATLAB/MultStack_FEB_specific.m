%%
% function res = MultStack(posFingerprint, ser_length, frame, path_seg, path_raw)

% num = int8(frame);
% posFingerprint = char(posFingerprint);
% path_seg = char(path_seg);
% path_raw  = char(path_raw);

% This is an automated file for MIJI integration with pipeline.
% This assumes that MIJI is already setup correctly. To do so, refer to "MIJ: Running ImageJ and Fiji within Matlab"
% http://bigwww.epfl.ch/sage/soft/mij/ or <https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab>
%%
cd ('D:\Manipulated\MATLAB_output\');
pos = "d0222r2p400200";
dir_n = sprintf('%s', pos);
mkdir(dir_n)
cd (dir_n);
path_seg = sprintf("E:/Microfluidics/RESULTS/%s/GFP_mKO_mKa", pos);
save_path_format= "%s_%s_%s.tiff";
save_path_temp = "temp_%s_%s_%s_%s.tiff";
ser_length = 71;
for num = 1:ser_length
    for p=["frame_pxMAT", "track_mask"]
        pad = num2str(num,'%05.f');
        file_name = sprintf('%s_%s.mat', p, pad);

        fullMatFileName = fullfile(path_seg,  file_name);

        if ~exist(fullMatFileName, 'file')
            message = sprintf('%s does not exist', fullMatFileName);
            % uiwait(warndlg(message));
        else
            if p == "frame_pxMAT"
                file = load(fullMatFileName).px_data;
            elseif p == "track_mask"
                file = load(fullMatFileName).data;
                %save_name_temp = sprintf(save_path_temp, 'orig' , pos, p, pad);
                %imwrite(file, save_name_temp);

                %file = permute(file,[2 1]);
                %save_name_temp = sprintf(save_path_temp,"permute", pos, p, pad);
                %imwrite(file, save_name_temp);

                file = uint16(file); %Convert to unint16. This will leave enough space for all the cells which could be present in a single frame

            end
            % file_copy = permute(file,[2 1]); # This should not be needed as it is converted to the ImagePlus object which is permuted
            % automatically
            save_name = sprintf(save_path_format,pos, p, pad);
            % imp = copytoImagePlus(file);
            imwrite(file, save_name)
        end
    end
end
%%
%username=getenv('USERNAME');
temp_loc = 'D:\Temp_images'; % This is the location for file outputs



temp = char(sprintf(user_path_format, version_num));
%%
addpath ('D:\Additional_Program Files\Fiji.app\scripts')
%addpath = char(Fiji_path)
ImageJ % Launch the imageJ kernel
%% This is the testing zone
% MATLAB uses differnt column and row indexing order

px_data_copy = permute(px_data,[2 1]); % Data from the expanded mask
data_copy = permute(data,[2 1]); % Data from the regular tracking mask

IJM.show('px_data_copy') % Good it is the permuted version which matches with the image target
IJM.show('data_copy');

%% 



IJM.createImage(px_data_copy)

% This is the version for MIJ, whereas the above uncommented is for IJM
user_path_format = 'C:Program Files/MATLAB/%s/java/%s';
version_num = version('-release');
javaaddpath("C:\Program Files\MATLAB\R2023a\java\ij.jar");
javaaddpath("C:\Program Files\MATLAB\R2023a\java\mij.jar");
% javaaddpath(genpath(char(sprintf(user_path_format, version_num, 'ij.jar'))));
MIJ.start()
MIJ.createImage(px_data)
% command = char(sprintf(save_path_format));
% MIJ.run(sprintf('saveAs("Tiff", "%s\%s.tif");', temp_loc,'test2'));
% MIJ.run(saveAs("C:/Users/pcnba/OneDrive/Desktop/test.tif"));


%Return that the function is complete
% res = sprintf('%s is done', posFingerprint);

% end

