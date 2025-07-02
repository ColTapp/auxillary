% SPLIT_ROIS_AND_EXPORT_TIMELAPSE
% Allows the user to split large microscopy timelapse images into ROIs, export each ROI as a TIF stack,
% and export a movie of each ROI. User can draw ROIs for each region/slice. Handles memory by loading slices only.
%
% USAGE:
%   1. Define a cell array of structs with parameters for each image set below.
%   2. Run the script. It will prompt for ROI selections and save crops and movies per ROI.


reload=1; %leave at 1 unless troubleshooting
if reload
    close all
    clear
    reload=1;
end


%% --------- USER IMAGE SET PARAMETERS (EDIT THIS BLOCK) -------------------
% Define your image sets here as an array of structs

% imageSets: Cell array of structs for each image set, ready to use in batch processing scripts.
% Each struct parameterizes a single image set according to your listed variables.

imageSets = {
    struct(...
        'maindir', 'E:\Microscopy\Julian\optotracer\finalruns\r1\Folder_20240417_20240415_ebba6well-top-JE-5_1149-5_agrA-5_mid-JE-6_1149-6_agrA-6_bot-icaA-4_icaA-5_icaA-6_1431', ...
        'spotrows', 3, ...
        'thr1', [0.001, 0.45], ...
        'thr2', [0.0015, 0.25], ...
        'thr3', [0.002, 0.1], ...
        'skipchan', [], ...
        'skiprow', [1,2], ...
        'maxfrrem', 1, ...
        'deltaT', 30, ...
        'addT', 0, ...
        'displayresize', 0.075, ...
        'rowoverlap', 0.06 ...
    ), ...
    struct(...
        'maindir', 'F:\Julian\optotracer\finalruns\r2\Folder_20240502_20240430_ebba6well-top-JE-6_JE2-6_JE2-6_mid-1149-6_1149-6_agrA-6_bot-icaA-6_icaA-6_agrA-6_1452', ...
        'spotrows', 3, ...
        'thr1', [0.001, 0.45], ...
        'thr2', [0.0015, 0.25], ...
        'thr3', [0.002, 0.1], ...
        'skipchan', [], ...
        'skiprow', [], ...
        'maxfrrem', 1, ...
        'deltaT', 30, ...
        'addT', 5*60, ...
        'displayresize', 0.075, ...
        'rowoverlap', 0.06 ...
    ), ...
    struct(...
        'maindir', 'F:\Julian\optotracer\finalruns\r3\Folder_20240506_20240504_ebba6well-top-JE-6_JE2-6_1149-6_mid-1149-6_agrA_agrA-6_bot-icaA-6_icaA-6_icaA-6_1545', ...
        'spotrows', 3, ...
        'thr1', [0.001, 0.45], ...
        'thr2', [0.0015, 0.25], ...
        'thr3', [0.002, 0.1], ...
        'skipchan', [], ...
        'skiprow', [], ...
        'maxfrrem', 1, ...
        'deltaT', 30, ...
        'addT', 0, ...
        'displayresize', 0.075, ...
        'rowoverlap', 0.06 ...
    ) ...
};



savedirname = 'savesV2'; % Name for savedir
bitdepth = 16;           % Image bitdepth
vidfr = 5;               % Video framerate
liveimages = 0.9;        % Display scaling/show images to user
chname1 = 'Ch1'; chname2 = 'Ch2'; chname3 = 'Ch3'; % Channel names
col1 = [1, 0, 0]; col2 = [1, 0, 0]; col3 = [0, 0, 1]; % Channel colors

%% --------- END OF USER PARAMETERS BLOCK ----------------------------------

for setIdx = 1:length(imageSets)
    % Load parameter struct for this image set
    p = imageSets{setIdx};
    % Intensity scaling
    maxin_val = (2^bitdepth)-1;
    p.thr1 = p.thr1 * maxin_val;
    p.thr2 = p.thr2 * maxin_val;
    p.thr3 = p.thr3 * maxin_val;

    % Registration optimizer/metric (default params, can be tuned)
    [optimizer, metric] = imregconfig('multimodal');
    optimizer.InitialRadius = 6.e-04;
    optimizer.MaximumIterations = 300;
    optimizer.GrowthFactor = 1.05;
    optimizer.Epsilon = 1.5e-06;

    % Directory setup
    if p.maindir(end) ~= filesep, p.maindir = [p.maindir, filesep]; end
    savedir    = [p.maindir, savedirname, filesep];
    vidsavedir = [savedir, 'movies', filesep];
    datsavedir = [savedir, 'crops', filesep];
    if ~isfolder(savedir), mkdir(savedir); end
    if ~isfolder(vidsavedir), mkdir(vidsavedir); end
    if ~isfolder(datsavedir), mkdir(datsavedir); end

    % List of all image position directories
    cd(p.maindir)
    alldir = dir('_*');
    posL = 1:numel(alldir);

    for posi = posL
        dird = [p.maindir, alldir(posi).name, filesep, 'stack1'];
        r = bfGetReader([dird, filesep, 'frame_t_0.ets']); % Use Bio-Formats to read only one image
        r = loci.formats.Memoizer(r);
        n_im = r.getSizeT - p.maxfrrem;
        n_channel = r.getSizeC;
        szx = r.getSizeX;
        szy = r.getSizeY;
        omeMeta = r.getMetadataStore();
        PxinUmX = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER));

        % Determine how to split image vertically into rows (for multi-row setups)
        if p.spotrows > 1
            splitpos = nan(p.spotrows,2);
            splitpos(1,1) = 1;
            splitpos(1,2) = round(szy/p.spotrows + p.rowoverlap*szy);
            for i = 2:p.spotrows
                splitpos(i,1) = round(splitpos(i-1,2) - 2*p.rowoverlap*szy);
                splitpos(i,2) = round(i*szy/p.spotrows + p.rowoverlap*szy);
            end
            splitpos(end,2) = szy;
        end

        % ROI selection
        if p.spotrows > 0
            roinames = {};
            roipos_all = {};
            for irow = 1:p.spotrows
                if isfield(p, 'skiprow') && ismember(irow, p.skiprow), continue, end
                msgbox('Please wait, image loading...');
                % Load and display one row at reduced size for ROI selection
                if p.spotrows > 1
                    bf = mat2gray(imresize(bfGetPlane(r, 1, 1, splitpos(irow,1), szx, splitpos(irow,2)-splitpos(irow,1)), ...
                        p.displayresize), p.thr1);
                else
                    bf = mat2gray(imresize(bfGetPlane(r, 1), p.displayresize), p.thr1);
                end
                closeMsgboxes();
                figure; imshow(bf)
                % Get number of ROIs for this row
                nroi0 = [];
                while isempty(nroi0)
                    answer = inputdlg('Number of ROIs:', 'ROI', [1], {'1'});
                    if isempty(answer), return; end
                    nroi0 = str2double(answer{1});
                    if isnan(nroi0) || nroi0 < 1, nroi0 = []; end
                end
                % Get ROI names and positions
                roinames_row = cell(nroi0,1);
                roipos_row = nan(nroi0,4);
                for i = 1:nroi0
                    roinames_row{i} = inputdlg(['Name for ROI ' num2str(i)], 'ROI names', 1, {['ROI',num2str(i)]});
                    title(['Draw rectangle for ', char(roinames_row{i})], 'Interpreter', 'none');
                    h = drawrectangle;
                    roipos_row(i,:) = customWait(h);
                end
                roinames = [roinames; roinames_row];
                roipos_all = [roipos_all; {roipos_row}];
            end
            % Combine ROI positions from all rows and correct for scale/offset
            roipos = [];
            for irow = 1:numel(roipos_all)
                roipos_row = floor(roipos_all{irow} * 1/p.displayresize);
                if p.spotrows > 1
                    roipos_row(:,2) = roipos_row(:,2) + round((irow-1)*szy/p.spotrows - (irow-1)*p.rowoverlap*szy);
                end
                roipos = [roipos; roipos_row];
            end
        else
            % Whole image as single ROI
            roinames = {''};
            roipos = [];
        end

        % Prepare folders for each ROI
        nroi = max(1, size(roipos,1));
        datsavedirR = cell(1, nroi);
        for i = 1:nroi
            datsavedirR{i} = [datsavedir, 'pos', num2str(posi, '%02.f'), 'roi', num2str(i, '%02.f'), '_', char(roinames{i}), filesep];
            if ~isfolder(datsavedirR{i}), mkdir(datsavedirR{i}); end
        end

        % --- MAIN LOOP: EXPORT TIFFS AND MOVIE FOR EACH ROI ---
        for i = 1:nroi
            % Define ROI cropping region
            if nroi > 1
                curroi = roipos(i,:);
                curroi(1:2) = max(curroi(1:2),1);
                curroi(3) = min(curroi(3), szx-curroi(1)-1);
                curroi(4) = min(curroi(4), szy-curroi(2)-1);
            end
            % Prepare video writer
            vidname = [vidsavedir, 'pos', num2str(posi, '%02.f'), 'roi', num2str(i, '%02.f'), '_', char(roinames{i}), '.avi'];
            vidOV = VideoWriter(vidname, 'Motion JPEG AVI');
            vidOV.FrameRate = vidfr; open(vidOV);

            for ti = 1:n_im
                % Channel indices for this timepoint
                baseIndex = (ti - 1) * n_channel;
                if n_channel == 2
                    bi = baseIndex + 1; gi = baseIndex + 2;
                elseif n_channel == 3
                    if ~isempty(p.skipchan)
                        if p.skipchan == 2
                            bi = baseIndex + 1; gi = baseIndex + 3;
                        elseif p.skipchan == 1
                            bi = baseIndex + 2; gi = baseIndex + 3;
                        end
                    else
                        bi = baseIndex + 1; gi = baseIndex + 2; ri = baseIndex + 3;
                    end
                end
                % Crop and scale images
                if nroi > 1
                    bcr = mat2gray(bfGetPlane(r, bi, curroi(1), curroi(2), curroi(3), curroi(4)), p.thr1);
                    gcr = mat2gray(bfGetPlane(r, gi, curroi(1), curroi(2), curroi(3), curroi(4)), p.thr2);
                else
                    bcr = mat2gray(bfGetPlane(r, bi), p.thr1);
                    gcr = mat2gray(bfGetPlane(r, gi), p.thr2);
                end
                if n_channel > 2 && isempty(p.skipchan)
                    if nroi > 1
                        rcr = mat2gray(bfGetPlane(r, ri, curroi(1), curroi(2), curroi(3), curroi(4)), p.thr3);
                    else
                        rcr = mat2gray(bfGetPlane(r, ri), p.thr3);
                    end
                else
                    rcr = zeros(size(bcr));
                end
                genname = [datsavedirR{i}, 'roi', num2str(i,'%02.f'), '_'];
                % Write TIFFs for each channel
                if bitdepth == 16
                    imwrite(uint16(bcr*maxin_val), [genname, chname1, 'Fr', num2str(ti,'%03.f'), '.tif']);
                    if n_channel > 1, imwrite(uint16(gcr*maxin_val), [genname, chname2, 'Fr', num2str(ti,'%03.f'), '.tif']); end
                    if n_channel > 2 && isempty(p.skipchan), imwrite(uint16(rcr*maxin_val), [genname, chname3, 'Fr', num2str(ti,'%03.f'), '.tif']); end
                else
                    imwrite(uint8(bcr*maxin_val), [genname, chname1, 'Fr', num2str(ti,'%03.f'), '.tif']);
                    if n_channel > 1, imwrite(uint8(gcr*maxin_val), [genname, chname2, 'Fr', num2str(ti,'%03.f'), '.tif']); end
                    if n_channel > 2 && isempty(p.skipchan), imwrite(uint8(rcr*maxin_val), [genname, chname3, 'Fr', num2str(ti,'%03.f'), '.tif']); end
                end
                % Create movie frame and write to video
                ov = movieframecreation(bcr, gcr, rcr, n_channel, col1, col2, col3, p.deltaT, ti, PxinUmX, p.displayresize, p.addT);
                if liveimages > 0
                    imshow(imresize(ov, liveimages)); pause(0.001);
                end
                writeVideo(vidOV, ov);
            end
            close(vidOV);
        end
    end
end

%% --- Helper function to close all msgboxes (optional utility)
function closeMsgboxes()
    h = findall(0,'Type','figure','Tag','Msgbox_Close');
    if ~isempty(h), close(h); end
end


function ov = movieframecreation(bf, gfp, rfp, n_channel,...
    col1, col2, col3, deltaT, ti, PxinUmX, resz, addT)
    %MOVIEFRAMECREATION Creates an RGB overlay frame from grayscale images for movie generation.
    %   bf, gfp, rfp: Grayscale images for the brightfield, GFP, and RFP channels
    %   n_channel: Number of channels (1 or more)
    %   col1, col2, col3: Color coefficients for channel mixing
    %   deltaT: Time between frames (minutes)
    %   ti: Current frame index
    %   PxinUmX: Pixel size in micrometers
    %   resz: Resize factor
    %   addT: Additional time in minutes
    %
    %   Returns:
    %       ov: RGB image frame (uint8)
    
    % Resize brightfield image
    bf = imresize(bf, resz);
    
    if n_channel > 1
        % Resize and mix in GFP and RFP channels
        gfp = imresize(gfp, resz);
        rfp = imresize(rfp, resz);
    
        % Red channel: BF + GFP*col1(1) + RFP*col2(1)
        rchan = bf + gfp*col1(1) + rfp*col2(1);
    
        % Green channel: BF + GFP*col1(2) + RFP*col2(2)
        gchan = bf + gfp*col1(2) + rfp*col2(2);
    
        % Blue channel: BF + GFP*col1(3)*1.5 - RFP*col2(2)*0.5
        bchan = bf + gfp*col1(3)*1.5 - rfp*col2(2)*0.5;
    
        % Clip channel values to [0, 0.999]
        rchan(rchan > 1) = 0.999;
        gchan(gchan > 1) = 0.999;
        bchan(bchan > 1) = 0.999;
    
        % Concatenate channels and scale to uint8
        ov = cat(3, rchan, gchan, bchan);
        ov = uint8(ov * 256);
    else
        % Single channel: replicate BF to all RGB channels
        ov = cat(3, bf, bf, bf);
        ov = uint8(ov * 256);
    end
    
    % Add scale bar and label
    if PxinUmX > 0.1
        % 1 mm scale bar for large pixel sizes
        barl = (1000 / PxinUmX) * resz;
        label = '1 mm';
        ypos = 0.94;
    else
        % 100 μm scale bar for small pixel sizes
        barl = (100 / PxinUmX) * resz;
        label = ['100 ', sprintf('%s', char(956), 'm')]; % μm
        ypos = 0.92;
    end
    
    % Calculate line and text positions
    linepos = [size(ov, 2)*0.97 - barl, size(ov, 1)*0.97, size(ov, 2)*0.97, size(ov, 1)*0.97];
    ov = insertShape(ov, 'Line', linepos, 'Color', 'white', 'LineWidth', 5);
    txtpos = [(linepos(3) + linepos(1)) / 2, linepos(2) * ypos];
    
    % Font size scales with image area
    txsize = round(6.5e-05 * (size(bf,1) * size(bf,2)));
    if PxinUmX <= 0.1
        txsize = round(0.5 * txsize);
    end
    
    ov = insertText(ov, txtpos, label, 'FontSize', txsize, 'TextColor', 'white', ...
        'AnchorPoint', 'CenterTop', 'BoxOpacity', 0);
    
    % Add time string (HH:MM) to corner
    time1 = (ti-1) * deltaT + addT * 60;   % total time in minutes
    tH = floor(time1 / 60);                % hours
    tM = time1 - tH * 60;                  % minutes
    
    % Pad with leading zero if needed
    tH = sprintf('%02d', tH);
    tM = sprintf('%02d', round(tM));
    tstring = [tH, ':', tM];
    
    ov = insertText(ov, [1, 1], tstring, 'FontSize', txsize, 'BoxColor', 'white');

end


function pos = customWait(hROI)
    % Listen for mouse clicks on the ROI
    l = addlistener(hROI,'ROIClicked',@clickCallback);
    
    % Block program execution
    uiwait;
    
    % Remove listener
    delete(l);
    
    % Return the current position
    pos = hROI.Position;
end

function clickCallback(~,evt)
    if strcmp(evt.SelectionType,'double')
        uiresume;
    end
end