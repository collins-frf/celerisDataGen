function main(sim_time, jj, wc_in)
%  This is a wrapper for Celeris Model
%    INPUTS:
%       forcast_date: a string for an hour to model format:
%           'YYYY-mm-ddTHH.MM.SS.Z' in zulu time
%       sim_time: wall time in seconds for simulation
%           Model is killed by wall time, calculate computational need by
%           domain and resolution
%
%   Written by Pat Lynett, USC and modified by Spicer Bak, USACE
%%  Start 
    w = warning ('off', 'all');
    clickableMap = false;
    %% INPUTS:
    forecast_count=1;       % only forecast for most recent data, can loop across this variable if more times desired
    num_frames=8000;         %num frames
    water_level_change=0;

    pause_int = .03;       
                            %virtual time is all that matters, however. we
                            %need 600 virtual seconds with approximately
                            %the same amount of virtual frames
                            
                            %frame rate will need adjusted for different
                            %wave conditions. Test is, where is elapsed at
                            %600 virtual seconds (we want 600 virtual
                            %seconds records, so record that amount of
                            %elapsed time
                            %official test on celerity for 10-01-15 WL0: 
                            %VT = 600,
                            %ET = 647,
    warmup_time = 400;
    cur_date = ['WC_', num2str(wc_in)];      %run name for write to file
    %cur_date = ['UAV_WC_str8_.5']
    avi_frame_rate = 10;     %since we are feeding FRAMES to the network and not video
                            %framerate doesnt matter except for
                            %visualization
    avi_quality = 100;
                            
    dx_target = 1;          % will set equal for y
    dataPrefix = "datafiles";        % this is where all forcing data will live
    stdOutput='D:';            % this is where image output lives 
    % used for folder labels (maybe add version prefix here)
    addTime  =  3600; % one day   % [seconds] added to back end to gather 'extra' data used for 
    %% Handle directories 
    %clean and create output directory
    frames_dir=[stdOutput];  % location for screen capture images
    
    %% Global Veriables
    % default values for all simulations:
    Courant=0.075;  % Courant number, controls time step
    Mannings_n=0.0035;  % Friction factor, quadratic law friction factor f/2*H*u^2
    min_depth=1e-5; % minimum depth allowable, in m, typically < 1/1000 of "offshore depth"
    
    % screen-size, for screen capture visualizations
    ss=get(0, 'ScreenSize');
    resol=[ss(3) ss(4)];       % automatically detect screen resolution
    
    % END TOP INPUT BLOCK
    %%
    % Database locations and corresponding "ID"
    numgrids=1;
    dir_names=cell(numgrids,1);
    wave_bc=zeros(numgrids,1);
%     grid_size=wave_bc;
    tidestatID=dir_names;
    
    % grid info -- all of this related to running multiple locations with
    % same code, could be streamlined
    ind_c=1;
    choice = ind_c;
    wave_bc(ind_c)=2;

    cd_home=cd;
    cur_timeZ=now+datenum(2010,1,1,7,0,0) - datenum(2010,1,1,0,0,0);  % simulation / forecast time

    % default values, will change only wave input boundary
    bc(1)=1; %westBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
    bc(2)=2; %eastBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
    bc(3)=2; %southBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
    bc(4)=2; %northBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
    sponge(1)=2; % westBoundary sponge width, only used if bc(1)=2
    sponge(2)=30; % eastBoundary sponge width, only used if bc(2)=2
    sponge(3)=30; % southBoundary sponge width, only used if bc(3)=2
    sponge(4)=30; % northBoundary sponge width, only used if bc(4)=2
    
    wave_boundary_str=num2str(wave_bc(choice));
    for ii=1:length(wave_boundary_str)
        
        wave_boundary_c=str2num(wave_boundary_str(ii));
        
        if wave_boundary_c==1 % waves through west boundary
            bc(1)=4; %westBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
        end
        
        if wave_boundary_c==2 % waves through east boundary
            bc(2)=4; %eastBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
        end
        
        if wave_boundary_c==3 % waves through south boundary
            bc(3)=4; %southBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
        end
        
        if wave_boundary_c==4 % waves through north boundary
            bc(4)=4; %northBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
        end
        
    end
    
    disp(['   Beginning New simulation: Loading Input Data...'])
    
    %load_FRFwave(array8m)
    load(['./WC/', cur_date, '.mat']);

    cd(cd_home)    %% get bathy
    H_toobig_factor=1;           % init for first bathy pass through (needed for below)
    load(['./bathy/celeris_gen_bathy' num2str(jj) '.mat'])
    %load(['z2017-02-28T00.15.00.Z.mat'])
    %disp(['loaded ./bathy/celeris_gen_bathy', num2str(jj), '.mat']) 
    load_bathy                % run bathy load script (loads/plots/shifts bathy with WL)
    
    % create spectrum file
    H=Hs1;
    T=1/FM1;
    theta=theta1;
    
    % determine the frequency cutoff by number of grid points in X and
    % Y, Current limitation of the model
    n_cutoff=min([1000,nx,ny]);
    
    spectrum_FRF_2D_interp(freq, direction, transpose(spec1), H, T ,n_cutoff)
    % this has a directional issue related to input shape of the spectra
    % Index in position 2 exceeds array bounds (must not exceed 72).
    % see load frf wave 

%% 
    load_bathy  % run bathy load script

    % find instruement node indices for output
    inst_ind=instruments*0;
    for i=1:length(instruments(:,1))
        inst_ind(i,1)=find(x>=instruments(i,1),1);
        inst_ind(i,2)=find(y>=instruments(i,2),1);
    end
    
    % strip range for output
    range_x=[min(inst_ind(:,1)) max(inst_ind(:,1))];
    range_y=median(inst_ind(:,2));
    
    % Determine time step
    max_depth=-min(min(B));
    ds=min(mean(diff(x)),mean(diff(y)));
    timestep=Courant*ds/sqrt(9.81*max_depth);
        
    % create colorbar for image overlay
    % cap color red limit at 2.5 m
    H = min(2.5/2,Hs1);
    
    % write the simulation control file
    file_name_cml='matlab_launch.cml';
   
    % set breaking parameters based on wave height.  These are VIZ
    % not pyhsical model parameters - they dont change simulation
    % results, these values are default determined in the model
    brk_limH=[0.65 1.5];
    brk_limthres=[0.2 0.45];
    plung_atten=0.02;
    %brk_limH = [0.1 1.0]
    %brk_limthres = [0.15 0.17]
    %plung_atten = .1
    
    if H>brk_limH(2) % energetic plunging
        brk1=brk_limthres(2);  % on/off threshold
    elseif H<brk_limH(2) % low energy spilling
        brk1=brk_limthres(1);
    else  % seomthing inbetween
        brk1=interp1(brk_limH,brk_limthres,H);
    end
    brk2=brk_limthres(2)-brk1+plung_atten;  % attenuation
       
    %------------
    % Type of input wave !not used for coastal spectrum-driven
    % waves!. just putting dummy values
    wave_type=2; %=1 for solitary wave, =2 for sine wave through boundary
    xco=0; % defines the initial x-location (m) of the crest of the solitary wave
    yco=0; % defines the initial y-location (m) of the crest of the solitary wave
    
    create_cml_photonadir(file_name_cml,timestep,Mannings_n,m_width,m_length,nx,ny,wave_type,H,T,theta,xco,yco,bc,sponge,min_depth,brk1,brk2,inst_ind,range_x,range_y,Courant)

    %insert pier image
    pier_image = imread('pier.jpg');
    [px,py,pz] = size(pier_image);
    
    delete array.txt
    delete time_axis.txt
%% set model to run         
    tic
    fprintf(' Model Running... please wait %d seconds\n', sim_time)

    import java.awt.*;
    import java.awt.event.*;
    rob=Robot;

    h = actxserver('WScript.Shell');
    h.Run('Celeris.exe');                     %Invokes
    pause(30)
    system('python C:/Users/acoll/Desktop/celerisDataGen/fix_celeris_rot.py')
    pause(warmup_time)                      % wait for sim warm up for 600 virtual seconds

    h.AppActivate('Celeris.exe');             %Brings to focus
    rob.keyPress(KeyEvent.VK_WINDOWS)
    rob.keyPress(KeyEvent.VK_UP)
    rob.keyRelease(KeyEvent.VK_WINDOWS)
    rob.keyRelease(KeyEvent.VK_UP)

    pause(10)

    %% model is running, Grab screen captures to 'visualize'
    cd_home=cd;

    for i=1:num_frames
        %imageData = screencapture(0, [350,0,resol(1)-1200,resol(2)]); % select a small desktop region, crop out user menu
        imageData = screencapture(0, [400,0,resol(1)-1200,resol(2)]); % select a small desktop region, crop out user menu
        [nxiD,nyiD,nziD]=size(imageData);
        if i==1
            imageData_anim=zeros(nxiD,nyiD,nziD,num_frames,'uint8');
            timedata=strings([1, num_frames]);
        end
        %imageData(661:660+px,101:500,:) = pier_image;
        if i==(num_frames-3)
            time_screenshot=screencapture(0, [0,0,resol(1),resol(2)]);%capture whole screen
        end
        imageData_anim(:,:,:,i)=imageData;  % straight screen capture stack
        toc
        timedata(i) = toc;
        pause(pause_int) %pause before grabbing next frame, currently .5 seconds for 1200 frames to make 10 minute timex gif
    end

    % kill celeris after vis loop is done so approximately 630 + 630
    % seconds
    ! taskkill /IM Celeris.exe > screen.txt   
    pause(30)  % sometimes celeris is slow to release file handles
    toc
    % write final surface to image
    disp([' DataGen: Model Complete!\nLoading output data and perfoming time series analysis'])

    %% write final surface to image
    cd(frames_dir)
    fname_c=['bathy_', num2str(jj), '_', cur_date, '.jpg'];
    imwrite(time_screenshot,fname_c); % save the captured image to file

    time_mat_name = ['bathy_', num2str(jj), '_', cur_date, '.mat'];
    filepath = [frames_dir, '\', num2str(time_mat_name), '.txt'];
    fid = fopen(filepath, 'wt');
    fprintf(fid, '%s\n', timedata{:});
    fclose(fid);
    save(time_mat_name, 'timedata');
    % write animat3ion to file   
    disp(['Finishing off animations and writing to file'])
    anim_name=['bathy_', num2str(jj), '_', cur_date];
    outputVideo = VideoWriter(fullfile(cd,anim_name));
    outputVideo.FrameRate = avi_frame_rate;
    outputVideo.Quality = avi_quality;
    open(outputVideo);
    writeVideo(outputVideo,imageData_anim);
    close(outputVideo);
    cd(cd_home)
    clear imageData_anim;
    clear timedata;

    disp('Done with ducksurf nowcast')
end