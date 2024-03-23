%% define parameters and load process config
param = Param();
cfg = ProcessConfig();
cfg.Cbv = Euler2DCM(cfg.installangle);
refpath = "dataset/REF_NAV.nav";
refdata = load(refpath);

%% load data
% imudata
imufid = fopen(cfg.imufilepath);
imudata = fread(imufid, [7, inf], 'double');
fclose(imufid);
imudata = imudata';
imustarttime = imudata(1, 1);
imuendtime = imudata(end, 1);

% gnss data
gnssdata = load(cfg.gnssfilepath);
gnssdata(:, 2:3) = gnssdata(:, 2:3) * param.D2R;
if (size(gnssdata, 2) < 13)
    cfg.usegnssvel = false;
end
gnssstarttime = gnssdata(1, 1);
gnssendtime = gnssdata(end, 1);

% odo data
if (cfg.useodonhc)
    odofid = fopen(cfg.odofilepath);
    ododata = fread(odofid, [2, inf], 'double');
    fclose(odofid);
    ododata = ododata';
end


%% save result
navpath = [cfg.outputfolder, '/NavResult'];
if cfg.usegnssvel
    navpath = [navpath, '_GNSSVEL'];
end
if cfg.useodonhc
    navpath = [navpath, '_ODONHC'];
end
navpath = [navpath, '.nav'];
navfp = fopen(navpath, 'wt');

imuerrpath = [cfg.outputfolder, '/ImuError.txt'];
imuerrfp = fopen(imuerrpath, 'wt');

stdpath = [cfg.outputfolder, '/NavSTD.txt'];
stdfp = fopen(stdpath, 'wt');


%% get process time
% start time and end time
if imustarttime > gnssstarttime
    starttime = imustarttime;
else
    starttime = gnssstarttime;
end
if imuendtime > gnssendtime
    endtime = gnssendtime;
else
    endtime = imuendtime;
end
if cfg.starttime < starttime
    cfg.starttime = starttime;
end
cfg.endtime = endtime;

if cfg.useodonhc
    odoupdatetime = ceil(cfg.starttime) + 0.5; % 第一个更新时刻在0.5s时
    num_to_getvel = 20; % 20个历元平均获取odovel
end

% data in process interval
imudata = imudata(imudata(:,1) >= cfg.starttime, :);
imudata = imudata(imudata(:,1) <= cfg.endtime, :);
gnssdata = gnssdata(gnssdata(:, 1) >= cfg.starttime, :);
gnssdata = gnssdata(gnssdata(:, 1) <= cfg.endtime, :);


%% for debug
disp("Start GNSS/INS Processing!");
lastprecent = 0;


%% initialization 
[kf, navstate] = Initialize(cfg);
laststate = navstate;

% data index preprocess
lastimu = imudata(1, :)';
thisimu = imudata(1, :)';
imudt = thisimu(1, 1) - lastimu(1, 1);
gnssindex = 1;
while gnssdata(gnssindex, 1) < thisimu(1, 1)
    gnssindex = gnssindex + 1;
end

refindex = 1;
while refdata(refindex,2) < thisimu(1, 1)
    refindex = refindex + 1;
end

zeroindex = 0;

if cfg.useodonhc
    odoindex = 1;
    while ododata(odoindex, 1) < thisimu(1, 1) && odoindex < size(ododata, 1)
        odoindex = odoindex + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MAIN PROCEDD PROCEDURE!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for imuindex = 2:size(imudata, 1)-1

    %% set value of last state
    lastimu = thisimu;
    laststate = navstate;
    thisimu = imudata(imuindex, :)';
    imudt = thisimu(1, 1) - lastimu(1, 1);


    %% compensate IMU error
    thisimu(2:4, 1) = (thisimu(2:4, 1) - imudt * navstate.gyrbias)./(ones(3, 1) + navstate.gyrscale);
    thisimu(5:7, 1) = (thisimu(5:7, 1) - imudt * navstate.accbias)./(ones(3, 1) + navstate.accscale);

    
    %% adjust GNSS index
    while (gnssindex <= size(gnssdata, 1) && gnssdata(gnssindex, 1) < lastimu(1, 1))
        gnssindex = gnssindex + 1;
    end
    % check whether gnss data is valid
    if (gnssindex > size(gnssdata, 1))
        disp('GNSS file END!');
        break;
    end

    %% adjust REF index
    while (refindex <= size(refdata, 1) && refdata(refindex, 2) < lastimu(1, 1))
        refindex = refindex + 1;
    end
    % check whether ref data is valid
    if (refindex > size(refdata, 1))
        disp('REF file END!');
        break;
    end

    %% determine whether gnss update is required
    if lastimu(1, 1) == gnssdata(gnssindex, 1)
        % do gnss update for the current state
        thisgnss = gnssdata(gnssindex, :)';
        imudt = thisimu(1, 1) - lastimu(1, 1);
        kf = GNSSUpdate(navstate, thisgnss, kf, cfg.antlever, cfg.usegnssvel, lastimu, imudt);
        [kf, navstate] = ErrorFeedback(kf, navstate);
        gnssindex = gnssindex + 1;
        laststate = navstate;
        
        % do propagation for current imu data
        imudt = thisimu(1, 1) - lastimu(1, 1);
        navstate = InsMech(laststate, lastimu, thisimu);
        kf = InsPropagate(navstate, thisimu, imudt, kf, cfg.corrtime);
    elseif (lastimu(1, 1) <= gnssdata(gnssindex, 1) && thisimu(1, 1) > gnssdata(gnssindex, 1))
        
        % ineterpolate imu to gnss time
        [firstimu, secondimu] = interpolate(lastimu, thisimu, gnssdata(gnssindex, 1));
        
        % do propagation for first imu
        imudt = firstimu(1, 1) - lastimu(1, 1);
        navstate = InsMech(laststate, lastimu, firstimu);
        kf = InsPropagate(navstate, firstimu, imudt, kf, cfg.corrtime);

        % do gnss update
        thisgnss = gnssdata(gnssindex, :)';
        kf = GNSSUpdate(navstate, thisgnss, kf, cfg.antlever, cfg.usegnssvel, firstimu, imudt);
        [kf, navstate] = ErrorFeedback(kf, navstate);
        gnssindex = gnssindex + 1;
        laststate = navstate;
        lastimu = firstimu;

        % do propagation for second imu
        imudt = secondimu(1, 1) - lastimu(1, 1);
        navstate = InsMech(laststate, lastimu, secondimu);
        kf = InsPropagate(navstate, secondimu, imudt, kf, cfg.corrtime);
    else
        %% only do propagation
        % INS mechanization
        navstate = InsMech(laststate, lastimu, thisimu);
        % error propagation
        kf = InsPropagate(navstate, thisimu, imudt, kf, cfg.corrtime);
    end

    %% ZUPT&&ZARU
    if (cfg.usezupt||cfg.usezaru)
        squared_sum = sqrt(refdata(refindex, 6).^2 + refdata(refindex, 7).^2 + refdata(refindex, 8).^2);
        if (squared_sum < 0.01)
            zeroindex = zeroindex + 1;
            if (mod(zeroindex, cfg.zeroupdatetimes) == 0)
                %% ZUPT
                if cfg.usezupt
                    kf = ZUPTUpdate(navstate, kf, cfg);
                    [kf, navstate] = ErrorFeedback(kf, navstate);
                end
    
                %% ZARU
                if cfg.usezaru
                    kf = ZARUUpdate(navstate, kf, thisimu, imudt, cfg);
                    [kf, navstate] = ErrorFeedback(kf, navstate);
                end
            end
        end
    end

    %% nhcodo
    if cfg.useodonhc
        %% update odo index
        while ododata(odoindex, 1) < thisimu(1, 1) && odoindex < size(ododata, 1)
            odoindex = odoindex + 1;
        end

        %% odonhc udpate
        if (cfg.useodonhc && lastimu(1, 1) <= odoupdatetime && odoupdatetime < thisimu(1, 1))
            startindex = odoindex - round(num_to_getvel / 2);
            endindex = odoindex + round(num_to_getvel / 2);
            if (startindex < 1)
                startindex = 1;
            end
            if (endindex > size(ododata, 1))
                endindex = size(ododata, 1);
            end
           
            % get odovel and update
            [odovel, valid] = getodovel(ododata(startindex:endindex, :), thisimu(1, 1));
            if valid
                odonhc_vel = [odovel; 0; 0];
                kf = ODONHCUpdate(navstate, odonhc_vel, kf, cfg, thisimu, imudt);
                [kf, navstate] = ErrorFeedback(kf, navstate);
            end
            odoupdatetime = odoupdatetime + 1 / cfg.odoupdaterate;
        end
    end

%     %% GNSS position update 
% 
%     if (thisimu(1, 1) >= gnssdata(gnssindex, 1) || abs(thisimu(1, 1) - gnssdata(gnssindex, 1)) < imudt / 3)    
%         t1 = clock;
%         thisgnss = gnssdata(gnssindex, :)';
%         kf = GNSSUpdate(navstate, thisgnss, kf, cfg.antlever, cfg.usegnssvel, thisimu, imudt);
%         [kf, navstate] = ErrorFeedback(kf, navstate);
%         gnssindex = gnssindex + 1;
%     end
    

    %% save data
    % write navresult to file
    nav = zeros(11, 1);
    nav(2, 1) = navstate.time;
    nav(3:5, 1) = [navstate.pos(1) * param.R2D; navstate.pos(2) * param.R2D; navstate.pos(3)];
    nav(6:8, 1) = navstate.vel;
    nav(9:11, 1) = navstate.att * param.R2D;
    fprintf(navfp, '%2d %12.6f %12.8f %12.8f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n', nav);

    % write imu error
    imuerror = zeros(13, 1);
    imuerror(1, 1) = navstate.time;
    imuerror(2:4, 1) = navstate.gyrbias * param.R2D * 3600;
    imuerror(5:7, 1) = navstate.accbias * 1e5;
    imuerror(8:10, 1) = navstate.gyrscale * 1e6;
    imuerror(11:13, 1) = navstate.accscale * 1e6;
    fprintf(imuerrfp, '%12.6f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n', imuerror);

    %write state std
    std = zeros(1, 22);
    std(1) = navstate.time;
    for idx=1:21
        std(idx + 1) = sqrt(kf.P(idx, idx));
    end
    std(8:10) = std(8:10) * param.R2D;
    std(11:13) = std(11:13) * param.R2D *3600;
    std(14:16) = std(14:16) * 1e5;
    std(17:22) = std(17:22) * 1e6;
    fprintf(stdfp, '%12.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \n', std);


    %% print processing information
    if (imuindex / size(imudata, 1) - lastprecent > 0.01) 
        disp("processing " + num2str(floor(imuindex * 100 / size(imudata, 1))) + " %!");
        lastprecent = imuindex / size(imudata, 1);
    end
end

% close file
fclose(imuerrfp);
fclose(navfp);
fclose(stdfp);

disp("GNSS/INS Integration Processing Finished!");

