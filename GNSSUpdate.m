function kf = GNSSUpdate(navstate, gnssdata, kf, antlever, usegnssvel, thisimu, dt)

    param = Param();

    %% GNSS position update
    % abandon gnss vel outlier 
    gnssposstd = gnssdata(5:7, 1);
    if gnssposstd(1, 1) > 5 || gnssposstd(2, 1) > 5 || gnssposstd(3, 1) > 5
        disp(['WARNING: Abandon gnss position measurement at: ', num2str(gnssdata(1, 1))]);
    else
        % measurement innovation
        DR = diag([navstate.Rm + navstate.pos(3), (navstate.Rn + navstate.pos(3))*cos(navstate.pos(1)), -1]);
        Z = DR*(navstate.pos - gnssdata(2:4, 1))+navstate.Cbn*antlever;%N系下的NED
        
        % measurement matrix and noise matrix
        R = diag(power(gnssdata(5:7, 1), 2));%m m m
        H = zeros(3, kf.RANK);
        H(1:3, 1:3) = eye(3);
        H(1:3, 7:9) = skew(navstate.Cbn * antlever);
        
        % update covariance and state vector
        K = kf.P * H' / (H * kf.P * H' + R);
        kf.x = kf.x + K*(Z - H*kf.x);
        kf.P=(eye(kf.RANK) - K*H) * kf.P * (eye(kf.RANK) - K*H)' + K * R * K';
    end

    %% GNSS velocity update
    if usegnssvel

        % abandon gnss vel outlier 
        gnssvelstd = gnssdata(11:13, 1);
        if gnssvelstd(1, 1) > 0.5 || gnssvelstd(2, 1) > 0.5 || gnssvelstd(3, 1) > 0.5
            disp(['WARNING: Abandon gnss velocity measurement at: ', num2str(gnssdata(1, 1))]);
        else

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% TODO: add your gnss velocity update code here!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % measurement innovation, noise, matrix
    
            % update covariance and state vector
            
        end
    end
end