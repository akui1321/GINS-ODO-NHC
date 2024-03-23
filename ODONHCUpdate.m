function kf = ODONHCUpdate(navstate, odonhc_vel, kf, cfg, thisimu, dt)

    param = Param();

    %% 适配代码做的参数准备
    pos = navstate.pos;
    vel = navstate.vel;
    Cbn = navstate.Cbn;
    Cnb = Cbn';
    Rm = navstate.Rm;
    Rn = navstate.Rn;
    e_omega = param.wie;
    omega_ib_b = thisimu(2:4, 1) / dt;
    Cbv = Euler2DCM(cfg.installangle);

    wie_n = [e_omega*cos(pos(1)); 0; -e_omega*sin(pos(1))];
    wen_n = [vel(2)/(Rn+pos(3)); -vel(1)/(Rm+pos(3)); -vel(2)*tan(pos(1))/(Rn+pos(3))];
    win_n = wie_n + wen_n;

    %% ODO/NHC update
    % measurement innovation
    Z = Cbv*Cnb*vel-Cbv*Cnb*skew(win_n)*Cbn*cfg.odolever-Cbv*skew(cfg.odolever)*omega_ib_b-odonhc_vel;

    if(~cfg.useodo)
        Z(1,:)=0;
    end

    if(~cfg.usenhc)
        Z(2:3,:)=0;
    end

    % measurement matrix and noise matrix
    R = diag(cfg.odonhc_measnoise);
    H = zeros(3, kf.RANK);
    H(1:3, 4:6) = Cbv*Cnb;
    H(1:3, 7:9) = -Cbv*Cnb*skew(vel);
    H(1:3, 10:12) = -Cbv*skew(cfg.odolever);
    H(1:3, 16:18) = -Cbv*skew(cfg.odolever)*diag(omega_ib_b);

    if(~cfg.useodo)
        H(1,:)=0;
    end

    if(~cfg.usenhc)
        H(2:3,:)=0;
    end

    % update covariance and state vector
    K = kf.P * H' / (H * kf.P * H' + R);
    kf.x = kf.x + K*(Z - H*kf.x);
    kf.P = (eye(kf.RANK) - K*H) * kf.P * (eye(kf.RANK) - K*H)' + K * R * K';

end