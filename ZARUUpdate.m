function kf = ZARUUpdate(navstate, kf, thisimu, dt, cfg)
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

    wie_n = [e_omega*cos(pos(1)); 0; -e_omega*sin(pos(1))];
    wen_n = [vel(2)/(Rn+pos(3)); -vel(1)/(Rm+pos(3)); -vel(2)*tan(pos(1))/(Rn+pos(3))];
    win_n = wie_n + wen_n;


    %% ZARU update
    % measurement innovation
    Z = omega_ib_b-Cnb*win_n-[0,0,0]';

    % measurement matrix and noise matrix
    R = diag(cfg.zaru_noise);
    H = zeros(3, kf.RANK);
    H(1:3, 10:12) = eye(3);
    H(1:3, 16:18) = diag(omega_ib_b);

    % update covariance and state vector
    K = kf.P * H' / (H * kf.P * H' + R);
    kf.x = kf.x + K*(Z - H*kf.x);
    kf.P = (eye(kf.RANK) - K*H) * kf.P * (eye(kf.RANK) - K*H)' + K * R * K';
end