function kf = ZUPTUpdate(navstate, kf, cfg)
    %% 适配代码做的参数准备
    vel = navstate.vel;

    %% ZUPT update
    % measurement innovation
    Z = vel-[0,0,0]';

    % measurement matrix and noise matrix
    R = diag(cfg.zupt_noise);
    H = zeros(3, kf.RANK);
    H(1:3, 4:6) = eye(3);

    % update covariance and state vector
    K = kf.P * H' / (H * kf.P * H' + R);
    kf.x = kf.x + K*(Z - H*kf.x);
    kf.P = (eye(kf.RANK) - K*H) * kf.P * (eye(kf.RANK) - K*H)' + K * R * K';
end