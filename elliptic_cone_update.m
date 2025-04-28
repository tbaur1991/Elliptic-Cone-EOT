function [X,P,meanInX,meanInY,varInX,varInY] = elliptic_cone_update(X,P,Y,sig_r,samples,numSamples,samples2,numSamples2,n_upd,meanInX,meanInY,varInX,varInY,tao,inside,filter,source)
    % dimensions
    nx = length(X);
    nMeas = size(Y,2);

    % iterate n times over n_upd measurements
    for j = 1:ceil(nMeas/n_upd)
        % get new measurement set
        if j <= floor(nMeas/n_upd)
            meas = Y(:,n_upd*j-n_upd+1:n_upd*j);
        else
            meas = Y(:,n_upd*(j-1)+1:end);
            % change number of measurements to be updated and UKF weights
            n_upd = mod(nMeas,n_upd);
            samples = samples2;
            numSamples = numSamples2;
        end

        % calculate samples for S2KF
        if strcmp(filter,'ERHM')
            Xu = [X; zeros(2*n_upd,1)]; 
            L = blkdiag(chol(P)', chol(eye(n_upd))', chol(sig_r^2*eye(n_upd))');
            z_predict = zeros(4*n_upd,numSamples);
        elseif strcmp(filter,'GAM')
            Xu = X; 
            L = chol(P)';
            z_predict = zeros(3*n_upd,numSamples);
        else
            error('Wrong filter setting!')
        end
        samples_upd = L*samples + Xu;
        if strcmp(filter,'ERHM')
            samples_upd(nx+1:nx+n_upd,:) = 1 - sqrt(1 - 0.5*(1 + erf(samples_upd(nx+1:nx+n_upd,:)/sqrt(2))));
        end
    
        % predict sigma points
        for i = 1:numSamples
            % extract sigma point parameters
            pos = samples_upd(1:3,i);
            or = samples_upd(4,i);
            a = c1(samples_upd(5,i),0,'lower');
            b = c1(samples_upd(6,i),0,'lower');
            h = c1(samples_upd(7,i),0,'lower');
            % measurements in local coordinates
            R_rot = [cos(or) -sin(or) 0; sin(or) cos(or) 0; 0 0 1];
            meas_loc = R_rot'*(meas - pos);
            % extrusion factors
            if strcmp(filter,'ERHM')
                us = samples_upd(8:7+n_upd,i);
                vs = samples_upd(8+n_upd:end,i);
            elseif strcmp(filter,'GAM')
                us = min(max(meas_loc(3,:)/h,0),1)';
                vs = zeros(n_upd,1);
            else
                us = 0; vs = 0;
            end
            if strcmp(source,'radial')
                % angle parameters
                thetas = atan2(meas_loc(2,:)*a, meas_loc(1,:)*b);
                % predict measurement
                zp = pos + R_rot*[(1 - us')*a.*cos(thetas); (1 - us')*b.*sin(thetas); us'*h + vs'];
                if strcmp(filter,'ERHM')
                    zpp = [zp; zp(3,:).^2] - [meas; meas(3,:).^2];
                elseif strcmp(filter,'GAM')
                    zpp = zp - meas;
                else
                    zpp = 0;
                end
                z_predict(:,i) = zpp(:);
            elseif strcmp(source,'projected')
                for n = 1:n_upd
                    % polygonal chain
                    ps = as_polygon((1 - us(n))*a,(1 - us(n))*b,25);
                    % project measurement to ellipse
                    zp = pos + R_rot*[project(ps, meas_loc(1:2,n)); us(n)*h + vs(n)];
                    if strcmp(filter,'ERHM')
                        zpp = [zp; zp(3)^2] - [meas(:,n); meas(3,n)^2];
                        z_predict(4*n-3:4*n,i) = zpp;
                    elseif strcmp(filter,'GAM')
                        zpp = zp - meas(:,n);
                        z_predict(3*n-2:3*n,i) = zpp;
                    end
                end
            else
                error('Wrong source setting!')
            end
        end
    
        % get predicted measurement
        z_pred = mean(z_predict, 2);
    
        % calculate measurement noise covariance matrix
        if inside
            % for asymmetric noise
            % first step measurements in local coordinates
            pos = X(1:3); or = X(4); 
            a = c1(X(5),0,'lower'); b = c1(X(6),0,'lower'); h = c1(X(7),0,'lower');
            R_rot = [cos(or) -sin(or) 0; sin(or) cos(or) 0; 0 0 1];
            meas_loc = R_rot'*(meas - pos); 
            % second step measurements inside or outside
            us = min(max(meas_loc(3,:)/h,0),1);
            inside_outside = meas_loc(1,:).^2./((1 - us)*(a - 3*sig_r)).^2 + meas_loc(2,:).^2./((1 - us)*(b - 3*sig_r)).^2 - 1;
            % indices of measurements inside
            idx = find(inside_outside < 0); 
            if strcmp(filter,'ERHM')
                % update estimate of inside measurement mean
                meanInX = 1/(1 + length(idx)/tao)*meanInX + 1/(tao + length(idx))*sum(abs(z_pred(4*idx-3)));
                meanInY = 1/(1 + length(idx)/tao)*meanInY + 1/(tao + length(idx))*sum(abs(z_pred(4*idx-2)));
                z_pred(4*idx-3) = z_pred(4*idx-3) + meanInX; z_pred(4*idx-2) = z_pred(4*idx-2) + meanInY;
                % update estimate of inside measurement variance
                varInX = 1/(1 + length(idx)/tao)*varInX + 1/(tao + length(idx))*sum((z_pred(4*idx-3) - meanInX).^2);
                varInY = 1/(1 + length(idx)/tao)*varInY + 1/(tao + length(idx))*sum((z_pred(4*idx-2) - meanInY).^2);
                % build covariance matrix
                R = repmat([sig_r^2*ones(1,2) 0 0], 1, n_upd);
                R(4*idx-3) = varInX; R(4*idx-2) = varInY; R = diag(R);
            elseif strcmp(filter,'GAM')
                % update estimate of inside measurement mean
                meanInX = 1/(1 + length(idx)/tao)*meanInX + 1/(tao + length(idx))*sum(abs(z_pred(3*idx-2)));
                meanInY = 1/(1 + length(idx)/tao)*meanInY + 1/(tao + length(idx))*sum(abs(z_pred(3*idx-1)));
                z_pred(3*idx-2) = z_pred(3*idx-2) + meanInX; z_pred(3*idx-1) = z_pred(3*idx-1) + meanInY;
                % update estimate of inside measurement variance
                varInX = 1/(1 + length(idx)/tao)*varInX + 1/(tao + length(idx))*sum((z_pred(3*idx-2) - meanInX).^2);
                varInY = 1/(1 + length(idx)/tao)*varInY + 1/(tao + length(idx))*sum((z_pred(3*idx-1) - meanInY).^2);
                % build covariance matrix
                R = sig_r^2*ones(1,3*n_upd);
                R(3*idx-2) = varInX; R(3*idx-1) = varInY; R = diag(R(:));
            else
                R = 0;
            end
        else
            if strcmp(filter,'ERHM')
                R = diag(repmat([sig_r^2*ones(1,2) 0 0], 1, n_upd));
            elseif strcmp(filter,'GAM')
                R = sig_r^2*eye(3*n_upd);
            else
                R = 0;
            end
        end
        
        % get state measurement cross-covariance matrix
        z_p = z_predict - z_pred;
        YY = z_p*z_p'/numSamples + R;
        C = (samples_upd(1:nx,:) - X)*z_p'/numSamples;
    
        % do S2KF measurement update
        K = C/YY;
        if strcmp(filter,'ERHM')
            X = X + K*(zeros(4*n_upd,1) - z_pred);
        elseif strcmp(filter,'GAM')
            X = X + K*(zeros(3*n_upd,1) - z_pred);
        end
        P = P - K*C';
        P = 0.5*(P+P');
    end
end