classdef GaussianWake < WakeModel
    % This wake model is developed based on
    % 1 - A new analytical model for wind-turbine wakes, M. Bastankhah, 2014
    % 2 - Experimental and theoretical study of yawed-turbine wakes, M. Bastankhah, 2016
    % 3 - A tutorial on the synthesis and validation of a closed-loop wind farm controller using a steady-state surrogate model, Bart M Doekemeijer, 2018
    % 4 - Wind plant power optimization through yaw control using a parametric model for wake effects - a cfd simulation study, P.M.O. Gebraad, 2016
    
    properties % Model Parameters
        %% Tuned for both
        def_alpha = 3.6429; % near wake parameter
        def_beta  = 0.089682; % near wake parameter
        def_ka    = 0.27694;% wake expansion parameter
        def_kb    = 0.0051879;% wake expansion parameter
        
        vel_alpha = 3.6429; % near wake parameter
        vel_beta  = 0.089682; % near wake parameter
        vel_ka    = 0.27694;% wake expansion parameter
        vel_kb    = 0.0051879;% wake expansion parameter
    end
    properties
        TI  = 0.05;
    end
    methods
        function Obj = GaussianWake()
            Obj.inflowAngle = 0;
        end
        function Obj = deflection(Obj,TI)
            rotorD   = Obj.turbine.bladeLength*2;
            Xi = Obj.turbine.farmPosX;
            Yi = Obj.turbine.farmPosY/rotorD;
            x_norm  = (Obj.axialDist-Xi)./rotorD;
            yawAngle = deg2rad(Obj.turbine.yawAngle);
            infAngle = deg2rad(Obj.inflowAngle);
            flowAngle = yawAngle + infAngle;
            
            CT =  Obj.turbine.rotorDat.CT;
            yawOffset = Obj.turbine.yawRotOffset;
            if exist('TI','var')
                I = TI;
            else
                I = Obj.TI;
            end
            
            % Eq. 6.12 [2] - Skew angle in near wake
            theta = yawOffset + (((0.3*flowAngle)./cos(flowAngle)).*(1-sqrt(1-(CT.*cos(flowAngle)))))+infAngle;
            
            % Eq 7 [3] - Wake expansion parameters
            ky = Obj.def_ka*I + Obj.def_kb;
            kz = Obj.def_ka*I + Obj.def_kb;
            
            % Eq 7.3 [2] - Boundary between near and far wake
            x0_norm = (cos(flowAngle).*(1+sqrt(1-CT)))./(sqrt(2).*(((Obj.def_alpha*4).*I)+((Obj.def_beta*2).*(1-sqrt(1-CT)))));
            
            % Eq 7.2 [2] - wake standard deviations
            sigma_y_norm = ky.*((x_norm-x0_norm))+((cos(flowAngle))./sqrt(8));
            sigma_z_norm = kz.*((x_norm-x0_norm))+(1./sqrt(8));
            
            % p.534 par.1 [2]
            nearWake_norm = (tan(theta).*x_norm);
            
            % Eq 7.4 [2]
            term1 = theta.*x0_norm;
            term2 = (theta./14.7).*(sqrt((cos(flowAngle))./(ky.*kz.*CT)));
            term3 = 2.9+(1.3.*sqrt(1-CT))-CT;
            innestTerm1 = 1.6.*sqrt((8.*sigma_y_norm.*sigma_z_norm)./(cos(flowAngle)));
            term4 = log(((1.6+sqrt(CT)).*(innestTerm1-sqrt(CT)))./((1.6-sqrt(CT)).*(innestTerm1+sqrt(CT))));
            farWake_norm = term1+(term2.*term3)*term4;
            
            % p.534 par.1 [2] - Combination of two wake sectors.
            deflectionYaw = farWake_norm.*(x_norm>x0_norm)+nearWake_norm.*(x_norm<=x0_norm);
            
            Obj.wakeCenterLoc = deflectionYaw+Yi;
        end
        function Obj = velocity(Obj)
            rotorD   = Obj.turbine.bladeLength*2;
            yawOffset = Obj.turbine.yawRotOffset;
            yawAngle = deg2rad(Obj.turbine.yawAngle);
            infAngle = deg2rad(Obj.inflowAngle);
            flowAngle = yawAngle + infAngle;
            Xi = Obj.turbine.farmPosX./rotorD;
            Yi = Obj.turbine.farmPosY./rotorD;
            Zi = Obj.turbine.farmPosZ./rotorD;
            x_norm = Obj.axialDist./rotorD-Xi;
            y_norm = Obj.lateralDist./rotorD-Yi;
            z_norm = Obj.verticalDist./rotorD-Zi;
            rotorD = 1;
            
            Ui = Obj.inflowVelocity./cos(infAngle);
            I = Obj.TI;
            CT = Obj.turbine.rotorDat.CT;
            
            % Eq. 6.4 [2] Wind Velocity at the rotor
            uR_norm = (CT.*cos(flowAngle))./(2.*(1-sqrt(1-(CT.*cos(flowAngle)))));
            % Eq. 6.7 [2] Wind velocity from x = 0 to x = x0
            u0_norm = sqrt(1-CT);
            % p.529 par last Normalized veloctiy deficit at x = x0
            C0_norm = 1-u0_norm;
            
            % p.530 par.2 [2] 3-D Wake area in the near wake region
            ellipseLat = rotorD.*cos(flowAngle).*sqrt(uR_norm./u0_norm)./2;
            ellipseVert = rotorD.*sqrt(uR_norm./u0_norm)./2;
            
            % Eq 6.16 [2]
            x0_norm = ((cos(flowAngle).*(1+sqrt(1-CT)))./(sqrt(2).*(((Obj.vel_alpha*4).*I) + (Obj.vel_beta*2).*(1-sqrt(1-CT)))));
            x0_norm = repmat(x0_norm,length(z_norm),1);
            
            % Eq 6.10 [2]
            sigma_z0_norm = sqrt((1+sqrt(1-(CT.*cos(flowAngle))))./(8.*(1+sqrt(1-CT))));
            sigma_y0_norm = sigma_z0_norm.*cos(flowAngle);
            
            
            % Near/Far Wake filters
            nearWakeInd = x_norm <= x0_norm & x_norm >= 0;
            nearWakeLen = sum(nearWakeInd,2);
            farWakeInd = x_norm > x0_norm;
            
            % p.532 par.4 [2] Characteristic velocity of the shear layer
            us = (Ui+u0_norm*Ui)./2;
            % p.532 par.4 [2] Characteristic relative velocity in the shear layer
            ue = (Ui-u0_norm*Ui)./2;
            
            % Eq 6.15 [2] Variation of the shear layer width
            s = (Ui./us).*(Obj.vel_alpha.*I + Obj.vel_beta.*(ue./Ui)).*x_norm;
            s = repmat(s,size(farWakeInd,1),1);
            s(farWakeInd) = 0;
            
            % Eq 7 [3] - Wake expansion parameters
            ky = Obj.vel_ka*I + Obj.vel_kb;
            kz = Obj.vel_ka*I + Obj.vel_kb;
            
            % Eq 7.2 [2] - wake standard deviations
            sigma_y_norm = ky.*((x_norm-x0_norm))+ (cos(flowAngle))/sqrt(8);
            sigma_z_norm = kz.*((x_norm-x0_norm))+ 1/sqrt(8);
            
            % Eq 7.1 [2] - inner term of the equation
            inTerm = sqrt(1-((CT.*cos(flowAngle))./(8.*sigma_y_norm.*sigma_z_norm)));
            inTerm(imag(inTerm)~=0) = 0;
            
            theta_C0 = yawOffset + (((0.3*flowAngle)./cos(flowAngle)).*(1-sqrt(1-(CT.*cos(flowAngle)))))+infAngle;
            M0 = C0_norm.*(2-C0_norm);
            C = 1 - sqrt(1-(((sigma_y0_norm.*sigma_z0_norm).*M0)./(sigma_y_norm.*sigma_z_norm)));
            E0 = (C0_norm.^2) - (3.*exp(1/12).*C0_norm) + (3.*exp(1/3));
            E  = (C.^2) - (3.*exp(1/12).*C) + 3.*exp(1/3);
            theta_C = (theta_C0.*(sigma_y0_norm.*sigma_z0_norm).*E0)./((sigma_y_norm.*sigma_z_norm).*E);
            theta_M = theta_C./(exp(-0.5));
            
            Obj.wakeVelocity = zeros(length(x_norm),length(y_norm),length(z_norm));
            Obj.wakeSkewAngle = zeros(length(x_norm),length(y_norm),length(z_norm));
            for iVert = 1:length(z_norm)
                % p.532 par.3 [2] - Potential core radius
                rpc = zeros(size(x_norm));
                rpcY0 = sqrt((1-((z_norm(iVert).^2)./(ellipseVert^2))).*(ellipseLat^2));
                rpcY0(imag(rpcY0)>0) = 0;
                z = z_norm(iVert);
                % Eq 7.1 [2] - 3rd term of the equation
                exp2 = exp(-0.5.*(((z)./(sigma_z_norm(iVert,:).*rotorD)).^2));
                for jLat = 1:length(y_norm)
                    y = y_norm(jLat);
                    % p.532 par.3 [2] - Lateral distance from the wake centre
                    rY = (y+Yi-Obj.wakeCenterLoc);
                    % rpc at x = 0 to 0 at x = x0
                    rpc(nearWakeInd(iVert,:)) = linspace(rpcY0,0,nearWakeLen(iVert));
                    
                    % Eq 6.13 [2] - Wind velocity at near wake region
                    nearWake = 1-(C0_norm.*exp(-(((abs(rY)-rpc).^2)./(2.*(s(iVert,:).^2)))));
                    nearWake(abs(rY)<rpc) = 1-C0_norm;
                    nearSkew = (nearWake.*0);
                    
                    % Eq 7.1 [2] - 2nd term of the equation
                    exp1 = exp(-0.5.*(((rY)./(sigma_y_norm(iVert,:).*rotorD)).^2));
                    % Eq 7.1 [2] - Wind velocity at far wake region
                    farWake = 1-((1-inTerm(iVert,:)).*exp1.*exp2);
                    
                    exp1_1 = exp(-(((rY+(sign(flowAngle)*sigma_y_norm(iVert,:).*rotorD)).^2)./(2*((sigma_y_norm(iVert,:).*rotorD).^2))));
                    farSkew = exp1_1.*exp2.*theta_M(iVert,:);
                    % Combining both regions
                    Obj.wakeVelocity(nearWakeInd(iVert,:),jLat,iVert) = nearWake(nearWakeInd(iVert,:));
                    Obj.wakeVelocity(farWakeInd(iVert,:),jLat,iVert) = farWake(farWakeInd(iVert,:));
                    Obj.wakeSkewAngle(nearWakeInd(iVert,:),jLat,iVert) = nearSkew(nearWakeInd(iVert,:));
                    Obj.wakeSkewAngle(farWakeInd(iVert,:),jLat,iVert) = farSkew(farWakeInd(iVert,:));
                    noWakeInd = x_norm < 0;
                    Obj.wakeVelocity(noWakeInd,jLat,iVert) = 1;
                    Obj.wakeSkewAngle(noWakeInd,jLat,iVert) = 0;
                end
            end
            Obj.wakeVelocity = Obj.wakeVelocity.*Ui.*cos(infAngle);
            Obj.wakeSkewAngle = rad2deg(Obj.wakeSkewAngle);
        end
    end
end
