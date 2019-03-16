classdef MultiZoneWake < WakeModel
    % Model Parameters
    properties %Deflection
        kd = 0.1410;
    end
    properties %Expansion
        ke  = 0.065;
        me = [-0.5 0.22 1];
    end
    properties %Velocity
        MU = [0.5 1 5.5];
        aU  = 5;
        bU  = 1.66;
    end
    methods
        function Obj = MultiZoneWake()
        end
        function Obj = deflection(Obj)
            x        = Obj.axialDist;
            Xi       = Obj.turbine.farmPosX;
            rotorD   = Obj.turbine.bladeLength*2;
            yawAngle = -Obj.turbine.yawAngle;
            yawOffset = Obj.turbine.yawRotOffset;
            % Jimenez Wake deflection model in combination with Multi-Zone
            % Wake model Eq. 20
            initSkewAngle       = yawOffset + ((cosd(yawAngle)^2)*sind(yawAngle)*Obj.turbine.rotorDat.CT)/2;
            % Eq. 21
            
            % Gebraad 2016 Eq. 10
            eqPart0 = 15.*((((2.*Obj.kd.*(x-Xi))./(rotorD))+1).^4);
            eqPart1 = initSkewAngle.*(eqPart0+(initSkewAngle.^2));
            eqPart2 = ((30.*Obj.kd)./rotorD).*((((2.*Obj.kd.*(x-Xi))./(rotorD))+1).^5);
            eqPart3 = (initSkewAngle.*rotorD.*(15+(initSkewAngle.^2)))./(30.*Obj.kd);
            % Deflection caused by the yaw misalignment Eq. 10
            deflectionYaw = (eqPart1./eqPart2)-eqPart3;
            Obj.wakeCenterLoc = (Obj.turbine.farmPosY + deflectionYaw)./rotorD;
        end
        function Obj = velocity(Obj)
            % Gebraad P., 2016
            Xi = Obj.turbine.farmPosX;
            Yi = Obj.turbine.farmPosY;
            xArray = Obj.axialDist;
            yArray  = Obj.lateralDist;
            yw = Obj.wakeCenterLoc-Yi;
            yawAngle = deg2rad(-Obj.turbine.yawAngle);
            rotorD = Obj.turbine.bladeLength*2;
            Ui = Obj.inflowVelocity;
            
            ci = zeros(length(xArray),length(yArray),length(yawAngle));
            for iYaw = 1:length(yawAngle)
                % Eq. 17
                mu = Obj.MU./(cosd(Obj.aU+(Obj.bU.*yawAngle(iYaw))));
                for jxDist = 1:length(xArray)
                    x = xArray(jxDist);
                    for kyDist = 1:length(yArray)
                        y = yArray(kyDist);
                        % Eq. 13
                        Dw = max((rotorD+(2.*Obj.ke.*Obj.me.*(x-Xi))),0);
                        % Eq. 15
                        r = abs(y-(Yi + yw(jxDist)));
                        % Eq. 17
                        nearWake = r <= Dw(1)./2;
                        farWake  = (Dw(1)./2 <= r) & (r <= Dw(2)./2);
                        mixWake  = (Dw(2)./2 <= r) & (r <= Dw(3)./2);
                        noWake   = r >= Dw(3)./2;
                        
                        if nearWake
                            wakeCase = 1;
                        elseif farWake
                            wakeCase = 2;
                        elseif mixWake
                            wakeCase = 3;
                        elseif noWake
                            wakeCase = 0;
                        else
                            error('Something went wrong in Wake Velocity calculation.')
                        end
                        
                        if wakeCase ~= 0 && x >= -(y-Yi)*tan(yawAngle)+Xi
                            % Eq. 16
                            ci(jxDist,kyDist,iYaw) = (rotorD./(rotorD+(2.*Obj.ke.*mu(wakeCase).*(x-Xi)))).^2;
                        end
                    end
                end
            end
            % Eq. 14
            f = @(a) ((4*a*(1-a)) - Obj.turbine.rotorDat.CT);
            axialInd = fzero(f,0);
            Obj.wakeVelocity = Ui.*(1-(2.*axialInd.*ci));
        end
    end
end