classdef WindFarm < matlab.mixin.SetGet
    properties
        wDir
        wSpd
        TI = 0.05;
        
        xPos
        yPos
        zPos
        yawAngle
        
        turbine
        wake
        
        flowField
        skewField
        TIField
        flowXField
        flowYField
        flowZField
        interpu
        interpv
        interpdir
        interpTI
        
        wakeCenterLoc
        cLineX
        
    end
    properties  % These are for method comparison
        wDirArray
        wSpdArray
        TIArray
    end
    properties  %(GetAccess = public,SetAccess=private ) % These are for method comparison
        addedTurbMethod = 'Niayifar'; % <'Crespo','Niayifar','off'>
        applyTurb       = ''; % <'VelOnly'>
        wakeFlowAngle   = 'freestream'; % <'freestream','wake','wake_max','wake_mean_angles'>
        wakeCombine     = 'on'; % <'on','off'>
        sampleWsp       = 'wake'; % <'freestream','wake'>
        samplePosition  = 0; % Where the averaging is done relative to the turbine [x/D,Numeric]
        sampleAngle     = 'noYaw'; % <'followRotor','noYaw'>
        shearOn         = true;
    end
    methods
        function Obj = WindFarm(varargin)
            % Creates a wind farm layout and does farm related calculations
            % EXAMPLE
            %
            % WindFarm(xPos,yPos,zPos,yawAngle,WindTurbineClass)
            %   - xPos: Array of normalized x coordinates of the position of turbines
            %   - yPos: Array of normalized y coordinates of the position of turbines
            %   - zPos: Array of normalized z coordinates of the position
            %   of turbines [Hub height = 0];
            %   - yawAngle: Array of yaw angle of each wind turbine in deg
            %   - WindTurbineClass: Wind turbine class in the farm [singular]
            
            if nargin == 6
                Obj.init_farm(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
            end
        end
        function init_farm(Obj,xPos,yPos,zPos,yawAngle,WindTurbineClass,WakeMod)
            assert(length(xPos)==length(yPos),'Length of x and y coordinate arrays do not match!')
            assert(length(xPos)==length(zPos),'Length of x,y and z coordinate arrays do not match!')
            assert(length(xPos)==length(yawAngle),'Length of x,y,z and yawAngle arrays do not match!')
            assert(isa(WindTurbineClass,'WindTurbine'),'The Wind Turbine input is not recognized.')
            assert(isa(WakeMod,'WakeModel'),'The Wake Model input is not recognized.')
            
            [~,posInd] = sort(xPos);
            xPos = xPos(posInd);
            yPos = yPos(posInd);
            zPos = zPos(posInd);
            yawAngle = yawAngle(posInd);
            
            Obj.xPos = xPos;
            Obj.yPos = yPos;
            Obj.zPos = zPos;
            Obj.yawAngle = yawAngle;
            rotorD = WindTurbineClass.bladeLength*2;
            Obj.turbine = eval(class(WindTurbineClass));
            Obj.wake = eval(class(WakeMod));
            
            for iTurb = 1:length(xPos)
                Obj.turbine(iTurb) = eval(class(WindTurbineClass));
                Obj.turbine(iTurb).farmPosX = xPos(iTurb).*rotorD;
                Obj.turbine(iTurb).farmPosY = yPos(iTurb).*rotorD;
                Obj.turbine(iTurb).farmPosZ = zPos(iTurb).*rotorD;
                Obj.turbine(iTurb).yawAngle = yawAngle(iTurb);
                Obj.turbine(iTurb).yawRotOffset = -0.009;
            end
            
        end
        function init_wake(Obj,axialDist,lateralDist,verticalDist,CTDat)
            rotorD = Obj.turbine(1).bladeLength*2;
            if ~exist('axialDist','var')
                axialDist   = (min(Obj.xPos(:))-3)*rotorD:rotorD/16:(max(Obj.xPos(:))+20)*rotorD;
                axialDist = axialDist';
            elseif isempty(axialDist)
                axialDist   = (min(Obj.xPos(:))-3)*rotorD:rotorD/16:(max(Obj.xPos(:))+20)*rotorD;
                axialDist = axialDist';
            end
            if size(axialDist,1)<size(axialDist,2)
                axialDist= axialDist';
            end
            if ~exist('lateralDist','var')
                lateralDist = floor((min(Obj.yPos(:))-3)*rotorD):rotorD/16:ceil((max(Obj.yPos(:))+3)*rotorD);
                lateralDist = lateralDist';
            elseif isempty(lateralDist)
                lateralDist = floor((min(Obj.yPos(:))-3)*rotorD):rotorD/16:ceil((max(Obj.yPos(:))+3)*rotorD);
                lateralDist = lateralDist';
            end
            if size(lateralDist,1)<size(lateralDist,2)
                lateralDist = lateralDist';
            end
            if ~exist('verticalDist','var')
                verticalDist = sort([20:20:200 102]);
                verticalDist = verticalDist';
            elseif isempty(verticalDist)
                verticalDist = 1:200;
                verticalDist = verticalDist';
            end
            if size(verticalDist,1)<size(verticalDist,2)
                verticalDist = verticalDist';
            end
            
            velocity_u = (zeros(length(axialDist),length(lateralDist),length(verticalDist))+1);
            
            if Obj.shearOn
                if Obj.wSpd==6
                    shearWSP = 0.152913982036/0.4*log((verticalDist+1.557061e-05)/1.557061e-05)./Obj.wSpd;
                else
                    shearWSP = 0.203885309382/0.4*log((verticalDist+1.557061e-05)/1.557061e-05)./Obj.wSpd;
                end
                for iShear = 1:length(shearWSP)
                    velocity_u(:,:,iShear) = velocity_u(:,:,iShear).*shearWSP(iShear);
                end
            end
            
            skewAng = zeros(length(axialDist),length(lateralDist),length(verticalDist));
            TI = zeros(length(axialDist),length(lateralDist),length(verticalDist))+Obj.TI;
            [X,Y,Z] = ndgrid(axialDist,lateralDist,verticalDist);
            
            
            for iTurb = 1:length(Obj.turbine)
                F = griddedInterpolant(X,Y,Z,velocity_u);
                [~,turbY] = WindTurbine.get_actuator_disc_position(0,Obj.wDir);
                turbY = linspace(turbY(1),turbY(2),40).*Obj.turbine(iTurb).bladeLength;
                xCoord = Obj.turbine(iTurb).farmPosX+(rotorD.*Obj.samplePosition);
                yCoord = Obj.turbine(iTurb).farmPosY+turbY;
                
                [xgrid, ygrid] = meshgrid(yCoord',verticalDist);
                areaMask = ((xgrid-Obj.yPos(iTurb).*Obj.turbine(iTurb).bladeLength*2).^2 + (ygrid-Obj.zPos(iTurb).*Obj.turbine(iTurb).bladeLength*2).^2) <= Obj.turbine(iTurb).bladeLength.^2;
                areaMask = areaMask';
                if strcmpi(Obj.sampleWsp,'freestream')
                    Obj.wSpdArray(iTurb) = Obj.wSpd;
                elseif strcmpi(Obj.sampleWsp,'wake')
                    dataMat = permute(F({xCoord yCoord verticalDist}),[2 3 1]);
                    Obj.wSpdArray(iTurb) = mean(dataMat(areaMask));
                end
                
                %% Effective wind flow angle at the rotor plane
                F = griddedInterpolant(X,Y,Z,skewAng);
                
                if strcmpi(Obj.wakeFlowAngle,'freestream')
                    Obj.wDirArray(iTurb) = 0;
                elseif strcmpi(Obj.wakeFlowAngle,'wake')
                    dataMat = permute(F({xCoord yCoord verticalDist}),[2 3 1]);
                    Obj.wDirArray(iTurb) = mean(dataMat(areaMask));
                elseif strcmpi(Obj.wakeFlowAngle,'wake_max')
                    dataMat = permute(F({xCoord yCoord verticalDist}),[2 3 1]);
                    val = dataMat(areaMask);
                    [~,inVal] = max(abs(val));
                    Obj.wDirArray(iTurb) = max(val(inVal));
                elseif strcmpi(Obj.wakeFlowAngle,'wake_mean_angles')
                    dataMat = permute(F({xCoord yCoord verticalDist}),[2 3 1]);
                    val = dataMat(areaMask);
                    Obj.wDirArray(iTurb) = atan2d(mean(sind(val)),mean(cosd(val)));
                end
                
                %% Effective TI at the rotor plane
                F = griddedInterpolant(X,Y,Z,TI);
                dataMat = permute(F({xCoord yCoord verticalDist}),[2 3 1]);
                Obj.TIArray(iTurb) = mean(dataMat(areaMask));
                %% Turbine Data Assignment
                Obj.wake(iTurb).turbine = Obj.turbine(iTurb);
                % Initial pitch angle [subject to change]
                Obj.wake(iTurb).turbine.pitchAngle = -2;
                % Yaw angle of the individual turbine
                Obj.wake(iTurb).turbine.yawAngle = Obj.yawAngle(iTurb);
                % Yaw deflection offset parameter caused by the rotation [This is for Controls]
                Obj.wake(iTurb).turbine.yawRotOffset = -0.009;
                % BEM calculation for the individual turbine
                if ~exist('CTDat','var')
                    Obj.wake(iTurb).turbine.set_BEM_data;
                    Obj.wake(iTurb).turbine.bem_solver(Obj.wSpdArray(iTurb)*Obj.wSpd);
                else
                    Obj.wake(iTurb).turbine.rotorDat.CT = CTDat(iTurb);
                end
                
                
                %% Wake Model run
                % Initial wind speed
                Obj.wake(iTurb).inflowVelocity  = Obj.wSpdArray(iTurb);
                % Initial wind speed
                Obj.wake(iTurb).inflowAngle  = Obj.wDirArray(iTurb);
                % Initial TI
                Obj.wake(iTurb).TI  = Obj.TIArray(iTurb);
                % Affected area masks
                latCond = lateralDist<=(Obj.turbine(iTurb).farmPosY+(rotorD*4)) & lateralDist>=(Obj.turbine(iTurb).farmPosY+(rotorD*-4));
                % axCond = axialDist>=min(xCoord);
                axCond = axialDist>=Obj.turbine(iTurb).farmPosX;
                
                % Only send the relevant area to the model
                Obj.wake(iTurb).lateralDist = lateralDist(latCond);
                Obj.wake(iTurb).axialDist   = axialDist(axCond);
                Obj.wake(iTurb).verticalDist = verticalDist;
                % Wake deflection and velocity model
                if strcmpi(Obj.applyTurb,'VelOnly')
                    Obj.wake(iTurb).deflection(0.05);
                else
                    Obj.wake(iTurb).deflection;
                end
                Obj.wake(iTurb).velocity;
                intVel_u = velocity_u(axCond,latCond,:);
                intSkew = skewAng(axCond,latCond,:);
                int_TI = TI(axCond,latCond,:);
                
                if strcmpi(Obj.wakeCombine,'off')
                    intVel_u = (Obj.wake(iTurb).wakeVelocity./Obj.wSpdArray(iTurb)).*(intVel_u./intVel_u);
                elseif strcmpi(Obj.wakeCombine,'on')
                    intVel_u = (Obj.wake(iTurb).wakeVelocity./Obj.wSpdArray(iTurb)).*intVel_u;
                end
                
                intSkew = Obj.wake(iTurb).wakeSkewAngle;
                
                int_aind = (1./(2.*cosd(Obj.yawAngle(iTurb))))*(1-sqrt(1-(Obj.wake(iTurb).turbine.rotorDat.CT.*cosd(Obj.yawAngle(iTurb))./cosd(Obj.wDirArray(iTurb)))));
                if strcmpi(Obj.addedTurbMethod,'Crespo')
                    int_TI = (int_TI.*0) + Obj.TI + 0.73.*(int_aind.^0.8325).*(0.05.^0.0325).*(((axialDist(axCond)-Obj.turbine(iTurb).farmPosX)./rotorD).^-0.32);
                elseif strcmpi(Obj.addedTurbMethod,'Niayifar')
                    int_TI = sqrt((((int_TI.*0) + Obj.TI).^2) + (0.73.*(int_aind.^0.8325).*(0.05.^0.0325).*(((axialDist(axCond)-Obj.turbine(iTurb).farmPosX)./rotorD).^-0.32)).^2);
                elseif strcmpi(Obj.addedTurbMethod,'off')
                    int_TI = (int_TI.*0) + Obj.TI;
                end
                
                velocity_u(axCond,latCond,:) = intVel_u;
                skewAng(axCond,latCond,:) = intSkew;
                TI(axCond,latCond,:) = int_TI;
                
            end
            Obj.flowField = velocity_u.*Obj.wSpd;
            Obj.skewField = skewAng;
            Obj.TIField = TI;
            [Obj.flowXField,Obj.flowYField,Obj.flowZField]  = ndgrid(axialDist./rotorD,lateralDist./rotorD,(verticalDist-102)./rotorD);
            Obj.interpStuff;
            
        end
        function interpStuff(Obj)
            
            Obj.interpu   = griddedInterpolant(Obj.flowXField,Obj.flowYField,Obj.flowZField,Obj.flowField,'linear','none');
            Obj.interpdir = griddedInterpolant(Obj.flowXField,Obj.flowYField,Obj.flowZField,Obj.skewField,'linear','none');
            Obj.interpv   = griddedInterpolant(Obj.flowXField,Obj.flowYField,Obj.flowZField,Obj.flowField.*tand(Obj.skewField),'linear','none');
            Obj.interpTI  = griddedInterpolant(Obj.flowXField,Obj.flowYField,Obj.flowZField,Obj.TIField,'linear','none');
            
        end
        function identify_wake_center(Obj)
            if isempty(Obj.wakeCenterLoc)
                axialArray = (Obj.flowXField(Obj.flowXField(:,1,1)>=3,1,1));
                lateralArray = Obj.flowYField(1,abs(Obj.flowYField(1,:,1)')<4,1)';
                [X,Y] = meshgrid(axialArray,lateralArray);
                if size(Obj.flowXField,3)>1
                    var = Obj.interpu(X,Y,Y.*0);
                else
                    var = Obj.interpu(X,Y);
                end
                ind95above = var<=Obj.wSpd*0.95;
                for iX = 1:length(axialArray)
                    if sum(ind95above(:,iX))>=6
                        indStart = find(ind95above(:,iX),1,'first')-1;
                        indEnd = find(ind95above(:,iX),1,'last')+1;
                        
                        profi = var(indStart:indEnd,iX);
                        arrayFit = lateralArray(indStart:indEnd);
                        [~,minInd] = min(profi);
                        Obj.wakeCenterLoc(iX) = arrayFit(minInd);
                        Obj.cLineX(iX) = axialArray(iX);
                    end
                end
                Obj.wakeCenterLoc = Obj.wakeCenterLoc';
                Obj.cLineX = Obj.cLineX';
            end
        end
        function plot_farm_layout(Obj,ax)
            if ~exist('ax','var')
                fig = figure;
                ax = axes(fig);
            end
            hold(ax,'on')
            for iTurb = 1:length(Obj.turbine)
                [turbx,turby,turbz] = WindTurbine.get_actuator_disc_position3(Obj.yawAngle(iTurb),Obj.wDir);
                plot3(ax,Obj.xPos(iTurb)+turbx,Obj.yPos(iTurb)+turby,turbz+(102/DTU_6MW.bladeLength/2),'-k','LineWidth',2)
            end
            xlim(ax,[min(Obj.xPos(:))-3 max(Obj.xPos(:))+20])
            ylim(ax,[min(Obj.yPos(:))-3 max(Obj.yPos(:))+3])
            view(ax,[0 90])
            daspect(ax,[1 1 1])
        end
        
        
    end
end