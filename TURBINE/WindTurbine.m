classdef WindTurbine < matlab.mixin.SetGet
    properties (Abstract=true,Constant)
        power
        numBlade
        bladeLength
        hubHeight
        cutIn
        cutOut
        ratedSpeed
        rotSpeed
    end
    properties (Abstract=true,Constant,Hidden)
        aeroFile
        bladeFile
    end
    % Blade
    properties
        aeroDat
        bladeDat
        rotorDat
    end
    % Farm
    properties
        farmPosX
        farmPosY
        farmPosZ
    end
    
    % Control
    properties
        yawAngle
        yawRotOffset
        tiltAngle
        pitchAngle
        %         omega
    end
    
    % Wake
    properties
        deflection
        isUpstream
    end
    
    % Ambient
    properties (Constant, Hidden)
        rho = 1.225;
        speedInt = 1;
    end
    
    methods
        function Obj = WindTurbine()
            Obj.aeroDat  = AeroData;
            Obj.bladeDat = BladeData;
        end
        function Obj = bem_solver(Obj,Vint)
            speedDiff = Obj.speedInt;
            if exist('Vint','var')
                Vo = Vint;
            else
                Vo = Obj.cutIn:speedDiff:Obj.cutOut;
            end

            B = Obj.numBlade;
            R = Obj.bladeLength;
            %             lambdaOpt = (Obj.rotSpeed(2)*R)./Obj.ratedSpeed;
            lambdaOpt = 7.7633;
            omega = (lambdaOpt.*Vo./R);
            omega(Vo>Obj.ratedSpeed) = Obj.rotSpeed(2);
            omega(omega<Obj.rotSpeed(1)) = Obj.rotSpeed(1);
            Obj.rotorDat.omega = omega;
           
            for iVo = 1:length(Vo)
                Obj.rotorDat.power(iVo) = 1E7;
                thetap = 0;
                for jR = 1:length(Obj.bladeDat.radius)-1
                    r    = Obj.bladeDat.radius(jR);
                    beta = Obj.bladeDat.beta(jR);
                    c    = Obj.bladeDat.chord(jR);
                    ToC  = Obj.bladeDat.toc;
                    Cd   = interp1([Obj.aeroDat.ToC]',[Obj.aeroDat.CD]',ToC(jR))';
                    Cl   = interp1([Obj.aeroDat.ToC]',[Obj.aeroDat.CL]',ToC(jR))';
                    AoA  = Obj.aeroDat(1).AoA;
                    
                    % Velocity component normal to the rotor
                    Vn = Vo(iVo);
                    % Velocity component in the rotational direction
                    Vi = r*omega(iVo);
                    
                    a       = 0;
                    aprime  = 0;
                    change  = 1;
                    while abs(change)>1E-3
                        % Flowangle
                        phi = atand(((1-a)*Vn)/((1+aprime)*Vi));
                        % Angle of the section
                        theta = beta + thetap;
                        % Effective AoA
                        alphanew = phi - theta;
                        % Load and drag coefficient extraction from the
                        % airfoil data
                        Cdnew = interp1(AoA,Cd,alphanew);
                        Clnew = interp1(AoA,Cl,alphanew);
                        % Normal coefficient
                        Cn = Clnew*cosd(phi)+Cdnew*sind(phi);
                        % Tangential coefficient
                        Ct = Clnew*sind(phi)-Cdnew*cosd(phi);
                        % Prandtl's tip loss correction factor
                        F = (2/pi)*(acos(exp(-((B/2)*((R-r)/(r*sind(phi)))))));
                        % Solidity
                        sigma = (c*B)/(2*pi*r);
                        % GLAUERT Correction
                        if a>1/3
                            % For simplification of induction the calculation
                            f1 = 0.25*(5-3*a);
                            % Thrust Coefficient
                            CT = ((1-a)^2*sigma*Cn)/(sind(phi)^2);
                            %Axial induction factor
                            anew = a-((4*a*(1-f1*a)-(CT/F))/((9*(a^2))-(10*a)+4));
                        else
                            %Axial induction factor
                            anew =  1/(((4*F*(sind(phi)^2))/(sigma*Cn))+1);
                        end
                        % Convergence check
                        change = (anew-a);
                        % Relaxation
                        anew = (0.1*anew) + ((1-0.1)*a);
                        % Tangential induction factor
                        aprimenew = 1/(((4*F*sind(phi)*cosd(phi))/(sigma*Ct))-1);
                        % Relaxation
                        aprimenew = (0.1*aprimenew) + ((1-0.1)*aprime);
                        a = anew;
                        aprime = aprimenew;
                    end
                    Obj.rotorDat.Vrelative(iVo,jR)=sqrt(((Vn*(1-a))^2)+((Vi*(1+aprime))^2));
                    Obj.rotorDat.axialInd(iVo,jR) = a;
                    %Relative wind speed
                    Obj.rotorDat.Pt(iVo,jR) = 0.5*Obj.rho*Obj.rotorDat.Vrelative(iVo,jR)^2*c*Ct; %Tangential Load
                    Obj.rotorDat.Pn(iVo,jR) = 0.5*Obj.rho*Obj.rotorDat.Vrelative(iVo,jR)^2*c*Cn; %Normal Load
                end
                Obj.rotorDat.Pt(iVo,jR+1) = 0;
                Obj.rotorDat.Pn(iVo,jR+1) = 0;
                
                % Aerodynamic torque
                Obj.rotorDat.torque(iVo) = B*(trapz(Obj.bladeDat.radius,Obj.bladeDat.radius.*Obj.rotorDat.Pt(iVo,:)'));
                % Aerodynamic power
                Obj.rotorDat.power(iVo) = omega(iVo)*Obj.rotorDat.torque(iVo);
                %                 end
                % Aerodynamic thrust
                Obj.rotorDat.thrust(iVo) = B*(trapz(Obj.bladeDat.radius,Obj.rotorDat.Pn(iVo,:)));
                % Power coefficient
                Obj.rotorDat.Cp(iVo) = Obj.rotorDat.power(iVo)/(0.5*Obj.rho*(Vn^3)*pi*(R^2));
                % Thrust coefficient
                Obj.rotorDat.CT(iVo) = Obj.rotorDat.thrust(iVo)/(0.5*Obj.rho*(Vn^2)*pi*(R^2));
                Obj.rotorDat.thetap(iVo) = thetap;
                Obj.rotorDat.lambda(iVo) = omega(iVo)*R/Vn;
            end
        end
        function plot_turbine_characteristics(Obj)
            Vo = Obj.cutIn:Obj.speedInt:Obj.cutOut;
            fig = figure;
            fig.Position = [286 321 635 575];
            ax1 = axes(fig);
            ax2 = axes(fig);
            ax3 = axes(fig);
            subplot(3,1,1,ax1);
            subplot(3,1,2,ax2);
            subplot(3,1,3,ax3);
            plot(ax1,Vo,Obj.rotorDat.lambda)
            ylabel(ax1,'Tip Speed Ratio (\lambda) [-]')
            plot(ax2,Vo,Obj.rotorDat.omega)
            ylabel(ax2,'Rotational Speed (\omega) [rad/s]')
            plot(ax3,Vo,Obj.rotorDat.thetap)
            ylabel(ax3,'Pitch Angle (\theta_{p}) [rad/s]')
            
            xlabel(ax3,'Wind Speed [m/s]')
            
            fig = figure;
            fig.Position = [858 291 699 495];
            ax1 = axes(fig);
            ax2 = axes(fig);
            ax3 = axes(fig);
            ax4 = axes(fig);
            subplot(2,2,1,ax1);
            subplot(2,2,2,ax2);
            subplot(2,2,3,ax3);
            subplot(2,2,4,ax4);
            plot(ax1,Vo,Obj.rotorDat.power./1E3)
            ylabel(ax1,'Aerodynamic Power [kW]')
            
            plot(ax2,Vo,Obj.rotorDat.Cp)
            ylabel(ax2,'Power Coefficient [-]')
            plot(ax3,Vo,Obj.rotorDat.thrust./1E3)
            ylabel(ax3,'Aerodynamic Thrust [kN]')
            xlabel(ax3,'Wind Speed [m/s]')
            plot(ax4,Vo,Obj.rotorDat.CT)
            ylabel(ax4,'Thrust Coefficient [-]')
            xlabel(ax4,'Wind Speed [m/s]')
            
            
        end
    end
    methods (Abstract = true)
        set_BEM_data(Obj)
    end
    methods (Static)
        % Calculates the tip positions of the actuator disc [Normalized]
        function varargout = get_actuator_disc_position(varargin)
            nargoutchk(0,3)
            meanDir = varargin{2};
            yawAngle = -50:50;
            tipPos = [cosd(meanDir+yawAngle);-sind(meanDir+yawAngle)];
            tipPos = [tipPos; -tipPos]';
            container = containers.Map(yawAngle,num2cell(tipPos,2));
            if nargin
                actDisc = container(varargin{1});
                xLocs = actDisc([1 3]);
                yLocs = actDisc([2 4]);
            else
                actDisc = container;
            end
            if nargout == 1
                varargout{1} = actDisc;
            elseif nargout == 2 && nargin
                varargout{1} = xLocs;
                varargout{2} = yLocs;
            elseif nargout == 3 && nargin
                varargout{1} = xLocs;
                varargout{2} = yLocs;
                varargout{3} = ones(size(yLocs));
            else
                disp(actDisc);
            end
        end
        % Calculates the tip positions of the actuator disc [Normalized]
        function [xLocs,yLocs,zLocs] = get_actuator_disc_position3(varargin)
            x = 0;
            y = 0;
            z = 0;
            r = 0.5;
            yawAngle = varargin{1};
            meanDir = varargin{2};
            
            th = 0:pi/50:2*pi;
            xLocs = cosd(meanDir+yawAngle)*(r * cos(th) + x);
            yLocs = -sind(meanDir+yawAngle)*(r * cos(th) + y);
            zLocs = r * sin(th) + z;
        end
        function scaleCoefs = spanwise_load_scaler(distArray)
            if size(distArray,2)>size(distArray,1)
                distArray = distArray';
            end
            turbine = DTU_6MW;
            turbine.set_BEM_data;
            turbine.bem_solver(8);
            steadyDat = turbine.rotorDat;
            
            bladeSpan = [-flipud(turbine.bladeDat.radius);0;turbine.bladeDat.radius]./(DTU_6MW.bladeLength*2);
            bladeLoads = [flipud(steadyDat.Pn');0;steadyDat.Pn'];
            interpLoads = interp1(bladeSpan,bladeLoads,distArray);
            scaleCoefs = interpLoads./sum(interpLoads);
            
        end
    end
end
