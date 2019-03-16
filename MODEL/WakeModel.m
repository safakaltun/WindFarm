classdef WakeModel < matlab.mixin.SetGet
    properties
        turbine
        inflowVelocity
        inflowAngle
        wakeCenterLoc
        wakeVelocity
        wakeSkewAngle
    end
    properties (SetObservable)
        axialDist
        lateralDist
        verticalDist
        axialDist_norm
        lateralDist_norm
        verticalDist_norm
    end
    methods
        function Obj = WakeModel()
            addlistener(Obj,'axialDist','PostSet',@Obj.standardize_dimensions);
            addlistener(Obj,'lateralDist','PostSet',@Obj.standardize_dimensions);
            addlistener(Obj,'verticalDist','PostSet',@Obj.standardize_dimensions);
            addlistener(Obj,'axialDist_norm','PostSet',@Obj.standardize_dimensions);
            addlistener(Obj,'lateralDist_norm','PostSet',@Obj.standardize_dimensions);
            addlistener(Obj,'verticalDist_norm','PostSet',@Obj.standardize_dimensions);
        end
    end
    methods (Abstract=true)
        deflection(Obj)
        velocity(Obj)
    end
    methods
        function rotate_2d(Obj)
            cornerPts = [Obj.axialDist(1) Obj.lateralDist(1);...
                Obj.axialDist(1) Obj.lateralDist(end);...
                Obj.axialDist(end) Obj.lateralDist(1);...
                Obj.axialDist(end) Obj.lateralDist(end);...
                ];
            R = [cosd(Obj.inflowAngle) -sind(Obj.inflowAngle); sind(Obj.inflowAngle) cosd(Obj.inflowAngle)];
            rotpoint = R*cornerPts';
            min(rotpoint(1,:))
            max(rotpoint(1,:))
            min(rotpoint(2,:))
            max(rotpoint(2,:))
        end
    end
    methods (Access = protected)
        function Obj = standardize_dimensions(Obj,varargin)
            propName = varargin{1}.Name;
            if size(Obj.(propName),1) > size(Obj.(propName),2)
                Obj.(propName) = Obj.(propName)';
            end
            
            if strcmpi(propName(end-3:end),'norm')
                Obj.(propName(1:end-5)) =  Obj.(propName).*Obj.turbine.bladeLength.*2;
            else
                Obj.([propName '_norm']) = Obj.(propName)./(Obj.turbine.bladeLength.*2);
            end
        end
    end
end