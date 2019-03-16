classdef DTU_6MW < WindTurbine
    properties (Constant)
        bladeLength = 77.017;
        power       = 6E6;
        cutIn       = 4;
        ratedSpeed  = 12;
        cutOut      = 25;
        rotSpeed    = [0.521 1.1520]; 
        numBlade    = 3;
        hubHeight   = 102;
    end
    properties (Constant,Hidden)
        aeroFile  = fullfile(regexprep(mfilename('fullpath'),'Codebase.*',''),'\Data\Turbine\DTU_6MW\DTU10MW_Rwt.pro');
        bladeFile = fullfile(regexprep(mfilename('fullpath'),'Codebase.*',''),'\Data\Turbine\DTU_6MW\DTU_6MW_Blade.dat'); 
    end
    % Turbine Props in RANS
    properties 
        posX     
        posY     
        posZ    
    end
    methods
        function Obj = DTU_6MW
        end
    end
    methods %BEM Methods
        function Obj = set_BEM_data(Obj)
            rawBlade = importdata(Obj.bladeFile);
            Obj.bladeDat.radius    = rawBlade.data(:,1); % Element radius
            Obj.bladeDat.chord     = rawBlade.data(:,2); % Chord length
            Obj.bladeDat.beta      = rawBlade.data(:,3); % Twist
            Obj.bladeDat.toc = rawBlade.data(:,4);       % Thickness/Chord
            
            fileAero = fopen(Obj.aeroFile,'rt');
            rawAero  = textscan(fileAero,'%s','Delimiter','\n');
            fclose(fileAero);
            
            airfoilName = regexprep(strtrim(rawAero{:}(~cellfun(@isempty,regexp(rawAero{:},'FFA')))),'-','_');
            startAirfoil = find(~cellfun(@isempty,regexp(rawAero{:},'FFA')))+1;
            endAirfoil = [find(~cellfun(@isempty,regexp(rawAero{:},'FFA')),length(startAirfoil)-1,'last')-1; length(rawAero{:})];
            for iFoil = 1:length(airfoilName)
                dat = str2num(char(rawAero{:}(startAirfoil(iFoil):endAirfoil(iFoil))));
                Obj.aeroDat(iFoil).name = airfoilName{iFoil};
                Obj.aeroDat(iFoil).AoA = dat(:,1);
                Obj.aeroDat(iFoil).CL  = dat(:,2);
                Obj.aeroDat(iFoil).CD  = dat(:,3);
                Obj.aeroDat(iFoil).ToC = str2num(strrep(airfoilName{iFoil},'FFA_W3_',''))/10;
                clear dat
            end
        end
    end
    
end
