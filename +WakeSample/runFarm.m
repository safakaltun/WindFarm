xPos = ([934.408191 2082.284962 3224.169811  934.408191 2082.284962 3224.169811]-934.408191)./154.03;
yPos = ([1292.412698 1292.412698 1292.412698 2220.513759 2250.165304 2250.165304]-1292.412698)./154.03;
zPos = [102 102 102 102 102 102]./154.03;
yawAngle = [-25 0 0 0 0 0];

range = 1:6;
xPos = xPos(range);
yPos = yPos(range);
zPos = zPos(range);
yawAngle = yawAngle(range);

Farm =  WindFarm(xPos,yPos,zPos,yawAngle,DTU_6MW,GaussianWake);
Farm.wDir = 270;
Farm.wSpd = 8;
Farm.TI = 0.05;
Farm.init_wake

Farm.plot_farm_layout;
surf(Farm.flowXField(:,:,1),Farm.flowYField(:,:,1),Farm.interpu({unique(Farm.flowXField)' unique(Farm.flowYField)' 0}),'EdgeAlpha',0)

ylabel('Lateral Distance [x/D]')
xlabel('Streamwise Distance [x/D]')
title([regexprep(class(Farm.wake),'Wake',' Wake') ' Model'])
cb = colorbar;
cb.Label.String = 'Streamwise Velocity [m/s]';
cb.Label.FontSize = 11;
set(gcf,'Position',[192 237 1226 583])
for i = 1:length(xPos)
    plot3(Farm.wake(i).axialDist_norm,Farm.wake(i).wakeCenterLoc,10+zeros(1,length(Farm.wake(i).wakeCenterLoc)),':k','LineWidth',1)
end
caxis([3 8])
