yawAngles =[ -25 -15 0];
saveDir = 'E:\Scientific\Dvl\Thesis\Data\Output\Studies\GaussianWake';
for i = 1:length(yawAngles)
    clear GObj
    GObj = GaussianWake;
    GObj.turbine = DTU_6MW;
    GObj.turbine.farmPosX = 0 ;
    GObj.turbine.farmPosY = 0;
    GObj.turbine.farmPosZ = 102;
    GObj.turbine.pitchAngle = -2;
    GObj.turbine.set_BEM_data;
    GObj.turbine.yawAngle = yawAngles(i);
    GObj.turbine.yawRotOffset = -0.0090;
    GObj.inflowVelocity  = 8;
    GObj.inflowAngle = -7;
    GObj.turbine.bem_solver(GObj.inflowVelocity);
    
    rotorD = GObj.turbine.bladeLength*2;
    GObj.lateralDist = -3*rotorD:rotorD/16:3*rotorD;
    GObj.axialDist   = (-3*rotorD:rotorD/16:20*rotorD)';
    GObj.verticalDist = sort([20:20:200 102]);
    GObj.deflection;
    GObj.velocity
    
    fig = figure;
    ax = axes(fig);
    dist = 102;
    ind = find(GObj.verticalDist>=dist,1,'first');
    surf(ax,GObj.axialDist./rotorD,GObj.lateralDist./rotorD,GObj.wakeVelocity(:,:,ind)','EdgeAlpha',0)
    view(ax,[0 90])
    daspect(ax,[1 1 1])
    xlabel(ax,'Axial Distance   [D]')
    ylabel(ax,'Lateral Distance [D]')
end
