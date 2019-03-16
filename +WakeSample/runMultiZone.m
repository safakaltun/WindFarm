yawAngles = [-30];
fig = findobj('Tag','MultiZone');
if isempty(fig)
    fig = figure('Tag','MultiZone');
else
    clf(fig)
end
fig.Position = [315 147 816 740];
kdNum = [0.02 0.17 0.5]
for i = 1:3
clear Obj
MObj = MultiZoneWake;
MObj.kd = kdNum(i);
MObj.turbine = DTU_6MW;
MObj.turbine.pitchAngle = -2;
MObj.turbine.set_BEM_data;
MObj.turbine.yawAngle = yawAngles;
MObj.turbine.yawRotOffset = 0;
MObj.inflowVelocity  = 8;
MObj.turbine.bem_solver(MObj.inflowVelocity);
rotorD = MObj.turbine.bladeLength*2;
MObj.lateralDist = -2*rotorD:2*rotorD;
MObj.axialDist   = 0*rotorD:15*rotorD;
MObj.turbine.farmPosX = 0;
MObj.turbine.farmPosY = 0;
MObj.turbine.farmPosZ = 102;
MObj.deflection;
MObj.velocity;


ax = axes(fig);
subplot(3,1,i,ax)
surf(ax,MObj.axialDist./rotorD,MObj.lateralDist./rotorD,MObj.wakeVelocity','EdgeAlpha',0)
hold(ax,'on');
[turbx,turby] = MObj.turbine.get_actuator_disc_position(yawAngles,270);
plot3(ax,turbx,turby,[10 10],'-k','LineWidth',1)
plot3(ax,MObj.axialDist./rotorD,MObj.wakeCenterLoc./rotorD,10.*ones(size(MObj.wakeCenterLoc)),'--k','LineWidth',1)
view(ax,[0 90])
daspect(ax,[1 1 1])
xlabel(ax,'Axial Distance   [D]')
ylabel(ax,'Lateral Distance [D]')
title(ax,['Yaw Angle: ' num2str(MObj.turbine.yawAngle) '\circ'])
xlim(ax,[-1 15])
end
