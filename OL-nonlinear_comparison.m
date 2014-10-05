% compare nonlinear system and ol linear system

% simForce (TimeSeries)
% [ttl y_sim_ol] (1 column time and 3 columns vector)

figure(1)
subplot(3,1,1)
plot(ttl,y_sim_ol(:,1),simForce.Time,simForce.data(:,1))
legend ('OL','Nonlinear')
title 'F_x'
subplot(3,1,2)
plot(ttl,y_sim_ol(:,2),simForce.Time,simForce.data(:,2))
legend ('OL','Nonlinear')
title 'F_y'
subplot(3,1,3)
plot(ttl,y_sim_ol(:,3),simPhiB.Time,simPhiB.data)
legend ('OL','Nonlinear')
title '\phi_B'

