function draw_rocket(state, m, M, L)
% x = [delta; q; alpha] 

%     x = state(1); % delta 
    theta = state(3);  % theta == alpha when the velocity lies on the x axis. 

    % dimensions
%     W = 1*sqrt(M/5);  % cart width
%     H = .5*sqrt(M/5); % cart height
%     wr = .2;          % wheel radius
%     mr = .3*sqrt(m);  % mass radius

    lr = 1; % rocket length
    % 
    % positions
    %%
    displaylim = [-0.5 1.5 -0.8 0.8];
    % rocket 
    rtail = [0 0];
    rhead = lr * [cos(theta) sin(theta)];
    
    % delta: 2/3 rocket length from tail. 
    
    
%     y = wr/2+H/2; % cart vertical position
%     pendx = x + L*sin(th);
%     pendy = y - L*cos(th);

    plot(displaylim(1 : 2), [0 0], 'm', 'LineWidth', 1)
%     plot([0 0], displaylim(3 : 4), 'm', 'LineWidth', 1)
    hold on
%     rectangle('Position', [x-W/2,y-H/2,W,H], 'Curvature', .1, 'FaceColor', [.5 0.5 1], 'LineWidth', 1.5); % Draw rocket
    plot([rtail(1) rhead(1)], [rtail(2) rhead(2)], '-*k', 'LineWidth', 30)
    
%     rectangle('Position',[x-.9*W/2,0,wr,wr],'Curvature',1,'FaceColor',[0 0 0],'LineWidth',1.5); % Draw wheel
%     rectangle('Position',[x+.9*W/2-wr,0,wr,wr],'Curvature',1,'FaceColor',[0 0 0],'LineWidth',1.5); % Draw wheel
%     plot([x pendx],[y pendy],'k','LineWidth',2); % Draw pendulum
%     rectangle('Position',[pendx-mr/2,pendy-mr/2,mr,mr],'Curvature',1,'FaceColor',[1 0.1 .1],'LineWidth',1.5);

    axis(displaylim);
%     axis equal
%     set(gcf,'Position',[100 100 1000 400])
    drawnow
    hold off

end