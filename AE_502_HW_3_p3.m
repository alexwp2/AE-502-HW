%% AE 502 HW 3 problem 3
e = 0.5;
a = 1;
i = 45 * pi/180;
w = [0.02, 0.1, 0.5];
omega = 10 * pi/180;

t = 0:1:100;

for j = 1:length(w)
    for m = 1:length(t)
        h(j,m) = e*sin(w(j)*t(m) + omega);
        k(j,m) = e*cos(w(j)*t(m) + omega);
        p(j,m) = tan(i/2)*sin(omega);
        q(j,m) = tan(i/2)*cos(omega);
    end
end


figure (1)
hold on;
axis equal;
plot(h(1,:),k(1,:));
% plot(h(1,:).*.8,k(1,:).*.8);
% plot(h(2,:).*.9,k(2,:).*.9);
% plot(h(3,:),k(3,:));
title('h vs k (w = 0.02)');
xlabel('h');
ylabel('k');

figure (2)
hold on;
axis equal;
plot(h(2,:),k(2,:));
title('h vs k (w = 0.1)');
xlabel('h');
ylabel('k');

figure (3)
hold on;
axis equal;
plot(h(3,:),k(3,:));
title('h vs k (w = 0.5)');
xlabel('h');
ylabel('k');

% figure (4)
% hold on;
% plot(p(1,:),q(1,:));
% plot(p(2,:),q(2,:));
% plot(p(3,:),q(3,:));




