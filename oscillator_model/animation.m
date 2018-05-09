function animation(x, t, slip)
% get number of masses
dim = size(x);
n = dim(2)/2;
% we`re onlz interested in the position values
x = x(:,1:n);

fps =30;
tanim = t(1):1/fps:t(end); % sampling points for the animation
xanim = interp1(t, x, tanim);
slipanim = interp1(t, slip, tanim);

% offset between the masses
offset = 0.5;
% size of 1 mass
a = 0.2;

slipdim = size(slipanim);
for i = 1:slipdim(1)
    for j = 1:slipdim(2)
        if slipanim(i,j) == 0
            color(i,j) = 'k';
        else
            color(i,j) = 'r';
        end
    end
end

massanim = zeros(n, 4);
massanim(:,2) = 0;  % y coordinates
massanim(:,3:4) = a; % fix side length
% set x coordinates
for i = 1:n
    massanim(i,1) = xanim(1,i) + (n-i)*offset - a/2;
end

% initialize figure
figure
axesHandle = gca;
xlim(axesHandle, [min(xanim(:,end)-a), max(xanim(:,1))+(n-1)*offset+a]);
ylim(axesHandle, [-3, 3] );
%axis equal;
for i = 1:n
    massobject(i) = rectangle('Position', massanim(i,:),...
                    'FaceColor', color(1,i), 'EdgeColor', color(1,i));
end
time = text(0,-2.5,num2str(tanim(1)));

%animation
for j = 1:length(tanim)
    for i = 1:n
    massanim(i,1) = xanim(j,i) + (n-i)*offset - a/2;
    set(massobject(i), 'Position', massanim(i,:), 'FaceColor', color(j,i),...
    'EdgeColor', color(j,i));
    end
    set(time,'String', num2str(tanim(j)));
    drawnow;
    pause(1/fps);
end


    

end