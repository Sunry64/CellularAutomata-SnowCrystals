clear
clc
close all

iterations = 100;
beta = 0.95;             %value for background level
gamma = 0.0035;           %constant added to receptive sites
crystal_size = 101;     %number of hex on top of center  (note only use odd numbers)

tic                     %just time stuff
        
A = beta * ones(crystal_size*2 +1, crystal_size*2 + 1 );        %Array to store states
A(crystal_size+1, crystal_size+1) = 1;                          %setting the center to ice
x = -crystal_size:1:crystal_size;                               %array to store x positions of hexagons
y = crystal_size:-1:-crystal_size;                              %array to store y positions of hexagons
[X, Y] = meshgrid(x,y);

%add offset to alternate rows and columns
for i = 1:2:crystal_size*2 +1
    X(i,:) = X(i,:) + 0.5;
end

Dis = sqrt(X.^2 + Y.^2);

% drawing(A,X,Y,Dis,crystal_size);
% drawnow
for i = 1 : iterations
    
    [odd_row, even_row] = neighbor_array(A, beta);

    receptive_odd = sum((odd_row  >= 1) , 3)>0;
    receptive_even = sum((even_row  >= 1) , 3)>0;
    receptive_logic = zeros(size(A));
    for k = 1:size(A,1)
        if mod(k,2) == 0
            receptive_logic(k,:) = receptive_even(k,:);
        else
            receptive_logic(k,:) = receptive_odd(k,:);
        end
    end
    unreceptive_logic = receptive_logic==0;
    receptive = receptive_logic .*  A;
    unreceptive = unreceptive_logic .*  A;

    receptive_update = (receptive_logic * gamma) + receptive;
    [odd_row, even_row] = neighbor_array(unreceptive, beta);

    unreceptive_ne = zeros(size(odd_row));
    for k = 1:size(A,1)
        if mod(k,2) == 0
            unreceptive_ne(k,:,:) = even_row(k,:,:);
        else
            unreceptive_ne(k,:,:) = odd_row(k,:,:);
        end
    end

    averaged_neighbors = (unreceptive_ne(:,:,1)*0.5) + (unreceptive_ne(:,:,2)*(1/12)) + (unreceptive_ne(:,:,3)*(1/12)) ...
        + (unreceptive_ne(:,:,4)*(1/12)) + (unreceptive_ne(:,:,5)*(1/12)) + (unreceptive_ne(:,:,6)*(1/12)) + (unreceptive_ne(:,:,7)*(1/12));

    final_update = averaged_neighbors + receptive_update;
    
%     if mod(i,10)==0
%         clf
%         drawing(final_update,beta)
%         drawnow
%     end
    clf
    drawing(final_update,beta)
    drawnow

    A = final_update;
    
end
toc

tic

clf
drawing(A,beta)
drawnow

toc



% ---------------------------------------------------------------
% Functions
% ---------------------------------------------------------------



function drawing(A, beta)

low = A < 1;
low_array = (low .* A) .* 0.5;
high = (low == 0);

temp_high = 1 - abs((high .* A) - 1);
logic_temp_high = (temp_high <= 0.5) .* high;
high_array = ((temp_high > 0.5).*temp_high) +  (0.5*logic_temp_high);

A = low_array + high_array;
A = kron(kron(A,ones(2,1)), ones(1,2));

back = beta*ones(size(A,1) , size(A,2)+1);
for ii = 1:size(A,1)
    choice = ceil(ii/2);
    if (mod(choice,2)==1) 
        back(ii, 2:end) = A(ii,:);
    else
        back(ii, 1:end-1) = A(ii,:);
    end
end
A = back;

imagesc(A)
colormap('gray')

end







function corners = hexagon(center, radius)

radius = radius / cosd(30);
angles = [30 90 150 210 270 330];

x = cosd(angles)*radius + center(1);
y = sind(angles)*radius + center(2);
corners = [x;y];

end







function [odd_row, even_row] = neighbor_array(A, beta)

top_left = [beta*ones(1,size(A,2)+1) ; [beta*ones(size(A,1),1) A]];
top_left = top_left(1:end-1 , 1:end-1);

top = [beta*ones(1, size(A,2)); A];
top = top(1:end-1 , :);

top_right = [beta*ones(1,size(A,2)+1) ; [A beta*ones(size(A,1),1)]];
top_right = top_right(1:end-1 , 2:end);

left = [beta*ones(size(A,1), 1) A];
left = left(: , 1:end-1);

right = [A beta*ones(size(A,1), 1)];
right = right(: , 2:end);

bottom_left = [[beta*ones(size(A,1),1) A] ; beta*ones(1,size(A,2)+1)];
bottom_left = bottom_left(2:end , 1:end-1);

bottom = [A ; beta*ones(1, size(A,2))];
bottom = bottom(2:end , :);

bottom_right = [[A beta*ones(size(A,1),1)] ; beta*ones(1,size(A,2)+1)];
bottom_right = bottom_right(2:end , 2:end);

odd_row = cat(3, A, top, top_right, left, right, bottom, bottom_right);
even_row = cat(3, A, top_left, top, left, right, bottom_left, bottom);

end






