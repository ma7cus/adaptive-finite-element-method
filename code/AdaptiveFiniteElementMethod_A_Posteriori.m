clc
close all
%% Setup for -u''(x) +mu'(x)+nu(x)= f(x)
%{
m = 100;
n = 250;

a = 0; %u(a) = A
b = 5; %u(b) = B

f = @(x) -6*x+3*m*x.^2+n*x.^3;

True_solution = @(x) x.^3;
dTrue = @(x) 3*x.^2;

A = True_solution(a);
B = True_solution(b);
%}


m = 10;
n = -20;


a = -0.05; %u(a) = A
b = 0.05; %u(b) = B

k=100; %Controls steepness

f = @(x) 2*k^2*tanh(k*x).*sech(k*x).^2+m*k*(sech(k*x)).^2+n*tanh(k*x);
True_solution = @(x) tanh(k*x);
dTrue = @(x) k*(sech(k*x)).^2;

A = True_solution(a);
B = True_solution(b);


%{
m = 10;
n = 20;

a = 0; %u(a) = A
b = 5; %u(b) = B

f = @(x) (1+n)*cos(x)-m*sin(x);

True_solution = @(x) cos(x);

A = True_solution(a);
B = True_solution(b);
%}
%% Initialising 
iteration = 1;
N = 5;
pts = linspace(a,b,N+1); %Set of evenly distributed xi points
starting_order = 1;
orders = starting_order*ones(1,N);

percentage_refined = 80; %The elements with top percentage_refined percent of worst errors will be refined

Global_TOL = 10^-10;
%% Starting iteration
[u,priori_error_vect_L2,priori_error_vect_H1,posteriori_error_vect] = IterateFiniteElementMethod(pts,orders,m,n,f,a,b,A,B,True_solution,dTrue,2);

global_error_L2 = norm(priori_error_vect_L2);
global_error_H1 = norm(priori_error_vect_H1);
global_error_posteriori= norm(posteriori_error_vect);
fprintf("Iteration %d \n",iteration)
fprintf("Global error L2: %e \n",global_error_L2)
fprintf("Global error H1: %e \n",global_error_H1)
fprintf("Global predicted error: %e \n",global_error_posteriori)
fprintf("Degrees of freedom: %d \n",length(u))
fprintf("\n")

%% Iterating
while global_error_L2 >= Global_TOL

N = length(pts)-1;
Percentile = prctile(posteriori_error_vect,percentage_refined);
% 
%% Flagging errors
flags = zeros(1,N);
pts_old = pts;
temp_pts = pts;
temp_orders = orders;
count = 0;
for i = 1:N
    if posteriori_error_vect(i) >= Percentile
        flags(i) = 1;
        new_interval = linspace(pts_old(i),pts_old(i+1),3); %Generates a new interval
        new_entry = new_interval(2); %Takes off the boundaries to avoid duplicates
        temp_pts = [temp_pts(1:i+count),new_entry,pts_old(i+1:length(pts))]; %Puts the new interval in the vector
        temp_orders = [temp_orders(1:i+count),orders(i),orders(i+1:length(orders))]; %Updates the list of orders for the new intervals
        count = count + 1;        
    end
end

flag_positions = find(flags);

N_temp = length(temp_pts)-1;

[~,~,~,posteriori_error_vect_interval] = IterateFiniteElementMethod(temp_pts,temp_orders,m,n,f,a,b,A,B,True_solution,dTrue,0);

e_vals_intervals=  zeros(1,N);

count = 0;
for i = 1:N
    if ismember(i,flag_positions)==1
        e_vals_intervals(i) = sqrt(posteriori_error_vect_interval(i+count)^2+posteriori_error_vect_interval(i+count+1)^2);
        count = count+1;
    else
        e_vals_intervals(i) = posteriori_error_vect_interval(i+count);
    end
end

%fprintf("r_vals_intervals: \n")
%disp(r_vals_intervals)

%% Updating orders
temp_orders = orders;
for i = 1:N
    if posteriori_error_vect(i) >= Percentile        
        temp_orders(i) = temp_orders(i)+1;
    end
end

[~,~,~,posteriori_error_vect_order] = IterateFiniteElementMethod(pts,temp_orders,m,n,f,a,b,A,B,True_solution,dTrue,0);

%fprintf("r_vals_orders: \n")
%disp(r_vals_orders)
%% Update
count = 0;
for i = 1:N
    if posteriori_error_vect(i) >= Percentile
        if  e_vals_intervals(i)<= posteriori_error_vect(i)&& e_vals_intervals(i)<=posteriori_error_vect_order(i)
            new_interval = linspace(pts(i+count),pts(i+count+1),3); %Generates a new interval
            new_entry = new_interval(2); %Takes off the boundaries to avoid duplicates
            pts = [pts(1:i+count),new_entry,pts(i+1+count:length(pts))]; %Puts the new interval in the vector
            orders = [orders(1:i+count),orders(i+count),orders(i+count+1:length(orders))]; %Updates the list of orders for the new intervals
            count = count + 1;  
        elseif posteriori_error_vect_order(i)<= posteriori_error_vect(i) && posteriori_error_vect_order(i)<=e_vals_intervals(i)
            orders(i+count) = orders(i+count)+1;
        end
    end
end

%% Run new model
[u,priori_error_vect_L2,priori_error_vect_H1,posteriori_error_vect] = IterateFiniteElementMethod(pts,orders,m,n,f,a,b,A,B,True_solution,dTrue,2);

iteration = iteration+1;

global_error_L2 = norm(priori_error_vect_L2);
global_error_H1 = norm(priori_error_vect_H1);
global_error_posteriori= norm(posteriori_error_vect);
fprintf("Iteration %d \n",iteration)
fprintf("Global error L2: %e \n",global_error_L2)
fprintf("Global error H1: %e \n",global_error_H1)
fprintf("Global predicted error: %e \n",global_error_posteriori)
fprintf("Degrees of freedom: %d \n",length(u))
fprintf("\n")

%fprintf("Using an order element pattern: ")
%disp(orders)
%fprintf("On points: ")
%disp(pts)
%fprintf("")
end