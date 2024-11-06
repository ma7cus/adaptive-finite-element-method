function [u,priori_error_vect_L2,priori_error_vect_H1,posteriori_error_vect] = IterateFiniteElementMethod(nodes,orders,m,n,f,a,b,A,B,True_solution,True_solution_derivative,Output)
    %IterateFiniteElementMethod evaluates the finite element method on the
    %PDE -u''(x)+m*u'(x)+n*u(x) = f(x) (u(a)=A,u(b) = B) over the elements
    %with boundaries 'nodes' at set of orders 'orders'.
    %[u,real_error,posteriori_error] = PerformFiniteElementMethod(pts,Order,m,n,f,a,b,A,B,True_solution,Output)
    %u is the vector of coefficients for the basis functions to build up
    %the approximation and real_error and posteriori_error are the vectors of real and approximated errors respectively of this
    %approximation over the intervals.

    %True solution is the correct solution to the PDE, used as a comparator
    %in the graphs. Output controls which graphs are outputted:
    % Output = 3 -> Approximation-real solution comparison graph with bases
    % and errors shown as well as matrix and vector of the system set up
    % Output = 2 -> Approximation-real solution comparison graph with bases and errors shown
    % Output = 1 -> Approximation-real solution comparison graph
    % Output = 0 -> No graphs outputted



       %% Calculating the set of bases on the interval [0,1] for the current order
    basis = cell(1,3); %Establishes the space for the bases (cells allows annonymous function entries)

    basis{1,1} = @(x) 1-x; %Linear downslope
    basis{2,1} = @(x) -1*ones(size(x)); %Derivative of the linear downslope
    basis{3,1} = @(x) zeros(size(x)); %Second derivative of the linear downslope

    basis{1,2} = @(x) x; %Linear upslope
    basis{2,2} = @(x) ones(size(x)); %Derivative of the linear upslope
    basis{3,2} = @(x) zeros(size(x)); %Second derivative of the linea upslope

    max_order = max(orders);

    N = length(nodes)-1; %This is the number of elements/intervals between the points we're evaluating

    %% General basis functions 0th derivative
    if max_order ~=1
        for pow = 2:(max_order)
            %Constructing the functions of higher order than linear
            funct = @(x) (-1)^(pow+1);
            for k = 0:(pow-1)
                funct = @(x) funct(x).*(x-k/(pow-1));
            end

            %Appending them to the basis matrix
            basis{1,pow+1} = @(x) funct(x); %This is the nth order basis function
        end
    end

    %% General basis functions 1st derivative
    if max_order ~=1
        for pow = 2:(max_order)
            %Constructing the functions of higher order than linear
            funct = @(x) zeros(size(x));
            for j = 0:(pow-1)
                funct_temp = @(x) (-1)^(pow+1)*ones(size(x));
                for k = 0:(pow-1)
                    if k ~=j
                        funct_temp = @(x) funct_temp(x).*(x-k/(pow-1));
                    end
                end
                funct = @(x) funct(x)+funct_temp(x);
            end

            %Appending them to the basis matrix
            basis{2,pow+1} = @(x) funct(x); %This is the nth order basis function 1st derivative
        end
    end

    %% General basis functions 2nd derivative
    if max_order ~=1
        for pow = 2:(max_order)
             %Constructing the functions of higher order than linear
            funct = @(x) zeros(size(x));
            for j = 0:(pow-1) %This runs through each term in the basis and omitts this in an evaluation
                for k = 0:(pow-1)
                    funct_temp = @(x) (-1)^(pow+1)*ones(size(x));
                    if k ~=j
                        for l = 0:(pow-1)
                            if l ~= k && l~=j
                                funct_temp = @(x) funct_temp(x).*(x-l/(pow-1));
                            end
                        end
                        funct = @(x) funct(x)+funct_temp(x);
                    end
                end
            end

            %Appending them to the basis matrix
            basis{3,pow+1} = @(x) funct(x); %This is the nth order basis function 2nd derivative
        end
    end

    %% Calculating the set of integrals of the bases and their derivatives on the interval [0,1]
    integral_matrix = zeros(max_order+1,max_order+1,3); %Establishes a space for the integrals of these basis functions

    for k = 1:(max_order+1)
        for j = 1:(max_order+1)
            integral_matrix(k,j,1) = integral(@(x) (basis{1,k}(x)).*(basis{1,j}(x)),0,1); %This is the integral of 0th derivatives of ith and jth basis functions
            integral_matrix(k,j,2) = integral(@(x) (basis{1,k}(x)).*(basis{2,j}(x)),0,1); %This is the integral of ith and (1st derivatives of jth) basis functions
            integral_matrix(k,j,3) = integral(@(x) (basis{2,k}(x)).*(basis{2,j}(x)),0,1); %This is the integral of 1st derivatives of (ith and jth basis functions)
        end
    end

    %% Calculating set of widths of elements 'h'
    h = zeros(1,N); %Set of gaps between points
    for k = 1:N
        h(k) = nodes(k+1)-nodes(k);
    end

    %% Storing the positions of all the required orders
    positions = zeros(max_order,N); %The ith row contains the ascending number of the elements with ith order contributions.

    count = 2;

    for i = 1:max_order
        for j = 1:N
            if orders(j)>=i
                if j == N && i==1
                else
                    positions(i,j) = count;
                    count=count+1;
                end

            end
        end
    end

    %% Initialising the matrix
    non_zero_max = ((N-1)+2*(N-2)+2)+(2*2*(max_order-1)*(N-1))+((max_order-1)^2*N); %This is the number of non-zero elements in the sparse matrix ((linear-linear interactions)+(linear-higher order interactions)+(higher order-higher order interactions)

    rows = zeros(non_zero_max,1); %This is the set of the row positions of elements in the matrix
    columns = zeros(1,non_zero_max); %This is the set of the column positions of elements in the matrix
    values = zeros(1,non_zero_max); %This is the set of values in the matrix (where a((rows(i)),(columns(i)) = values(i))

    row_total = count;

    %% Appending the equation u0 = A
    rows(1) = 1;
    columns(1) = 1;
    values(1) = 1;

    start = 2; %We set up 'start' to be the position in rows,columns and values after the last element appended
    %% Appending the contributions of basis functions 1-Order*N

    for k = 1:N %This iterates over the N intervals in the system
        Order = orders(k);
        %% Linear down-down interactions (k,k)
        %This appends all the interaction of the linear downslope with itself on interval k
        if  k~=1 %There's no downslope for the first interval

            %The last thing to be put in the loop is the down-down interaction in position (k+1,k+1) for interval k-1. This overlaps with (k,k) for interval k so we recall the same element with start

            values(start) = values(start) + 1/h(k)*integral_matrix(1,1,3)+m*integral_matrix(1,1,2)+n*h(k)*integral_matrix(1,1,1); %This is the integral from the interaction on the kth interval of the two downs of theta(k)

            %(The integral is rescaled from the kth interval to interval [0,1] so we need only refer to the bases we set up initially)
            start = start+1;

        end

        %% Linear up-down interactions (k,k+1)
        %This appends all the interaction of the linear downslope with linear upslope on interval k
        if  k~=1 && k~=N %There's no upslope for the last interval or downslope for the first interval

            rows(start) = k;
            columns(start) = k+1;
            values(start) = (1/h(k)*integral_matrix(1,2,3)+m*(integral_matrix(1,2,2))+n*(h(k)*integral_matrix(1,2,1))); %This is the integral from the interaction on the kth interval of the down of theta(k) and the up of theta(k+1)

            rows(start+1) = k+1;
            columns(start+1) = k;
            values(start+1) = (1/h(k)*integral_matrix(2,1,3)+m*(integral_matrix(2,1,2))+n*(h(k)*integral_matrix(2,1,1)));  %This is the integral from the interaction on the kth interval of the up of theta(k+1) and the down on theta(k)
            start = start+2;

        end

        %% Linear- higher order interactions (k+1,k+(j-1)N)
        if Order>1
            for j = 2:(Order)
                %Down-Higher order interactions (k+1, - )
                %This appends all the higher order interactions of the linear downslope on the kth interval
                if k~=N
                    rows(start) = k+1;
                    columns(start) = positions(j,k);
                    values(start) = (1/h(k)*integral_matrix(2,j+1,3)+m*(integral_matrix(2,j+1,2))+n*(h(k)*integral_matrix(2,j+1,1))); %This is the integral from the interaction on the kth interval of the down of theta(k) and the order j basis function (in position j+1 as there are 2 linear functions)

                    rows(start+1) = positions(j,k);
                    columns(start+1) = k+1;
                    values(start+1) = (1/h(k)*integral_matrix(j+1,2,3)+m*(integral_matrix(j+1,2,2))+n*(h(k)*integral_matrix(j+1,2,1)));  %This is the integral from the interaction on the kth interval of the down of theta(k) and the order j basis function

                    start = start+2;
                end

                %Up-Higher order interactions (k, - )
                %This appends all the higher order interactions of the linear upslope on the kth interval
                if k~=1
                    rows(start) = k;
                    columns(start) = positions(j,k);
                    values(start) = 1/h(k)*integral_matrix(1,j+1,3)+m*(integral_matrix(1,j+1,2))+n*(h(k)*integral_matrix(1,j+1,1)); %This is the integral from the interaction on the kth interval of the up of theta(k) and the order j basis function

                    rows(start+1) = positions(j,k);
                    columns(start+1) = k;
                    values(start+1) = 1/h(k)*integral_matrix(j+1,1,3)+m*(integral_matrix(j+1,1,2))+n*(h(k)*integral_matrix(j+1,1,1));  %This is the integral from the interaction on the kth interval of the down of theta(k) and the order j basis function

                    start = start+2;
                end
            end
        end

        %% Higher order-higher order interactions

        for Order_row = 2:Order
            for Order_col = 2:Order
                rows(start) = positions(Order_row,k);
                columns(start) = positions(Order_col,k);
                values(start) = 1/h(k)*integral_matrix(Order_row+1,Order_col+1,3)+m*(integral_matrix(Order_row+1,Order_col+1,2))+n*(h(k)*integral_matrix(Order_row+1,Order_col+1,1));%This is the integral from the interaction on the kth interal of the Order_row and Order_col basis functions

                start = start+1;
            end
        end

        %% Linear down-down interactions (k+1,k+1)
        if k~=N
            rows(start) = k+1;
            columns(start) = k+1;
            values(start) = 1/h(k)*integral_matrix(2,2,3)+m*(integral_matrix(2,2,2))+n*(h(k)*integral_matrix(2,2,1)); %This is the integral from the interaction on the kth interval of the two linear downs of theta(k)
        end

    end

    %% Appending the equation u(2N) = B
    rows(start) = row_total;
    columns(start) = row_total;
    values(start) = 1;

    rows = rows(1:start);
    columns = columns(1:start);
    values = values(1:start);

    %% Generating the matrix
    S = sparse(rows,columns,values,row_total,row_total);

    if Output == 3
        fprintf("Matrix = \n")
        disp(full(S))
    end
    %% Computing the corresponding vector

    F = zeros(row_total,1);

    F(1) = A;
    F(row_total) = B;

    for k = 1:N %This iterates over the intervals
        Order = orders(k);
        %% Upslopes
        if k ~= (N)
            F(k+1) = F(k+1) + h(k)*integral(@(x)f(h(k)*x+nodes(k)).*(basis{1,2}(x)),0,1); %Upslope interactions on the kth interval (theta(k+1)
        end

        %% Downslopes
        if k ~= 1
            F(k) = F(k) + h(k)*integral(@(x)f(h(k)*x+nodes(k)).*(basis{1,1}(x)),0,1); %Downslope interactions on the kth interval (theta(k))
        end

        %% Higher orders

        if Order ~=1
            for j = 2:(Order)
                F(positions(j,k)) = F(positions(j,k)) + h(k)*integral(@(x)f(h(k)*x+nodes(k)).*(basis{1,j+1}(x)),0,1); %Higher order interactions on the kth interval (theta(k+(j-1)N))
            end
        end
    end

    %% Removing interactions with theta(0) (down)
    F(2) = F(2)-A*(1/h(1)*integral_matrix(1,2,3)+m*(integral_matrix(2,1,2))+n*(h(1)*integral_matrix(1,2,1))); %Removes the interaction of theta(0) (down) with theta(1) (up)

    Order = orders(1);
    if Order~=1
        for k = 1:Order-1
            F(positions(k+1,1)) = F(positions(k+1,1))-A*(1/h(1)*integral_matrix(1,k+2,3)+m*(integral_matrix(k+2,1,2))+n*(h(1)*integral_matrix(1,k+2,1))); %This removes the interaction of theta(0) (down) with any other bases
        end
    end


    %% Removing interactions on the last interval (up)
    F(N) = F(N)-B*(1/h(N)*integral_matrix(2,1,3)+m*(integral_matrix(1,2,2))+n*(h(N)*integral_matrix(2,1,1))); %Removes the interaction of theta(end) (up) with theta(N-1) (down)

    Order = orders(N);
    if Order~=1
        for k = 1:Order-1
            F(positions(k+1,N)) = F(positions(k+1,N))-B*(1/h(N)*integral_matrix(2,k+2,3)+m*(integral_matrix(k+2,2,2))+n*(h(N)*integral_matrix(2,k+2,1))); %This removes the interaction of theta(N) (up) with any other bases
        end
    end

    if Output == 3
        fprintf("Vector = \n")
        disp(F')
    end

    %% Solving the sysytem
    u = linsolve(full(S),F);

    if Output == 3
        fprintf("Values = \n")
        disp(u')
    end

    %% Evaluating modelled points and true error
    element_pts = 1000; %Number of points evaluated on each interval
    x_vals = zeros(1,element_pts*N);
    y_vals = zeros(1,element_pts*N);

    temppts = linspace(0,1,element_pts);

    priori_error_vect_L2 = zeros(1,N);
    priori_error_vect_H1 = zeros(1,N);

    for k = 1:N
        interval = (element_pts*(k-1)+1):(element_pts*k);
        x_vals(interval) = h(k)*temppts+nodes(k);
        Order = orders(k);
        interval_model = @(x) zeros(size(x));
        interval_model_derivative = @(x) zeros(size(x));

        % Linear contributions
        if k ~=N
            interval_model = @(x) interval_model(x) + u(k)*basis{1,1}(x)+u(k+1)*basis{1,2}(x);
            interval_model_derivative = @(x) interval_model_derivative(x) + u(k)*basis{2,1}(x)+u(k+1)*basis{2,2}(x);
        else
            interval_model = @(x) interval_model(x) + u(k)*basis{1,1}(x)+u(row_total)*basis{1,2}(x);
            interval_model_derivative = @(x) interval_model_derivative(x) + u(k)*basis{2,1}(x)+u(row_total)*basis{2,2}(x);
        end

        % Higher-order contributions
        if Order~=1
            for j = 1:Order-1
                interval_model = @(x) interval_model(x) + u(positions(j+1,k)).*basis{1,j+2}(x);
                interval_model_derivative = @(x) interval_model_derivative(x) + u(positions(j+1,k)).*basis{2,j+2}(x);
            end
        end
        y_vals(interval) = interval_model(temppts);

        %L2 norm
        L2 = h(k)*(integral(@(x)(True_solution(h(k)*x+nodes(k))-interval_model(x)).^2,0,1));

        %L2' norm
        L2_derivative = integral(@(x)(True_solution_derivative(h(k)*x+nodes(k))-interval_model_derivative(x)).^2,0,1);

        priori_error_vect_L2(k) = sqrt(L2);
        priori_error_vect_H1(k) = sqrt(L2+L2_derivative);
    end

    %% Calculating aposteriori error
    posteriori_error_vect = zeros(1,N);

    for k = 1:N
        Order = orders(k);
        res = @(x) f(h(k)*x+nodes(k));

        % Linear contributions
        if k ~=N
            res = @(x) res(x) - (u(k).*((m/h(k)).*basis{2,1}(x)+n.*basis{1,1}(x))+u(k+1).*((m/h(k)).*basis{2,2}(x)+n.*basis{1,2}(x)));
        else
            res = @(x) res(x) - (u(k).*((m/h(k)).*basis{2,1}(x)+n.*basis{1,1}(x))+u(row_total).*((m/h(k)).*basis{2,2}(x)+n.*basis{1,2}(x)));
        end

        % Higher-order contributions
        if Order~=1
            for j = 1:Order-1
                res = @(x) res(x)-u(positions(j+1,k)).*(-1*(1/h(k)).^2.*basis{3,j+2}(x)+(m/h(k)).*basis{2,j+2}(x)+n.*basis{1,j+2}(x));
            end
        end

        res = @(x) res(x).^2;

        res = @(x) res(x).*abs((1-x).*x);

        %1/(p_j(p_{j+1}))||R(u_h)(x_{j+1}-x)(x-x_j)||_{L2}
        RL2 = (1/(Order*(Order+1)))*h(k)*integral(@(x)res(x),0,1);
        %R3 = (h(k)/Order)^2*h(k)*integral(@(x)res(x),0,1);
        posteriori_error_vect(k) = sqrt(RL2);
    end

    %% Creating and plotting functions

    if Output~=0
        if Output == 1
            %% Plotting model
            figure
            hold on
            p1 = plot(x_vals,y_vals,'b');
            yline(0)
            plot(nodes,zeros(length(nodes)),'k|')

            %% Plotting true solution for comparison
            pts2 = linspace(a,b,100); %Set of points the true function is evaluated over
            p2 = plot(pts2,True_solution(pts2),'r');
            legend([p1 p2],'Approximation','True Solution','Location','northeast')
            title(sprintf("Approximation to the function with %d elements at max order %d",N,max_order))

            hold off
        else
            %% Plotting model
            figure
            subplot(2,1,1)
            hold on
            p1 = plot(x_vals,y_vals,'b');

            %yline(0) %Matlab code
            plot(xlim, [0, 0], 'k--'); % Updated Octave code

            plot(nodes,zeros(length(nodes)),'k|')

            %% Plotting true solution for comparison
            pts2 = linspace(a,b,100); %Set of points the true function is evaluated over
            p2 = plot(pts2,True_solution(pts2),'r');
            legend([p1 p2],'Approximation','True Solution','Location','northeast')
            title(sprintf("Approximation to the solution with %d elements at max order %d",N,max_order))
            hold off

            %% Plotting errors over the intervals
            subplot(2,1,2)
            hold on

            for i = 1:N
                interval = (element_pts*(i-1)+1):(element_pts*i);
                plot(x_vals(interval),log(sqrt(priori_error_vect_L2(i)))*ones(1,element_pts),'b')
                plot(x_vals(interval),log(sqrt(priori_error_vect_H1(i)))*ones(1,element_pts),'r')
                plot(x_vals(interval),log(sqrt(posteriori_error_vect(i)))*ones(1,element_pts),'g')
            end

            legend('Actual norm error (L2)','Actual norm error (H1)','Predicted norm error','Location','southoutside')
            title("Log of error norm over the intervals")
            %axis([a b -20 0])
            hold off
        end
    end
end
