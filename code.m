clc;
clear all;

n = 21;          % Number of points
h = 0.1;         % Step size
x = linspace(0,2,n);   % Generate evenly spaced points between 0 and 2
T = zeros(19,1);   % Initialize T as a column vector of zeros

% Loop to iterate and update T
for iter = 1:20
    J = jac();          % Call the jac() function to calculate Jacobian matrix
    R = r(T);           % Call the r() function to calculate the residual vector
    T1 = T + gauss(J,-R);   % Solve the system of equations using Gauss elimination method
    error = norm2(T1,T);    % Calculate the 2-norm error of the difference between T1 and T
    T = T1;             % Update T with T1
end

Temperature = ones(n,1);    % Initialize Temperature as a column vector of ones
Temperature(1) = 5;        % Set the Temperature of point A to 5

% Loop to populate Temperature based on T
for i = 2:20
    Temperature(i) = T(i-1);
end

Temperature(21) = 2.38144;   % Set the Temperature at point B to 2.38144

plot(x, Temperature, '-o','MarkerEdgeColor','r')   % Plot Temperature against x
xlabel('x');
ylabel('Temperature');


% Define the 'a' function
function y = a(x)
    y = -(x+3)/(x+1);
end

% Define the 'b' function
function y = b(x)
    y = (x+3)/(x+1)^2;
end

% Define the 'r' function to calculate the residual vector
function R = r(y)
    h = 0.1;
    R = zeros(19,1);    % Initialize R as a column vector of zeros
    
    R(1) = y(2)(1/h^2 + a(0.1)/(2*h)) + y(1)(b(0.1) - 2/h^2) + 5*(1/h^2 - a(0.1)/(2*h)) - 2.2 - 3*b(0.1);
    
    for i = 2:18
        R(i) = y(i+1)(1/h^2 + a(0.1*i)/(2*h)) + y(i)(b(0.1*i) - 2/h^2) + y(i-1)(1/h^2 - a(0.1*i)/(2*h)) - 2(i*h+1) - 3*b(i*h);
    end
    
    R(19) = 2.38144*(1/h^2 + a(1.9)/(2*h)) + y(19)(b(1.9) - 2/h^2) + y(18)(1/h^2 - a(1.9)/(2*h)) - 5.8 - 3*b(1.9);
end

% Define the 'jac' function to calculate the Jacobian matrix
function J = jac()
    h = 0.1;
    J = zeros(19,19);    % Initialize J as a 19x19 matrix of zeros
    
    J(1,1) = b(0.1) - 2/h^2;
    J(1,2) = 1/h^2 + a(0.1)/(2*h);
    
    for i = 2:18
        J(i,i-1) = 1/h^2 - a(0.1*i)/(2*h);
        J(i,i) = b(0.1*i) - 2/h^2;
        J(i,i+1) = 1/h^2 + a(0.1*i)/(2*h);
    end
    
    J(19,18) = 1/h^2 - a(1.9)/(2*h);
    J(19,19) = b(1.9) - 2/h^2;
end

% Define the 'gauss' function for solving a system of equations using Gauss elimination method
function u = gauss(v, w)
    n = length(v);
    
    % Making the matrix Diagnol Dominant
    for k = 1:n
        s = k;
        m = abs(v(k,k));
        
        % Find the maximum element in the column below v(k,k)
        for j = k+1:n
            if abs(v(j,k)) > m
                s = j;
                m = abs(v(j,k));
            end
        end
        
        % Swap rows in v and w
        temp = v(s,:);
        v(s,:) = v(k,:);
        v(k,:) = temp;
        
        t = w(s);
        w(s) = w(k);
        w(k) = t;
        
        % Check if v(k,k) is still zero
        if v(k,k) == 0
            % Checking any rows above for a non-zero value for the kth
            % column
            for j = 1:n
                if v(j,k) ~= 0
                    x = j;
                    break
                end
            end
            % Applying Row operation to make the pivot non-zero
            w(k) = w(k) + w(x);
            for i = 1:n
                v(k,i) = v(k,i) + v(x,i);
            end
        end
    end
    
    % Gauss Elimination
    for k = 1:n-1
        for i = k+1:n
            m = v(i,k)/v(k,k);
            for j = 1:n
                v(i,j) = v(i,j) - m*v(k,j);
            end
            w(i) = w(i) - m*w(k);
        end
    end
    
    %Backward Substitution
    u = zeros(n,1);
    u(n) = w(n)/v(n,n);
    
    for i = n-1:-1:1
        s = w(i);
        for j = i+1:n
            s = s - v(i,j)*u(j);
        end
        u(i) = s/v(i,i);
    end
end

% Define the 'norm2' function to calculate the 2-norm of the difference between two vectors
function y = norm2(i, j)
    l = 0;
    for x = 1:4
        l = l + (i(x) - j(x))^2;
    end
    y = abs(l^0.5);
end
