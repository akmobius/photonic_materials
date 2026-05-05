% Function to find the roots of an arbitrary function g(x) using Newton's
% method.
% 
% ~~~~~~~~~~~~~~~~~~
% ----- INPUTS -----
% ~~~~~~~~~~~~~~~~~~
% g = arbitrary function, passed as a symbolic object.
% x_0 = initial guess.
% tolerance = percent error at which the equation is considered solved.
% n = maximum number of iterations.
% 
% --------- ---------
% ~~~~~ OUTPUTS ~~~~~
% --------- ---------
% solution = the solution of the equation.
% error = percent error as a function of iteration number.
% 
function [solution, error] = NewtonsMethod(g, x_0, tolerance, n)
    % === TYPE YOUR CODE BELOW === %
    
    % Initialize variables
    current_x = x_0;
    error = zeros(n, 1);
    solution = x_0;
    
    % Calculate the derivative of g symbolically
    dg = diff(g);
    
    for i = 1:n
        % Evaluate g(x) and g'(x) at the current x
        g_val = double(subs(g, current_x));
        dg_val = double(subs(dg, current_x));
        
        % Avoid division by zero if the derivative is too small
        if abs(dg_val) < 1e-12
            fprintf('Warning: Derivative is near zero. Newton''s method may fail.\n');
            error = error(1:i-1);
            return;
        end
        
        % Newton-Raphson update formula
        next_x = current_x - (g_val / dg_val);
        
        % Calculate percent relative error
        if next_x ~= 0
            current_error = abs((next_x - current_x) / next_x) * 100;
        else
            current_error = abs(next_x - current_x) * 100;
        end
        
        % Store error and update values
        error(i) = current_error;
        solution = next_x;
        current_x = next_x;
        
        % Check for convergence
        if current_error < tolerance
            error = error(1:i);
            return;
        end
    end
    
    if i == n
        fprintf('Warning: Maximum iterations reached without meeting tolerance.\n');
    end
end