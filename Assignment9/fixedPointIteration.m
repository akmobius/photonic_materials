function [solution, error] = fixedPointIteration(f, x_0, tolerance, n)
    % === TYPE YOUR CODE BELOW === %
    
    % Initialize variables
    current_x = x_0;
    error = zeros(n, 1); % Preallocate error array for efficiency
    solution = x_0;      % Default solution to initial guess
    
    for i = 1:n
        % Evaluate the function at the current x
        % f is a symbolic object, so we substitute current_x into it
        next_x = double(subs(f, current_x));
        
        % Calculate percent relative error
        if next_x ~= 0
            current_error = abs((next_x - current_x) / next_x) * 100;
        else
            current_error = abs(next_x - current_x) * 100;
        end
        
        % Store error
        error(i) = current_error;
        
        % Update solution and current x
        solution = next_x;
        current_x = next_x;
        
        % Check for convergence
        if current_error < tolerance
            % Trim the error array to the number of iterations actually performed
            error = error(1:i);
            return;
        end
    end
    
    % Warning if the maximum number of iterations is reached without convergence
    if i == n
        fprintf('Warning: Maximum iterations reached without meeting tolerance.\n');
    end
end