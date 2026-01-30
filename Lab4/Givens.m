A = [0, 20, -14;
     3, 27, -4;
     4, 11, -2];
[m, n] = size(A);
R = A;
Q = eye(m);

% Iterate through columns (k) and rows (j) below the diagonal
for k = 1:min(m-1,n)
    for j = k+1:m
        if R(j,k) ~= 0
            % Calculate Givens parameters c and s
            % We annihilate the element at (j,k) using (k,k)
            a = R(k,k);
            b = R(j,k);
            r = sqrt(a^2 + b^2);
            c = a / r;
            s = -b / r; 

            % Construct the Givens Rotation Matrix G
            G = eye(m);
            G(k,k) = c;
            G(j,j) = c;
            G(k,j) = -s;
            G(j,k) = s;

            % Update R (Apply rotation to A)
            R = G * R;
            
            % Update Q (Accumulate the transpose of the rotations)
            % Since G_n...G_1 A = R  =>  A = G_1' ... G_n' R
            Q = Q * G';
            
            % Display the Givens matrix for this step
            disp(G);
        end
    end
end

disp(Q);
disp(R); 