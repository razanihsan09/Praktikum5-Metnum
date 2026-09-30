function [x, L, U, y] = dekomposisi_lu(A, b, tampil)
    if nargin < 2
        A = [ 2,  1, -1;
              4,  3,  1;
             -2,  1,  2 ];
        b = [ 3; 9; 4 ];
    end
    if nargin < 3, tampil = true; end

    b = b(:);
    n = size(A, 1);
    U = A;
    L = eye(n);

    % --- Pembentukan L dan U (catat semua pengali m) ---
    for i = 1:(n - 1)
        if abs(U(i, i)) < eps
            error('Pivot nol pada baris %d, diperlukan pivoting.', i);
        end
        m = U(i+1:n, i) / U(i, i);
        L(i+1:n, i) = m;
        U(i+1:n, :) = U(i+1:n, :) - m * U(i, :);
        U(i+1:n, i) = 0;
    end

    % --- Forward Substitution: L*y = b ---
    y = zeros(n, 1);
    for i = 1:n
        y(i) = (b(i) - L(i, 1:i-1) * y(1:i-1)) / L(i, i);
    end

    % --- Backward Substitution: U*x = y ---
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end

    if tampil
        fprintf('=== METODE DEKOMPOSISI LU ===\n');
        fprintf('Matriks L =\n'); disp(L);
        fprintf('Matriks U =\n'); disp(U);
        fprintf('Bukti A = L*U:\n'); disp(L * U);
        fprintf('Selisih maksimum |A - L*U| = %.3e\n', max(abs(A(:) - reshape(L*U, [], 1))));
        fprintf('Hasil y (L*y = b) =\n'); disp(y);
        fprintf('Solusi (U*x = y):\n');
        for i = 1:n, fprintf('  x%d = %.10f\n', i, x(i)); end
    end
end
