function [x, U, c] = eliminasi_gauss(A, b, tampil)
    if nargin < 2
        A = [ 2,  1, -1;
              4,  3,  1;
             -2,  1,  2 ];
        b = [ 3; 9; 4 ];
    end
    if nargin < 3, tampil = true; end

    n   = size(A, 1);
    Aug = [A, b(:)];

    if tampil
        fprintf('=== METODE ELIMINASI GAUSS ===\n');
        fprintf('Matriks augmented awal [A | b]:\n'); disp(Aug);
    end

    % --- Forward Elimination: ubah A menjadi segitiga atas ---
    for i = 1:(n - 1)
        if abs(Aug(i, i)) < eps
            error('Pivot nol pada baris %d, diperlukan pivoting.', i);
        end
        m = Aug(i+1:n, i) / Aug(i, i);
        Aug(i+1:n, :) = Aug(i+1:n, :) - m * Aug(i, :);
        Aug(i+1:n, i) = 0;
        if tampil
            fprintf('Pengali m (kolom %d): %s\n', i, mat2str(m.', 6));
            fprintf('Setelah eliminasi kolom %d:\n', i); disp(Aug);
        end
    end

    U = Aug(:, 1:n);
    c = Aug(:, end);

    % --- Backward Substitution: U*x = c ---
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (c(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end

    if tampil
        fprintf('Matriks U =\n'); disp(U);
        fprintf('Vektor c (b setelah eliminasi) =\n'); disp(c);
        fprintf('Solusi:\n');
        for i = 1:n, fprintf('  x%d = %.10f\n', i, x(i)); end
    end
end
