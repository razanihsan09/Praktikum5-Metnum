function [x, D, c] = eliminasi_gauss_jordan(A, b, tampil)
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
        fprintf('=== METODE GAUSS-JORDAN ===\n');
        fprintf('Matriks augmented awal [A | b]:\n'); disp(Aug);
    end

    % --- Forward Elimination: nolkan elemen di bawah diagonal ---
    for i = 1:(n - 1)
        if abs(Aug(i, i)) < eps
            error('Pivot nol pada baris %d, diperlukan pivoting.', i);
        end
        m = Aug(i+1:n, i) / Aug(i, i);
        Aug(i+1:n, :) = Aug(i+1:n, :) - m * Aug(i, :);
        Aug(i+1:n, i) = 0;
    end
    if tampil
        fprintf('Setelah forward elimination (segitiga atas):\n'); disp(Aug);
    end

    % --- Backward Elimination: nolkan elemen di atas diagonal ---
    for i = n:-1:2
        m = Aug(1:i-1, i) / Aug(i, i);
        Aug(1:i-1, :) = Aug(1:i-1, :) - m * Aug(i, :);
        Aug(1:i-1, i) = 0;
        if tampil
            fprintf('Pengali m (kolom %d, ke atas): %s\n', i, mat2str(m.', 6));
            fprintf('Setelah backward elimination kolom %d:\n', i); disp(Aug);
        end
    end

    D = Aug(:, 1:n);
    c = Aug(:, end);

    % --- Solusi langsung: x_i = c_i / D_ii (tanpa substitusi mundur) ---
    x = c ./ diag(D);

    if tampil
        fprintf('Matriks diagonal D =\n'); disp(D);
        fprintf('Vektor c akhir =\n'); disp(c);
        fprintf('Solusi:\n');
        for i = 1:n, fprintf('  x%d = %.10f\n', i, x(i)); end
    end
end
