%вспомогательная функция, возвращающая матрицу степеней
function [X] = matrix_of_deg(n)
matr = 0:n-1;
X = flip(matr);
end
