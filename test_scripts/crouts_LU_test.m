A = [1 3 4 8; 2 1 2 3; 4 3 5 8; 9 2 7 4];
L = zeros(length(A));
U = zeros(length(A));

L(1:length(A),1) = A(1:length(A),1);
for i = 1:length(A)
    U(i,i) = 1;
end
U(1,2:length(A)) = A(1,2:length(A))/L(1,1);

for j = 2 : length(A)
    for k = j : length(A)
        L(k,j) = A(k,j) - L(k,1:j-1)*U(1:j-1,j);
    end
    for k = j+1 : length(A)
        U(j,k) = (A(j,k) - (L(j,1:j-1)*U(1:j-1,k)))/L(j,j);
    end
end

% forward substitution
b = ones(length(A),1);
y = zeros(length(A),1);
for i = 1:length(A)
        if i == 1
            y(1) = b(1)/L(1,1);
        else
            y(i) = (b(i)-L(i,1:i-1)*y(1:i-1))/L(i,i);
        end
end

% backward substitution
x = zeros(length(A),1);
for i = length(A):-1:1   
    if i == length(A)
        x(length(A)) = y(length(A))/U(length(A),length(A));
    else
        x(i) = (y(i)-U(i,i+1:length(A))*x(i+1:length(A)))/U(i,i);
    end
end