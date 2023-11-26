function y = interp1d(A, idx, xi)

    if xi < idx(1)
        r = 1;
    elseif xi < idx(end)
        r = max(find(idx <= xi));
    else
        r = length(idx) - 1;
    end

    DA = (xi - idx(r)) / (idx(r + 1) - idx(r));
    y = A(r) + (A(r + 1) - A(r)) * DA;
end
