function d = dot2(A,B)
    % dot product along 2nd dimension with singleton extension
    d = sum(bsxfun(@times,A,B),2);
end