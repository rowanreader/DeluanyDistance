function logical = ind2logical(indices,len)
    logical = false(len,1);
    logical(indices(~isnan(indices))) = true;
end
