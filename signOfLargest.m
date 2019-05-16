function sgn = signOfLargest(coeff)
    [~,I] = max(abs(coeff));
    sgn = sign(coeff(I));
    if sgn==0, sgn=1; end
end
