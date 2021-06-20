function relerr = relative_error_D2n(xref, xest)
% Relative error between xref (reference) and xest, up to an element of
% D_2n

xest1 = align_to_reference(xest, xref);
relerr1 = norm(xref-xest1) / norm(xref);

xest2 = align_to_reference(reverse(xest), xref);
relerr2 = norm(xref-xest2) / norm(xref);

relerr = min(relerr2, relerr1);

end