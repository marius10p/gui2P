function lam = normalize_clusters(M, iclust)
Nk = max(M(:));

lam = M;
for i = 1:Nk
    ix = find(iclust==i);
    nT0 = numel(ix);
    if nT0>0
        vM = lam(ix);
        vM = vM/sum(vM.^2)^.5;
        lam(ix) = vM;
    end
end