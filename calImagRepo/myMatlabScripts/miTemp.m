function mi = miTemp(n_t,i_t,Npix)
M = round(max(n_t));

PL = zeros(M+1,Npix);
PL_relSE = PL;
for j=1:Npix
    inds = find(i_t==j);%conditional ensemble for emissions when mouse is in pixel m
    vec = hist(n_t(inds),0:M)';
    PL_relSE(:,j) = 1./sqrt(vec);%relative standard error of the m-th column
    PL(:,j) = vec/length(inds);%vec/sum(vec);
end

%calculate marginal distribution of spike counts
SpkCountDist = hist(n_t,0:M)' / length(n_t);%normalized distrubution of counts.    


Pi = hist(i_t,1:Npix)/length(i_t);
mi = sum( negplogp(SpkCountDist)  - negplogp(PL)*Pi' );

end

function f = negplogp(p)
 f = zeros(size(p));
 inds = find(p>0);
 f(inds) = - p(inds) .* log(p(inds))/log(2);
end
