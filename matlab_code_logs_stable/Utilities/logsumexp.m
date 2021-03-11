function ret = logsumexp(l, opt_d)
if nargin > 1
    m = max(l, [], opt_d);
    m(isinf(m)) = 0;
    ret = m + log(sum(exp(l - m), opt_d));
else
    m = max(l);
    m(isinf(m)) = 0;
    ret = m + log(sum(exp(l - m)));
end