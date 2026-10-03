function A = profile_tail_area(pi, a, limitName)
if nargin < 3 || isempty(limitName)
    limitName = 'frontier';
end
pi = pi(:)';
L = numel(pi)-1;
A = 0;
switch lower(char(limitName))
    case 'frontier'
        for l = 0:L
            s = 0;
            for k = 0:l
                s = s + (l-k+1)*a^k/factorial(k);
            end
            A = A + pi(l+1)*exp(-a)*s;
        end
    case 'parallel'
        for l = 0:L
            n = l+1;
            s = 0;
            for j = 1:n
                s = s + (-1)^(j+1)*nchoosek(n,j)*exp(-j*a)/j;
            end
            A = A + pi(l+1)*s;
        end
    otherwise
        error('Unknown repair limit: %s',limitName);
end
end
