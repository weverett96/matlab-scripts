function S = ewma(dat, alpha)
% EWMA dynamic conditional beta following EWMA variance process.
%
% Arguments:
% dat: TxK matrix containing time series data.
% alpha: smoothing parameter

    dat(isnan(dat)) = 0;
    T = size(dat,1);
    S = alpha*(dat(T,:)'*dat(T,:)) + (1-alpha)*ewma_recursion(dat,alpha,T-1);

    function S = ewma_recursion(dat, alpha, T)
        if T > 0
            S = alpha*(dat(T,:)'*dat(T,:)) + (1-alpha)*ewma_recursion(dat,alpha,T-1);
        else
            N = size(dat,1);
            % Use starting tenth of data for initialization.
            dat_init = dat(1:ceil(N/10),:);
            S = (1/(N-1))*(dat_init - ones(size(dat_init,1),1)*mean(dat_init))'*...
                (dat_init - ones(size(dat_init,1),1)*mean(dat_init));
        end
    end
end