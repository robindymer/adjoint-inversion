SP = fullfile(fileparts(mfilename('fullpath')), 'results');
files = dir(fullfile(SP, 'res_*.mat'));
fid = fopen(fullfile(SP, 'summary.csv'), 'w');
fprintf(fid, 'method,iter,t,f,err,nFwd,nAdj,nHv,cost\n');
fprintf('%-24s %5s %7s %10s %7s %5s %5s %5s | cost to f<: %6s %6s %6s %6s\n', 'method','iters','cost','f_final','err','fwd','adj','Hv','1','1e-1','1e-2','1e-3');
for i = 1:numel(files)
    S = load(fullfile(SP, files(i).name));
    h = S.hist;
    cost = (4.3*h.nFwd + 1.3*h.nAdj + 4.1*h.nHv)/5.6;
    for k = 1:numel(h.f)
        fprintf(fid, '%s,%d,%g,%g,%g,%d,%d,%d,%g\n', S.method, k-1, h.t(k), h.f(k), h.err(k), h.nFwd(k), h.nAdj(k), h.nHv(k), cost(k));
    end
    thr = [1 1e-1 1e-2 1e-3]; c = nan(size(thr));
    for j = 1:numel(thr)
        k = find(h.f < thr(j), 1); if ~isempty(k), c(j) = cost(k); end
    end
    fprintf('%-24s %5d %7.1f %10.3e %7.4f %5d %5d %5d | %6.1f %6.1f %6.1f %6.1f\n', erase(erase(files(i).name, 'res_'), '.mat'), numel(h.f)-1, cost(end), h.f(end), h.err(end), h.nFwd(end), h.nAdj(end), h.nHv(end), c);
end
fclose(fid);
