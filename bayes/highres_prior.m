% Prior of the high-res runs and the tag of its result files (run_highres_*.m):
%
%   gammaFactor   the increment std is gammaFactor times the scaled default of smoothness_prior.m;
%                 1 gives the old file names, otherwise the tag is _g<gammaFactor>
%   priorType     'gaussian' (default) or 'tv' (tv_prior.m with the matching scale b, times
%                 gammaFactor), which adds _tv to the tag
function [prior, tag] = highres_prior(m_p, gammaFactor, priorType)
    if nargin < 2, gammaFactor = 1; end
    if nargin < 3, priorType = 'gaussian'; end
    gamma = gammaFactor*1e-3*sqrt(25/(m_p - 1));
    tag = '';
    if gammaFactor ~= 1
        tag = sprintf('_g%g', gammaFactor);
    end
    switch priorType
        case 'gaussian'
            prior = smoothness_prior(m_p, gamma);
        case 'tv'
            prior = tv_prior(m_p, gamma*sqrt(2/pi));
            tag = [tag, '_tv'];
        otherwise
            error('Unknown priorType: %s', priorType);
    end
end
