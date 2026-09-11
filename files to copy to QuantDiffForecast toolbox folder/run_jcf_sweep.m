% =========================================================================
% RUN JCF LAMBDA SWEEP
% =========================================================================
% Runs create_ensemble_jcf() across a range of lambda values in one go.
% Each run writes its own set of tagged files, so nothing overwrites
% anything else - including the calibration-only baseline produced by
% create_ensemble_simple().
%
% lambda = 1 reproduces the calibration-only baseline. Keep it in the sweep:
% it is a free correctness check. If ensemble_results_JCF-lam100.csv does not
% match ensemble_results.csv to numerical precision, something is wrong.
%
% The sweep accepts a numeric vector or a cell array. A cell entry may be
% the string 'length', which weights each period by its number of time
% points: lambda = C / (C + H). That is a non-arbitrary reference point, but
% it treats an in-sample and an out-of-sample point as equally informative,
% and with a long calibration period it sits close to calibration-only.
%
% Usage:
%   run_jcf_sweep()                              % default sweep
%   run_jcf_sweep([0 0.5 1])                     % numeric lambdas
%   run_jcf_sweep({0, 0.5, 'length', 1})         % include length-weighting
%   run_jcf_sweep({0, 'length'}, overrides)      % also override config
%
% =========================================================================

function run_jcf_sweep(lambdas, overrides)

    if nargin < 1 || isempty(lambdas)
        lambdas = {0, 0.25, 0.5, 0.75, 'length', 1};
    end
    if isnumeric(lambdas)
        lambdas = num2cell(lambdas);
    end
    if nargin < 2
        overrides = struct();
    end

    fprintf('\n########################################################\n');
    fprintf('   JCF LAMBDA SWEEP: %s\n', strjoin(cellfun(@describe_lambda, ...
        lambdas, 'UniformOutput', false), ', '));
    fprintf('########################################################\n');

    failed = {};

    for k = 1:numel(lambdas)

        cfg = overrides;
        cfg.lambda = lambdas{k};

        lab = describe_lambda(lambdas{k});

        fprintf('\n--- [%d/%d] lambda = %s ---\n', k, numel(lambdas), lab);

        try
            create_ensemble_jcf(cfg);
        catch err
            warning('lambda = %s failed: %s', lab, err.message);
            failed{end+1} = lab; %#ok<AGROW>
        end
    end

    fprintf('\n########################################################\n');
    if isempty(failed)
        fprintf('   SWEEP COMPLETE - all %d lambda values succeeded\n', numel(lambdas));
    else
        fprintf('   SWEEP COMPLETE - failed for lambda = %s\n', strjoin(failed, ', '));
    end
    fprintf('   Next: compare_ensembles_jcf()\n');
    fprintf('########################################################\n\n');
end

function s = describe_lambda(x)
    if ischar(x) || isstring(x)
        s = char(x);
    else
        s = sprintf('%.2f', x);
    end
end
