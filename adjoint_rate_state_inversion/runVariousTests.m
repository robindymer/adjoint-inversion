% Gradient comparison, one block lambda
opts.m = {[5,6]};
adjopt.compareWithFD(@Elastic2DMultiblockOpt, @pars.lambda1Block, {'lambda'}, opts, []);

% Gradient comparison, two blocks lambda
opts.m = {[5,6], [5,6]};
opts.nRec = 10;
adjopt.compareWithFD(@Elastic2DMultiblockOpt, @pars.lambda2Blocks, {'lambda'}, opts, []);

% FD convergence test, two blocks mu
opts.m = {[5,6], [5,6]};
opts.nRec = 10;
adjopt.runFDConvergence(@Elastic2DMultiblockOpt, @pars.mu2Blocks, {'mu'}, opts, [0:-2:-8])

%--- Acoustic-elastic ----
% optFunc = @ElasticAcoustic2DMultiblockOpt;
optFunc = @ElasticAcoustic2DMultiblockOptNew;

% mu, 1 block on either side
adjopt.compareWithFD(optFunc, @pars.muElAc, {'mu'}, [], []);

% mu, 2 blocks on either side
adjopt.compareWithFD(optFunc, @pars.muElAc2Blocks, {'mu'}, [], []);

% lambda, 2 blocks on either side
adjopt.compareWithFD(optFunc, @pars.lambdaElAc2Blocks, {'lambda'}, [], []);

% ------------------------