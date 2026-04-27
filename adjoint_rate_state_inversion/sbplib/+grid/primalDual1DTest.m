function tests = primalDual1DTest()
    tests = functiontests(localfunctions);
end


function testErrorInvalidLimits(testCase)
     in  = {
        {10,{1}},
        {10,[0,1]},
        {10,{1,0}},
    };

    for i = 1:length(in)
        testCase.verifyError(@()grid.primalDual1D(in{i}{:}),'grid:primalDual1D:InvalidLimits',sprintf('in(%d) = %s',i,toString(in{i})));
    end
end

function testCompiles(testCase)
    in  = {
        {5, {0,1}},
    };

    out = {
        {[0; 0.25; 0.5; 0.75; 1], [0; 0.125; 0.375; 0.625; 0.875; 1]},
    };

    for i = 1:length(in)
        [gp, gd] = grid.primalDual1D(in{i}{:});
        testCase.verifyEqual(gp.points(),out{i}{1});
        testCase.verifyEqual(gd.points(),out{i}{2});
    end
end
