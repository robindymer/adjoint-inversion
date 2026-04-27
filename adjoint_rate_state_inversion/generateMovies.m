% Generate relevant movies
for DS = [10,4]
	marmousiForwardTests(DS,false,'CmapAcBott',@pars.marmousi2AcousticSources,'bottom');
	marmousiMoviesElasticAcoustic('CmapAcBott',DS,true, 6e-4, 5e-10, 5e-10);

	marmousiForwardTests(DS,false,'CmapAcSurf',@pars.marmousi2AcousticSources,'surface');
	marmousiMoviesElasticAcoustic('CmapAcSurf',DS,true,  6e-4, 5e-10, 5e-10);

	marmousiForwardTests(DS,false,'CmapElSeafloor10p7');
	marmousiMoviesElasticAcoustic('CmapElSeafloor10p7',DS,true,  8e-4, 1e-9, 1e-9);

end

DS = 4;
marmousiSnapshots('CmapElSeafloor10p7',DS, true,  8e-4, 1e-9, 1e-9);
marmousiSpaceTime('CmapElSeafloor10p7',DS, true,  1e-3, 2e-9, 2e-9);

