function checkSymmetryAndEigenvalues(M)

    Lmax = eigs((M+M')/2, 100, 'largestreal');
    if sum(isnan(Lmax)) > 0
        disp('Will have to compute eigenvalues using eig, which may be very slow.')
        answer = input('Continue(y) or abort(n)?', 's');
        if strcmp(answer, 'y')
            Lmax = max(eig(full(M)));
        else
            return
        end
    else
        Lmax = max(Lmax);
    end

    s = max(max(abs(M-M')));
    Mmax = max(max(abs(M)));

    disp(['Largest eigenvalue: ' num2str(Lmax)]);
    disp(['Largest symmetry violation: ' num2str(s)]);
    disp(['Largest relative symmetry violation: ' num2str(s/Mmax)]);

end

