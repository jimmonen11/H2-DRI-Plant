function [winnerbox] = minimax(Nsamp, Ndim, Nmatch)

%outputs the winnerbox to run GP on given number of samples (Nsamp) that
%want to run for the GP, the number of dimensions the GP has (Ndim) and the
%number of matches want to run (Nmatch)

%Low discrepency sampling to find points to run Simulink model off of
%Finds the maximum minimum distance between any two points for a single
%Nmatch then loops through the number of Nmatches to find the maximum
%minimum distance for all the uniform random samples

for i = 1:Nmatch
    r = 2; %initialize r - no difference can be bigger than 1 on unit box
    box = rand(Nsamp, Ndim);

    %sum over the number of dimensions squaring each
    %r = sqrt(x1^2 + x2^2 + x3^2 + ...) just leave out sqrt
    r_matrix = sum(eu_dist(box, box).^2, 3);
    
    %Find the minimum of the lower triangual of r_matrix
    for ii = 2:Nsamp
        for jj = 1:ii-1
            r = min(r, r_matrix(ii,jj));
        end
    end
    
    if i == 1 | r > winner_r
        winner_r = r;
        winnerbox = box;
    end
    
end
        
        
        
        
        