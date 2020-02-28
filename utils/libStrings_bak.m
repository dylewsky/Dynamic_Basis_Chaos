function stringLib = libStrings(yin,nVars,polyorder,usesine)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

n = size(yin,1);
yout = zeros(n,1+nVars+(nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))+11);

ind = 1;
% poly order 0
stringLib(ind) = '1';
ind = ind+1;

% poly order 1
for i=1:nVars
    stringLib = [stringLib,sprintf('x%d',i)];
    ind = ind+1;
end

if(polyorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            stringLib = [stringLib, sprintf('x%d * x%d',i,j)];
            ind = ind+1;
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                stringLib = [stringLib, sprintf('x%d * x%d * x%d',i,j,k)];
                ind = ind+1;
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    stringLib = [stringLib, sprintf('x%d * x%d * x%d * x%d',i,j,k,l)];
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        stringLib = [stringLib, sprintf('x%d * x%d * x%d * x%d * x%d',i,j,k,l)];
                        ind = ind+1;
                    end
                end
            end
        end
    end
end

if(usesine)
    for k=1:10
        stringLib = [stringLib, sprintf('sin(%d*x)',k), sprintf('cos(%d*x)',k)];
        ind = ind+2;
    end
end
disp('');