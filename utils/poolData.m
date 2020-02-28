function yout = poolData(yin,nVars,polyorder,usesine)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

n = size(yin,1);
% yout = zeros(n,1+nVars+(nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))+11);

% If polyorder is scalar, use all orders up to that
% If it's a vector, just use the orders it contains
if (length(polyorder) == 1)
    allPolyOrders = 0:polyorder;
else
    allPolyOrders = polyorder;
end

ind = 1;

% poly order 0
if ismember(0, allPolyOrders)
    yout(:,ind) = ones(n,1);
    ind = ind+1;
end
    
% poly order 1
if ismember(1, allPolyOrders)
    for i=1:nVars
        yout(:,ind) = yin(:,i);
        ind = ind+1;
    end
end

% poly order 2
if ismember(2, allPolyOrders)
    for i=1:nVars
        for j=i:nVars
            yout(:,ind) = yin(:,i).*yin(:,j);
            ind = ind+1;
        end
    end
end

% poly order 3
if ismember(3, allPolyOrders)
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
                ind = ind+1;
            end
        end
    end
end

% poly order 4
if ismember(4, allPolyOrders)
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l);
                    ind = ind+1;
                end
            end
        end
    end
end

% poly order 5
if ismember(5, allPolyOrders)
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*yin(:,m);
                        ind = ind+1;
                    end
                end
            end
        end
    end
end

if(usesine)
    for k=1:10
        yout = [yout sin(k*yin) cos(k*yin)];
        ind = ind+2*nVars; %added to original code
    end
end