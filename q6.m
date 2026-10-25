
[x, it, xall, fall, range] = golden(@molecule, 0.1, 1.5, 10e-6, 75);
fprintf('%6.2e \n', x);

err = zeros(1, it);
for i = 1:it
    err(i) = abs(range(i, 2) - range(i, 1));
end

fprintf('x f(x) error\n');
for i = 1:it
    fprintf('%3d %15.12f %15.12f %10.3e\n', i-2, xall(i), fall(i), err(i));
end
fprintf("\n");


[x, it, xall, fall, range] = succparab(@molecule, 0.2, 1.4, 10e-6, 100);
fprintf('%6.2e \n', x);

succerr = zeros(1, length(range));
for i = 1:length(range)
    succerr(i) = abs(range(i, 2) - range(i, 1));
end

fprintf('x f(x) error\n');
for i = 1:it
    fprintf('%3d %15.12f %15.12f %10.3e\n', i-2, xall(i), fall(i), succerr(i));
end
fprintf("\n");


[x, it, xall, resall] = newton(@molecule_der, 0.8, 10e-6, 80);
fprintf('%6.2e \n', x);

disp(xall);



newtonfall = zeros(1, length(xall));
for i = 1:length(xall)
    newtonfall(i) = molecule(xall(i));
end

fprintf('x f(x) error\n');
for i = 1:it
    fprintf('%3d %15.12f %15.12f %10.3e\n', i-1, xall(i), newtonfall(i), resall(i));
end
fprintf("\n");

newerr = zeros(1, it);
for i = 1:it
    newerr(i) = abs(resall(i));
end



function y = test_func1(x)
    y = -(1/(x^2 + 1));
end

% golden section search


function [x, it, xall, fall, range] = golden(fun, a, b, tol, maxit)
    it = 0;
    range = [a b];

    tau = (-1 + sqrt(5))/2;
    x1 = a + (1 - tau) * (b - a);
    f1 = fun(x1);
    x2 = a + tau * (b - a);
    f2 = fun(x2);
    
    xall = x2;
    fall = f2;


    for i = 1:maxit
        if (b - a < tol)    
            break;
        end

        it = it + 1; 
        if f1 >= f2
            a = x1;
            x1 = x2;
            f1 = f2;

            x2 = a + tau * (b - a);
            f2 = fun(x2);

            xall = [xall; x2]; 
            fall = [fall; f2];
            range = [range; a b];
        else 
            b = x2;
            x2 = x1;
            f2 = f1;

            x1 = a + (1-tau) * (b - a);
            f1 = fun(x1);

            xall = [xall; x1]; 
            fall = [fall; f1];
            range = [range; a b];
        end

    end
    x = (a + b)/2;
end

% successive parabolic interpolation

function [x, it, xall, fall, range] = succparab(fun, a, b, tol, maxit)
    it = 0;
    range = [a b];


    % get initial three points
    u = a;
    fu = fun(u);
    v = (a+b)/2;
    fv = fun(v);
    w = b;
    fw = fun(w);
    
    xall = [u; v; w];
    fall = [fu; fv; fw];

    for i = 1:maxit
        it = it + 1;

        % calculate the p and q
        p = (v-u)^2*(fv-fw) - (v-w)^2*(fv-fu);
        q = -2 * ((v-u) * (fv-fw) - (v-w) * (fv-fu));

        u = w;
        fu = fw;
        w = v;
        fw = fv;
        v = v + p/q;
        fv = fun(v);
        
        
        M = max([u v w]);
        m = min([u v w]);

        xall = [xall; v];
        fall = [fall; fv];
        range = [range; m M];

        if((M-m) < tol)
            break;
        end
    end
    
    x = v;

end


% newton's 
function [x, it, xall, resall] = newton(fun, a, tol, maxit)
    it = 0;
    x = a;
    xall = x;
    [first_der, second_der] = fun(x);
    resall = first_der;

    for i = 1:maxit
        it = it + 1;
        if(abs(first_der) < tol)
            break;
        end
        x = x - first_der/second_der;
        
        xall = [xall; x];

        [first_der, second_der] = fun(x);
        resall = [resall; first_der];
    end
end


function y = molecule(x)
    y = 1/(x^12) - 2/(x^6);
end

function [first_der, second_der] = molecule_der(x)
    first_der = 12/(x^7) - 12/(x^13);
    second_der = 156/(x^14) - 84/(x^8);
end
















