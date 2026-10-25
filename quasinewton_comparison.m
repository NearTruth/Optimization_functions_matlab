
disp_output()
%get_graphs()

get_graphs_hess_mod()

[x, it, xall, resall, fall, spdall] = newton_mod(@fun, [0.9, 0.7], 10e-8, 50);

[new_order, new_order_sq] = get_order(xall, it);
disp(new_order);
disp(new_order_sq);

[x, it, xall, resall, fall, stall, ninc, ndec] = steepest_parab(@fun, [0.9, 0.7], 10e-8, 300);

[parab_order, parab_order_sq] = get_order(xall, it);
disp(parab_order);
disp(parab_order_sq);

[x, it, xall, resall, fall, stall, ninc, ndec, total_iter] = steepest_backtrace(@fun, [0.9, 0.7], 10e-8, 300);

[back_order, back_order_sq] = get_order(xall, it);
disp(back_order);
disp(back_order_sq);

function [order, order_sq] = get_order(xall, it)
    x_min = [1, 1/sqrt(2)];
    order = [];
    order_sq = [];
    for i = 1:it
        order = [order; norm([xall(i+1, 1), xall(i+1, 2)] - x_min) / norm([xall(i, 1), xall(i, 2)] - x_min)];
        order_sq = [order_sq; norm([xall(i+1, 1), xall(i+1, 2)] - x_min) / norm([xall(i, 1), xall(i+1, 2)] - x_min)^2];
    end
    
end

function get_graphs()
    [f_min, ~, ~] = fun([1, 1/sqrt(2)]);
    [x, new_it, new_xall, new_resall, fall, ~] = newton_mod(@fun, [0.9, 0.7], 10e-8, 50);
    new_error = fall - f_min;
    
    [x, para_it, para_xall, para_resall, fall, stall, ninc, ndec] = steepest_parab(@fun, [0.9, 0.7], 10e-8, 300);
    para_error = fall - f_min;


    [x, back_it, back_xall, back_resall, fall, stall, ninc, ndec, total_iter] = steepest_backtrace(@fun, [0.9, 0.7], 10e-8, 300);
    back_error = fall - f_min;
    
    sz = 14;
    h = semilogy(1:new_it+1, new_error, 'r-o', 1:para_it+1, para_error, 'g-x', 1:back_it+1, back_error, 'b-*');
    %axis([6 320 1e-10 1e-1]);
    %axis tight; ax = axis; axis([ax(1)*0.8 ax(2)*1.2 ax(3)*0.8 ax(4)*1.2]);
    title("Function error vs iteration")
    hax = gca;
    hlx = xlabel('iteration'); set(hlx, 'FontSize', sz);
    %px = get(hlx, 'Position'); set(hlx, 'Position', [px(1) 2*px(2) px(3)])
    hly = ylabel('function error -- log scale'); set(hly, 'FontSize', sz);
    %ht = title('convergence study'); set(ht, 'FontSize', sz);
    hl = legend('Newton', 'SD no backtrack', 'SD with backtrack', ...
                'Location', 'SouthEast');
    set(hl, 'FontSize', sz);
    
    figure();
    sz = 14;
    h2 = semilogy(1:new_it+1, new_resall, 'r-o', 1:para_it+1, para_resall, 'g-x', 1:back_it+1, back_resall, 'b-*');
    %axis([6 320 1e-10 1e-1]);
    %axis tight; ax = axis; axis([ax(1)*0.8 ax(2)*1.2 ax(3)*0.8 ax(4)*1.2]);
    title("Gradient norm vs iteration")
    hax = gca;
    hlx = xlabel('iteration'); set(hlx, 'FontSize', sz);
    %px = get(hlx, 'Position'); set(hlx, 'Position', [px(1) 2*px(2) px(3)])
    hly = ylabel('gradient norm -- log scale'); set(hly, 'FontSize', sz);
    %ht = title('convergence study'); set(ht, 'FontSize', sz);
    hl = legend('Newton', 'SD no backtrack', 'SD with backtrack', ...
                'Location', 'SouthEast');
    set(hl, 'FontSize', sz);

end

function get_graphs_hess_mod()
    [f_min, ~, ~] = fun([1, 1/sqrt(2)]);


    [x, new_it, xall, new_resall, fall, spdall] = newton_nomod(@fun, [10, 1.5], 10e-8, 50);
    new_error = fall - f_min;

    [x, para_it, xall, para_resall, fall, spdall] = newton_mod(@fun, [10, 1.5], 10e-8, 50);
    para_error = fall - f_min;

    [x, back_it, xall, back_resall, fall, spdall] = newton_mod_delta(@fun, [10, 1.5], 10e-8, 50);
    back_error = fall - f_min;

    sz = 14;
    h = semilogy(1:new_it+1, new_error, 'r-o', 1:para_it+1, para_error, 'g-x', 1:back_it+1, back_error, 'b-*');
    %axis([6 320 1e-10 1e-1]);
    %axis tight; ax = axis; axis([ax(1)*0.8 ax(2)*1.2 ax(3)*0.8 ax(4)*1.2]);
    title("Function error vs iteration")
    hax = gca;
    hlx = xlabel('iteration'); set(hlx, 'FontSize', sz);
    %px = get(hlx, 'Position'); set(hlx, 'Position', [px(1) 2*px(2) px(3)])
    hly = ylabel('function error -- log scale'); set(hly, 'FontSize', sz);
    %ht = title('convergence study'); set(ht, 'FontSize', sz);
    hl = legend('Newton no mod', 'Newton delta = 1', 'Newton small delta', ...
                'Location', 'SouthEast');
    set(hl, 'FontSize', sz);
    
    figure();
    sz = 14;
    h2 = semilogy(1:new_it+1, new_resall, 'r-o', 1:para_it+1, para_resall, 'g-x', 1:back_it+1, back_resall, 'b-*');
    %axis([6 320 1e-10 1e-1]);
    %axis tight; ax = axis; axis([ax(1)*0.8 ax(2)*1.2 ax(3)*0.8 ax(4)*1.2]);
    title("Gradient norm vs iteration")
    hax = gca;
    hlx = xlabel('iteration'); set(hlx, 'FontSize', sz);
    %px = get(hlx, 'Position'); set(hlx, 'Position', [px(1) 2*px(2) px(3)])
    hly = ylabel('gradient norm -- log scale'); set(hly, 'FontSize', sz);
    %ht = title('convergence study'); set(ht, 'FontSize', sz);
    hl = legend('Newton no mod', 'Newton delta = 1', 'Newton small delta', ...
                'Location', 'SouthEast');
    set(hl, 'FontSize', sz);

end



function disp_output()

    [x, it, xall, resall, fall, spdall] = newton_nomod(@fun, [-10, 10], 10e-8, 50);
    fprintf('iter, minimizer, function, gradient_norms, spd? \n');
    for i = 1:it+1
        fprintf('%3d, [%6.2e %6.2e], %6.2e, %6.2e, %d\n', i, xall(i, 1), xall(i, 2), fall(i), resall(i), spdall(i));
    end
    
    disp("===============================");
    [x, it, xall, resall, fall, spdall] = newton_nomod(@fun, [10, 1.5], 10e-8, 50);
     fprintf('iter, minimizer, function, gradient_norms, spd? \n');
    for i = 1:it+1
        fprintf('%3d, [%6.2e %6.2e], %6.2e, %6.2e, %d \n', i, xall(i, 1), xall(i, 2), fall(i), resall(i), spdall(i));
    end
    disp("mod+++++++++++++++++++")
    [x, it, xall, resall, fall, spdall] = newton_mod(@fun, [10, 1.5], 10e-8, 50);
    fprintf('iter, minimizer, function, gradient_norms, spd? \n');
    for i = 1:it+1
        fprintf('%3d, [%6.2e %6.2e], %6.2e, %6.2e, %d\n', i, xall(i, 1), xall(i, 2), fall(i), resall(i), spdall(i));
    end
    
    disp("mod+++++++++++++++++++")
    [x, it, xall, resall, fall, spdall] = newton_mod_delta(@fun, [10, 1.5], 10e-8, 50);
    fprintf('iter, minimizer, function, gradient_norms, spd? \n');
    for i = 1:it+1
        fprintf('%3d, [%6.2e %6.2e], %6.2e, %6.2e, %d\n', i, xall(i, 1), xall(i, 2), fall(i), resall(i), spdall(i));
    end
    
    disp("===============================");
    [x, it, xall, resall, fall, spdall] = newton_mod(@fun, [-1, -1], 10e-8, 50);
     fprintf('iter, minimizer, function, gradient_norms, spd? \n');
    for i = 1:it+1
        fprintf('%3d, [%6.2e %6.2e], %6.2e, %6.2e, %d \n', i, xall(i, 1), xall(i, 2), fall(i), resall(i), spdall(i));
    end
    
    disp("===============================");
    [x, it, xall, resall, fall, spdall] = newton_mod(@fun, [0.9, 0.7], 10e-8, 50);
     fprintf('iter, minimizer, function, gradient_norms, spd? \n');
    for i = 1:it+1
        fprintf('%3d, [%6.2e %6.2e], %6.2e, %6.2e, %d \n', i, xall(i, 1), xall(i, 2), fall(i), resall(i), spdall(i));
    end
    
    
    disp("______________________________");
    [x, it, xall, resall, fall, stall, ninc, ndec] = steepest_parab(@fun, [0.9, 0.7], 10e-8, 300);
    
    fprintf('iter, minimizer, function, gradient_norms, ninc, ndec \n');
    for i = 1:it+1
        if i <= 10 || mod(i, 10) == 0 || i == 33
            fprintf('%3d, [%6.2e %6.2e], %6.2e, %6.2e, %3d, %3d \n', i, xall(i, 1), xall(i, 2), fall(i), resall(i), ninc(i), ndec(i));
        end
    end
    
    fprintf("%4d\n", it);

    disp("+++++++++++++++++++++++++++++++");
    [x, it, xall, resall, fall, stall, ninc, ndec, total_iter] = steepest_backtrace(@fun, [0.9, 0.7], 10e-8, 300);
    fprintf('iter, minimizer, function, gradient_norms, ninc, ndec \n');
    for i = 1:it+1
        if i <= 10 || mod(i, 10) == 0 || i == 64
            fprintf('%3d, [%6.2e %6.2e], %6.2e, %6.2e, %3d, %3d \n', i, xall(i, 1), xall(i, 2), fall(i), resall(i), ninc(i), ndec(i));
        end
    end
    fprintf("%4d\n", total_iter);
    fprintf("%4d\n", it);
    fprintf("%6.2e\n", (total_iter - it)/it);
end


function [inv] = get_inverse(hess)
    inv = [hess(2,2), -hess(1,2); -hess(2,1), hess(1,1)] ./ (hess(2,2) * hess(1,1) - hess(1,2) * hess(2,1));
end

function [val, grad, hess] = fun(x)
    val = (x(1) - 1)^2 + 2 * (2 * x(2)^2 - x(1))^2;
    grad = transpose([(2 * (x(1) - 1) - 4 * (2 * x(2)^2 - x(1))), (16 * x(2) * (2*x(2)^2 - x(1)))]);
    hess = [(6), (-16 * x(2)); (-16 * x(2)), (96 * x(2)^2 - 16 * x(1))];
end


function [new_hess] = update_hess(hess)
    [V, D] = eig(hess);
    new_diag_vec = zeros(1,2);
    for i = 1:2
        if D(i,i) < 0
            new_diag_vec(i) = 1 - D(i,i); %sqrt(eps) - D(i,i);
        end
    end
   
    new_diag = diag(new_diag_vec); 
    update_diag = V * new_diag * V.';
    new_hess = hess + update_diag;
end



% newton's method no mod ==================================================
function [x, it, xall, resall, fall, spdall] = newton_nomod(fun, a, tol, maxit)
    it = -1;
    x = a;
        
    [val, grad, hess] = fun(x);

    [~, flag] = chol(hess);
    if flag == 0
        spdall = 1;
    else
        spdall = 0;
    end

    hess_inv = get_inverse(hess);
    
    n_grad = norm(grad);
    
   
    xall = x;
    
    resall = n_grad;
    fall = val;

    for i = 1:maxit
        it = it + 1;
        x_new = x - (hess_inv * grad).';

        if(n_grad < tol)
            return;
        end
        
        xall = [xall; x_new];
        x = x_new;

        [val, grad, hess] = fun(x);
        [~, flag] = chol(hess);
        if flag == 0
            spdall = [spdall; 1];
        else
            spdall = [spdall; 0];
        end

        hess_inv = get_inverse(hess);
        n_grad = norm(grad);

        resall = [resall; n_grad];
        fall = [fall; val];
    end
end

% newton's method with mod ================================================
function [x, it, xall, resall, fall, spdall] = newton_mod(fun, a, tol, maxit)
    it = -1;
    x = a;
        
    [val, grad, hess] = fun(x);

    [~, flag] = chol(hess);
    if flag == 0
        spdall = 1;
    else
        spdall = 0;

        % modify the hessian to be spd
        [new_hess] = update_hess(hess);
        hess = new_hess;
    end

    hess_inv = get_inverse(hess);
    n_grad = norm(grad);

    xall = x;
    
    resall = n_grad;
    fall = val;

    for i = 1:maxit
        it = it + 1;
        update = (hess_inv * grad).';
        if flag == 0
            stepsize = 1;
        else
            % find a suitable stepsize
            [stepsize, ~, ~, ~] = find_step(fun, x, 1, val, grad, -update, 100, 0);
        end
        x_new = x - stepsize * update;

        if(n_grad < tol)
            return;
        end
        
        xall = [xall; x_new];
        x = x_new;

        [val, grad, hess] = fun(x);
        [~, flag] = chol(hess);
        if flag == 0
            spdall = [spdall; 1];
        else
            spdall = [spdall; 0];
    
            % modify the hessian to be spd
            [new_hess] = update_hess(hess);
            hess = new_hess;
        end

        hess_inv = get_inverse(hess);
        n_grad = norm(grad);

        resall = [resall; n_grad];
        fall = [fall; val];
    end
end


function [new_hess] = update_hess_small(hess)
    [V, D] = eig(hess);
    new_diag_vec = zeros(1,2);
    for i = 1:2
        if D(i,i) < 0
            new_diag_vec(i) = sqrt(eps) - D(i,i);
        end
    end
   
    new_diag = diag(new_diag_vec); 
    update_diag = V * new_diag * V.';
    new_hess = hess + update_diag;
end

% newton's method with mod small delta================================================
function [x, it, xall, resall, fall, spdall] = newton_mod_delta(fun, a, tol, maxit)
    it = -1;
    x = a;
        
    [val, grad, hess] = fun(x);

    [~, flag] = chol(hess);
    if flag == 0
        spdall = 1;
    else
        spdall = 0;

        % modify the hessian to be spd
        [new_hess] = update_hess_small(hess);
        hess = new_hess;
    end

    hess_inv = get_inverse(hess);
    n_grad = norm(grad);

    xall = x;
    
    resall = n_grad;
    fall = val;

    for i = 1:maxit
        it = it + 1;
        
        update = (hess_inv * grad).';

        if flag == 0
            stepsize = 1;
        else
            % find a suitable stepsize
            [stepsize, ~, ~, ~] = find_step(fun, x, 1, val, grad, -update, 100, 0);
        end


        x_new = x - stepsize * update;

        if(n_grad < tol)
            return;
        end
        
        xall = [xall; x_new];
        x = x_new;

        [val, grad, hess] = fun(x);
        [~, flag] = chol(hess);
        if flag == 0
            spdall = [spdall; 1];
        else
            spdall = [spdall; 0];
    
            % modify the hessian to be spd
            [new_hess] = update_hess_small(hess);
            hess = new_hess;
        end

        hess_inv = get_inverse(hess);
        n_grad = norm(grad);

        resall = [resall; n_grad];
        fall = [fall; val];
    end
end

% steepest descent no backtrack ===========================================
function [x, it, xall, resall, fall, stall, ninc, ndec] = ...
steepest_parab(fun, a, tol, maxit)
    it = -1;
    x = a;

    [val, grad, hess] = fun(x);
    n_grad = norm(grad);

    step = n_grad^2 / dot(grad, (hess * grad));

    xall = x;
    resall = n_grad;
    fall = val;
    stall = step;
    ninc = 0;
    ndec = 0;

    for i = 1:maxit
        it = it + 1;
        x_new = x - step .* grad.';

        if(n_grad < tol)
            return;
        end

        x = x_new;

        [val, grad, hess] = fun(x);
        n_grad = norm(grad);
        
        step = n_grad^2 / dot(grad, (hess * grad));
        
        xall = [xall; x];
        resall = [resall; n_grad];
        fall = [fall; val];
        stall = [stall; step];
        ninc = [ninc; 0];
        ndec = [ndec; 0];
    end
end

% backtrack ==============================================================
function [step, total_iter, num_inc, num_dec] = find_step(fun, x, init_step, val, grad, descent, maxit, total_iter)
    c1 = 0.001;
    c2 = 0.9;
    if init_step == 0
        step = 1;
    else
        step = init_step;
    end

    d = dot(grad, descent);
    a = 0;
    b = Inf;
    num_inc = 0;
    num_dec = 0;
    
    for i = 1:maxit
        total_iter = total_iter + 1;
        [new_val, new_grad, ~] = fun(x + (step .* descent));
        if new_val > val + c1 * step * d  % no sufficient decrease
            b = step;
            step = (a + b) / 2;
            num_dec = num_dec + 1;
        elseif dot(new_grad, descent) < c2 * d
            a = step;
            if b == Inf
                step = 2 * a;
            else
                step = (a + b) / 2;
            end
            num_inc = num_inc + 1;
        else
            break;
        end
    end
end

% steepest descent with backtrack =========================================
function [x, it, xall, resall, fall, stall, ninc, ndec, total_iter] = ...
steepest_backtrace(fun, a, tol, maxit)
    it = -1;
    total_iter = -1;
    x = a;

    [val, grad, ~] = fun(x);
    n_grad = norm(grad);
   

    [step, total_iter, num_inc, num_dec] = find_step(fun, x, 0, val, grad, -grad.', 100, total_iter);
    
    %disp(step);

    prev_step = step;
    xall = x;
    resall = n_grad;
    fall = val;
    stall = step;
    ninc = num_inc;
    ndec = num_dec;

    for i = 1:maxit
        it = it + 1;
        total_iter = total_iter + 1;

        x_new = x - step .* grad.';
        if(n_grad < tol)
            return;
        end
        x = x_new;

        prev_norm_sq = n_grad^2;

        [val, grad, ~] = fun(x);
        n_grad = norm(grad);

        init_step = prev_step * (prev_norm_sq / (n_grad ^ 2));
        
        [step, total_iter, num_inc, num_dec] = find_step(fun, x, init_step, val, grad, -grad.', 100, total_iter);
        
        prev_step = step;
        xall = [xall; x];
        resall = [resall; n_grad];
        fall = [fall; val];
        stall = [stall; step];
        ninc = [ninc; num_inc];
        ndec = [ndec; num_dec];
    end
end




