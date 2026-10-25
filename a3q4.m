
display_graphs();

function display_graphs()
    get_graphs([-10, 10]);

    get_graphs([10, 1.5]);
    
    get_graphs([-1, -1]);
    
    get_graphs([0.9, 0.7]);

end


function get_graphs(x_init)
    [f_min, ~, ~] = fun([1, 1/sqrt(2)]);
    [x, trust_it, xall, trust_resall, fall, Dall, rhoall, dlall, tauall] = ...
        trustreg(@fun, x_init, 1e-7, 80, 0.1, 1, 0.1);
    trust_error = fall - f_min;
    
    fprintf('it x y fun norm area rho dl tau\n');
    for i = 1:trust_it
        fprintf('%3d %11.7f %11.7f %10.2e %9.2e %9.2e %7.4f %2d %7.4f', ...
            i, xall(i, 1), xall(i, 2), fall(i), trust_resall(i), Dall(i), rhoall(i), dlall(i), tauall(i));
        fprintf("\n");
    end


    [x, conj_it, xall, conj_resall, fall, stall, ball] = ...
        cgnl(@fun, x_init, 1e-7, 80);
    conj_error = fall - f_min;
    
    fprintf('it x y fun norm alpha beta\n');
    for i = 1:conj_it
        fprintf('%3d %11.7f %11.7f %10.2e %9.2e %9.5f %10.2e %9.2e\n', ...
            i, xall(i, 1), xall(i, 2), fall(i), conj_resall(i), stall(i), ball(i));
        fprintf("\n");
    end


    disp(length(trust_error));
    disp(trust_it);
    disp(length(conj_error));
    disp(conj_it);


    figure();
    sz = 14;
    h = semilogy(1:trust_it, trust_error, 'r-o', 1:conj_it+1, conj_error, 'g-x');
    %axis([6 320 1e-10 1e-1]);
    %axis tight; ax = axis; axis([ax(1)*0.8 ax(2)*1.2 ax(3)*0.8 ax(4)*1.2]);
    title("Function error vs iteration")
    hax = gca;
    hlx = xlabel('iteration'); set(hlx, 'FontSize', sz);
    %px = get(hlx, 'Position'); set(hlx, 'Position', [px(1) 2*px(2) px(3)])
    hly = ylabel('function error -- log scale'); set(hly, 'FontSize', sz);
    %ht = title('convergence study'); set(ht, 'FontSize', sz);
    hl = legend('Trust', 'Conjugate', ...
                'Location', 'SouthEast');
    set(hl, 'FontSize', sz);
    
    figure();
    sz = 14;
    h = semilogy(1:trust_it, trust_resall, 'r-o', 1:conj_it+1, conj_resall, 'g-x');
    %axis([6 320 1e-10 1e-1]);
    %axis tight; ax = axis; axis([ax(1)*0.8 ax(2)*1.2 ax(3)*0.8 ax(4)*1.2]);
    title("Gradient norm vs iteration")
    hax = gca;
    hlx = xlabel('iteration'); set(hlx, 'FontSize', sz);
    %px = get(hlx, 'Position'); set(hlx, 'Position', [px(1) 2*px(2) px(3)])
    hly = ylabel('gradient norm -- log scale'); set(hly, 'FontSize', sz);
    %ht = title('convergence study'); set(ht, 'FontSize', sz);
    hl = legend('Trust', 'Conjugate', ...
                'Location', 'SouthEast');
    set(hl, 'FontSize', sz);
    
end




function [val, grad, hess] = fun(x)
    val = (x(1) - 1)^2 + 2 * (2 * x(2)^2 - x(1))^2;
    grad = transpose([(2 * (x(1) - 1) - 4 * (2 * x(2)^2 - x(1))), (16 * x(2) * (2*x(2)^2 - x(1)))]);
    hess = [(6), (-16 * x(2)); (-16 * x(2)), (96 * x(2)^2 - 16 * x(1))];
end



function [val_m] = model(fun_val, fun_grad, fun_hess, dir)
    val_m = fun_val + dot(fun_grad, dir) + 0.5 * dot(dir, fun_hess * dir);
end


function [rho] = check_model(f_base, f_dir, m_base, m_dir)
    rho = (f_base - f_dir) / (m_base - m_dir);
end

function [root1, root2] = quad(a ,b, c)
    disc = b^2 - 4 * a * c; 
    if disc < 0
        root1 = Inf;
        root2 = Inf;
        return
    end
    sq_disc = sqrt(b^2 - 4 * a * c );
    root1 = (-b - sq_disc)/(2*a);
    root2 = (-b + sq_disc)/(2*a);

end

% pray that this is not wrong
function [root1, root2, val_max] = dist_quad(p_U, p_B, Delta_sq)
    U_sq = dot(p_U, p_U);
    B_U_diff = p_B - p_U;
    U_BU = dot(B_U_diff, p_U);
    diff_sq = dot(B_U_diff, B_U_diff);
    coeffs = [(U_sq - 2 * U_BU + diff_sq - Delta_sq), (2 * (U_BU - diff_sq)), diff_sq];
    [root1, root2] = quad(coeffs(1), coeffs(2), coeffs(3));
end

function [val] = dogleg_quad(p_U, p_B, tau)
    val = p_U + (tau-1) * (p_B - p_U); 
end

function [val] = dogleg_quad_test(p_U, p_B, tau)
    U_sq = dot(p_U, p_U);
    B_U_diff = p_B - p_U;
    U_BU = dot(B_U_diff, p_U);
    diff_sq = dot(B_U_diff, B_U_diff);
    val = (U_sq - 2 * U_BU + diff_sq) + tau * (2 * (U_BU - diff_sq)) + tau^2 * diff_sq;
end

function [p_cauchy, tau_cauchy, dl_cauchy] = cauchy_point(g, g_norm, g_hess_dot, Delta)
  
    if(g_hess_dot < 0)
        tau_cauchy = 1;
        p_cauchy = -Delta/g_norm * g;
        dl_cauchy = -3;
        return;
    end
    strange_ratio = g_norm^3/(Delta * g_hess_dot);


    if (strange_ratio >= 1)
        tau_cauchy = 1;
        p_cauchy = -Delta/g_norm * g;
        dl_cauchy = -2;
        return;
    else
        tau_cauchy = strange_ratio;
        p_cauchy = -tau_cauchy * Delta/g_norm * g;
        dl_cauchy = -1;
        return;
    end

end

function [p, tau, dl] = dogleg(fun, x, Delta, fx, grad_f, hess_f)
    g_norm = norm(grad_f);
    g_dot = g_norm^2;
    g_hess_dot = dot(grad_f, hess_f * grad_f);
    p_U = -g_dot/g_hess_dot * grad_f;
    p_B = -inv(hess_f) * grad_f;
    
    % check if p_B is in the trust region
    if norm(p_B) <= Delta
        p = p_B;
        tau = 2;
        dl = 1;
    else
        % calculate the quadratic
        [root1, root2] = dist_quad(p_U, p_B, Delta^2);
        if(max(root1, root2) == Inf)
            [p, tau, dl] = cauchy_point(grad_f, g_norm, g_hess_dot, Delta);
        elseif(root1 <= 2 && root1 >= 1 && root2 <= 2 && root2 >= 1)
            % both are in the region
            if(root2 > root1)
                p = dogleg_quad(p_U, p_B, root2);
                tau = root2;
                dl = 2;
            else
                p = dogleg_quad(p_U, p_B, root1);
                tau = root1;
                dl = 3;
            end

        elseif(root1 <= 2 && root1 >= 1 && (root2 > 2 || root2 < 1))
            p = dogleg_quad(p_U, p_B, root1);
            tau = root1;
            dl = 3;
        elseif((root1 > 2 || root1 < 1) && root2 <= 2 && root2 >= 1)
            p = dogleg_quad(p_U, p_B, root2);
            tau = root2;
            dl = 2;
        else
            % neither are in region get cauchy point
            [p, tau, dl] = cauchy_point(grad_f, g_norm, g_hess_dot, Delta);
        end
    end

    
end

% it, xall, resall, fall, Dall, rhoall, dlall, tauall
function [x, it, xall, resall, fall, Dall, rhoall, dlall, tauall] = ...
    trustreg(fun, x0, tol, maxit, Delta0, Dmax, eta)
    
    delta_k = Delta0;
    x_k = x0;
    xall = x0;
    [val, grad, hess] = fun(x_k);
    n_grad = norm(grad);
    resall = n_grad;
    fall = val;
    Dall = Delta0;


    [p_k, tau_k, dl_k] = dogleg(fun, x_k, delta_k, val, grad, hess);
    [val_dir,~,~] = fun(x_k + p_k.');

    [m_dir] = model(val, grad, hess, p_k);

    [rho] = check_model(val, val_dir, val, m_dir);
    
    norm_dir = norm(p_k);

    rhoall = rho;
    dlall = dl_k;
    tauall = tau_k;

    for i = 1:maxit

        if(n_grad < tol)
            x = x_k;
            it = i;
            return;
        end
        
        if rho > 0.75 && norm_dir == delta_k
            delta_k = min(2 * delta_k, Dmax);
            x_k = x_k + p_k.';

        elseif 0.25 <= rho
            %delta_k stays constant
            x_k = x_k + p_k.';

        elseif (eta < rho) && (rho < 0.25)
            delta_k = 0.25 * delta_k;
            x_k = x_k + p_k.';
        else
            delta_k = 0.25 * delta_k;
            %x_k stays constant
        end

        % find min direction in specified area
        % calculuate cauchy point
        xall = [xall; x_k];
        [val, grad, hess] = fun(x_k);
        n_grad = norm(grad);
        resall = [resall; n_grad];
        fall = [fall; val];
        Dall = [Dall; delta_k];
        
        
        [p_k, tau_k, dl_k] = dogleg(fun, x_k, delta_k, val, grad, hess);
        [val_dir,~,~] = fun(x_k + p_k.');
        
        [m_dir] = model(val, grad, hess, p_k);
        
        [rho] = check_model(val, val_dir, val, m_dir);
        
        norm_dir = norm(p_k);
        
        rhoall = [rhoall; rho];
        dlall = [dlall; dl_k];
        tauall = [tauall; tau_k];

    end
    x = x_k;
    it = i;
end


%_______________________________________________________ conjugate gradient

function [x, it, xall, resall, fall, stall, ball] = ...
    cgnl(fun, x0, tol, maxit)
    x_k = x0.';
    xall = x_k.';
    [val, grad, hess] = fun(x_k.');
    g_dot = dot(grad, grad);
    prod = g_dot;
    resall = sqrt(prod);
    fall = val;
    ball = [];
    stall = [];
    
    if(sqrt(prod) < tol)
        x = x_k.';
        it = 0;
        return;
    end
    p_k = -grad;
    prod0 = prod;

    for i = 0:maxit
        
        g_hess_dot = dot(grad, hess * grad);
        alpha_k = g_dot/g_hess_dot;
        
        stall = [stall; alpha_k];

        x_k = x_k + alpha_k*p_k;
        [val, new_grad, hess] = fun(x_k.');

        prodnew = dot(new_grad, new_grad);
        
        if(sqrt(prodnew) <= tol)
            disp(sqrt(prodnew));

            x = x_k.';
            it = i;
            return;
        end

        
        PR_dot = -dot(new_grad, (new_grad - grad));
        beta_PR = PR_dot/g_dot;
        new_beta = min(0.75, max(beta_PR, 0));  

        ball = [ball; new_beta];

        grad = new_grad;
        g_dot = prodnew;

        p_k = -new_grad + new_beta * p_k;

        % find min direction in specified area
        % calculuate cauchy point
        xall = [xall; x_k.'];
        resall = [resall; sqrt(g_dot)];
        fall = [fall; val];   

    end
    x = x_k.';
    it = i+1;
end













