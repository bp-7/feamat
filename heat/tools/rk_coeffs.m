function coeffs = rk_coeffs(method)

if (strcmp('BE',method))
    alpha = zeros(1);
    alpha(1,1) = 1;
    
    b = zeros(1,1);
    b(1) = 1;  
    
    coeffs.alpha = alpha;
    coeffs.b = b;
    coeffs.nstages = size(alpha,1);
    
elseif (strcmp('rk2_mid',method))
    alpha = zeros(1);
    alpha(1,1) = 0.5;
    
    b = zeros(1,1);
    b(1) = 1;  
    
    coeffs.alpha = alpha;
    coeffs.b = b;
    coeffs.nstages = size(alpha,1);

elseif (strcmp('rk2_trap',method))
    alpha = zeros(2);
    alpha(2,1) = 0.5;
    alpha(2,2) = 0.5;
    
    b = zeros(2,1);
    b(1) = 0.5;
    b(2) = 0.5;
    
    coeffs.alpha = alpha;
    coeffs.b = b;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('sdirk4',method))
    alpha = zeros(5);
    alpha(1,1) = 0.25;
    
    alpha(2,1) = 0.5;
    alpha(2,2) = 0.25;
    
    alpha(3,1) = 0.34;
    alpha(3,2) = -0.04;
    alpha(3,3) = 0.25;
    
    alpha(4,1) = 0.2727941176470588;
    alpha(4,2) = -0.0503676470588235;
    alpha(4,3) = 0.0275735294117647;
    alpha(4,4) = 0.25;
    
    alpha(5,1) = 1.0416666666666667;
    alpha(5,2) = -1.0208333333333333;
    alpha(5,3) = 7.8125;
    alpha(5,4) = -7.0833333333333333;
    alpha(5,5) = 0.25;

    
    b = zeros(5,1);
    b(1) = 1.0416666666666667;
    b(2) = -1.0208333333333333;
    b(3) = 7.8125;
    b(4) = -7.0833333333333333;
    b(5) = 0.25;
    
    coeffs.alpha = alpha;
    coeffs.b = b;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('skdir4_s4',method))
    alpha = zeros(4);
    alpha(1,1) = 0.097961082941;
   
    alpha(2,1) = 0.262318069183;
    alpha(2,2) = 0.097961082941;
    
    alpha(3,1) = 0.230169419019;
    alpha(3,2) = 0.294466719347;
    alpha(3,3) = 0.097961082941;
    
    alpha(4,1) = 0.210562684389;
    alpha(4,2) = 0.269382888280;
    alpha(4,3) = 0.307008634881;
    alpha(4,4) = 0.097961082941;
    
    b = zeros(4,1);
    b(1) = 0.222119403264;
    b(2) = 0.282060762166;
    b(3) = 0.236881213175;
    b(4) = 0.258938621395;
    
    coeffs.alpha = alpha;
    coeffs.b = b;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('sdirk3',method))
    p = (3-sqrt(3))/6;
    alpha = zeros(2);
    alpha(1,1) = p;
    
    alpha(2,1) = 1-2*p;
    alpha(2,2) = p;
    
    b = zeros(2,1);
    b(1) = 0.5;
    b(2) = 0.5;
    
    coeffs.alpha = alpha;
    coeffs.b = b;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('rk3_crouzeix',method))
    alpha = zeros(2);
    alpha(1,1) = 0.7886751345948128;
    
    alpha(2,1) = -0.5773502691896257;
    alpha(2,2) = 0.7886751345948128;
    
    b = zeros(2,1);
    b(1) = 0.5;
    b(2) = 0.5;
    
    coeffs.alpha = alpha;
    coeffs.b = b;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('rk4_crouzeix',method))
    x = 1.13715804260325770;
    alpha = zeros(3);
    alpha(1,1) = 1.0685790213016288;
    
    alpha(2,1) = -0.5685790213016288;
    alpha(2,2) = 1.0685790213016288;
    
    alpha(3,1) = 2.1371580426032577;
    alpha(3,2) = -3.2743160852065154;
    alpha(3,3) = 1.0685790213016288;

    b = zeros(3,1);
    b(1) = 0.128886400515720395;
    b(2) = 0.742227198968559154;
    b(3) = 0.128886400515720395;
    
    coeffs.alpha = alpha;
    coeffs.b = b;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('sdirk3_2',method))
    x = 0.4358665215;
    alpha = zeros(3);
    alpha(1,1) = x;
    
    alpha(2,1) = 0.5*(1-x);
    alpha(2,2) = x;
    
    alpha(3,1) = -1.5*x*x + 4*x - 0.25;
    alpha(3,2) = 1.5*x*x - 5*x + 1.25;
    alpha(3,3) = x;
    
    b = zeros(3,1);
    b(1) = -1.5*x*x + 4*x - 0.25;
    b(2) = 1.5*x*x - 5*x + 1.25;
    b(3) = x;
    
    coeffs.alpha = alpha;
    coeffs.b = b;
    coeffs.nstages = size(alpha,1);

    
elseif (strcmp('sdirk3_3',method))
    alpha = zeros(3);
    
    alpha(1,1) = 1;

    alpha(2,1) = 1;
    alpha(2,2) = 1;
    
    alpha(3,1) = 0.5;
    alpha(3,2) = -0.5;
    alpha(3,3) = 1;
    
    b = zeros(3,1);
    b(1) = 0.5;
    b(2) = -0.5;
    b(3) = 1;
    
    coeffs.alpha = alpha;
    coeffs.b = b;
    coeffs.nstages = size(alpha,1);
else
    
    error('Runge-Kutta method not implemented');
end

return