*Aufgabe 1 a)


Set
   
    i /1*5/          ! Indizes für x Variablen
    j /1*3/          ! Indizes für mu
    iter /1*1000/ ;  ! Iterationszähler
    
Parameter
    mu(j) /0,0,0/,          ! Lagrange-Multiplier für k=1
    delta /1e-2/,           ! Stoppkriterium
    C(j) /-100, 100/,       ! Beschränkung für mu
    f_mu,                   ! Wert der Funktion f für mu
    delta_k,                ! Differenz für Stoppkriterium
    f_hat_k                 ! Wert für f_hat
    f_hat_k_minus_1 /1e20/ ;! initialer (unendlicher) Wert für f_hat (Hier wird Wert der vorherigen It. gespeichert)
    

Variable
    x(i),                   ! Enstscheidungsvaraiblen x1 bis x5
    z;                      ! Zielfuntionswert der Lagrangerelaxierung
    
Positive Variable
    delta_x;
    
Equation
    lag_function,                   ! Lagrange Relaxierung
    con1, con2, con3, con4, con5;   ! Ungleichungsrestriktionen
    

*Lagrange relaxierung Definition
lag_function.. z=e= -x('1') + x('2') + 10*x('3') + 5*x('4') + 6*x('5') 
            + mu('1') * (-x('1') + x('2') + 5*x('3') + 2*x('4') + 4*x('5') - 5)
            + mu('2') * (5*x('1') - 2*x('2') + 9*x('3') + x('4') - 2*x('5') + 3)
            + mu('3') * (2*x('1') + 2*x('2') + x('3') + 5*x('4') - 3*x('5') - 9);
            
* Constraints
con1.. x('1') =l= 10;
con2.. x('2') =l= 10;
con3.. x('3') =l= 10;
con4.. x('4') =l= 10;
con5.. x('5') =l= 10;
x.lo(i) = 0;
    
Model RelaxedProblem /obj, con1, con2, con3, con4, con5/;

*Schnittebenenverfahren

loop(iter,

*Minimiere
solve RelaxedProblem minimizing z using lp;
f_mu = z.l;






)