Set
  i Indizes für x Variablen
    /i1, i2, i3, i4, i5/,         
  j Indizes für mu
    /j1, j2, j3/,
  iter Iterationszähler
    /1*1000/;

Parameter
    A(j, i) /j1.i1 -1, j1.i2 1, j1.i3 5, j1.i4 2, j1.i5 4,
            j2.i1 5, j2.i2 -2, j2.i3 9, j2.i4 1, j2.i5 -2,
            j3.i1 2, j3.i2 2, j3.i3 1, j3.i4 5, j3.i5 -3/,
    b(j) /j1 5, j2 -3, j3 9/,
    c(i) /i1 -1, i2 1, i3 10, i4 5, i5 6/,
    delta /10e-2/,
    f_mu             * Wert der Funktion f für mu
    delta_k          * Differenz für Stoppkriterium
    f_hat_k          * Wert für f_hat
    f_hat_k_minus_1 /1e20/
    subparam(j)      * Parameter to store subgradient values;
* initialer (unendlicher) Wert für f_hat (Hier wird Wert der vorherigen It. gespeichert)

Variable
    mu(j)    Lagrange-Multiplier für k=1
    x(i)     Entscheidungsvariablen x1 bis x5
    z        Zielfunktionswert der Lagrangerelaxierung
    g(j)     Wert des Gradienten in Zeile j
    sub(j)   Wert des Subgradienten in Zeile j
    f_hat_k_var      * Variable for f_hat_k;

Positive Variable
    delta_x, y(j);

Equation
    lag_function    Lagrange Relaxierung
    grad_mu(j)      Gradienten (zeilenweise)
    subgr(j)        Subgradienten (zeilenweise)
    con1
    con2
    con3
    con4
    con5  Ungleichungsrestriktionen
    objfunc
    bounds;

* Lagrange Relaxierung Definition
lag_function.. z =e= -x('i1') + x('i2') + 10*x('i3') + 5*x('i4') + 6*x('i5')
            + mu('j1') * (-b('j1'))
            + mu('j2') * (-b('j2'))
            + mu('j3') * (-b('j3'))
            + mu('j1') * (-x('i1'))
            + mu('j1') * x('i2')
            + mu('j1') * 5 * x('i3')
            + mu('j1') * 2 * x('i4')
            + mu('j1') * 4 * x('i5')
            + mu('j2') * 5 * x('i1')
            + mu('j2') * ((-2) * x('i2'))
            + mu('j2') * 9 * x('i3')
            + mu('j2') * x('i4')
            + mu('j2') * (-2) * x('i5')
            + mu('j3') * 2 * x('i1')
            + mu('j3') * 2 * x('i2')
            + mu('j3') * x('i3')
            + mu('j3') * 5 * x('i4')
            + mu('j3') * (-3) * x('i5');

* Constraints
con1.. x('i1') =l= 10;
con2.. x('i2') =l= 10;
con3.. x('i3') =l= 10;
con4.. x('i4') =l= 10;
con5.. x('i5') =l= 10;
x.lo(i) = 0;

* Definition der Gradienten (zeilenweise)
*grad_mu(j).. g(j) =e= sum(i, A(j,i) * x(i)) - b(j);

* Definition der Subgradienten (zeilenweise)
subgr(j).. sub(j) =e= g(j);

* Beschränkungen für mu
mu.lo(j) = -100;
mu.up(j) = 100;

Model RelaxedProblem / lag_function, con1, con2, con3, con4, con5 /;

* Initialize mu
mu.l(j) = 0;

* Initialize f_hat_k to a very high number
f_hat_k_minus_1 = 1e20;

* Define Optimization Problem to find mu_k+1
Variable y(j);


* Set bounds for y
y.lo(j) = -100;
y.up(j) = 100;

* Define Optimization Problem to find mu_k+1
* Equation for the objective function
objfunc.. f_hat_k_var =e= f_mu + sum(j, subparam(j) * (y(j) - mu.l(j)));

Model FindMu /objfunc/;

mu.l(j) = 0;
* Calculate the subgradients using the optimal x(mu_k)
grad_mu(j).. g(j) =e= sum(i, A(j,i) * x.l(i)) - b(j);
*Schnittebenenverfahren

loop(iter,

* Calculate f_mu
    Solve RelaxedProblem using nlp minimizing z;
    f_mu = z.l;
* Calculate the subgradients using the optimal x(mu_k)
    g.l(j) = sum(i, A(j,i) * x.l(i)) - b(j);
    sub.l(j) = g.l(j);
    subparam(j) = g.l(j);
* Define delta_k
    delta_k = f_hat_k_minus_1 - f_mu;
* Check Stopping Criterion
    if (delta_k < delta, break);
* Update f_hat
    f_hat_k = min(f_hat_k_minus_1, f_mu + sum(j, g.l(j) * (mu.l(j) - mu.l(j))));
    
* Solve FindMu to find mu_k+1
    Solve FindMu using nlp maximizing f_hat_k_var;

* Update mu_k+1
    mu.l(j) = y.l(j);

* Increment k and update f_hat_k_minus_1
    f_hat_k_minus_1 = f_hat_k;
);

* Final Solution
Display mu.l, x.l, z.l;
