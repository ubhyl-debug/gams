*Aufgabe 1 a)


Set
  i Indizes für x Variablen
    /i1, i2, i3, i4, i5/,         
  j Indizes für mu
    /j1, j2, j3/,
  iter Iterationszähler
    /1*1000/ ; 
    
Parameter
    A(j, i) /j1.i1 -1, j1.i2 1, j1.i3 5, j1.i4 2, j1.i5 4,
            j2.i1 5, j2.i2 -2, j2.i3 9, j2.i4 1, j2.i5 -2,
            j3.i1 2, j3.i2 2, j3.i3 1, j3.i4 5, j3.i5 -3/,
    b(j) /j1 5, j2 -3, j3 9/,
    c(i) /i1 -1, i2 1, i3 10, i4 5, i5 6/,
    delta   Stoppkriterium
    /1e-2/,
    f_mu      Wert der Funktion f für mu,
    delta_k   Differenz für Stoppkriterium,
    f_hat_k   Wert für f_hat,
    f_hat_k_minus_1  initialer (unendlicher) Wert für f_hat (Hier wird Wert der vorherigen It. gespeichert)
    /1e20/;

Variable
    mu(j)    Lagrange-Multiplier für k=1,
    x(i)     Enstscheidungsvaraiblen x1 bis x5,
    z        Zielfuntionswert der Lagrangerelaxierung,
    g(j)     Wert des Gradienten in Zeile j
    sub(j)   Wert des Subgradienten in Zeilej;
    

Positive Variable
    delta_x;
    
Equation
    lag_function, Lagrange Relaxierung,
    grad_mu(j)    Gradienten(zeilenweise),
    subgr(j)        Subgradienten(zeilenweise),
    con1, con2, con3, con4, con5  Ungleichungsrestriktionen;
    

*Lagrange relaxierung Definition
lag_function.. z=e= -x('i1') + x('i2') + 10*x('i3') + 5*x('i4') + 6*x('i5') 
            + mu('j1') * (-x('i1') + x('i2') + 5*x('i3') + 2*x('i4') + 4*x('i5') - b('j1'))
            + mu('j2') * (5*x('i1') - 2*x('i2') + 9*x('i3') + x('i4') - 2*x('i5') - b('j2'))
            + mu('j3') * (2*x('i1') + 2*x('i2') + x('i3') + 5*x('i4') - 3*x('i5') - b('j3'));
            
* Constraints
con1.. x('i1') =l= 10;
con2.. x('i2') =l= 10;
con3.. x('i3') =l= 10;
con4.. x('i4') =l= 10;
con5.. x('i5') =l= 10;
x.lo(i) = 0;

*Definition der Gradienten (zeilenweise)
grad_mu(j).. g(j)=e= sum(i, A(j,i)* x(i))-b(j);

*Definition der Subgradienten (zeilenweise)
*subgr(j)..sub(j)=e=      





*Beschränkungen für mu
mu.lo(j) = -100;
mu.up(j) = 100;


Model RelaxedProblem / lag_function, con1, con2, con3, con4, con5/;


*Schnittebenenverfahren

loop(iter,
*f_hat_k_minus_1

*Definiere delta_k
*def. stimmt, aber f_mu muss noch Wert zugwiesen werden deshalb rot
delta_k = f_hat_k_minus_1 - f_mu;   

*Prüfe Stoppkriterium
if (delta_k < delta, break);


)