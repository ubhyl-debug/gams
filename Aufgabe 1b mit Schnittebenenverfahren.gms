Sets
    s /1*10/
    iter /1*50/;
    
* Maximale Anzahl an Iterationen

Parameters
    ps(s) /
        1 0.008879787280035216
        2 0.2531272729529044
        3 0.05557756667548357
        4 0.12094171132547692
        5 0.035887726627956466
        6 0.024573168299676793
        7 0.1574063690984937
        8 0.20204437133785533
        9 0.00393791885131994
        10 0.1376241075507977 /;

Parameters
    h1(s) /
        1 -87
        2 60
        3 15
        4 -37
        5 -86
        6 46
        7 73
        8 58
        9 -59
        10 1 /;

Parameters
    h2(s) /
        1 94
        2 45
        3 -92
        4 0
        5 87
        6 78
        7 45
        8 -62
        9 34
        10 41 /;

Variables
    x1, x2, x3, x4
    y1(s), y2(s), y3(s), y4(s)
    obj
    alpha;
    
*Schnittebenenwert

Positive Variables
    x1, x2, x3, x4
    y1(s), y2(s), y3(s), y4(s);

Parameters
    lambda1, lambda2, lambda3, mu1(s), mu2(s)
    LB, UB;
    
* Untere und obere Schranke

Equations
    obj_def
    constr1
    constr2
    constr3
    constr4(s)
    constr5(s)
    cut_def;

obj_def ..
    obj =e= 2*x1 - 7*x2 + 4*x3 - 2*x4 + sum(s, ps(s)*(2*y1(s) + 24*y2(s) + y3(s) + 8*y4(s)));

constr1 ..
    x1 + x2 + 2*x3 + 5*x4 =e= 20;

constr2 ..
    2*x1 + 5*x2 - 7*x3 + 11*x4 =e= 5;

constr3 ..
    -x1 + 6*x2 + x3 - 5*x4 =e= -4;

constr4(s) ..
    -5*x1 + 2*x2 - 4*x3 + 3*x4 - 2*y1(s) + 4*y2(s) + y3(s) + 2*y4(s) =e= h1(s);

constr5(s) ..
    3*x1 + x2 + 8*x3 + 8*x4 + y1(s) + y2(s) + y3(s) - 5*y4(s) =e= h2(s);

cut_def ..
    alpha =g= obj + lambda1 * (20 - x1 - x2 - 2*x3 - 5*x4)
                 + lambda2 * (5 - 2*x1 - 5*x2 + 7*x3 - 11*x4)
                 + lambda3 * (-4 + x1 - 6*x2 - x3 + 5*x4)
                 + sum(s, mu1(s) * (h1(s) + 5*x1 - 2*x2 + 4*x3 - 3*x4 + 2*y1(s) - 4*y2(s) - y3(s) - 2*y4(s)))
                 + sum(s, mu2(s) * (h2(s) - 3*x1 - x2 - 8*x3 - 8*x4 - y1(s) - y2(s) - y3(s) + 5*y4(s)));

Model SOPmodel /obj_def, constr1, constr2, constr3, constr4, constr5/;
Model CutModel /cut_def/;

* Initialisierung
LB = -inf;
UB = +inf;

* Lagrange-Multiplikatoren initialisieren
lambda1 = 0;
lambda2 = 0;
lambda3 = 0;
mu1(s) = 0;
mu2(s) = 0;

alpha.l = +inf;

loop(iter,
    
* Lösen des aktuellen Problems
    Solve SOPmodel using lp minimizing obj;

    
* Aktualisierung der oberen und unteren Schranken
    if (obj.l < UB,
        UB = obj.l;
    );

    if (alpha.l > LB,
        LB = alpha.l;
    );

    
* Abbruchbedingung
    if (UB - LB < 1e-6,
        break;
    );

    
* Schnittebene aktualisieren
    cut_def.m = obj.l + lambda1 * (20 - x1.l - x2.l - 2*x3.l - 5*x4.l)
                 + lambda2 * (5 - 2*x1.l - 5*x2.l + 7*x3.l - 11*x4.l)
                 + lambda3 * (-4 + x1.l - 6*x2.l - x3.l + 5*x4.l)
                 + sum(s, mu1(s) * (h1(s) + 5*x1.l - 2*x2.l + 4*x3.l - 3*x4.l + 2*y1.l(s) - 4*y2.l(s) - y3.l(s) - 2*y4.l(s)))
                 + sum(s, mu2(s) * (h2(s) - 3*x1.l - x2.l - 8*x3.l - 8*x4.l - y1.l(s) - y2.l(s) - y3.l(s) + 5*y4.l(s)));

    Solve CutModel using lp minimizing alpha;

    
* Aktualisierung der Lagrange-Multiplikatoren (Subgradientenmethode)
    lambda1 = lambda1 + 0.1 * (20 - x1.l - x2.l - 2*x3.l - 5*x4.l);
    lambda2 = lambda2 + 0.1 * (5 - 2*x1.l - 5*x2.l + 7*x3.l - 11*x4.l);
    lambda3 = lambda3 + 0.1 * (-4 + x1.l - 6*x2.l - x3.l + 5*x4.l);
    mu1(s) = mu1(s) + 0.1 * (h1(s) + 5*x1.l - 2*x2.l + 4*x3.l - 3*x4.l + 2*y1.l(s) - 4*y2.l(s) - y3.l(s) - 2*y4.l(s));
    mu2(s) = mu2(s) + 0.1 * (h2(s) - 3*x1.l - x2.l - 8*x3.l - 8*x4.l - y1.l(s) - y2.l(s) - y3.l(s) + 5*y4.l(s));

    alpha.l = UB;
);

Display x1.l, x2.l, x3.l, x4.l, y1.l, y2.l, y3.l, y4.l, obj.l, LB, UB;
