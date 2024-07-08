Sets
    s /1*10/;

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
    x1(s), x2(s), x3(s), x4(s)
    y1(s), y2(s), y3(s), y4(s)
    obj;

Positive Variables
    x1(s), x2(s), x3(s), x4(s)
    y1(s), y2(s), y3(s), y4(s);

Parameters
    lambda1(s), lambda2(s), lambda3(s), mu1(s), mu2(s);

* Lagrange Multiplikatoren initialisieren
lambda1(s) = 0;
lambda2(s) = 0;
lambda3(s) = 0;
mu1(s) = 0;
mu2(s) = 0;

Equations
    obj_def
    constr1(s)
    constr2(s)
    constr3(s)
    constr4(s)
    constr5(s);

obj_def ..
    obj =e= sum(s, ps(s)*(2*x1(s) + 2*y1(s) - 7*x2(s) + 24*y2(s) + 4*x3(s) + y3(s) - 2*x4(s) + 8*y4(s))
                + lambda1(s) * (20 - x1(s) - x2(s) - 2*x3(s) - 5*x4(s))
                + lambda2(s) * (5 - 2*x1(s) - 5*x2(s) + 7*x3(s) - 11*x4(s))
                + lambda3(s) * (-4 + x1(s) - 6*x2(s) - x3(s) + 5*x4(s))
                + mu1(s) * (h1(s) + 5*x1(s) - 2*x2(s) + 4*x3(s) - 3*x4(s) + 2*y1(s) - 4*y2(s) - y3(s) - 2*y4(s))
                + mu2(s) * (h2(s) - 3*x1(s) - x2(s) - 8*x3(s) - 8*x4(s) - y1(s) - y2(s) - y3(s) + 5*y4(s)));

constr1(s) ..
    x1(s) + x2(s) + 2*x3(s) + 5*x4(s) =e= 20;

constr2(s) ..
    2*x1(s) + 5*x2(s) - 7*x3(s) + 11*x4(s) =e= 5;

constr3(s) ..
    -x1(s) + 6*x2(s) + x3(s) - 5*x4(s) =e= -4;

constr4(s) ..
    -5*x1(s) + 2*x2(s) - 4*x3(s) + 3*x4(s) - 2*y1(s) + 4*y2(s) + y3(s) + 2*y4(s) =e= h1(s);

constr5(s) ..
    3*x1(s) + x2(s) + 8*x3(s) + 8*x4(s) + y1(s) + y2(s) + y3(s) - 5*y4(s) =e= h2(s);

Model SOPmodelLagrange /all/;

Solve SOPmodelLagrange using lp minimizing obj;

Display x1.l, x2.l, x3.l, x4.l, y1.l, y2.l, y3.l, y4.l, obj.l;
