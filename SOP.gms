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
    x1, x2, x3, x4
    y1(s), y2(s), y3(s), y4(s)
    obj;
    
Positive Variables
    x1, x2, x3, x4
    y1(s), y2(s), y3(s), y4(s);

Equations
    obj_def
    constr1
    constr2
    constr3
    constr4(s)
    constr5(s);

obj_def ..
    obj =e= 2*x1 - 7*x2 + 4*x3 - 2*x4 + sum(s, ps(s)*(2*y1(s) + 4*y2(s) + y3(s) + 8*y4(s)));

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

* Non-negativity constraints are implicit in GAMS for positive variables
*Positive Variables x1, x2, x3, x4, y1(s), y2(s), y3(s), y4(s);

Model SOPmodel /all/;

Solve SOPmodel using lp minimizing obj;

Display x1.l, x2.l, x3.l, x4.l, y1.l, y2.l, y3.l, y4.l, obj.l;
