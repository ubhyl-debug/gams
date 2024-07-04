* Define the optimization problem
Set
   i /1*5/;

* Define the variables
Variables
   x(i)        'decision variables'
   obj         'objective function';

* Set the objective function direction (minimize)
Positive Variables x;

* Define the objective function
Equation
   objfunc     'objective function'
   eq1         'constraint 1'
   eq2         'constraint 2'
   eq3         'constraint 3'
   ub          'upper bounds on x';

objfunc ..   obj =e= -x('1') + x('2') + 10*x('3') + 5*x('4') + 6*x('5');

eq1 ..      -x('1') + x('2') + 5*x('3') + 2*x('4') + 4*x('5') =e= 5;

eq2 ..       5*x('1') - 2*x('2') + 9*x('3') + x('4') - 2*x('5') =e= -3;

eq3 ..       2*x('1') + 2*x('2') + x('3') + 5*x('4') - 3*x('5') =e= 9;

ub(i) ..     x(i) =l= 10;

* Specify the model
Model lpModel /all/;

* Solve the model using CPLEX
Solve lpModel using lp minimizing obj;

* Display the results
Display x.l, obj.l;
