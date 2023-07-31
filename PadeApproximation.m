## Pade approximation of a function in neighbourhood of 0

## After first run of the program we have to wait for a few 
## seconds until symbolic pagkage loads

## test-case 1: (example from lesson)
## beginning: -2
## end: 2
## function: exp(x)
## degree of numerator: 2
## degree of denominator: 2

## test-case 2:
## beginning: -2*pi
## end: 2*pi
## function: sin(x)
## degree of numerator: 3
## degree of denominator: 4

clear;
pkg load symbolic ## using symbolic package for representation of functions
syms x ## variable 

prompt = "Enter beginning of the interval: ";
try
  a = input(prompt);
catch
  a = 0;  
end_try_catch

prompt = "Enter end of the interval: ";
try
  b = input(prompt);
catch
  b = 0;  
end_try_catch

while(a>=b || a>=0 || b<=0 || not(isreal(a)) || not(isreal(b)))
  disp("Try again! Beginning of the interval must be negative number and end must be positive!")
  prompt = "Enter beginning of the interval: ";
  try
    a = input(prompt);
  catch
    a = 0;  
  end_try_catch

  prompt = "Enter end of the interval: ";
  try
    b = input(prompt);
  catch
    b = 0;  
  end_try_catch
endwhile

prompt = "Enter function to approximate as a function of variable x: ";

## needed for catching majority of wrong input function
## adds @(x) in front of the string, makes it a proper input for "str2func"
function f = to_function_handle(f1)
  f = strjoin({"@(x)",f1}, ' ');
endfunction

## please enter e^x as exp(x)
while(1)
  input_function = input(prompt, 's');
  try
    my_f = str2func(to_function_handle(input_function))(0);
    if (isnan(my_f) || isinf(my_f)) ## throw error if function is udefined at 0
      error();
    endif
    my_function = sym(input_function);
    break;
  catch
    disp("Invalid function, try again: (function must be defined and finite at 0)")
  end_try_catch
endwhile

prompt = "Input m: (degree of numerator polynomial) ";
m = input(prompt);

prompt = "Input k: (degree of denominator polynomial) ";
k = input(prompt);

N = m+k; ## index

## create McLaurin polynomial of input function
my_taylor = taylor(my_function, 'expansionPoint', 0, 'order', N+1);

disp("")
disp("Taylor expansion of input function: ")
disp(my_taylor)
disp("")

[coef, t] = coeffs(my_taylor);

## coeffitients of Taylor expansion sorted by ascending power of x
j=0;
M = length(coef);
c = sym(zeros(1,N+1));
for i=1:(N+1)
  if(M-j>0 && t(M-j)==x^(i-1))
    c(i) = coef(M-j);
    j = j+1;
  else
    c(i) = 0;
  endif
endfor 

disp("Coefficients of Taylor expansion: ")
disp(c)
disp("")
disp("I'm computing..")
disp("")

#{ 
difference from lessons -> indexing begins with 1

Principle of code:

Block matrix looking like this (for m = 1 and k = 2 -> N = 3) 
represents our linear equations when we assume b(1) = 1.

          k+2=N+1           m+2..N+1       
             |               |    |
  b(2) b(3) b(4)  a(1) a(2) a(3) a(4)  right side
|------------------------------------|  |-----|
|  0    0    0  | -1    0    0    0  |  |-c(1)|
| c(1)  0    0  |  0   -1    0    0  |  |-c(2)|
| c(2) c(1)  0  |  0    0   -1    0  |  |-c(3)|     
| c(3) c(2) c(1)|  0    0    0   -1  |  |-c(4)|
|------------------------------------|  |-----|

We also suppose that b(k+2)..b(N+1) and a(m+2)..a(N+1) 
are 0 so the coresponding columns are omitted and 
the matrix can be simplified to the following form:

  b(2) b(3)  a(1) a(2)  right side
|---------------------|  |-----|
|  0    0  | -1    0  |  |-c(1)|
| c(1)  0  |  0   -1  |  |-c(2)|
| c(2) c(1)|  0    0  |  |-c(3)|
| c(3) c(2)|  0    0  |  |-c(4)|
|---------------------|  |-----|

We might notice this matrix can have some problems with
regularity, so we need to check it in our code.

Matrix is representing system of N+1 linear equations.
Unknown variables are b(2)...b(k+1) and a(1)..a(m+1). 

Block matrix is composed of 2 matrices with following 
dimensions:

Left Matrix:  L_m - (N+1)x(k)    
                                --->  [L_m, R_m] - (N+1)x(N+1)
Right Matrix: R_m - (N+1)x(m+1)
#}

L_m = sym(zeros(N+1, k));
for i=1:(N+1)
  for j=1:min((i-1), k)
    L_m(i,j) = c(i-j);
  endfor
endfor

R_m = sym(zeros(N+1, m+1));
for i=1:(m+1)
  R_m(i,i) = -1 ;
endfor

my_matrix = [L_m, R_m];

## checking regularity of matrix
if det(my_matrix) == 0
  disp("")
  disp("Sorry, matrix is singular, I can't find solution. I recommend changing m or k.")
  disp("")
  return
endif

if N>4
  disp("I almost got it..")
  disp("")
endif


right_side = sym(zeros(N+1, 1)); 
for i=1:N+1
  right_side(i,1) = -c(i);
endfor

## solving the system
ans_coeffs = my_matrix\right_side;

P_m = sym(0); ## polynomial in numerator
P_k = sym(1); ## polynomial in denominator

for i=1:k
  P_k = P_k + ans_coeffs(i)*x^i;
endfor

for i=(k+1):(N+1)
  P_m = P_m + ans_coeffs(i)*x^(i-(k+1));
endfor

disp("Result function: ")
R_mk = P_m/P_k;
disp(R_mk)

clf('reset');
ezplot(R_mk,[a, b]);
hold on;
ezplot(my_function, [a, b]); 