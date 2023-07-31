## discrete orthogonalization of real functions (type="dis")
## continuous orthogonalization of real functions (type="con")

## After first run of the program we have to wait for a few 
## seconds until symbolic pagkage loads

## test-case 1: 
## type: con
## Number of functions: 4
## 1
## x
## x^2
## x^3
## Begining of the interval: -1
## End of the interval: 1
## Weight function: 1-x^2

## test-case 2:
## type: dis
## Number of functions: 4
## 1
## x
## x^2
## x^3
## Points: 0 1 2 3 (press enter after each point)
## Weight function: 1

## test-case 2:
## type: con
## Number of functions: 5
## 1
## sin(x)
## cos(x)
## sin(2*x)
## cos(2*x)
## Begining of the interval: -pi
## End of the interval: pi
## Weight function: 1

clear;
clear -global;
pkg load symbolic

disp("Input type of ortogonalization: (""con"" for continuous or ""dis"" for discrete)")
prompt="";

type = input(prompt, 's'); 

while not(strcmp(type, "con") || strcmp(type, "dis"))
  disp("Wrong input, please input ""con"" or ""dis"": ")
  prompt="";
  type = input(prompt, 's');
endwhile

prompt = "Number of functions: ";
global n = input(prompt);

while n<2
  prompt="Please input integer greater than 1: ";
  n = input(prompt);
endwhile

my_functions = cell(n,1);
my_ort_functions = cell(n,1);

disp("Enter all functions as functions of variable x. Functions must be linearly independent!");

## adds @(x) in front of the string, makes it a proper input for "str2func"
function f = to_function_handle(f1)
  f = strjoin({"@(x)",f1}, ' ');
endfunction

## I tried to find out whether input function is in the right format
## this is the best what I was able to do..
for i=1:n
  prompt = strjoin({"Enter function ", int2str(i), ": "},'');
  
  while(1)
    my_functions{i} = input(prompt, 's');
    try
      my_f = str2func(to_function_handle(my_functions{i}))(1); ## this line didn't throw error for "fsgrgihg" and similar "functions"
      if isnan(my_f) ## so I did it like this
        error();
      endif
      break;
    catch
      disp("Invalid function, try again: ")
    end_try_catch
  endwhile
  
endfor

global i_a;
global i_b;
global my_points;

if strcmp(type, "con")
  prompt = "Begining of interval: ";
  i_a = input(prompt);
  prompt = "End of interval: ";
  i_b = input(prompt);
else
  disp(strjoin({"Enter", int2str(n), "points: (after each point press enter)"}, ' '))
  for i=1:n
    prompt = "";
    point = input(prompt);
    my_points(i) = point;
  endfor
endif

if strcmp(type, "con")
  disp("Weight function must be continuous, real and positive on the interval (a,b) and it's integral over the interval must converge!")
else
  disp("Weight function must be positive in all points!")
endif

prompt = "Enter weight function: ";
global w_str = input(prompt, 's');

## takes 2 strings and returns string representing sum of input strings 
function f = add_functions(f1, f2) 
  f=strjoin({f1, "+(", f2, ")"}, ''); 
endfunction

## takes 2 strings and returns string representing difference of input strings
function f = sub_functions(f1, f2)
  f = strjoin({f1, "-(", f2, ")"}, '');
endfunction

## takes 2 strings and returns string representing multiple of input strings
function f = mult_functions(f1, f2)
  f = strjoin({"(", f1, ")*(", f2, ")"}, '');
endfunction

## takes constant and string and returns string representing multiple of string by constant
function f = mult(c, f1)
  f = strjoin({sprintf('%.5f',c), "*(", f1, ")"}, '');
endfunction

## calculates dot product of two functions in continuous orthogonalization
function f = cont_scalar(f1, f2)
  global i_a;
  global i_b;
  global w_str;
  
  to_integrate = str2func(to_function_handle(mult_functions(mult_functions(f1,f2),w_str)));
  f = quad(to_integrate, i_a, i_b); ##integrate over the interval
endfunction

## calculates dot product of two functions in discrete orthogonalization
function f = disc_scalar(f1, f2)
  global n;
  global w_str;
  global my_points;
  
  ## make functions from strings so we can evaluate them in points
  f1_func = str2func(to_function_handle(f1)); 
  f2_func = str2func(to_function_handle(f2));
  w = str2func(to_function_handle(w_str));
  
  f=0;
  for i=1:n ##sum
    index = my_points(i);
    f = f+f1_func(index)*f2_func(index)*w(index);
  endfor
endfunction

## first function stays the same
my_ort_functions{1} = my_functions{1};

## formula for orthogonalization
for i=2:n
  sum = "0";
  
  for j=1:(i-1)  
    ## use dot product correctly depending on type of orthogonalization
    if strcmp(type, "dis")
      scal_up = disc_scalar(my_functions{i}, my_ort_functions{j});
      scal_down = disc_scalar(my_ort_functions{j}, my_ort_functions{j});
    else
      scal_up = cont_scalar(my_functions{i}, my_ort_functions{j});
      scal_down = cont_scalar(my_ort_functions{j}, my_ort_functions{j});
    endif
    
    ## using my functions to calculate sum and multiplication by constant
    to_add = mult((scal_up/scal_down), my_ort_functions{j});
    sum = add_functions(sum, to_add);
    
  endfor
  ## using my function to calculate difference
  my_ort_functions{i} = sub_functions(my_functions{i}, sum);
end;

## using the "symbolic" package to simplify functions 
syms x

## converting my string representation to symbolic expressions
func = cell2sym(my_ort_functions); 

## displaying result
disp("")
disp("Result of calculations:")
disp("")
for i=1:n
  disp(strjoin({"Function ", int2str(i), ":"},''));
  disp(func(i))
end;

## checking for errors
for i=1:n
  for j=1:(i-1)
    ## for each pair of functions calculate dot product
    ## if the dot product is less than some treshold we declare it as a success
    ## if the dot product is more than the treshold we display it for the manual check
    if strcmp(type, "dis")
      error = disc_scalar(my_ort_functions{i}, my_ort_functions{j});
    else
      error = cont_scalar(my_ort_functions{i}, my_ort_functions{j});
    endif
      
    if error<0.00001 ## adjust treshold here
      disp(strjoin({"Function", int2str(i), "and function", int2str(j), "are OK"}, ' ')); 
    else
      disp(strjoin({"Function", int2str(i), "and function", int2str(j), "have dot product", sprintf('%.5f',error)},' '));
    endif
  endfor
endfor

