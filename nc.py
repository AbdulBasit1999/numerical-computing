from Equation import Expression
import os
from sympy import*
from numpy import*
from math import *
from mpmath.calculus.odes import ODEMethods
from sympy.solvers.ode.ode import ode_1st_exact
FIXEDDIGITS = 7
def equation(streq,var,val):
    fn = Expression(streq,[var])
    return fn(float(val))

def bisection():
    eq = input("ENTER EQUATION: ")
    a = float(input("ENTER LOWER BOUND : "))
    b = float(input("ENTER UPPER BOUND : "))
    tol = float(input("ENTER TOLERENCE VALUE : "))
    fd = int(input("ENTER NUMBER OF DECIMAL PLACES : "))
    fa = round(equation(eq,"x",a),fd)
    i = 0
    if(equation(eq,"x",a) * equation(eq,"x",b) > 0):
        print("\nBISECTION CAN NOT START\n")
        return
    else:
        print("a\t\t\tb\t\t\tc\t\t\tf(c)\n")
        while(i <= 100):
            c = round((a+b)/2,fd)
            fc = round(equation(eq,"x",c),fd)
            print(str(a) + "\t\t\t" + str(b) + "\t\t\t" + str(c) + "\t\t\t" + str(fc) + "\n")
            if(fc == 0 or (b-a)/2 < tol):
                print("\nTHE PROCEDURE HAS COMPLETED THE ROOT IS : " + str(c) + "\n")
                break
            i = i+ 1
            if(fa * fc > 0):
                a = c
                fa = fc
            else:
                b = c
        if(i == 100):
            print("\nTHE ROOT FINDING IS IMPOSSIBLE\n")
    n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))


def deri(streq,var,val):
    x = Symbol(var)
    dev = Derivative(streq,x,evaluate = True)
    return dev.subs(x,val)

def newtonmethod():
    eq = input("ENTER EQUATION: ")
    a = float(input("ENTER INITIAL APPROXIMATION : "))
    tol = float(input("ENTER TOLERENCE VALUE : "))
    fd = int(input("ENTER NUMBER OF DECIMAL PLACES : "))
    print("n\t\t\t\tp")
    i = 0
    while(i <= 100):
        print(str(i) + "\t\t\t\t" + str(a))
        p = round(a - (equation(eq,"x",a)/deri(eq,"x",a)),fd)
        if( abs(p - a) < tol):
            print(str(i+1) + "\t\t\t\t" + str(p))
            print("THE PROCEDURE WAS SUCCESSFUL ROOT IS : " + str(p))
            break
        i = i + 1
        a = p
    if(i == 100):
        print("THE PROCEDURE WAS UNSUCCESSFUL")

def secant():
    eq = input("ENTER EQUATION: ")
    p0 = float(input("ENTER INITIAL APPROXIMATION 1 : "))
    p1 = float(input("ENTER INITIAL APPROXIMATION 2 : "))
    tol = float(input("ENTER TOLERENCE VALUE : "))
    fd = int(input("ENTER NUMBER OF DECIMAL PLACES : "))
    i = 0
    print("n\t\t\t\tp\n")
    print(str(i) + "\t\t\t\t" + str(p0) + "\n")
    i = i + 1
    print(str(i) + "\t\t\t\t" + str(p1) + "\n")
    i = i + 1
    while(i<= 100):
        fp0 = equation(eq,"x",p0)
        fp1 = equation(eq,"x",p1)
        p2 = round(((p0 * fp1)-(p1 * fp0))/(fp1 - fp0) , fd)
        print(str(i) + "\t\t\t\t" + str(p2) + "\n")
        if(abs(p2 - p1) < tol):
            print("\nTHE PROCEDURE COMPLETED ROOT IS : " + str(p2))
            break
        i = i + 1
        p0 = p1
        p1 = p2
    if(i == 100):
        print("\nPROCEDURE IS NOT COMPLETED\n")

def falseposition():
    eq = input("ENTER EQUATION: ")
    a = float(input("ENTER LOWER BOUND : "))
    b = float(input("ENTER UPPER BOUND : "))
    tol = float(input("ENTER TOLERENCE VALUE : "))
    fd = int(input("ENTER NUMBER OF DECIMAL PLACES : "))
    i = 1
    fa = round(equation(eq,"x",a) , fd)
    if(equation(eq,"x",a) * equation(eq,"x",b) > 0):
        print("\nREGULA FALSI CAN NOT START\n")
        return
    else:
        print("a\t\t\tb\t\t\tc\t\t\tf(c)\n")
        while(i <= 100):
            fb = round(equation(eq,"x",b) , fd)
            c = round(((a * fb)-(b * fa))/(fb - fa),fd)
            fc = round(equation(eq,"x",c),fd)
            print(str(a) + "\t\t\t" + str(b) + "\t\t\t" + str(c) + "\t\t\t" + str(fc) + "\n")
            if(fc == 0 or abs(fc) < tol):
                print("\nTHE PROCEDURE HAS COMPLETED THE ROOT IS : " + str(c) + "\n")
                break
            i = i+ 1
            if(fa * fc > 0):
                a = c
                fa = fc
            else:
                b = c
        if(i == 100):
            print("\nTHE ROOT FINDING IS IMPOSSIBLE\n")   

def lagrange_interpol():
    eq = input("ENTER EQUATION: ")
    n = int(input("ENTER VALUE OF n : "))
    fd = int(input("ENTER DECIAML PLACES : "))
    x = []
    for i in range(n+1):
        y = float(input("\nENTER x" + str(i) + " : "))
        x.append(y)
    y = float(input("\nENTER VALUE OF X FOR APPROXIMATION : "))
    pn = 0
    for i in range(n+1):
        l0 = 1
        for j in range(n+1):
            if(i != j):
                l0 = l0 * (y - x[j])/(x[i] - x[j])
        pn = round(pn + (l0 * equation(eq,"x",x[i])),fd)  
    print(pn)

def divideddiff():
    n = int(input("ENTER TOTAL NO OF VALUES: "))
    table = [[ 0 for i in range(n + 1)]for j in range(n)]
    print("ENTER VALUE OF x : \n")
    for i in range(n):
        table[i][0] = float(input("ENTER x" + str(i) + " : "))
    for i in range(n):
        table[i][1] = float(input("ENTER fx" + str(i) + " : "))

    fd = int(input("ENTER DECIMAL PLACES : "))
    for i in range(2,n + 1): 
        for j in range(n + 1 - i):
            table[j][i] = round(((table[j][i-1] - table[j+1][i-1]) / (table[j][0] - table[i+j - 1 ][0])),fd)
    for i in range(n): 
        for j in range(n + 1 - i): 
            print(table[i][j], "\t",end = " ");
        print("")
    xn = float(input("ENTER VALUE YOU WANT TO APPROXIMATE: "))
    sum= 0
    m = 0
    for i in range(2,n+1):
        pro = 1
        for j in range(m):
            pro = pro * (xn - table[j][0])
        sum = sum + (pro * table[0][i])
        m = m + 1
    print(round(sum,fd))

def multiply(x,n):
    pro = 1
    for i in range(n):
        pro = pro * (x - i)
    return pro

def stirling():
    n = int(input("ENTER TOTAL NO OF VALUES: "))
    table = [[ 0 for i in range(n + 1)]for j in range(n)]
    print("ENTER VALUE OF x : \n")
    for i in range(n):
        table[i][0] = float(input("ENTER x" + str(i) + " : "))
    for i in range(n):
        table[i][1] = float(input("ENTER fx" + str(i) + " : "))

    fd = int(input("ENTER DECIMAL PLACES : "))
    for i in range(2, n + 1):
        for j in range(n + 1 - i):
            table[j][i] = round((table[j + 1][i - 1] -  table[j][i - 1]),fd)

    for i in range(n): 
        for j in range(n + 1 - i): 
            print(table[i][j], "\t",end = " ")
        print("")
    xn = float(input("ENTER VALUE YOU WANT TO APPROXIMATE: "))
    mp = floor(n/2)
    if(table[mp - 1][0] > xn or xn > table[mp + 1][0]):
        print("\nSTRILING'S FORMULA NOT APPLICABLE\n")
        return
    h = table[1][0] - table[0][0]
    s = (xn - table[mp][0]) / h
    p = 0
    pro = 1
    sum = 0
    for i in range(1,n + 1):
        if(i % 2 != 0):
            if(p < 3):
                sum = sum + ((table[mp][i] * pow(s,p))/factorial(i-1))
            else:
                sum = sum + ((multiply(pow(s,2),p - 1) * table[mp][i])/factorial(i - 1))
        else:
            mp -= 1
            if(p < 3):
                sum = sum + ((((table[mp + 1][i] + table[mp][i])/2) * pow(s,p))/factorial(i-1))
            else:
                sum = sum + ((multiply(pow(s,2),p - 1) /factorial(i - 1)) * ((table[mp + 1][i] + table[mp - 1][i])/2))
        p+=1 
    print("\nTHE APPROXIMATION IS : " + str(round(sum,fd))) 

def diff():
    print("\n\n1 : FORWARD DIFFERENTIATION BY VALUES OF H")
    print("\n2 : DIFFERENTIATION BY VALUES OF F(X)")
    ch = int(input("\nENTER YOUR CHOICE : "))
    if(ch == 1):
        n = int(input("ENTER NUMBER OF VALUES OF H : "))
        x = []
        for i in range(n):
            y = float(input("ENTER VALUE OF h : "))
            x.append(y)
        eq = str(input("ENTER EQUATION : "))
        x0 = float(input("ENTER VALUE OF x0 YOU WANT TO APPROXIMATE : "))
        fd = int(input("ENTER NUMBER OF DECIMAL PLACES : "))
        for i in range(n):
            a = 1
            b = ((a - equation(eq,"x",x0))/x[i])
            c = abs(x[i])/(2 * pow(x[i],2))
            print(str(round(a,fd)) + "\t\t" + str(round(b,fd)) + "\t\t" + str(round(c,fd)) + "\n")
    else:
        n = int(input("ENTER NUMBER OF VALUES OF x : "))
        x = []
        for i in range(n):
            y = float(input("ENTER VALUE OF x : "))
            x.append(y)
        fx = []
        for i in range(n):
            y = float(input("ENTER VALUE OF fx : "))
            fx.append(y)
        fd = int(input("ENTER NUMBER OF DECIMAL PLACES : "))
        h = x[1] - x[0]
        for i in range(n):
            if(i != (n-1)):
                a = ((fx[i + 1] - fx[i])/h)
                print(str(x[i]) + "\t\t" + str(fx[i]) + "\t\t" + str(round(a,fd)) + "\n")
            else:
                a = ((fx[i] - fx[i - 1])/h)
                print(str(x[i]) + "\t\t" + str(fx[i]) + "\t\t" + str(round(a,fd)) + "\n")

def threepoint():
    print("\n1 : THREE POINT DIFFERENTIATION WITH EQUATION")
    print("\n2 : THREE POINT DIFFERENTIATION WITHOUT EQUATION")
    ch = int(input("\nENTER YOUR CHOICEC : "))
    n = int(input("ENTER NUMBER OF VALUES : "))
    x = []
    for i in range(n):
        y = float(input("ENTER VALUE OF X : "))
        x.append(y)
    h = x[1] - x[0]
    if(ch == 1):
        eq = str(input("\nENTER EQUATION : "))
        fd = int(input("ENTER NUMBER OF DECIMAL PLACES : "))
        print("\n\t\tx\t\t\tf(x)\t\t\tf'(x)\n")
        for i in range(n):
            if(i < (floor(n/2) - 1)):
                ans = round((1/(2*h))*((-3 * equation(eq,"x",x[i])) + (4 * equation(eq,"x",x[i + 1])) - equation(eq,"x",x[i + 2])),fd)
            elif(i > (floor(n/2) + 1)):
                ans = round((1/(2*h))*((-3 * equation(eq,"x",x[n])) + (4 * equation(eq,"x",x[n - 1])) - equation(eq,"x",x[n - 2])),fd)
            else:
                ans = round((1/(2*h))*( equation(eq,"x",x[1]) - equation(eq,"x",x[2])),fd)
                print("\n\t\t" + str(x[i]) + "\t\t\t" + str(round(equation(eq,"x",x[i]),fd)) + "\t\t\t" + str(ans))
    else:
        fx = []
        for i in range(n):
            y = float(input("ENTER VALUE OF f(x) : "))
            fx.append(y)
        fd = int(input("ENTER NUMBER OF DECIMAL PLACES : "))
        print("\n\t\tx\t\t\tf(x)\t\t\t\tf'(x)\n")
        for i in range(n):
            if(i < (floor(n/2) - 1)):
                ans = (1/(2*h))*((-3 * fx[i]) + (4 * fx[i + 1]) - fx[i + 2])
            elif(i >= (floor(n/2) + 1)):
                ans = (1/(2*h))*((-3 * fx[n - 1]) + (4 * fx[n - 2]) - x[n - 3])
            else:
                ans = (1/(2*h))*( fx[i + 1] - fx[i - 1])
            print("\n\t\t" + str(x[i]) + "\t\t\t" + str(fx[i]) + "\t\t\t" + str(round(ans,fd)))

def fivepoint():
    print("\n1 : FIVE POINT DIFFERENTIATION WITH EQUATION")
    print("\n2 : FIVE POINT DIFFERENTIATION WITHOUT EQUATION")
    ch = int(input("\nENTER YOUR CHOICEC : "))
    n = int(input("ENTER NUMBER OF VALUES : "))
    x = []
    for i in range(n):
        y = float(input("ENTER VALUE OF X : "))
        x.append(y)
    h = x[1] - x[0]
    if(ch == 1):
        eq = str(input("\nENTER EQUATION : "))
        fd = int(input("ENTER NUMBER OF DECIMAL PLACES : "))
        print("\n\t\tx\t\t\tf(x)\t\t\tf'(x)\n")
        for i in range(n):
            if(i < (floor(n/2) - 2)):
                ans = round((1/(12*h))*((-25 * equation(eq,"x",x[i])) + (48 * equation(eq,"x",x[i + 1])) - (36 * equation(eq,"x",x[i + 2])) + (16 * equation(eq,"x",x[i + 3])) - (3 * equation(eq,"x",x[i + 4]))),fd)
            elif(i > (floor(n/2) + 2)):
                ans = round((1/(12*h))*((-25 * equation(eq,"x",x[n - 1])) + (48 * equation(eq,"x",x[n - 2])) - (36 * equation(eq,"x",x[n - 3])) + (16 * equation(eq,"x",x[n - 4])) - (3 * equation(eq,"x",x[n - 5]))),fd)
            else:
                ans = round((1/(12*h))*( equation(eq,"x",x[i - 2]) - (8 * equation(eq,"x",x[i - 1])) + (8 * equation(eq,"x",x[i + 1])) - equation(eq,"x",x[i + 2])),fd)
                print("\n\t\t" + str(x[i]) + "\t\t\t" + str(round(equation(eq,"x",x[i]),fd)) + "\t\t\t" + str(ans))
    else:
        fx = []
        for i in range(n):
            y = float(input("ENTER VALUE OF f(x) : "))
            fx.append(y)
        fd = int(input("ENTER NUMBER OF DECIMAL PLACES : "))
        print("\n\t\tx\t\t\tf(x)\t\t\t\tf'(x)\n")
        for i in range(n):
            if(i < (floor(n/2) - 1)):
                ans = (1/(12*h))*((-25 * fx[i]) + (48 * fx[i + 1]) - (36 * fx[i + 2]) + (16 * fx[i + 3]) - (3 * fx[i + 4]))
            elif(i >= (floor(n/2) + 1)):
                ans = (1/(12*h))*((-25 * fx[n - 1]) + (48 * fx[n - 2]) - (36 * fx[n - 3]) + (16 * fx[n - 4]) - (3 * fx[n - 5]))
            else:
                ans = (1/(12*h))*( fx[i - 2] - (8 * fx[i - 1]) + (8 * fx[i + 1]) - fx[i + 2])
            print("\n\t\t" + str(x[i]) + "\t\t\t" + str(fx[i]) + "\t\t\t" + str(round(ans,fd)))

def traprule():
    n = int(input("ENTER VALUE OF N : "))
    eq = str(input("ENTER EQUATION : "))
    a = float(input("ENTER LOWER BOUND A : "))
    b = float(input("ENTER UPPER BOUND B : "))
    fd = int(input("ENTER NUMBER OF DECIMAL PLACES : "))
    h = (b - a)/n
    x = []
    x.append(a)
    for i in range(n):
        y = x[i] + h
        x.append(y)
    sum = 0
    for i in range(n + 1):
        if(i == 0 or i == (n - 1)):
            sum = sum + equation(eq,"x",x[i])
        else:
            sum = sum + (2 * equation(eq,"x",x[i]))
    ans = (h/2) * sum
    print("\nAPROXIMATED VALUE IS : " + str(round(ans,fd)))

def simpsonrule():
    n = int(input("ENTER VALUE OF N : "))
    eq = str(input("ENTER EQUATION : "))
    a = float(input("ENTER LOWER BOUND A : "))
    b = float(input("ENTER UPPER BOUND B : "))
    fd = int(input("ENTER NUMBER OF DECIMAL PLACES : "))
    h = (b - a)/n
    x = []
    x.append(a)
    for i in range(n):
        y = x[i] + h
        x.append(y)
    sum = 0
    for i in range(n + 1):
        if(i == 0 or i == (n)):
            sum = sum + equation(eq,"x",x[i])
            print(" 1 " + str(equation(eq,"x",x[i])) + " , ") 
        elif(i % 2 != 0):
            sum = sum + (4 * equation(eq,"x",x[i]))
            print(" 2 " + str(equation(eq,"x",x[i])) + " , ") 
        else:
            sum = sum + (2 * equation(eq,"x",x[i]))
            print(" 3 " + str(equation(eq,"x",x[i])) + " , ")
        print("i : " + str(i) + " , " + str(sum) + "\n")
    ans = (h/3) * sum
    print("\nAPROXIMATED VALUE IS : " + str(round(ans,fd)))

def compsitmidpoint():
    n = int(input("ENTER VALUE OF N : "))
    eq = str(input("ENTER EQUATION : "))
    a = float(input("ENTER LOWER BOUND A : "))
    b = float(input("ENTER UPPER BOUND B : "))
    fd = int(input("ENTER NUMBER OF DECIMAL PLACES : "))
    h = (b - a)/(n + 2)
    x = []
    for i in range(-1,n + 1):
        y = a + ((i + 1)*h)
        x.append(y)
    sum = 0
    for i in range(n / 2):
        if(i % 2 == 0):
            sum = sum + (4 * equation(eq,"x",x[i]))
    ans = (2 * h) * sum
    print("\nAPROXIMATED VALUE IS : " + str(round(ans,fd)))

def chap2():
    ch = 0
    while(ch != 5):
        os.system("cls")
        print("1 : BISECTION METHOD\n")
        print("2 : NEWTON METHOD\n")
        print("3 : SECANT METHOD\n")
        print("4 : REGULAR FALSE POSITION\n")
        print("5 : EXIT\n")
        ch = int(input("ENTER YOUR CHOICE : "))
        if(ch == 1):
            os.system("cls")
            bisection()
        elif(ch == 2):
            os.system("cls")
            newtonmethod()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 3):
            os.system("cls")
            secant()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 4):
            os.system("cls")
            falseposition()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 5):
            break

class CALCULATOR1:
	def __init__(self):
		pass;
	def calc(self , eq , var , a , b):
		eq=eq.replace("^" , "**");
		s = var.split(" ");
		x,y = symbols(var);
		eq = eq.replace(s[0] , "x");
		eq = eq.replace(s[1] , "y");
		
		expr = eval(eq).evalf();
		val = expr.subs([(x,a) , (y,b)]).evalf();
		return val;
		
class ODE:
	
	def format_table(self , table , method):
		print("x" + "  				" +   method)
		
		for entry in self.table:
			gaps= ""
			for i in range(0, 30-len(str(entry[0]))):
					gaps = gaps + " ";
			print(str(entry[0]) + gaps + str(entry[1]));
		
	
	def __init__(self):
		self.eq = "";
		self.vars = "";
		self.table = [];
		self.h=0.0;
		self.at = 0.0;
		self.x0 = 0.0;
		self.y0 = 0.0;
		self.equation= CALCULATOR1();
		
		
		pass
	
	def setvals(self):
		self.table = [];
		if(self.eq ==""):
			
			self.vars = input("variables = ");
			self.eq = input("equation = ");
			self.h = float(eval(input("h = ")));
			self.at =	float(eval(input("value at = ")));
			self.x0 = float(eval(input("x0 = ")));
			self.y0 = float(eval(input("y0 = ")));
		else:
			x = input("do you want to use last equation : ");
			if(x == "no" or x == "NO" or x == "n" or x =="N" ):
				self.vars = input("variables = ");
				self.eq = input("equation = ");
				self.h = float(eval(input("h = ")));
				self.at =	float(eval(input("value at = ")));
				self.x0 = float(eval(input("x0 = ")));
				self.y0 = float(eval(input("y0 = ")));
			else:
				self.h = float(eval(input("h = ")));
				self.at =	float(eval(input("value at = ")));
				self.x0 = float(eval(input("x0 = ")));
				self.y0 = float(eval(input("y0 = ")));
				pass;	
	def midpoint(self):
		self.setvals();
		print("func = " + self.eq)
		i = 0;
		self.table = []
		while(round(self.x0 ,10) != round(self.at , 10)):
			arr= [];
			val = self.equation.calc(self.eq , self.vars , self.x0 +(self.h/2) , self.y0 + ((self.h*(self.equation.calc(self.eq , self.vars , self.x0 , self.y0)))/2 ) 		);
			self.y0 = self.y0 + (self.h * val);
			self.x0 = self.x0+self.h;
			arr.append(round(self.x0 , FIXEDDIGITS));
			arr.append(round(self.y0, FIXEDDIGITS) );
			self.table.append(arr);
			
			pass;
		
		self.format_table(self.table , "midpoint");	
		print("\nat   " + str(self.table[len(self.table)-1][0] )   + "   y` is   "  +str(self.table[len(self.table)-1][1] ))
		
				
	def euler(self):
		
		self.setvals();
		
		print("func = " + self.eq);
		i = 0;
		
		while(round(self.x0 ,10) != round(self.at , 10)):
			arr= [];
			val = self.equation.calc(self.eq , self.vars , self.x0 , self.y0);
			self.y0 = self.y0 + (self.h * val);
			self.x0 = self.x0+self.h;
			arr.append(round(self.x0 , FIXEDDIGITS ));
			arr.append(round(self.y0 , FIXEDDIGITS));
			self.table.append(arr);
			pass;
		
		self.format_table(self.table , "euler");	
		print("\nat   " + str(self.table[len(self.table)-1][0] )   + "   y` is   "  +str(self.table[len(self.table)-1][1] ))
		#truev = TrueValue();
		#err = self.table[len(self.table)-1][1] - truev.calc_diff2(self.eq , self.vars , self.table[len(self.table)-1][0])
		
		
		pass;
	
	def modifiedeuler(self):
		self.setvals();
		
		print("func = " + self.eq);
		i = 0;
		
		while(round(self.x0 ,10) != round(self.at , 10)):
			arr= [];
			val = self.equation.calc(self.eq , self.vars , self.x0 , self.y0);
			self.y0 = self.y0 + ( (self.h/2) *( val + self.equation.calc(self.eq , self.vars , self.x0 + self.h , self.y0+  (self.h*	val			) 	)  					)  );
			self.x0 = self.x0+self.h;
			arr.append(round(self.x0 , FIXEDDIGITS ));
			arr.append(round(self.y0 , FIXEDDIGITS));
			self.table.append(arr);
			
			pass;
		
		self.format_table(self.table , "modified euler");	
		
		print("\nat   " + str(self.table[len(self.table)-1][0] )   + "   y` is   "  +str(self.table[len(self.table)-1][1] ))
		
		
		pass;			

def chap3():
    ch = 0
    while(ch != 4):
        os.system("cls")
        print("1 : LAGRANGE INTERPOLATION\n")
        print("2 : DIVIDED DIFFERENCE\n")
        print("3 : STIRLING FORMULA\n")
        print("4 : EXIT\n")
        ch = int(input("ENTER YOUR CHOICE : "))
        if(ch == 1):
            os.system("cls")
            lagrange_interpol()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 2):
            os.system("cls")
            divideddiff()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 3):
            os.system("cls")
            stirling()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 4):
            break

def chap4():
    ch = 0
    while(ch != 7):
        os.system("cls")
        print("1 : FORWARD BACKWARD DIFFERENTIATION\n")
        print("2 : 3 POINTS DIFFERENTIATION\n")
        print("3 : 5 POINTS DIFFERENTIATION\n")
        print("4 : TRAP'S RULE INTEGRATION\n")
        print("5 : SIMPSON'S RULE INTEGRATION\n")
        print("6 : COMPOSITE MIDPOINT RULE INTEGRATION\n")
        print("7 : EXIT\n")
        ch = int(input("ENTER YOUR CHOICE : "))
        if(ch == 1):
            os.system("cls")
            diff()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 2):
            os.system("cls")
            threepoint()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 3):
            os.system("cls")
            fivepoint()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 4):
            os.system("cls")
            traprule()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 5):
            os.system("cls")
            simpsonrule()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 6):
            os.system("cls")
            compsitmidpoint()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 7):
            break

def chap5():
    o1 = ODE()
    ch = 0
    while(ch != 4):
        os.system("cls")
        print("1 : EULER METHOD\n")
        print("2 : MODIFIED EULER\n")
        print("3 : MIDPOINT\n")
        print("4 : EXIT\n")
        ch = int(input("ENTER YOUR CHOICE : "))
        if(ch == 1):
            os.system("cls")
            o1.euler()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 2):
            os.system("cls")
            o1.modifiedeuler()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 3):
            os.system("cls")
            o1.midpoint()
            n = int(input("\n\n\t\t\t\tPRESS 1 TO EXIT"))
        elif(ch == 4):
            break

def mainlayout():
    ch = 0
    while(ch != 5):
        os.system("cls")
        print("\n\n*********************WELCOME TO NUMERICAL COMPUTING PROGRAM*************************\n\n")
        print("1 : CHAPTER 2\n")
        print("2 : CHAPTER 3\n")
        print("3 : CHAPTER 4\n")
        print("4 : CHAPTER 5\n")
        print("5 : EXIT\n")
        ch = int(input("ENTER YOUR CHOICE : "))
        if(ch == 1):
            chap2()
        elif(ch == 2):
            chap3()
        elif(ch == 3):
            chap4()
        elif(ch == 4):
            chap5()
        elif(ch == 5):
            print("\n\n...........PROGRAM DEVELOPED BY\nABDUL BASIT\nJAWWAD MATEEN\nHAMZA SIDDIQUI\nMAAZ...........")
            break

mainlayout()