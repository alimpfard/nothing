# An Introduction to Programming in Java \-\- whatever the rest was
### Ali Mohammad Pur Fard -- anothertest
## Assignments 2.1.\*
### 2.1.4  Write a static method eq()
that takes two int arrays as arguments and returns true if the arrays have the same length and all corresponding pairs of of elements are equal, and false otherwise.

```java
public static boolean eq(int[] ar0, int[] ar1){
  /*Honest-to-god answer : *///return (ar0.length==ar1.length)&&(!Collections.disjoint(ar0,ar1));
  if(ar0.length!=ar1.length) return false;
  for(int i=0; i<ar0.length; i++)
    if(ar0[i]!=ar1[i]) return false;
  return true;
}
```

### 2.1.5 Write a static method areTriangular()
that takes three double arguments and returns true if they could be the sides of a triangle (none of them is greater than or equal to the sum of the other two).

```java
public static boolean areTriangular(double a, double b, double c){
  return a<b+c&&b<a+c&&c<a+b;
}
```

### 2.1.8 + Give the function-call trace for `java Harmonic 3 5`

**Why even ask for this?????**

```java
3 :
        at Harmonic.harmonic(Harmonic.java:9)
        at Harmonic.main(Harmonic.java:16)
5 :
        at Harmonic.harmonic(Harmonic.java:9)
        at Harmonic.main(Harmonic.java:16)
```

### 2.1.12 Consider the static method duplicate() below.

```java
public static String duplicate(String s)
{
   String t = s + s;
   return t;
}
```

What does the following code fragment do?
```java
String s = "Hello";
s = duplicate(s);
String t = "Bye";
t = duplicate(duplicate(duplicate(t)));
StdOut.println(s + t);
```

+ #### *prints -> "HelloHelloByeByeByeByeByeByeByeBye\\n"*

###  2.1.15 Given two stars with angles of declination and right ascension (d1, a1) and (d2, a2),
the angle they subtend is given by the formula

2 arcsin((sin2(d/2) + cos (d1)cos(d2)sin2(a/2))1/2)

where a1 and a2 are angles between − 180 and 180 degrees, d1 and d2 are angles between − 90 and 90 degrees, a = a2 − a1, and d = d2 − d1. Write a program to take the declination and right ascension of two stars as command-line arguments and print the angle they subtend. Hint : Be careful about converting from degrees to radians.

```java
class IDEK{
  public static double deg2rad(double deg){
    return deg*Math.PI/180;
  }
  public static double calculateDatFormula(double d0, double d1, double a0, double a1){
    return 2*Math.arcsin((Math.sin(deg2rad(d1-d0)/2))*(Math.sin(deg2rad(d1-d0)/2))
            + Math.cos(deg2rad(d0))*Math.cos(deg2rad(d1))*Math.sin(deg2rad(a1-a0)/2)*Math.sin(deg2rad(a1-a0)/2));
  }
  public static void Main(String[] args){
    if(!args.length>3) return;
    System.out.println(calculateDatFormula(Double.parseDouble(arg[0]),Double.parseDouble(args[2]), Double.parseDouble(args[1]),Double.parseDouble(args[3])));
  }
}
```

### 2.1.19 Write a static method histogram()
that takes an int array a[] and an integer m as arguments and returns an array of length m whose ith element is the number of times the integer i appeared in a[]. Assuming the values in a[] are all between 0 and m-1, the sum of the values in the returned array should equal a.length.

```java
public static int occ_num(int[] arr, int i){
  int k=0;
  for(int j : arr) k+=(j==i)?1:0;
  return k;
}
public static int[] histogram(int[] arr, int m){
  int[] ret = new int[m];
  for(int i=0; i<m; i++)
    ret[i] = occ_num(arr,i);
    return ret;
}
public static void main(String[] args) {
  int[] a = new int[10]{1,2,1,3,4,4,5,3,2,4};
  assert histogram(a,6).stream().mapToInt(i->i).sum().equals(a.length) : "Well that's embarrasing...my code is wrong.";
}
```

### 2.1.21 Write a static method multiply()
that takes two square matrices of the same dimension as arguments and produces their product (another square matrix of that same dimension). Extra credit : Make your program work whenever the number of columns in the first matrix is equal to the number of rows in the second matrix.

```java
public static Double[][] multiply(Double[][] m0, Double[][] m1){
  int m0_xc = m0.length, m0_yc = m0[0].length, m1_xc = m1.length, m1_yc = m1[0].length;
  if(m0_yc!=m1_xc) throw new IllegalArgumentException("Can't multiply a "+m0_xc+"x"+m0_yc+" matrix with a "+m1_xc+"x"+m1_yc+"matrix.");
  Double[][] m2 = new Double[m0_xc][m1_yc];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                m2[i][j] = 0.00000;
        for (int i = 0; i < m0_xc; i++)
            for (int j = 0; j < m1_yc; j++)
                for (int k = 0; k < m0_yc; k++)
                    m2[i][j] += m0[i][k] * m1[k][j];
        return m2;
}
```

### 2.1.27 Harmonic numbers.
Write a program Harmonic that contains three static methods harmoinc(), harmoincSmall(), and harmonicLarge() for computing the harmonic numbers. The harmonicSmall() method should just compute the sum (as in Program 1.3.5), the harmonicLarge() method should use the approximation Hn = loge(n ) + γ + 1/(2n ) − 1/(12n + 1/(120n4) (the number γ = 0.577215664901532... is known as Euler’ constant), and the harmonic() method should call harmonicSmall() for n < 100 and harmonicLarge() otherwise.

```java
public static double datCalc(int n){
	  final double gamma = 0.577215664901532 ;
	  return Math.log(n)+gamma+1/(2*n)-1/(12*n+1/(120*Math.pow(n,4)));
	}
	public static double harmonic(int Count){
	  double x=0;
	  if(Count<100)
	  {x = harmonicSmall(Count);}
	  else
	  {x = harmonicLarge(Count);}
	  return x;
	}
	public static double harmonicSmall(int Count){
	  double sum = 0.0;
	        for (int i = 1; i <= Count; i++) {
	            sum += 1.0 / i;
	        }
	        return sum;
	}
	public static double harmonicLarge(int Count){
	  double ret = 0;
	  for (; Count>0; Count--) {
	    ret += datCalc(Count);
	  }
	  return ret;
	}

```

### 2.1.28 Black–Scholes option valuation.
The Black–Scholes formula supplies the theoretical value of a European call option on a stock that pays no dividends, given the current stock price s, the exercise price x, the continuously compounded risk-free interest rate r, the volatility σ, and the time (in years) to maturity t. The Black–Scholes value is given by the formula s Φ(a) − xe−rtΦ(b),where Φ(z) is the Gaussian cumulative distribution function, Image, and Image. Write a program that takes s, r, σ, and t from the command line and prints the Black–Scholes value.

```java
//Bleh
public static double BlackScholes(double S, double X, double T, double r, double v) //no respect for camelCase
{
    double d1, d2;
    d1=(Math.log(S/X)+(r+v*v/2)*T)/(v*Math.sqrt(T));
    d2=d1-v*Math.sqrt(T);
    return S*CND(d1)-X*Math.exp(-r*T)*CND(d2);
}

// The cumulative normal distribution function
public static double CND(double X) //maybe some PascalCase will fix me right up?
{
    double L, K, w ;
    double a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937, a4 = -1.821255978, a5 = 1.330274429;
    L = Math.abs(X);
    K = 1.0 / (1.0 + 0.2316419 * L);
    w = 1.0 - 1.0 / Math.sqrt(2.0 * Math.PI) * Math.exp(-L *L / 2) * (a1 * K + a2 * K *K + a3
            * Math.pow(K,3) + a4 * Math.pow(K,4) + a5 * Math.pow(K,5));
    if (X < 0.0)
        w= 1.0 - w;
    return w;
}
public static void main(String[] args) {
  System.out.println(BlackScholes(Double.parseDouble(args[0]),Double.parseDouble(args[1]),Double.parseDouble(args[2]),Double.parseDouble(args[3]),Double.parseDouble(args[4])));
}
```

### 2.1.29 Fourier spikes.
Write a program that takes a command-line argument n and plots the function

(cos(t) + cos(2t) + cos(3t) + ... + cos(nt)) / n

for 500 equally spaced samples of t from − 10 to 10 (in radians).

```java
  public static double da_function(double t, int n){ double x=0;
    final int nb=n;
    for(;n>-1;n--)
    x+=Math.cos(t*n);
    return x/n;}
  public static void main(String[] args) {
    StdDraw.setXScale(-10,10);
    int n = Integer.parseInteger(args[1]);
    for(double x=-10; x<10; x+=20/n))
    StdDraw.point(x, da_function(x,n));
  }
```

### 2.1.32 Chords.
Develop a version of PlayThatTune that can handle songs with chords (including harmonics). Develop an input format that allows you to specify different durations for each chord and different amplitude weights for each note within a chord. Create test files that exercise your program with various chords and harmonics, and create a version of Für Elise that uses them.

+ Do what now? what is this 'PlayThatTune'?

### 2.1.34 Binomial distribution.
Write a function
```java
public static double binomial(int n, int k, double p)
```
to compute the probability of obtaining exactly k heads in n biased coin flips (heads with probability p) using the formula

f (n, k, p) = pk(1-p)n-k n! / (k!(n-k)!)

Hint : To stave off overflow, compute x = ln f (n, k, p) and then return ex. In main(), take n and p from the command line and check that the sum over all values of k between 0 and n is (approximately) 1. Also, compare every value computed with the normal approximation

f (n, k, p) ≈ ϕ (np, np(1-p))

```java
public static double lnfactorial(int n){
  double out=1;
  while((n--)>0) out*=n;
  return Math.log(out);
}
public static double binomial(int n, int k, double p){
  double x = Math.log(Math.pow(p,k)*Math.pow(1-p,n-k)*Math.exp(lnfactorial(n)-lnfactorial(k)-lnfactorial(n-k)));
  return Math.exp(x);
}
public static void main(String[] args) {
  double x=0, p=Double.parseDouble(args[1]);
  int n=Integer.parseInt(args[0]);
  for(int k=0;k<=n;k++)
  x+=binomial(n,k,p);
  assert x==1 : "Well that went wrong...";
}
```
### There is no \#39
