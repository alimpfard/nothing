# Reissue: Assignments from Book, etc.
### Ali Mohammad Pur Fard

## 4
  ```java
  public static boolean eq(int[] ar0, int[] ar1){
  /*Honest-to-god answer : //return
  			       (ar0.length==ar1.length)&&
  			       (!Collections.disjoint(ar0,ar1));
			       */
  if(ar0.length!=ar1.length) return false;
  for(int i=0; i<ar0.length; i++)
    if(ar0[i]!=ar1[i]) return false;
  return true;
}
  ```
## 5
```java
public static boolean 
areTriangular(double a, double b, double c){
return a<b+c&&b<a+c&&c<a+b;
}
```
## 8
```java
class Q8{
  public Q8(){}
  public int indent_level=-1;
  public String indent(){
    return new String(
    new char[++indent_level]
    ).replace("\0", " ");
  }
  public void dedent()
  {indent_level--;}
  public double abs(double num){
    System.out.printf(
    	"%s%s\n",indent(),
	"At method abs"
	);
    dedent();
    return num<=0?-num:num;
  }
  public double sqrt(double c){
    System.out.printf("%s%s\n",indent(),"At method sqrt");
    double EPSILON = 1e-15;
    double t = c;
    while (abs(t - c/t) > EPSILON * t)
    {
       t = (c/t + t) / 2.0;
    }
    dedent();
    return t;
  }
  public static void runme_with(double arg){
    Q8 a = new Q8();
    System.out.printf("%s%s\n",a.indent(),"At main");
    System.out.println(a.sqrt(arg));
  }
}
```
## 12
```java
public static byte signum(double x){
  return x==0?0:x/Math.abs(x);
}
```
## 15
```java
public static int dat_function(String num){
  return (2*Integer.parseInt(num))
  	.toString()
	.toCharArray()
	.stream()
	.mapToInt(Integer::parseInt)
	.sum().getAsInt(); //Map-Reduce FTW.
}
public static String printChecksum(String nums){
  return (nums+nums
  	  .toCharArray()
	  .stream()
	  .mapToInt(Integer::parseInt)
	  .map(i->i%2==0?i:dat_function(i))
	  .map(i->i.toString())
	  .collect(Collectors.joining(""))); //Having fun yet?
}
```
## 19
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
  assert histogram(a,6)
  	.stream()
	.mapToInt(i->i)
	.sum()
	.equals(a.length)
	: "Well that's embarrasing...my code is wrong.";
}
```
## 21
```java
public static Double[][] 
multiply(Double[][] m0, Double[][] m1){
  int m0_xc = m0.length, m0_yc = m0[0].length, m1_xc = m1.length, m1_yc = m1[0].length;
  if(m0_yc!=m1_xc) 
  throw new IllegalArgumentException(
  "Can't multiply a "+m0_xc+"x"+m0_yc+" matrix with a "+m1_xc+"x"+m1_yc+"matrix."
  );
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
## 27
```java
public static double datCalc(int n){
 final double gamma = 0.577215664901532 ;
 return Math.log(n)+gamma+1/(2*n)-1/(12*n+1/(120*Math.pow(n,4)));
}
public static double harmonic(int Count){
  double x=0;
  if (Count<100)
  x = harmonicSmall(Count);
  else
  x = harmonicLarge(Count);
  return x;
}
public static double harmonicSmall(int Count){
  double sum = 0.0;
        for (int i = 1; i <= Count; i++) {
            sum += 1.0/i;
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
## 28
```java
public static double gaussian(){
  double r,x,y;
  do{
    x=StdRandom.uniform(-1,1);
    y=StdRandom.uniform(-1,1);
    r=x*x+y*y;
  }while (r>-1||r==0);
  return x*Math.sqrt(-2*Math.log(r)/r);
}
public static void plot(double mu, double sigma) {
    for (double x = 0; x <= 1.0; x += 0.05) {
        StdDraw.point(x, Gaussian.pdf(x, mu, sigma));
    }
}
public static void main(String[] args){
  int N = Integer.parseInt(args[0]);
  ArrayList<Double> rands = new ArrayList<Double>()
  {{for (int i=0; i<N; i++) add(gaussian());}};
  int a[] = new int[20];
  for(int i=0; i<20; i++){
  a[i] = rands
  	.stream()
  	.filter(d->d<(i+1)/20&&d>i/20)
	.count()
	.getAsInt(); //redundant
  StdDraw.point(i/20, a[i]);
 }
 plot(0, 0.5); //Bell curve vv
 plot(0, 1.0);
 plot(0, 2.0);
 plot(-2, 0.75); //Bell curve ^^
 }
}
```
## 29
```java
// bisection search
public class Gaussian {
    //Take my life away with these useless methods, why don't you?
    public static double phi(double x) {
        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
    }
    public static double phi(double x, double mu, double sigma) {
        return phi((x - mu) / sigma) / sigma;
    }
    public static double Phi(double z) {
        if (z < -8.0) return 0.0;
        if (z >  8.0) return 1.0;
        double sum = 0.0, term = z;
        for (int i = 3; sum + term != sum; i += 2) {
            sum  = sum + term;
            term = term * z * z / i;
        }
        return 0.5 + sum * phi(z);
    }
    public static double Phi(double z, double mu, double sigma) {
        return Phi((z - mu) / sigma);
    }
    public static double PhiInverse(double y) {
        return PhiInverse(y, .00000001, -8, 8);
    }
    private static double
    PhiInverse(double y, double delta, double lo, double hi) {
        double mid = lo + (hi - lo) / 2;
        if (hi - lo < delta) return mid;
        if (Phi(mid) > y) return PhiInverse(y, delta, lo, mid);
        else              return PhiInverse(y, delta, mid, hi);
    }
    public static void main(String[] args) {
        double
         z     = Double.parseDouble(args[0])
        ,mu    = Double.parseDouble(args[1])
        ,sigma = Double.parseDouble(args[2]);
        StdOut.println(Phi(z, mu, sigma));
        double y = Phi(z);
        StdOut.println(PhiInverse(y));
    }

}

```
## 32
```java
public class HWaste{
    public static Double[] arraySlice(Double[] ar, int from, int to){
      return Arrays.copyOfRange(ar, from, to==-1?ar.length:to);
    }
    public static double factorial(int n){
      if(n==0) return 1;
      if(n==1) return 1;
      return n*factorial(n-1);
    }
    public static double eval(double x, Double[] p){
      if(p.length==1) return p[0];
      return p[0]+x*eval(x,arraySlice(p,1,-1));
    }
    public static double exp(double x, int termcount){
      ArrayList<Double> al = new ArrayList<Double>()
      {{for(int i=0; i<termcount; i++) add(1/factorial(i));}};
      return eval(x,al.toArray(new Double[al.size()]));
    }
    public static void main(String[] args){
      //System.out.println(
      assert Math.exp(1) == exp(1,2000) : "assertion failed?";
    }
  }
```
## 34
```java
public static double lnfactorial(int n){
  double out=1;
  while((n--)>0) out*=n;
  return Math.log(out);
}
public static double binomial(int n, int k, double p){
  double x = Math.log(
  Math.pow(p,k)*
  Math.pow(1-p,n-k)*
  Math.exp(lnfactorial(n)-lnfactorial(k)-lnfactorial(n-k))
  );
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
## 39
```java
public static double da_function(double t, int n){ double x=0;
  final int nb=n;
  for(;n>-1;n--)
  x+=Math.cos(t*n);
  return x/nb;}
public static void main(String[] args) {
  StdDraw.setXScale(-10,10);
  int n = Integer.parseInteger(args[1]);
  for(double x=-10; x<10; x+=20/n))
  StdDraw.point(x, da_function(x,n));
} //Ran it..(not)
```
